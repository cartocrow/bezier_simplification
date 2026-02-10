#include "frontend.h"

#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QPushButton>
#include <QCheckBox>
#include <QSpinBox>
#include <QImageReader>
#include <QLabel>
#include <QScrollArea>
#include <QFileDialog>
#include <QProgressDialog>

#include <cartocrow/widgets/double_slider.h>

#include "library/read_graph_gdal.h"
#include "read_ipe_bezier_spline.h"

#include <cartocrow/renderer/voronoi_drawer.h>
#include <cartocrow/core/transform_helpers.h>

void saveGraphIntoTopoSet(const BaseGraph& graph, TopoSet<Inexact>& topoSet) {
    std::unordered_set<const BaseGraph::Edge*> visited;

    for (auto eit = graph.edges_begin(); eit != graph.edges_end(); ++eit) {
        if (visited.contains(&*eit)) {
            continue;
        }

        auto initial = eit;
        // find the start, if it exists.
        auto current = initial;
        while (current->source()->degree() == 2 && current->prev() != initial) {
            current = current->prev();
        }
        auto start = current;
        BaseGraph::Edge_const_handle end;
        if (current->prev() == initial) {
            end = initial;
        }
        else {
            current = initial;
            while (current->target()->degree() == 2 && current->next() != initial) {
                current = current->next();
            }
            end = current;
        }
        auto arcIndex = start->data().index;

        int arcSize = 1;

        visited.insert(&*start);

        if (start != end) {
            int iterations = 0;
            for (auto arcEit = start->next(); iterations < graph.number_of_edges(); arcEit = arcEit->next()) {
                visited.insert(&*arcEit);
                ++arcSize;
                ++iterations;
                if (arcEit == end) break;
            }
            if (iterations >= graph.number_of_edges()) {
                throw std::runtime_error("Problem in traversing graph to save it into TopoSet");
            }
        }

        int nPoints = std::max(std::min(1000.0, 31998.0 / arcSize), 5.0);

        Polyline<Inexact> arc;
        start->curve().samplePoints(isStraight(start->curve()) ? 2 : nPoints, std::back_inserter(arc));

        if (start != end) {
            for (auto arcEit = start->next(); true; arcEit = arcEit->next()) {
                arc.pop_back();
                arcEit->curve().samplePoints(isStraight(arcEit->curve()) ? 2 : nPoints, std::back_inserter(arc));
                if (arcEit == end) break;
            }
        }

        topoSet.arcs[arcIndex] = arc;
    }
}

void BezierSimplificationDemo::loadInput(const std::filesystem::path& path) {
    double prevScale = getScale();

    m_baseGraph.clear();
    m_debugEdge = std::nullopt;
    m_forcer.m_withinDistEdgeComponents.clear();
    m_forcer.m_delaunay.clear();
    m_approxGraph.clear();
    m_voronoiPainting = PaintingRenderer{};

    std::unordered_map<Point<Inexact>, BaseGraph::Vertex_handle> pToV;

    auto getVertex = [&pToV, this](const Point<Inexact>& p) {
        if (pToV.contains(p)) {
            return pToV.at(p);
        } else {
            auto newV = m_baseGraph.insert_vertex(p);
            pToV[p] = newV;
            return newV;
        }
    };

    if (path.extension() == ".ipe") {
        auto splines = ipeSplinesToIsolines(path);

        for (const auto& spline : splines) {
            if (spline.empty()) continue;
            auto lastV = getVertex(spline.curves()[0].source());
            for (const auto &curve: spline.curves()) {
                auto targetV = getVertex(curve.target());
                m_baseGraph.add_edge(lastV, targetV, curve);
                lastV = targetV;
            }
        }
    } else {
        auto [regionSet, spatialRef] = readRegionSetUsingGDAL(path);
        m_spatialRef = spatialRef;
        m_toposet = TopoSet<Inexact>(regionSet);

        std::vector<Point<Inexact>> points;
        for (const auto& arc : m_toposet.arcs) {
            std::copy(arc.vertices_begin(), arc.vertices_end(), std::back_inserter(points));
        }

        auto bbox = CGAL::bbox_2(points.begin(), points.end());

        m_transform = fitInto(bbox, Rectangle<Inexact>(0, 0, 1000, 1000));

        for (int i = 0; i < m_toposet.arcs.size(); ++i) {
            auto& arc = m_toposet.arcs[i];
            for (auto eit = arc.edges_begin(); eit != arc.edges_end(); ++eit) {
                auto sourceV = getVertex(m_transform.transform(eit->source()));
                auto targetV = getVertex(m_transform.transform(eit->target()));
                auto curve = CubicBezierCurve(sourceV->point(), targetV->point());

                // There should be no duplicates but there are somehow a couple, so we do this for now...
                bool oppositeExists = false;
                bool alreadyExists = false;
                for (auto ieit = targetV->incident_edges_begin(); ieit != targetV->incident_edges_end(); ++ieit) {
                    if ((*ieit)->target() == sourceV && (*ieit)->curve() == curve.reversed()) {
                        oppositeExists = true;
                        break;
                    }
                    if ((*ieit)->source() == sourceV && (*ieit)->curve() == curve) {
                        alreadyExists = true;
                        break;
                    }
                }
                if (!oppositeExists && !alreadyExists) {
                    auto eh = m_baseGraph.add_edge(sourceV, targetV, curve);
                    eh->data().index = i;
                }
            }
        }
    }

    m_baseGraph.orient();
    m_original = m_baseGraph;
    m_graph.reset();

    auto bbox = m_baseGraph.bbox();
    m_referencePolygon = bbox;

    m_renderer->fitInView(transform(bbox, m_transform.inverse()));

    complexity->setMaximum(m_graph.number_of_edges());
    complexity->setMinimum(m_graph.number_of_edges());
    complexity->setValue(m_graph.number_of_edges());
    complexityLog->setMaximum(log(m_graph.number_of_edges() + -1.5));
    complexityLog->setMinimum(log(m_graph.number_of_edges() + -1.5));
    complexityLog->setValue(log(m_graph.number_of_edges() + -1.5));
    desiredComplexity->setMinimum(1);
    desiredComplexity->setMaximum(m_graph.number_of_edges());
    m_minDist->setMaximum(getScale() * 30);
    m_minAdjDist->setMaximum(getScale() * 30);
    m_minComponentLength->setMaximum(getScale() * 100);
    m_minDist->setValue(getScale() / prevScale * m_minDist->value());
    m_minAdjDist->setValue(getScale() / prevScale * m_minAdjDist->value());
    m_minComponentLength->setValue(getScale() / prevScale * m_minComponentLength->value());
    updateComplexityInfo();
}

CGAL::Object transform(const CGAL::Object& o, const CGAL::Aff_transformation_2<Exact>& trans) {
    using Exact_SDG_traits = CGAL::Segment_Delaunay_graph_traits_2<Exact>;
    Segment<Exact> s;
    Line<Exact> l;
    Ray<Exact> r;
    CGAL::Parabola_segment_2<Exact_SDG_traits> parab;

    std::variant<Segment<Exact>, Line<Exact>, Ray<Exact>, CGAL::Parabola_segment_2<Exact_SDG_traits>> var;
    if (CGAL::assign(s, o)) {
        var = s.transform(trans);
    }
    if (CGAL::assign(l, o)) {
        var = l.transform(trans);
    }
    if (CGAL::assign(r, o)) {
        var = r.transform(trans);
    }
    if (CGAL::assign(parab, o)) {
        var = cartocrow::transform(parab, trans);
    }

    return var;
}

Forcer::VoronoiEdge transform(const Forcer::VoronoiEdge& vEdge, const CGAL::Aff_transformation_2<Inexact>& trans) {
    if (auto sP = std::get_if<Segment<Inexact>>(&vEdge)) {
        auto& s = *sP;
        return s.transform(trans);
    }
    if (auto psP = std::get_if<CGAL::Parabola_segment_2<Forcer::Gt>>(&vEdge)) {
        auto& ps = *psP;
        return cartocrow::transform(ps, trans);
    }
    if (auto lP = std::get_if<Line<Inexact>>(&vEdge)) {
        auto& l = *lP;
        return l.transform(trans);
    }
    if (auto rP = std::get_if<Ray<Inexact>>(&vEdge)) {
        auto& r = *rP;
        return r.transform(trans);
    }
}

void BezierSimplificationDemo::repaintVoronoi() {
    m_voronoiPainting = PaintingRenderer();
    m_voronoiPainting.setMode(GeometryRenderer::stroke);
    auto voronoiDrawer = VoronoiDrawer<MinimumDistanceForcer<std::monostate, std::monostate>::Gt>(&m_voronoiPainting);
    for (auto eit = m_forcer.m_delaunay.finite_edges_begin(); eit != m_forcer.m_delaunay.finite_edges_end(); ++eit) {
        if (m_forcer.filterVoronoiEdge(*eit)) continue;
        using Exact_SDG_traits = CGAL::Segment_Delaunay_graph_traits_2<Exact>;
        CGAL::Object o = transform(exact_primal(*eit, m_forcer.m_delaunay), pretendExact(m_transform.inverse()));

        m_voronoiPainting.setStroke(Color{ 210, 210, 210 }, 1.0);

        Segment<Exact> s;
        Line<Exact> l;
        Ray<Exact> r;
        CGAL::Parabola_segment_2<Exact_SDG_traits> parab;
        if (CGAL::assign(s, o)) {
            voronoiDrawer << s;
        }
        if (CGAL::assign(l, o)) {
            voronoiDrawer << l;
        }
        if (CGAL::assign(r, o)) {
            voronoiDrawer << r;

        }
        if (CGAL::assign(parab, o)) {
            voronoiDrawer << parab;
        }
    }

    for (const auto& comp : m_forcer.m_withinDistEdgeComponents) {
        if (m_forcer.filterComponent(comp)) continue;

        for (const auto &[vEdge, dEdge]: comp) {
            std::visit([&](const auto &geom) {
                Color red(255, 50, 50);
                m_voronoiPainting.setStroke(red, 2.0);
                voronoiDrawer << geom;
            }, transform(vEdge, m_transform.inverse()));
        }

    }
}

void BezierSimplificationDemo::updateComplexityInfo() {
    m_backup.clear();
    complexity->setMinimum(std::min((int) m_graph.number_of_edges(), complexity->minimum()));
    complexity->setValue(m_graph.number_of_edges());
    complexityLabel->setText(QString::fromStdString("#Edges: " + std::to_string(m_graph.number_of_edges())));
    complexityLog->setMinimum(std::min(log(m_graph.number_of_edges() + 0.5), complexityLog->minimum()));
    complexityLog->setValue(log(m_graph.number_of_edges() + 0.5));
}

void BezierSimplificationDemo::addIOTab() {
    auto* ioSettings = new QWidget();
    m_tabs->addTab(ioSettings, "IO");
    auto* vLayout = new QVBoxLayout(ioSettings);
    vLayout->setAlignment(Qt::AlignTop);
    auto* loadFileButton = new QPushButton("Load file");
    vLayout->addWidget(loadFileButton);

    auto* loadReferenceDataButton = new QPushButton("Load reference data");
    vLayout->addWidget(loadReferenceDataButton);

    auto* loadReferencePolygonButton = new QPushButton("Load reference polygon for image");
    vLayout->addWidget(loadReferencePolygonButton);

    auto* clearReferenceDataButton = new QPushButton("Clear reference data");
    vLayout->addWidget(clearReferenceDataButton);

    auto* exportButton = new QPushButton("Export");
    vLayout->addWidget(exportButton);

    auto* stackPolygons = new QCheckBox("Export: stacked polygons instead of holes");
    stackPolygons->setChecked(false);
    vLayout->addWidget(stackPolygons);

    m_editControlPoints = new QCheckBox("Edit control points");
    m_editControlPoints->setChecked(false);
    vLayout->addWidget(m_editControlPoints);

    auto* editAlignTangents = new QCheckBox("Edit: align tangents");
    editAlignTangents->setChecked(true);
    vLayout->addWidget(editAlignTangents);

    auto* resetEditsButton = new QPushButton("Edit: reset edits");
    vLayout->addWidget(resetEditsButton);

    connect(loadFileButton, &QPushButton::clicked, [this]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getOpenFileName(this, tr("Select input file (.ipe or .shp)"), startDir).toStdString();
        if (filePath == "") return;
        loadInput(filePath);
    });

    connect(loadReferenceDataButton, &QPushButton::clicked, [this]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getOpenFileName(this, tr("Select reference file (.tif or .shp)"), startDir).toStdString();
        if (filePath == "") return;
        if (filePath.extension() == ".tif") {
            m_referenceData.push_back(QImage(filePath.string().c_str()));
        } else {
            auto geometrySet = readGeometrySetUsingGDAL(filePath);
            m_referenceData.push_back(geometrySet);
        }
    });

    connect(clearReferenceDataButton, &QPushButton::clicked, [this]() {
        m_referenceData.clear();
    });

    connect(loadReferencePolygonButton, &QPushButton::clicked, [this]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getOpenFileName(this, tr("Select reference polygon for image (.shp)"), startDir).toStdString();
        if (filePath == "") return;
        auto regionSet = readRegionSetUsingGDAL(filePath);
        std::vector<PolygonWithHoles<Inexact>> pgns;
        regionSet.first[0].geometry.polygons_with_holes(std::back_inserter(pgns));
        m_referencePolygon = pgns[0].bbox();
    });

    connect(exportButton, &QPushButton::clicked, [this, stackPolygons]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getSaveFileName(this, tr("Select isolines"), startDir).toStdString();
        if (filePath.empty()) return;
        saveGraphIntoTopoSet(m_baseGraph, m_toposet);
        exportTopoSetUsingGDAL(filePath, m_toposet, m_transform.inverse(), m_spatialRef, stackPolygons->isChecked());
    });

    connect(m_editControlPoints, &QCheckBox::stateChanged, [this]() {
        if (m_backup.number_of_vertices() == 0) {
            m_backup = m_baseGraph;
        }
        m_editables.clear();
        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            auto sourceControlEditable = std::make_shared<ControlPoint>(eit->curve().sourceControl(), std::pair(eit, false));
            m_editables.push_back(sourceControlEditable);
            auto targetControlEditable = std::make_shared<ControlPoint>(eit->curve().targetControl(), std::pair(eit, true));
            m_editables.push_back(targetControlEditable);
        }
        for (auto vit = m_baseGraph.vertices_begin(); vit != m_baseGraph.vertices_end(); ++vit) {
            auto endpointEditable = std::make_shared<ControlPoint>(vit->point(), vit);
            m_editables.push_back(endpointEditable);
        }
        m_renderer->repaint();
    });

    connect(resetEditsButton, &QPushButton::clicked, [this]() {
        resetEdits();
        m_renderer->repaint();
    });

    connect(m_renderer, &GeometryWidget::dragStarted, [this](const Point<Inexact>& p) {
        m_dragging = nullptr;

        std::optional<std::shared_ptr<ControlPoint>> closest;
        double minDist = std::numeric_limits<double>::infinity();

        for (auto editable : m_editables) {
            auto diff = m_renderer->convertPoint(editable->point) - m_renderer->convertPoint(p.transform(m_transform));
            auto d2 = diff.x() * diff.x() + diff.y() * diff.y();
            if (d2 < 400 && d2 < minDist) {
                minDist = d2;
                closest = editable;
            }
        }

        if (closest.has_value()) {
            m_dragging = *closest;
        }
        m_renderer->repaint();
    });

    connect(m_renderer, &GeometryWidget::dragMoved, [this, editAlignTangents](const Point<Inexact>& px) {
        auto p = px.transform(m_transform);
        if (m_dragging != nullptr) {
            m_dragging->point = p;
            if (auto vhP = std::get_if<Vertex_handle>(&m_dragging->type)) {
                auto vh = *vhP;
                vh->point() = p; // update
                for (auto ieit = vh->incident_edges_begin(); ieit != vh->incident_edges_end(); ++ieit) {
                    auto eh = *ieit;
                    if (eh->source() == vh) {
                        auto& c = eh->curve();
                        auto diff = c.sourceControl() - c.source();
                        c = CubicBezierCurve(p, p + diff, c.targetControl(), c.target());

                        for (const auto& editable : m_editables) {
                            if (auto otherEdgeEditableP = std::get_if<std::pair<Edge_handle, bool>>(&editable->type)) {
                                if (otherEdgeEditableP->first == eh && !otherEdgeEditableP->second) {
                                    editable->point = p + diff;
                                    break;
                                }
                            }
                        }
                    } else {
                        assert(eh->target() == vh);
                        auto& c = eh->curve();
                        auto diff = c.targetControl() - c.target();
                        c = CubicBezierCurve(c.source(), c.sourceControl(), p + diff, p);

                        for (const auto& editable : m_editables) {
                            if (auto otherEdgeEditableP = std::get_if<std::pair<Edge_handle, bool>>(&editable->type)) {
                                if (otherEdgeEditableP->first == eh && otherEdgeEditableP->second) {
                                    editable->point = p + diff;
                                    break;
                                }
                            }
                        }
                    }
                }

            } else if (auto eiP = std::get_if<std::pair<Edge_handle, bool>>(&m_dragging->type)) {
                auto [eh, b] = *eiP;
                auto& c = eh->curve();
                c = CubicBezierCurve(c.source(), !b ? p : c.sourceControl(), b ? p : c.targetControl(), c.target());

                if (editAlignTangents->isChecked()) {
                    if (b && eh->target()->degree() == 2) {
                        auto next = eh->next();
                        auto& nc = next->curve();
                        auto diff = nc.sourceControl() - nc.source();
                        auto dist = sqrt(diff.squared_length());
                        auto vec = c.target() - c.targetControl();
                        vec /= sqrt(vec.squared_length());
                        auto newSourceControl = nc.source() + vec * dist;
                        nc = CubicBezierCurve(nc.source(), newSourceControl, nc.targetControl(), nc.target());
                        for (const auto& editable : m_editables) {
                            if (auto otherEdgeEditableP = std::get_if<std::pair<Edge_handle, bool>>(&editable->type)) {
                                if (otherEdgeEditableP->first == next && !otherEdgeEditableP->second) {
                                    editable->point = newSourceControl;
                                }
                            }
                        }
                    } else if (!b && eh->source()->degree() == 2) {
                        auto prev = eh->prev();
                        auto& pc = prev->curve();
                        auto diff = pc.target() - pc.targetControl();
                        auto dist = sqrt(diff.squared_length());
                        auto vec = c.sourceControl() - c.source();
                        vec /= sqrt(vec.squared_length());
                        auto newTargetControl = pc.target() - vec * dist;
                        pc = CubicBezierCurve(pc.source(), pc.sourceControl(), newTargetControl, pc.target());
                        for (const auto& editable : m_editables) {
                            if (auto otherEdgeEditableP = std::get_if<std::pair<Edge_handle, bool>>(&editable->type)) {
                                if (otherEdgeEditableP->first == prev && otherEdgeEditableP->second) {
                                    editable->point = newTargetControl;
                                }
                            }
                        }
                    }
                }
            }
            m_renderer->repaint();
        }
    });

    connect(m_renderer, &GeometryWidget::dragEnded, [this](const Point<Inexact>& p) {
        m_dragging = nullptr;
        m_renderer->repaint();
    });
}

void BezierSimplificationDemo::resetEdits() {
    m_editControlPoints->setChecked(false);

    if (m_backup.number_of_vertices() > 0) {
        m_baseGraph = m_backup;
    }
}

void BezierSimplificationDemo::addSimplificationTab() {
    auto* simplificationSettings = new QWidget();
    m_tabs->addTab(simplificationSettings, "Simplification");
    auto* vLayout = new QVBoxLayout(simplificationSettings);
    vLayout->setAlignment(Qt::AlignTop);

    auto* initializeButton = new QPushButton("Initialize");
    vLayout->addWidget(initializeButton);

    auto* stepButton = new QPushButton("Step");
    vLayout->addWidget(stepButton);

    auto* step10Button = new QPushButton("Step (x10)");
    vLayout->addWidget(step10Button);

    auto* step100Button = new QPushButton("Step (x100)");
    vLayout->addWidget(step100Button);

    auto* runToComplexityButton = new QPushButton("Run to specified complexity");
    vLayout->addWidget(runToComplexityButton);
    desiredComplexity = new QSpinBox();
    vLayout->addWidget(desiredComplexity);

    auto* undoButton = new QPushButton("Undo");
    vLayout->addWidget(undoButton);

    auto* undo10Button = new QPushButton("Undo (x10)");
    vLayout->addWidget(undo10Button);

    auto* undo100Button = new QPushButton("Undo (x100)");
    vLayout->addWidget(undo100Button);

    auto* redoButton = new QPushButton("Redo");
    vLayout->addWidget(redoButton);

    auto* redo10Button = new QPushButton("Redo (x10)");
    vLayout->addWidget(redo10Button);

    auto* redo100Button = new QPushButton("Redo (x100)");
    vLayout->addWidget(redo100Button);

    complexityLabel = new QLabel("#Edges: ");
    complexity = new QSlider();
    complexityLog = new DoubleSlider();
    complexity->setMaximum(m_graph.number_of_edges());
    complexity->setMinimum(m_graph.number_of_edges());
    complexity->setValue(m_graph.number_of_edges());
    complexityLog->setMaximum(log(m_graph.number_of_edges() + 0.5));
    complexityLog->setMinimum(log(m_graph.number_of_edges() + 0.5));
    complexityLog->setValue(log(m_graph.number_of_edges() + 0.5));
    desiredComplexity->setMinimum(1);
    desiredComplexity->setMaximum(m_graph.number_of_edges());
    desiredComplexity->setValue(1);

    complexity->setOrientation(Qt::Horizontal);
    vLayout->addWidget(complexityLabel);
    vLayout->addWidget(complexity);
    complexityLog->setOrientation(Qt::Horizontal);
    complexityLog->setPrecision(10000);
    vLayout->addWidget(complexityLog);

    connect(complexity, &QSlider::valueChanged, [this](int value) {
        resetEdits();
        m_graph.recallComplexity(value);
        m_renderer->repaint();
        updateComplexityInfo();
    });

    connect(complexityLog, &DoubleSlider::valueChanged, [this](double value) {
        resetEdits();
        m_graph.recallComplexity(std::exp(value));
        m_renderer->repaint();
        updateComplexityInfo();
    });

    connect(initializeButton, &QPushButton::clicked, [this]() {
        m_collapse.initialize();
    });

    connect(stepButton, &QPushButton::clicked, [this]() {
        resetEdits();
        m_collapse.step();
        m_renderer->repaint();
        updateComplexityInfo();
    });

    connect(step10Button, &QPushButton::clicked, [this]() {
        resetEdits();
        for (int i = 0; i < 10; ++i) m_collapse.step();
        m_renderer->repaint();
        updateComplexityInfo();
    });

    connect(step100Button, &QPushButton::clicked, [this]() {
        resetEdits();
        for (int i = 0; i < 100; ++i) m_collapse.step();
        m_renderer->repaint();
        updateComplexityInfo();
    });

    connect(undoButton, &QPushButton::clicked, [this]() {
        resetEdits();
        m_graph.backInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(undo10Button, &QPushButton::clicked, [this]() {
        resetEdits();
        for (int i = 0; i < 10; ++i)
            m_graph.backInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(undo100Button, &QPushButton::clicked, [this]() {
        resetEdits();
        for (int i = 0; i < 100; ++i)
            m_graph.backInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(redoButton, &QPushButton::clicked, [this]() {
        resetEdits();
        m_graph.forwardInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(redo10Button, &QPushButton::clicked, [this]() {
        resetEdits();
        for (int i = 0; i < 10; ++i)
            m_graph.forwardInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(redo100Button, &QPushButton::clicked, [this]() {
        resetEdits();
        for (int i = 0; i < 100; ++i)
            m_graph.forwardInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(runToComplexityButton, &QPushButton::clicked, [this]() {
        auto startComplexity = m_graph.number_of_edges();
        auto targetComplexity = desiredComplexity->value();
        QProgressDialog progress("Simplifying...", "Abort", 0, startComplexity - targetComplexity, this);
        progress.setWindowModality(Qt::WindowModal);
        progress.setMinimumDuration(1000);
        m_collapse.runToComplexity(desiredComplexity->value(), [&progress, startComplexity, this](int i) {
                                       progress.setValue(startComplexity - i);
                                       if (i % 100 == 0)
                                           m_renderer->repaint();
                                   },
                                   [&progress]() { return progress.wasCanceled();
                                   });
        updateComplexityInfo();
        m_renderer->repaint();
    });
}

void BezierSimplificationDemo::addMinimumDistanceTab() {
    auto* minimumDistanceSettings = new QWidget();
    m_tabs->addTab(minimumDistanceSettings, "Minimum distance");
    auto* vLayout = new QVBoxLayout(minimumDistanceSettings);
    vLayout->setAlignment(Qt::AlignTop);
    minimumDistanceSettings->setLayout(vLayout);

    auto* ignoreBbox = new QCheckBox("Ignore bounding box");
    ignoreBbox->setChecked(true);
    m_forcer.m_ignoreBbox = ignoreBbox->isChecked();
    vLayout->addWidget(ignoreBbox);

    auto* minAngleLabel = new QLabel("Minimum angle between adjacent sites");
    vLayout->addWidget(minAngleLabel);
    m_minAngle = new DoubleSlider(Qt::Horizontal);
    m_minAngle->setMinimum(0);
    m_minAngle->setMaximum(std::numbers::pi);
    m_minAngle->setValue(1);
    m_forcer.m_minAngle = m_minAngle->value();
    vLayout->addWidget(m_minAngle);

    auto* minAdjDistLabel = new QLabel("Minimum adjacency distance");
    vLayout->addWidget(minAdjDistLabel);
    auto* minAdjDist = new DoubleSlider(Qt::Horizontal);
    vLayout->addWidget(minAdjDist);

    auto* minAdjDistSpinBox = new QDoubleSpinBox();
    minAdjDistSpinBox->setSuffix(" m");
    vLayout->addWidget(minAdjDistSpinBox);

    m_minAdjDist = std::make_unique<DoubleSliderSpinBox>(minAdjDist, minAdjDistSpinBox);
    m_minAdjDist->setMinimum(0);
    m_minAdjDist->setMaximum(getScale() * 30);
    m_minAdjDist->setValue(getScale() * 30);
    m_forcer.m_minAdjDist = m_minAdjDist->value() / getScale();

    auto* minDistLabel = new QLabel("Minimum distance between lines");
    vLayout->addWidget(minDistLabel);
    auto* minDist = new DoubleSlider(Qt::Horizontal);
    vLayout->addWidget(minDist);

    auto* minDistSpinBox = new QDoubleSpinBox();
    minDistSpinBox->setSuffix(" m");
    vLayout->addWidget(minDistSpinBox);

    m_minDist = std::make_unique<DoubleSliderSpinBox>(minDist, minDistSpinBox);
    m_minDist->setMinimum(0);
    m_minDist->setMaximum(getScale() * 30);
    m_minDist->setValue(getScale() * 6);
    m_forcer.m_requiredMinDist = minDist->value();

    auto* minComponentLengthLabel = new QLabel("Minimum component length");
    vLayout->addWidget(minComponentLengthLabel);
    auto* minComponentLength = new DoubleSlider(Qt::Horizontal);
    
    vLayout->addWidget(minComponentLength);

    auto* minCompLengthSpinBox = new QDoubleSpinBox();
    minCompLengthSpinBox->setSuffix(" m");
    vLayout->addWidget(minCompLengthSpinBox);

    m_minComponentLength = std::make_unique<DoubleSliderSpinBox>(minComponentLength, minCompLengthSpinBox);
    m_minComponentLength->setMinimum(0);
    m_minComponentLength->setMaximum(getScale() * 100);
    m_minComponentLength->setValue(0);
    m_forcer.m_requiredLength = m_minComponentLength->value();

    auto* mdInitializeButton = new QPushButton("Initialize / reset");
    vLayout->addWidget(mdInitializeButton);
    auto* mdStepButton = new QPushButton("Step");
    vLayout->addWidget(mdStepButton);
    auto* mdStep100Button = new QPushButton("Step (x100)");
    vLayout->addWidget(mdStep100Button);
    auto* mdReconstructButton = new QPushButton("Reconstruct");
    vLayout->addWidget(mdReconstructButton);
    auto* undoReconstructButton = new QPushButton("Undo reconstruct");
    vLayout->addWidget(undoReconstructButton);

    connect(ignoreBbox, &QCheckBox::stateChanged, [this, ignoreBbox] {
        m_forcer.m_ignoreBbox = ignoreBbox->isChecked();
        m_forcer.recomputeDelaunay();
        m_forcer.recomputeAuxiliary();
        repaintVoronoi();
        m_renderer->repaint();
    });

    connect(&*m_minAdjDist, &DoubleSliderSpinBox::valueChanged, [this](double v) {
        m_forcer.m_minAdjDist = v / getScale();
        m_forcer.recomputeDelaunay();
        m_forcer.recomputeAuxiliary();
        repaintVoronoi();
        m_renderer->repaint();
    });

    connect(&*m_minDist, &DoubleSliderSpinBox::valueChanged, [this](double v) {
        m_forcer.m_requiredMinDist = v / getScale();
        m_forcer.recomputeAuxiliary();
        repaintVoronoi();
        m_renderer->repaint();
    });

    connect(m_minAngle, &DoubleSlider::valueChanged, [this](double v) {
        m_forcer.m_minAngle = v;
        m_forcer.recomputeAuxiliary();
        repaintVoronoi();
        m_renderer->repaint();
    });

    connect(&*m_minComponentLength, &DoubleSliderSpinBox::valueChanged, [this](double v) {
        m_forcer.m_requiredLength = v / getScale();
        repaintVoronoi();
        m_renderer->repaint();
    });

    connect(mdInitializeButton, &QPushButton::clicked, [this]() {
        m_forcer.m_g = approximateBezierGraph(m_baseGraph, std::min(20.0, 1000000.0 / m_baseGraph.number_of_edges()));
        m_forcer.initialize();

        repaintVoronoi();

        m_renderer->repaint();
    });

    connect(mdStepButton, &QPushButton::clicked, [this]() {
        m_forcer.step();
        repaintVoronoi();
        m_renderer->repaint();
    });

    connect(mdStep100Button, &QPushButton::clicked, [this]() {
        QProgressDialog progress("Pushing isolines...", "Abort", 0, 100, this);
        progress.setWindowModality(Qt::WindowModal);
        progress.setMinimumDuration(1000);
        for (int i = 0; i < 100; ++i) {
            progress.setValue(i);
            if (progress.wasCanceled()) break;
            if (!m_forcer.step()) break;
            m_renderer->repaint();
        }
        repaintVoronoi();
        m_renderer->repaint();
    });

    connect(mdReconstructButton, &QPushButton::clicked, [this]() {
        m_beforeReconstruct = m_baseGraph;
        m_baseGraph = reconstructBezierGraph(m_forcer.m_g, m_minDist->value() / getScale() * m_minDist->value() / getScale() / 16);
        m_forcer.m_g.clear();
        m_forcer.m_withinDistEdgeComponents.clear();
        m_forcer.m_delaunay.clear();
        m_voronoiPainting = {};
        m_renderer->repaint();
    });

    connect(undoReconstructButton, &QPushButton::clicked, [this]() {
        m_baseGraph = m_beforeReconstruct;
        m_renderer->repaint();
    });
}

void BezierSimplificationDemo::addDrawingTab() {
    auto* drawSettings = new QLabel();
    m_tabs->addTab(drawSettings, "Drawing");
    auto* vLayout = new QVBoxLayout(drawSettings);
    vLayout->setAlignment(Qt::AlignTop);

    showEdgeDirection = new QCheckBox("Show edge direction");
    vLayout->addWidget(showEdgeDirection);
    showOldVertices = new QCheckBox("Show old vertices");
    vLayout->addWidget(showOldVertices);
    showNewVertices = new QCheckBox("Show new vertices");
    vLayout->addWidget(showNewVertices);
    showNewControlPoints = new QCheckBox("Show new control points");
    vLayout->addWidget(showNewControlPoints);
    showDebugInfo = new QCheckBox("Show debug info");
    vLayout->addWidget(showDebugInfo);

    connect(showEdgeDirection, &QCheckBox::stateChanged, [this]() {
        m_renderer->repaint();
    });
    connect(showOldVertices, &QCheckBox::stateChanged, [this]() {
        m_renderer->repaint();
    });
    connect(showNewVertices, &QCheckBox::stateChanged, [this]() {
        m_renderer->repaint();
    });
    connect(showNewControlPoints, &QCheckBox::stateChanged, [this]() {
        m_renderer->repaint();
    });
    connect(showDebugInfo, &QCheckBox::stateChanged, [this]() {
        m_renderer->repaint();
    });
}

double BezierSimplificationDemo::getScale() const {
    auto unitXT = Vector<Inexact>(1, 0).transform(m_transform.inverse());
    return unitXT.x();
}

void BezierSimplificationDemo::addPaintings() {
    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        int i = 0;
        for (const auto& refData : m_referenceData) {
            renderer.setStroke(m_colors.at(i), 2.0);
            if (auto qImageP = std::get_if<QImage>(&refData)) {
                auto& qImage = *qImageP;
                const auto formats = QImageReader::supportedImageFormats();
                if (!formats.contains("tiff")) {
                    std::cerr << "No support for tiff" << std::endl;
                }
                if (auto gw = dynamic_cast<GeometryWidget*>(&renderer)) {
                    gw->drawImage(m_referencePolygon, qImage);
                }
            } else if (auto geometrySetP = std::get_if<GeometrySet<Inexact>>(&refData)) {
                auto& geometrySet = *geometrySetP;
                for (const auto& geom : geometrySet.geometries) {
                    std::visit([&](const auto& g) { renderer.draw(g); }, geom);
                }
            }
            ++i;
        }
    }, "Reference data");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setStroke(Color(200, 200, 200), 2.0);
        for (auto eit = m_original.edges_begin(); eit != m_original.edges_end(); ++eit) {
            renderer.draw(eit->curve().transform(m_transform.inverse()));
        }
        if (showOldVertices->isChecked()) {
            for (auto vit = m_original.vertices_begin(); vit != m_original.vertices_end(); ++vit) {
                renderer.draw(vit->point().transform(m_transform.inverse()));
            }
        }
    }, "Original");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::fill);
        renderer.setFill(Color(230, 230, 230));
        renderer.draw(Circle<Inexact>(m_renderer->mousePosition(), m_minDist->value() * m_minDist->value()));
    }, "Min. dist. disk");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        auto inv = m_transform.inverse();
        renderer.setMode(GeometryRenderer::stroke);
        renderer.setStroke(Color(0, 0, 0), 2.0);
        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            if (showDebugInfo->isChecked() && eit->data().collapse.has_value()) {
                renderer.setStroke(Color(50, 200, 50), 2.0);
            } else {
                renderer.setStroke(Color(0, 0, 0), 2.0);
            }
            renderer.draw(eit->curve().transform(inv));
            if (showEdgeDirection->isChecked()) {
                renderer.setStroke(Color(0, 0, 0), 6.0);
                renderer.draw(eit->curve().split(0.25).first.transform(inv));
            }
        }
        // Control polylines
        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            if (showNewControlPoints->isChecked()) {
                renderer.setStroke(Color(0, 255, 0), 1.0);
                Polyline<Inexact> pl;
                for (int c = 0; c < 4; ++c) pl.push_back(eit->curve().control(c));
                renderer.draw(pl.transform(inv));
            }
        }
        // Control points
        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            if (showNewControlPoints->isChecked()) {
                renderer.setStroke(Color(0, 0, 255), 2.0);
                renderer.draw(eit->curve().sourceControl().transform(inv));
                renderer.draw(eit->curve().targetControl().transform(inv));
                renderer.setStroke(Color(255, 0, 255), 2.0);
                renderer.draw(eit->curve().source().transform(inv));
                renderer.draw(eit->curve().target().transform(inv));
            }
        }
        if (showNewVertices->isChecked()) {
            for (auto vit = m_baseGraph.vertices_begin(); vit != m_baseGraph.vertices_end(); ++vit) {
                renderer.draw(vit->point().transform(inv));
            }
        }
    }, "Simplification");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        m_voronoiPainting.paint(renderer);
    }, "Segment Voronoi diagram");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        auto inv = m_transform.inverse();
        renderer.setStroke(Color{0, 0, 0}, 2.0);
        for (auto eit = m_forcer.m_g.edges_begin(); eit != m_forcer.m_g.edges_end(); ++eit) {
            renderer.draw(eit->curve().transform(inv));
        }
        for (auto vit = m_forcer.m_g.vertices_begin(); vit != m_forcer.m_g.vertices_end(); ++vit) {
            renderer.draw(vit->point().transform(inv));
        }
    }, "Linearized");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        auto inv = m_transform.inverse();
        if (!m_debugEdge.has_value()) return;
        renderer.setStroke(Color(255, 0, 0), 2.0);
        auto eh = *m_debugEdge;
        renderer.draw(eh->curve());
        const auto& d = eh->data();
        renderer.setStroke(Color(50, 200, 50), 2.0);
        if (d.collapse.has_value()) {
            const auto& col = *d.collapse;
            renderer.draw(col.before.transform(inv));
            renderer.draw(col.after.transform(inv));
            auto bb = col.before.transform(inv).bbox();
            auto bbC = Point<Inexact>((bb.xmin() + bb.xmax()) / 2, (bb.ymin() + bb.ymax())  / 2);
            renderer.drawText(bbC, std::to_string(col.cost));
        }
        if (d.hist == nullptr) return;
        auto op = std::dynamic_pointer_cast<curved_simplification::detail::CollapseEdgeOperation<BaseGraph>>(d.hist);
        renderer.setStroke(Color(200, 50, 200), 2.0);
        renderer.draw(op->m_c0.transform(inv));
        renderer.draw(op->m_c1.transform(inv));
        renderer.draw(op->m_c2.transform(inv));
    }, "Debug edge");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        auto inv = m_transform.inverse();
        if (!m_editControlPoints->isChecked()) return;

        // Control polylines
        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            renderer.setStroke(Color(0, 255, 0), 1.0);
            Polyline<Inexact> pl;
            for (int c = 0; c < 4; ++c) pl.push_back(eit->curve().control(c));
            renderer.draw(pl.transform(inv));
        }

        for (auto editable : m_editables) {
            if (auto vhP = std::get_if<Vertex_handle>(&editable->type)) {
                renderer.setStroke(Color(255, 0, 255), 2.0);
                renderer.draw(editable->point.transform(inv));
            } else if (auto eiP = std::get_if<std::pair<Edge_handle, bool>>(&editable->type)) {
                renderer.setStroke(Color(0, 0, 255), 2.0);
                auto [eh, b] = *eiP;
                if (b) {
                    renderer.draw(editable->point.transform(inv));
                } else {
                    renderer.draw(editable->point.transform(inv));
                }
            }
        }
    }, "Control points");
}

BezierSimplificationDemo::BezierSimplificationDemo() : m_graph(m_baseGraph), m_collapse(m_graph, Traits()), m_forcer(m_approxGraph, 1.0) {
	setWindowTitle("CartoCrow: Bézier simplification");
	m_renderer = new GeometryWidget();
	m_renderer->setDrawAxes(false);
    m_renderer->setMaxZoom(1000000000);
    m_renderer->setMinZoom(0.00000001);
	setCentralWidget(m_renderer);

	auto* dockWidget = new QDockWidget();
	addDockWidget(Qt::RightDockWidgetArea, dockWidget);
    m_tabs = new QTabWidget();
	dockWidget->setWidget(m_tabs);

    addIOTab();
    addSimplificationTab();
    addMinimumDistanceTab();
    addDrawingTab();

	loadInput("splines.ipe");

    connect(m_renderer, &GeometryWidget::clicked, [this](Point<Inexact> pt) {
        double minDist2Edge = std::numeric_limits<double>::infinity();
        std::optional<Graph::Edge_handle> closest;
        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            double minDist2 = CGAL::squared_distance(pt.transform(m_transform), eit->curve().nearest(pt).point);
            if (minDist2 < minDist2Edge) {
                minDist2Edge = minDist2;
                closest = eit;
            }
        }
        if (m_debugEdge == closest) {
            m_debugEdge = std::nullopt;
            repaint();
            return;
        }
        m_debugEdge = closest;
        Traits debugTraits(true);
        debugTraits.determineCollapse(*m_debugEdge);
        repaint();
    });

    addPaintings();
}

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
	BezierSimplificationDemo demo;
	demo.show();
	app.exec();
}
