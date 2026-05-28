#include "frontend.h"

#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QPushButton>
#include <QCheckBox>
#include <QListWidget>
#include <QSpinBox>
#include <QImageReader>
#include <QLabel>
#include <QScrollArea>
#include <QFileDialog>
#include <QProgressDialog>
#include <QShortcut>

#include <cartocrow/widgets/double_slider.h>

#include "library/read_graph_gdal.h"
#include "library/read_ipe_bezier_spline.h"
#include "library/vertex_snap.h"
#include "library/topojson_export.h"

#include <cartocrow/renderer/voronoi_drawer.h>
#include <cartocrow/renderer/painting_renderer.h>
#include <cartocrow/core/transform_helpers.h>

#define DEBUG 0

void saveGraphIntoTopoSet(const Graph& graph, BezierTopoSet& topoSet) {
    std::unordered_set<const Graph::Edge*> visited;

    for (auto eit : graph.edges()) {
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
        Graph::Edge_const_handle end;
        if (current->source()->degree() == 2) {
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

        CubicBezierSpline arc;
        arc.appendCurve(start->curve());

        if (start != end) {
            for (auto arcEit = start->next(); true; arcEit = arcEit->next()) {
                arc.appendCurve(arcEit->curve());
                if (arcEit == end) break;
            }
        }

        topoSet.arcs[arcIndex] = arc;
    }
}

bool liesOnBbox(const Point<Inexact>& p, const Box& bbox) {
    return abs(p.x() - bbox.xmin()) < M_EPSILON ||
        abs(p.x() - bbox.xmax()) < M_EPSILON ||
        abs(p.y() - bbox.ymin()) < M_EPSILON ||
        abs(p.y() - bbox.ymax()) < M_EPSILON;
}

void insertSplinesIntoGraph(const std::vector<CubicBezierSpline>& splines, Graph& g) {
    std::unordered_map<Point<Inexact>, Graph::Vertex_handle> pToV;

    auto getVertex = [&pToV, &g](const Point<Inexact>& p) {
        if (pToV.contains(p)) {
            return pToV.at(p);
        }
        else {
            auto newV = g.add_vertex(p);
            pToV[p] = newV;
            return newV;
        }
    };

    for (const auto& spline : splines) {
        if (spline.empty()) continue;
        auto lastV = getVertex(spline.curves()[0].source());
        for (const auto& curve : spline.curves()) {
            auto targetV = getVertex(curve.target());
            g.add_edge(lastV, targetV, curve);
            lastV = targetV;
        }
    }
}

void insertTopoSetIntoGraph(const BezierTopoSet& ts, Graph& g, bool ignoreBbox) {
    std::unordered_map<Point<Inexact>, Graph::Vertex_handle> pToV;

    auto getVertex = [&pToV, &g](const Point<Inexact>& p) {
        if (pToV.contains(p)) {
            return pToV.at(p);
        }
        else {
            auto newV = g.add_vertex(p);
            pToV[p] = newV;
            return newV;
        }
    };

    auto bboxT = ts.bbox();

    for (int i = 0; i < ts.arcs.size(); ++i) {
        const auto& arc = ts.arcs[i];
        for (auto eit = arc.curves_begin(); eit != arc.curves_end(); ++eit) {
            auto sourceV = getVertex(eit->source());
            auto targetV = getVertex(eit->target());
            auto curve = *eit;

            // There should be no duplicates but there are somehow a couple, so we do this for now...
            bool oppositeExists = false;
            bool alreadyExists = false;
            for (auto ieit = targetV->incident_edges_begin(); ieit != targetV->incident_edges_end(); ++ieit) {
                if ((*ieit)->target() == sourceV && (*ieit)->curve() == curve.reversed()) {
                    oppositeExists = true;
                    std::cout << "Found an opposite edge!" << std::endl;
                    break;
                }
                if ((*ieit)->source() == sourceV && (*ieit)->curve() == curve) {
                    alreadyExists = true;
                    std::cout << "Found a duplicate edge!" << std::endl;
                    break;
                }
            }

            if (!oppositeExists && !alreadyExists) {
                auto eh = g.add_edge(sourceV, targetV, curve);
                eh->data().index = i;
                bool edgeOnBbox = ignoreBbox && liesOnBbox(sourceV->point(), bboxT) && liesOnBbox(targetV->point(), bboxT);
                eh->data().collapse_allowed = !edgeOnBbox;
            }
        }
    }
}

void BezierSimplificationDemo::loadInput(const std::filesystem::path& path) {
    baseModified(false);
    double prevScale = getScale();

    m_baseGraph.clear();
    m_debugEdge = std::nullopt;
    m_forcer.m_withinDistEdgeComponents.clear();
    m_forcer.m_delaunay.clear();
    m_approxGraph.clear();
    m_voronoiPainting = PaintingRenderer{};

    if (path.extension() == ".ipe") {
        auto splines = ipeSplinesToIsolines(path);
        insertSplinesIntoGraph(splines, m_baseGraph);
    }
    else if (path.extension() == ".json") {
        auto [topoSet, spatialRef] = parseBezierTopoJson(path);
        m_spatialRef = spatialRef;

        m_transform = fitInto(topoSet.bbox(), Rectangle<Inexact>(0, 0, 1000, 1000));
        m_graphPainting->m_drawSettings.m_trans = m_transform.inverse();
        m_editGraphPainting->m_drawSettings.m_trans = m_transform.inverse();

        m_toposet = topoSet.transform(m_transform);
        insertTopoSetIntoGraph(m_toposet, m_baseGraph, m_ignoreBbox->isChecked());
    }
    else {
        auto [regionSet, spatialRef] = readRegionSetUsingGDAL(path);
        m_spatialRef = spatialRef;

        m_transform = fitInto(regionSet.bbox(), Rectangle<Inexact>(0, 0, 1000, 1000));
        m_graphPainting->m_drawSettings.m_trans = m_transform.inverse();
        m_editGraphPainting->m_drawSettings.m_trans = m_transform.inverse();

        regionSet = regionSet.transform(m_transform);

        snapVertices(regionSet);
        m_toposet = bezierTopoSetFromStraightTopoSet(regionSetToTopoSet(regionSet));
        insertTopoSetIntoGraph(m_toposet, m_baseGraph, m_ignoreBbox->isChecked());
    }

    m_baseGraph.orient();
    m_original = m_baseGraph;
    m_graph.reset();

    auto bbox = m_baseGraph.bbox();
    m_referencePolygon = bbox;

    m_renderer->fitInView(transform(bbox, m_transform.inverse()));

    m_complexitySliders->setMaximum(m_graph.number_of_edges());
    m_complexitySliders->setMinimum(m_graph.number_of_edges());
    m_complexitySliders->setValue(m_graph.number_of_edges());

    m_desiredComplexity->setMinimum(1);
    m_desiredComplexity->setMaximum(m_graph.number_of_edges());
    m_minDist->setMaximum(getScale() * 30);
    m_minAdjDist->setMaximum(getScale() * 30);
    m_minComponentLength->setMaximum(getScale() * 100);
    m_minCrumbArea->setMaximum(getScale() * 100);
    m_minDist->setValue(getScale() / prevScale * m_minDist->value());
    m_minAdjDist->setValue(getScale() / prevScale * m_minAdjDist->value());
    m_minComponentLength->setValue(getScale() / prevScale * m_minComponentLength->value());
    m_minCrumbArea->setValue(getScale() / prevScale * m_minCrumbArea->value());
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

// Rewrite the below functions so that they may accept a graph and store the movement in an operation batch.
void moveVertex(Vertex_handle vh, Point<Inexact> p, Graph* g = nullptr) {
    if (g == nullptr) {
        vh->point() = p;
    } else {
        g->moveVertex(vh, p);
    }
    for (auto ieit = vh->incident_edges_begin(); ieit != vh->incident_edges_end(); ++ieit) {
        auto eh = *ieit;
        if (eh->source() == vh) {
            auto& c = eh->curve();
            auto diff = c.sourceControl() - c.source();
            auto newCurve = CubicBezierCurve(p, p + diff, c.targetControl(), c.target());
            if (g == nullptr) {
                c = newCurve;
            }
            else {
                g->replaceCurve(eh, newCurve);
            }
        }
        else {
            assert(eh->target() == vh);
            auto& c = eh->curve();
            auto diff = c.targetControl() - c.target();
            auto newCurve = CubicBezierCurve(c.source(), c.sourceControl(), p + diff, p);
            if (g == nullptr) {
                c = newCurve;
            }
            else {
                g->replaceCurve(eh, newCurve);
            }
        }
    }
}

void moveInnerControlPoint(InnerControlPoint cp, Point<Inexact> p, bool alignTangents, Graph* g = nullptr) {
    auto [eh, b] = cp;
    auto& c = eh->curve();
    auto newCurve = CubicBezierCurve(c.source(), !b ? p : c.sourceControl(), b ? p : c.targetControl(), c.target());
    if (g == nullptr) {
        c = newCurve;
    }
    else {
        g->replaceCurve(eh, newCurve);
    }

    if (alignTangents) {
        if (b && eh->target()->degree() == 2) {
            auto next = eh->next();
            auto& nc = next->curve();
            auto diff = nc.sourceControl() - nc.source();
            auto dist = sqrt(diff.squared_length());
            auto vec = c.target() - c.targetControl();
            vec /= sqrt(vec.squared_length());
            auto newSourceControl = nc.source() + vec * dist;
            auto newCurve = CubicBezierCurve(nc.source(), newSourceControl, nc.targetControl(), nc.target());
            if (g == nullptr) {
                nc = newCurve;
            }
            else {
                g->replaceCurve(next, newCurve);
            }
        }
        else if (!b && eh->source()->degree() == 2) {
            auto prev = eh->prev();
            auto& pc = prev->curve();
            auto diff = pc.target() - pc.targetControl();
            auto dist = sqrt(diff.squared_length());
            auto vec = c.sourceControl() - c.source();
            vec /= sqrt(vec.squared_length());
            auto newTargetControl = pc.target() - vec * dist;
            auto newCurve = CubicBezierCurve(pc.source(), pc.sourceControl(), newTargetControl, pc.target());
            if (g == nullptr) {
                pc = newCurve;
            }
            else {
                g->replaceCurve(prev, newCurve);
            }
        }
    }
}

void moveControlPoint(ControlPoint cp, Point<Inexact> p, bool alignTangents, Graph* g = nullptr) {
    if (auto vhP = std::get_if<Vertex_handle>(&cp.type)) {
        moveVertex(*vhP, p, g);
    }
    else if (auto eiP = std::get_if<InnerControlPoint>(&cp.type)) {
        moveInnerControlPoint(*eiP, p, alignTangents, g);
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
    m_complexitySliders->setMinimum(std::min((int)m_graph.number_of_edges(), m_complexitySliders->minimum()));
    m_complexitySliders->setValue(m_graph.number_of_edges());
}

void BezierSimplificationDemo::baseModified(bool modified) {
    m_renderer->setVisibility(m_graphPainting, !modified);
    m_renderer->setVisibility(m_editGraphPainting, modified);
    m_tabs->setTabEnabled(1, !modified);
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

    auto* referenceDataList = new QListWidget();
    vLayout->addWidget(referenceDataList);
    referenceDataList->setMaximumHeight(200);

    auto* clearReferenceDataButton = new QPushButton("Clear reference data");
    vLayout->addWidget(clearReferenceDataButton);

    auto* checkIntersectionsButton = new QPushButton("Check intersections");
    vLayout->addWidget(checkIntersectionsButton);

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

    auto shortcut = new QShortcut(QKeySequence(tr("e", "Edit control points")), this);
    connect(shortcut, &QShortcut::activated, [this]() {
        m_editControlPoints->toggle();
    });

    auto* resetEditsButton = new QPushButton("Edit: reset edits");
    vLayout->addWidget(resetEditsButton);

    auto undoShortcut = new QShortcut(QKeySequence(tr("Ctrl+z", "Undo")), this);
    connect(undoShortcut, &QShortcut::activated, [this]() {
        m_editGraph.backInTime();
        updateEditables();
        m_renderer->repaint();
    });
    auto redoShortcut = new QShortcut(QKeySequence(tr("Ctrl+Shift+z", "Redo")), this);
    connect(redoShortcut, &QShortcut::activated, [this]() {
        m_editGraph.forwardInTime();
        updateEditables();
        m_renderer->repaint();
    });

    connect(loadFileButton, &QPushButton::clicked, [this, loadFileButton]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getOpenFileName(this, tr("Select input file (.ipe or .shp)"), startDir).toStdString();
        if (filePath == "") return;
        loadInput(filePath);

        std::stringstream ss;
        ss << "Load file (" << filePath.filename().string() << ")";
        loadFileButton->setText(QString::fromStdString(ss.str()));
    });

    connect(loadReferenceDataButton, &QPushButton::clicked, [this, referenceDataList]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getOpenFileName(this, tr("Select reference file (.tif or .shp)"), startDir).toStdString();
        if (filePath == "") return;

        auto* listItem = new QListWidgetItem(QString::fromStdString(filePath.filename().string()), referenceDataList);
        listItem->setCheckState(Qt::Checked);

        if (filePath.extension() == ".tif") {
            m_referenceData.emplace_back(QImage(filePath.string().c_str()), filePath, listItem);
        } else {
            auto geometrySet = readGeometrySetUsingGDAL(filePath);
            m_referenceData.emplace_back(geometrySet, filePath, listItem);
            PaintingRenderer& pr = m_referenceData.back().painting;
            GeometryRenderer& r = pr;
            for (const auto& geom : geometrySet.geometries) {
                std::visit([&](const auto& g) { r.draw(g); }, geom);
            }
        }
        
        listItem->setCheckState(Qt::Checked);
    });

    connect(clearReferenceDataButton, &QPushButton::clicked, [this, referenceDataList]() {
        m_referenceData.clear();
        referenceDataList->clear();
    });

    connect(loadReferencePolygonButton, &QPushButton::clicked, [this]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getOpenFileName(this, tr("Select reference polygon for image (.shp)"), startDir).toStdString();
        if (filePath == "") return;
        auto regionSet = readRegionSetUsingGDAL(filePath);
        std::vector<PolygonWithHoles<Inexact>> pgns;
        m_referencePolygon = regionSet.first.regions[0].geometry.polygons_with_holes[0].bbox();
    });

    connect(checkIntersectionsButton, &QPushButton::clicked, [this]() {
        checkIntersections();
        m_renderer->repaint();
    });

    connect(exportButton, &QPushButton::clicked, [this, stackPolygons]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getSaveFileName(this, tr("Select folder to create a folder with .shp and auxiliary files."), startDir).toStdString();
        if (filePath.empty()) return;

        std::ofstream out(filePath / "beziers.txt");

        auto& theBaseGraph = m_editBaseGraph.number_of_edges() > 0 ? m_editBaseGraph : m_baseGraph;

        saveGraphIntoTopoSet(theBaseGraph, m_toposet);
        std::string jsonFileName = filePath.stem().string() + ".json";
        exportTopoSetToJson(filePath / jsonFileName, m_toposet.transform(m_transform.inverse()), m_spatialRef);
        exportTopoSetUsingGDAL(filePath, approximate(m_toposet), m_transform.inverse(), m_spatialRef, stackPolygons->isChecked());
    });

    connect(m_editControlPoints, &QCheckBox::stateChanged, [this]() {
        if (!m_editControlPoints->isChecked()) {
            m_renderer->repaint();
            return;
        }
        if (m_editBaseGraph.number_of_vertices() == 0) {
            m_editBaseGraph = m_baseGraph;
        }
        updateEditables();
        m_renderer->repaint();
    });

    connect(resetEditsButton, &QPushButton::clicked, [this]() {
        resetEdits();
        m_renderer->repaint();
    });

    connect(m_renderer, &GeometryWidget::dragStarted, [this](const Point<Inexact>& p) {
        if (!m_editControlPoints->isChecked()) return;

        m_dragging = nullptr;

        std::optional<std::shared_ptr<ControlPoint>> closest;
        double minDist = std::numeric_limits<double>::infinity();

        for (const auto& editable : m_editables) {
            auto diff = m_renderer->convertPoint(editable->point().transform(m_transform.inverse())) - m_renderer->convertPoint(p);
            auto d2 = diff.x() * diff.x() + diff.y() * diff.y();
            if (d2 < 400 && d2 < minDist) {
                minDist = d2;
                closest = editable;
            }
        }

        if (closest.has_value()) {
            m_dragging = std::make_shared<Dragging>(**closest, (*closest)->point());
        }
        m_renderer->repaint();
    });

    connect(m_renderer, &GeometryWidget::dragMoved, [this, editAlignTangents](const Point<Inexact>& px) {
        if (!m_editControlPoints->isChecked()) return;

        auto p = px.transform(m_transform);
        if (m_dragging != nullptr) {
            baseModified(true);
            moveControlPoint(m_dragging->controlPoint, p, editAlignTangents->isChecked());
            m_renderer->repaint();
        }
    });

    connect(m_renderer, &GeometryWidget::dragEnded, [this, editAlignTangents](const Point<Inexact>& p) {
        if (!m_editControlPoints->isChecked()) return;
        if (m_dragging == nullptr) return;

        auto finalPosition = m_dragging->controlPoint.point();
        // Move control point back
        moveControlPoint(m_dragging->controlPoint, m_dragging->startPoint, editAlignTangents->isChecked());
        // Move it in one go
        m_editGraph.startBatch();
        moveControlPoint(m_dragging->controlPoint, finalPosition, editAlignTangents->isChecked(), &m_editGraph);
        m_editGraph.endBatch();
        m_dragging = nullptr;
        m_renderer->repaint();
    });

    connect(referenceDataList, &QListWidget::itemChanged, [this]() {
        m_renderer->repaint();
    });
}

void BezierSimplificationDemo::resetEdits() {
    baseModified(false);
    m_editControlPoints->setChecked(false);

    if (m_editBaseGraph.number_of_vertices() > 0) {
        m_editBaseGraph = {};
    }
}

void BezierSimplificationDemo::addSimplificationTab() {
    auto* simplificationSettings = new QWidget();
    m_tabs->addTab(simplificationSettings, "Simplification");
    auto* vLayout = new QVBoxLayout(simplificationSettings);
    vLayout->setAlignment(Qt::AlignTop);

    m_ignoreBbox = new QCheckBox("Ignore bounding box");
    m_ignoreBbox->setChecked(true);
    vLayout->addWidget(m_ignoreBbox);

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
    m_desiredComplexity = new QSpinBox();
    vLayout->addWidget(m_desiredComplexity);

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

    m_complexitySliders = std::make_unique<ComplexitySliders>(vLayout);
    m_complexitySliders->setMaximum(m_graph.number_of_edges());
    m_complexitySliders->setMinimum(m_graph.number_of_edges());
    m_complexitySliders->setValue(m_graph.number_of_edges());

    connect(&*m_complexitySliders, &ComplexitySliders::valueChanged, [this](double v) {
        resetEdits();
        m_graph.recallComplexity(v);
        m_renderer->repaint();
//        m_backup.clear();
    });

    m_desiredComplexity->setMinimum(1);
    m_desiredComplexity->setMaximum(m_graph.number_of_edges());
    m_desiredComplexity->setValue(1);

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
        m_complexitySliders->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(undo10Button, &QPushButton::clicked, [this]() {
        resetEdits();
        for (int i = 0; i < 10; ++i)
            m_graph.backInTime();
        m_complexitySliders->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(undo100Button, &QPushButton::clicked, [this]() {
        resetEdits();
        for (int i = 0; i < 100; ++i)
            m_graph.backInTime();
        m_complexitySliders->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(redoButton, &QPushButton::clicked, [this]() {
        resetEdits();
        m_graph.forwardInTime();
        m_complexitySliders->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(redo10Button, &QPushButton::clicked, [this]() {
        resetEdits();
        for (int i = 0; i < 10; ++i)
            m_graph.forwardInTime();
        m_complexitySliders->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(redo100Button, &QPushButton::clicked, [this]() {
        resetEdits();
        for (int i = 0; i < 100; ++i)
            m_graph.forwardInTime();
        m_complexitySliders->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(runToComplexityButton, &QPushButton::clicked, [this]() {
        auto startComplexity = m_graph.number_of_edges();
        auto targetComplexity = m_desiredComplexity->value();
        QProgressDialog progress("Simplifying...", "Abort", 0, startComplexity - targetComplexity, this);
        progress.setWindowModality(Qt::WindowModal);
        progress.setMinimumDuration(1000);
        m_collapse.runToComplexity(m_desiredComplexity->value(), [&progress, startComplexity, this](int i) {
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
    m_minAngle->setMaximum(std::numbers::pi / 2);
    m_minAngle->setValue(std::numbers::pi / 4);
    m_forcer.m_minAngle = m_minAngle->value();
    vLayout->addWidget(m_minAngle);

    auto* minAdjDistLabel = new QLabel("Minimum adjacency distance");
    vLayout->addWidget(minAdjDistLabel);

    auto* minAdjDistSpinBox = new QDoubleSpinBox();
    minAdjDistSpinBox->setSuffix(" m");
    vLayout->addWidget(minAdjDistSpinBox);

    auto* minAdjDist = new DoubleSlider(Qt::Horizontal);
    vLayout->addWidget(minAdjDist);

    m_minAdjDist = std::make_unique<DoubleSliderSpinBox>(minAdjDist, minAdjDistSpinBox);
    m_minAdjDist->setMinimum(0);
    m_minAdjDist->setMaximum(getScale() * 30);
    m_minAdjDist->setValue(getScale() * 30);
    m_forcer.m_minAdjDist = m_minAdjDist->value() / getScale();

    auto* minDistLabel = new QLabel("Minimum distance between lines");
    vLayout->addWidget(minDistLabel);

    auto* minDistSpinBox = new QDoubleSpinBox();
    minDistSpinBox->setSuffix(" m");
    vLayout->addWidget(minDistSpinBox);

    auto* minDist = new DoubleSlider(Qt::Horizontal);
    vLayout->addWidget(minDist);

    m_minDist = std::make_unique<DoubleSliderSpinBox>(minDist, minDistSpinBox);
    m_minDist->setMinimum(0);
    m_minDist->setMaximum(getScale() * 30);
    m_minDist->setValue(getScale() * 6);
    m_forcer.m_requiredMinDist = minDist->value();

    auto* minComponentLengthLabel = new QLabel("Minimum component length");
    vLayout->addWidget(minComponentLengthLabel);
    auto* minComponentLength = new DoubleSlider(Qt::Horizontal);

    auto* minCompLengthSpinBox = new QDoubleSpinBox();
    minCompLengthSpinBox->setSuffix(" m");
    vLayout->addWidget(minCompLengthSpinBox);

    vLayout->addWidget(minComponentLength);

    m_minComponentLength = std::make_unique<DoubleSliderSpinBox>(minComponentLength, minCompLengthSpinBox);
    m_minComponentLength->setMinimum(0);
    m_minComponentLength->setMaximum(getScale() * 100);
    m_minComponentLength->setValue(0);
    m_forcer.m_requiredLength = m_minComponentLength->value();

    auto* minCrumbAreaLabel = new QLabel("Minimum crumb area");
    vLayout->addWidget(minCrumbAreaLabel);
    auto* minCrumbArea = new DoubleSlider(Qt::Horizontal);
    auto* minCrumbAreaSpinBox = new QDoubleSpinBox();
    minCrumbAreaSpinBox->setSuffix(" km2");
    vLayout->addWidget(minCrumbAreaSpinBox);
    vLayout->addWidget(minCrumbArea);

    m_minCrumbArea = std::make_unique<DoubleSliderSpinBox>(minCrumbArea, minCrumbAreaSpinBox);
    m_minCrumbArea->setMinimum(0);
    m_minCrumbArea->setMaximum(getScale() * 100);
    m_minCrumbArea->setValue(0);

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

    connect(&*m_minCrumbArea, &DoubleSliderSpinBox::valueChanged, [this](double v) {
        m_forcer.m_minLoopArea = v / getScale() * 1E6 / getScale();
        m_forcer.initialize();
        repaintVoronoi();
        m_renderer->repaint();
    });

    connect(mdInitializeButton, &QPushButton::clicked, [this]() {
        auto& theBaseGraph = m_editBaseGraph.number_of_edges() > 0 ? m_editBaseGraph : m_baseGraph;
        m_forcer.m_g = approximateBezierGraph(theBaseGraph, std::min(20.0, 1000000.0 / theBaseGraph.number_of_edges()));
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
        auto& theBaseGraph = m_editBaseGraph.number_of_edges() > 0 ? m_editBaseGraph : m_baseGraph;
        m_beforeReconstruct = theBaseGraph;
        theBaseGraph = reconstructBezierGraph(m_forcer.m_g, m_minDist->value() / getScale() * m_minDist->value() / getScale() / 16);
        m_forcer.m_g.clear();
        m_forcer.m_withinDistEdgeComponents.clear();
        m_forcer.m_delaunay.clear();
        m_voronoiPainting = {};
        m_renderer->repaint();
    });

    connect(undoReconstructButton, &QPushButton::clicked, [this]() {
        auto& theBaseGraph = m_editBaseGraph.number_of_edges() > 0 ? m_editBaseGraph : m_baseGraph;
        theBaseGraph = m_beforeReconstruct;
        m_renderer->repaint();
    });
}

void BezierSimplificationDemo::addDrawingTab() {
    auto* drawSettings = new QLabel();
    m_tabs->addTab(drawSettings, "Drawing");
    auto* vLayout = new QVBoxLayout(drawSettings);
    vLayout->setAlignment(Qt::AlignTop);

    m_editHighlight = new QCheckBox("Edit highglighting");
    vLayout->addWidget(m_editHighlight);
    auto showEdgeDirection = new QCheckBox("Show edge direction");
    vLayout->addWidget(showEdgeDirection);
    m_showOldVertices = new QCheckBox("Show old vertices");
    vLayout->addWidget(m_showOldVertices);
    auto showNewVertices = new QCheckBox("Show new vertices");
    vLayout->addWidget(showNewVertices);
    auto showNewControlPoints = new QCheckBox("Show new control points");
    vLayout->addWidget(showNewControlPoints);
    auto showArcIndex = new QCheckBox("Show arc index");
    vLayout->addWidget(showArcIndex);
    m_showDebugInfo = new QCheckBox("Show debug info");
    vLayout->addWidget(m_showDebugInfo);
    auto* doTransform = new QCheckBox("Transform to 1000 x 1000");
    vLayout->addWidget(doTransform);

    m_showCurvature = new QCheckBox("Show curvature");
    vLayout->addWidget(m_showCurvature);
    m_curvatureScale = new DoubleSlider();
    m_curvatureScale->setOrientation(Qt::Horizontal);
    m_curvatureScale->setMinimum(0);
    m_curvatureScale->setValue(100);
    m_curvatureScale->setMaximum(1000);
    vLayout->addWidget(m_curvatureScale);

    connect(showEdgeDirection, &QCheckBox::stateChanged, [this](bool v) {
        m_graphPainting->m_drawSettings.m_showEdgeDirection = v;
        m_editGraphPainting->m_drawSettings.m_showEdgeDirection = v;
        m_renderer->repaint();
    });
    connect(m_showOldVertices, &QCheckBox::stateChanged, [this]() {
        m_renderer->repaint();
    });
    connect(showNewVertices, &QCheckBox::stateChanged, [this](bool v) {
        m_graphPainting->m_drawSettings.m_showVertices = v;
        m_editGraphPainting->m_drawSettings.m_showVertices = v;
        m_renderer->repaint();
    });
    connect(showNewControlPoints, &QCheckBox::stateChanged, [this](bool v) {
        m_graphPainting->m_drawSettings.m_showControlPoints = v;
        m_editGraphPainting->m_drawSettings.m_showControlPoints = v;
        m_renderer->repaint();
    });
    connect(showArcIndex, &QCheckBox::stateChanged, [this](bool v) {
        m_graphPainting->m_drawSettings.m_showArcIndex = v;
        m_editGraphPainting->m_drawSettings.m_showArcIndex = v;
        m_renderer->repaint();
    });
    connect(m_showDebugInfo, &QCheckBox::stateChanged, [this]() {
        m_renderer->repaint();
    });
    connect(doTransform, &QCheckBox::stateChanged, [this, doTransform]() {
        auto prevScale = getScale();
        auto minDistOld = m_minDist->value();
        auto minAdjDistOld = m_minAdjDist->value();
        auto minComponentLengthOld = m_minComponentLength->value();
        auto minCrumbAreaOld = m_minCrumbArea->value();

        if (doTransform->isChecked()) {
            m_backupTransform = m_transform;
            m_transform = CGAL::IDENTITY;
        } else {
            m_transform = m_backupTransform;
        }
        m_graphPainting->m_drawSettings.m_trans = m_transform.inverse();
        m_editGraphPainting->m_drawSettings.m_trans = m_transform.inverse();

        m_minDist->setMaximum(getScale() * 30);
        m_minAdjDist->setMaximum(getScale() * 30);
        m_minComponentLength->setMaximum(getScale() * 100);
        m_minCrumbArea->setMaximum(getScale() * 100);
        m_minDist->setValue(getScale() / prevScale * minDistOld);
        m_minAdjDist->setValue(getScale() / prevScale * minAdjDistOld);
        m_minComponentLength->setValue(getScale() / prevScale * minComponentLengthOld);
        m_minCrumbArea->setValue(getScale() * prevScale * minCrumbAreaOld);
        repaintVoronoi();
        Rectangle<Inexact> rect = m_baseGraph.bbox();
        auto rectT = rect.transform(m_transform.inverse());
        Box boxT(rectT.xmin(), rectT.ymin(), rectT.xmax(), rectT.ymax());
        m_renderer->fitInView(boxT);
        m_renderer->repaint();
    });
    connect(m_showCurvature, &QCheckBox::stateChanged, [this]() {
        m_renderer->repaint();
    });
    connect(m_curvatureScale, &DoubleSlider::valueChanged, [this]() {
        m_renderer->repaint();
    });
}

double BezierSimplificationDemo::getScale() const {
    auto unitXT = Vector<Inexact>(1, 0).transform(m_transform.inverse());
    return unitXT.x();
}

double BezierSimplificationDemo::screenDistanceToMouse2(const Point<Inexact>& p) const {
    auto mp = m_renderer->mousePosition();
    auto diff = m_renderer->convertPoint(p) - m_renderer->convertPoint(mp);
    return diff.x() * diff.x() + diff.y() * diff.y();
}

void BezierSimplificationDemo::addPaintings() {
    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        int i = 0;
        for (const auto& refData : m_referenceData) {
            if (refData.listItem->checkState() == Qt::Unchecked) {
                ++i;
                continue;
            }

            renderer.setStroke(m_colors.at(i), 2.0);

            if (auto qImageP = std::get_if<QImage>(&refData.data)) {
                const auto formats = QImageReader::supportedImageFormats();
                if (!formats.contains("tiff")) {
                    std::cerr << "No support for tiff" << std::endl;
                }
                if (auto gw = dynamic_cast<GeometryWidget*>(&renderer)) {
                    gw->drawImage(m_referencePolygon, *qImageP);
                }
            }
            else {
                refData.painting.paint(renderer);
            }
            ++i;
        }
    }, "Reference data");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setStroke(Color(200, 200, 200), 2.0);
        for (auto eit : m_original.edges()) {
            renderer.draw(eit->curve().transform(m_transform.inverse()));
        }
        if (m_showOldVertices->isChecked()) {
            for (auto vit : m_original.vertices()) {
                renderer.draw(vit->point().transform(m_transform.inverse()));
            }
        }
    }, "Original");

    auto minDistPainting = m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::fill);
        renderer.setFill(Color(230, 230, 230));
        renderer.draw(Circle<Inexact>(m_renderer->mousePosition(), m_minDist->value() * m_minDist->value()));
    }, "Min. dist. disk");
    m_renderer->setVisibility(minDistPainting, false);

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        if (!m_showCurvature->isChecked()) return;
        renderer.setMode(GeometryRenderer::stroke);
        renderer.setStroke(Color(200, 0, 200), 2.0);
        renderer.setStrokeOpacity(10);
        for (auto eit : m_baseGraph.edges()) {
            int tSteps = 1000;
            auto curve = eit->curve().transform(m_transform.inverse());
            for (int tStep = 0; tStep <= tSteps; ++tStep) {
                double t = static_cast<double>(tStep) / tSteps;
                auto c = curve.curvature(t);
                auto n = curve.normal(t);
                n /= sqrt(n.squared_length());
                auto p = curve.position(t);
                renderer.draw(Segment<Inexact>(p, p + m_curvatureScale->value() * c * n));
            }
        }
    }, "Curvature");

    m_graphPainting = std::make_shared<GraphPainting>(m_baseGraph);
    m_renderer->addPainting(m_graphPainting, "Simplification");

    m_editGraphPainting = std::make_shared<GraphPainting>(m_editBaseGraph);
    m_renderer->addPainting(m_editGraphPainting, "Edited simplification");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        m_voronoiPainting.paint(renderer);
    }, "Segment Voronoi diagram");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        auto inv = m_transform.inverse();
        renderer.setStroke(Color{0, 0, 0}, 2.0);
        for (auto eit : m_forcer.m_g.edges()) {
            renderer.draw(eit->curve().transform(inv));
        }
        for (auto vit : m_forcer.m_g.vertices()) {
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
            const auto& col = d.collapse->result;
            if (auto* col3P = std::get_if<curved_simplification::detail::Collapse3>(&col)) {
                auto& col3 = *col3P;
                renderer.draw(col3.before.transform(inv));
                renderer.draw(col3.after.transform(inv));
                auto bb = col3.before.transform(inv).bbox();
                auto bbC = Point<Inexact>((bb.xmin() + bb.xmax()) / 2, (bb.ymin() + bb.ymax()) / 2);
                renderer.drawText(bbC, std::to_string(d.collapse->cost));
            } else if (auto* col2P = std::get_if<curved_simplification::detail::Collapse2>(&col)) {
                auto& col2 = *col2P;
                renderer.draw(col2.replacement.transform(inv));
                auto bb = col2.replacement.transform(inv).bbox();
                auto bbC = Point<Inexact>((bb.xmin() + bb.xmax()) / 2, (bb.ymin() + bb.ymax()) / 2);
                renderer.drawText(bbC, std::to_string(d.collapse->cost));
            }
        }
        if (d.hist == nullptr) return;
        auto opC = std::dynamic_pointer_cast<curved_simplification::detail::CollapseEdgeOperation<BaseGraph>>(d.hist);
        if (opC != nullptr) {
            renderer.setStroke(Color(200, 50, 200), 2.0);
            renderer.draw(opC->m_c0.transform(inv));
            renderer.draw(opC->m_c1.transform(inv));
            renderer.draw(opC->m_c2.transform(inv));
        }
        auto opM = std::dynamic_pointer_cast<curved_simplification::detail::MergeEdgeOperation<BaseGraph>>(d.hist);
        if (opM != nullptr) {
            renderer.setStroke(Color(200, 50, 200), 2.0);
            renderer.draw(opM->m_c.transform(inv));
            renderer.draw(opM->m_cPrev.transform(inv));
        }
    }, "Debug edge");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        auto inv = m_transform.inverse();
        if (!m_editControlPoints->isChecked()) return;

        auto mp = m_renderer->mousePosition();

        int nearestArc = -1;
        double minD2 = std::numeric_limits<double>::infinity();

        if (m_editHighlight->isChecked()) {
            // Find nearest control polyline
            for (auto eit : m_editBaseGraph.edges()) {
                auto mpT = mp.transform(m_transform);
                Polyline<Inexact> pl;
                for (int c = 0; c < 4; ++c) pl.push_back(eit->curve().control(c));
                auto n = pl.nearest(mpT);
                auto d2 = CGAL::squared_distance(n, mpT);
                if (d2 < minD2) {
                    nearestArc = eit->data().index;
                    minD2 = d2;
                }
            }

            renderer.setStrokeOpacity(255);
            renderer.setStroke(Color(0, 0, 0), 4);
            renderer.setMode(GeometryRenderer::stroke);
            for (auto eit : m_editBaseGraph.edges()) {
                if (eit->data().index == nearestArc) {
                    renderer.draw(eit->curve().transform(inv));
                }
            }
        }

        // Control polylines
        for (auto eit : m_editBaseGraph.edges()) {
            renderer.setStroke(Color(0, 255, 0), 1.0);
            Polyline<Inexact> pl;
            for (int c = 0; c < 4; ++c) pl.push_back(eit->curve().control(c));
            auto opacity = nearestArc == -1 || eit->data().index == nearestArc ? 255 : 50;
            renderer.setStrokeOpacity(opacity);
            renderer.draw(pl.transform(inv));
        }

        for (const auto& editable : m_editables) {
            auto point = editable->point().transform(inv);

            int opacity = 50;
            if (auto vhP = std::get_if<Vertex_handle>(&editable->type)) {
                renderer.setStroke(Color(255, 0, 255), 2.0);
                auto vh = *vhP;
                if (vh->degree() == 2 && (nearestArc == -1 || vh->outgoing()->data().index == nearestArc)) {
                    opacity = 255;
                }
            } else if (auto eiP = std::get_if<std::pair<Edge_handle, bool>>(&editable->type)) {
                renderer.setStroke(Color(0, 0, 255), 2.0);
                auto [eh, b] = *eiP;
                if (nearestArc == -1 || eh->data().index == nearestArc) {
                    opacity = 255;
                }
            }
            renderer.setStrokeOpacity(opacity);
            renderer.draw(point);
        }
        renderer.setStrokeOpacity(255);
    }, "Control points");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        if (!m_shiftDown && !m_altDown) return;
        if (!m_editControlPoints->isChecked()) return;
        auto mp = m_renderer->mousePosition();
        
        auto inv = m_transform.inverse();

        if (m_shiftDown) {
            // naive...
            std::optional<BaseGraph::Edge_handle> closest;
            std::optional<CubicBezierCurve::CurvePoint> nearest;
            double minDist = std::numeric_limits<double>::infinity();
            for (auto eit : m_editBaseGraph.edges()) {
                auto n = eit->curve().nearest(mp.transform(m_transform));
                auto d2 = screenDistanceToMouse2(n.point.transform(inv));
                
                if (d2 < minDist && d2 < 400) {
                    minDist = d2;
                    nearest = n;
                    closest = eit;
                }
            }
            if (!closest.has_value()) return;

            m_shiftNearest = {*closest, *nearest};

            renderer.setStroke(Color(255, 0, 0), 1.0);
            renderer.draw(nearest->point.transform(m_transform.inverse()));
        }
        if (m_altDown) {
            // naive...
            std::optional<BaseGraph::Vertex_handle> closest;
            double minDist = std::numeric_limits<double>::infinity();
            for (auto vit : m_editBaseGraph.vertices()) {
                if (vit->degree() != 2) continue;
                auto d2 = screenDistanceToMouse2(vit->point().transform(inv));
                if (d2 < minDist && d2 < 400) {
                    minDist = d2;
                    closest = vit;
                }
            }
            if (!closest.has_value()) return;

            m_altNearest = *closest;

            renderer.setStroke(Color(255, 0, 0), 1.0);
            renderer.draw((*closest)->point().transform(m_transform.inverse()));
        }
    }, "Cut point");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        m_intersectionsPainting.paint(renderer);
    }, "Intersections");
}

void BezierSimplificationDemo::updateEditables() {
    m_editables.clear();
    for (auto eit : m_editBaseGraph.edges()) {
        auto sourceControlEditable = std::make_shared<ControlPoint>(std::pair(eit, false));
        m_editables.push_back(sourceControlEditable);
        auto targetControlEditable = std::make_shared<ControlPoint>(std::pair(eit, true));
        m_editables.push_back(targetControlEditable);
    }
    for (auto vit : m_editBaseGraph.vertices()) {
        auto endpointEditable = std::make_shared<ControlPoint>(vit);
        m_editables.push_back(endpointEditable);
    }
}

BezierSimplificationDemo::BezierSimplificationDemo() : m_graph(m_baseGraph), m_collapse(m_graph, Traits()), m_forcer(m_approxGraph, 1.0), m_editGraph(m_editBaseGraph) {
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

    connect(m_renderer, &GeometryWidget::clicked, [this](Point<Inexact> pt) {
        if (m_shiftDown && m_shiftNearest.has_value()) {
            auto& [e, cp] = *m_shiftNearest;
            auto [toPoint, fromPoint] = e->curve().split(cp.t);
            auto edgeData = e->data();
            m_editGraph.startBatch();
            auto v = m_editGraph.splitEdge(e, toPoint, fromPoint);
            v->incoming()->data().index = edgeData.index;
            v->outgoing()->data().index = edgeData.index;
            m_editGraph.endBatch();

            m_shiftNearest = std::nullopt;

            baseModified(true);
            updateEditables();
            m_renderer->repaint();
            return;
        }
        if (m_altDown && m_altNearest.has_value()) {
            auto vh = *m_altNearest;
            auto inc = vh->incoming()->curve();
            auto out = vh->outgoing()->curve();
            CubicBezierSpline spline;
            spline.appendCurve(inc);
            spline.appendCurve(out);
            auto pl = spline.polyline(100);
            std::vector<Point<Inexact>> points(pl.vertices_begin(), pl.vertices_end());
            auto newCurve = fitCurve(points, inc.tangent(0), -out.tangent(1));
            auto edgeData = vh->outgoing()->data();
            m_editGraph.startBatch();
            auto eh = m_editGraph.merge_edge_with_prev(vh->outgoing(), newCurve);
            eh->data().index = edgeData.index;
            m_editGraph.endBatch();

            m_altNearest = std::nullopt;

            baseModified(true);
            updateEditables();
            m_renderer->repaint();

            return;
        }
#if DEBUG
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
#endif
    });

    addPaintings();
}

void BezierSimplificationDemo::checkIntersections() {
    Rectangle<Inexact> rect(0, 0, 1000, 1000);
    BezierCurveQuadTree<BaseGraph> bcqt(rect, 10, 0.05);

    m_intersectionsPainting = {};
    GeometryRenderer& r = m_intersectionsPainting;
    r.setMode(GeometryRenderer::stroke);
    r.setStroke(Color(255, 0, 0), 2.0);

    auto inv = m_transform.inverse();

    auto& theBaseGraph = m_editBaseGraph.number_of_edges() > 0 ? m_editBaseGraph : m_baseGraph;

    for (auto eit : theBaseGraph.edges()) {
        if (eit->curve().selfIntersects()) {
            std::cout << "Found self-intersection!" << std::endl;
            r.draw(eit->curve().transform(inv));
        }
        bcqt.insert(eit);
    }

    for (auto eit : theBaseGraph.edges()) {
        auto box = eit->curve().bbox();
        Rectangle<Inexact> rect(box.xmin(), box.ymin(), box.xmax(), box.ymax());
        bcqt.findOverlapped(rect, [&](const BaseGraph::Edge_handle& other) {
            if (other == eit) return false;

            if (other->curve().sub(0.01, 0.99).intersects(eit->curve())) {
                std::cout << "Found intersection! " << other->curve().source().transform(inv) << " -> " << other->curve().target().transform(inv) << " and " << eit->curve().source().transform(inv) << " -> " << eit->curve().target().transform(inv) << std::endl;
                r.draw(eit->curve().transform(inv));
                r.draw(other->curve().transform(inv));
                return true;
            }
            return false;
            });
    }
}

void BezierSimplificationDemo::keyPressEvent(QKeyEvent *event) {
    if (event->modifiers().testFlag(Qt::ShiftModifier)) {
        m_shiftDown = true;
        m_renderer->repaint();
    }
    if (event->modifiers().testFlag(Qt::AltModifier)) {
        m_altDown = true;
        m_renderer->repaint();
    }
}

void BezierSimplificationDemo::keyReleaseEvent(QKeyEvent *event) {
    if (event->text() == "") {
        m_shiftDown = false;
        m_altDown = false;
        m_renderer->repaint();
    }
}

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
	BezierSimplificationDemo demo;
	demo.show();
	app.exec();
}
