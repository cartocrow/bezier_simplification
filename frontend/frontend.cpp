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

#include "double_slider.h"

#include "library/read_graph_gdal.h"
#include "read_ipe_bezier_spline.h"

#include <cartocrow/renderer/voronoi_drawer.h>
#include <cartocrow/core/transform_helpers.h>

void saveGraphIntoTopoSet(const BaseGraph& graph, TopoSet<Inexact>& topoSet) {
    std::unordered_set<const BaseGraph::Edge*> visited;

    //int nSegs = std::max(5.0, (double)31998 / graph.number_of_edges());

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
                //arc.pop_back();
                //arcEit->curve().samplePoints(nSegs, std::back_inserter(arc));
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
        start->curve().samplePoints(nPoints, std::back_inserter(arc));

        if (start != end) {
            for (auto arcEit = start->next(); true; arcEit = arcEit->next()) {
                arc.pop_back();
                arcEit->curve().samplePoints(nPoints, std::back_inserter(arc));
                if (arcEit == end) break;
            }
        }

        topoSet.arcs[arcIndex] = arc;
    }
}

void BezierSimplificationDemo::loadInput(const std::filesystem::path& path) {
    m_baseGraph.clear();
    m_debugEdge = std::nullopt;

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
//                    eh->data().

                }
            }
        }
    }

    m_baseGraph.orient();
    m_original = m_baseGraph;
    m_collapse.initialize();

    auto bbox = m_baseGraph.bbox();
    m_referencePolygon = bbox;

    m_renderer->fitInView(bbox);
}

BezierSimplificationDemo::BezierSimplificationDemo() : m_graph(m_baseGraph), m_collapse(m_graph, Traits()), m_forcer(m_approxGraph, 1.0) {
	setWindowTitle("Bézier simplification");
	m_renderer = new GeometryWidget();
	m_renderer->setDrawAxes(false);
    m_renderer->setMaxZoom(1000000000);
    m_renderer->setMinZoom(0.00000001);
	setCentralWidget(m_renderer);

	auto* dockWidget = new QDockWidget();
	addDockWidget(Qt::RightDockWidgetArea, dockWidget);
	auto* vWidget = new QWidget();
	auto* vLayout = new QVBoxLayout(vWidget);
	vLayout->setAlignment(Qt::AlignTop);
	dockWidget->setWidget(vWidget);

    auto* ioSettings = new QLabel("<h3>Input / output</h3>");
    vLayout->addWidget(ioSettings);

    auto* loadFileButton = new QPushButton("Load file");
    vLayout->addWidget(loadFileButton);

    auto* exportButton = new QPushButton("Export");
    vLayout->addWidget(exportButton);

    auto* stackPolygons = new QCheckBox("Export: stacked polygons instead of holes");
    stackPolygons->setChecked(false);
    vLayout->addWidget(stackPolygons);

    auto* editControlPoints = new QCheckBox("Edit control points");
    editControlPoints->setChecked(false);
    vLayout->addWidget(editControlPoints);

    auto* editAlignTangents = new QCheckBox("Edit: align tangents");
    editAlignTangents->setChecked(true);
    vLayout->addWidget(editAlignTangents);

    auto* loadReferenceDataButton = new QPushButton("Load reference data");
    vLayout->addWidget(loadReferenceDataButton);

    auto* loadReferencePolygonButton = new QPushButton("Load reference polygon for image");
    vLayout->addWidget(loadReferencePolygonButton);

    auto* simplificationSettings = new QLabel("<h3>Simplification</h3>");
    vLayout->addWidget(simplificationSettings);

	auto* stepButton = new QPushButton("Step");
	vLayout->addWidget(stepButton);

    auto* step10Button = new QPushButton("Step (x10)");
    vLayout->addWidget(step10Button);

    auto* step100Button = new QPushButton("Step (x100)");
    vLayout->addWidget(step100Button);

    auto* runToComplexityButton = new QPushButton("Run to specified complexity");
    vLayout->addWidget(runToComplexityButton);
    auto* desiredComplexity = new QSpinBox();
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

	loadInput("splines.ipe");

    auto* complexityLabel = new QLabel("#Edges: ");
    auto* complexity = new QSlider();
    auto* complexityLog = new DoubleSlider();
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

    auto* minimumDistanceSettings = new QLabel("<h3>Minimum distance</h3>");
    vLayout->addWidget(minimumDistanceSettings);

    m_minDist = new DoubleSlider(Qt::Horizontal);
    m_minDist->setMaximum(30);
    m_minDist->setMinimum(0);
    m_minDist->setValue(6);

    m_forcer.m_requiredMinDist = m_minDist->value();
    vLayout->addWidget(m_minDist);

    auto* mdInitializeButton = new QPushButton("Initialize / reset");
    vLayout->addWidget(mdInitializeButton);
    auto* mdStepButton = new QPushButton("Step");
    vLayout->addWidget(mdStepButton);
    auto* mdStep100Button = new QPushButton("Step (x100)");
    vLayout->addWidget(mdStep100Button);
    auto* mdReconstructButton = new QPushButton("Reconstruct");
    vLayout->addWidget(mdReconstructButton);

    auto* drawSettings = new QLabel("<h3>Drawing</h3>");
    vLayout->addWidget(drawSettings);

    auto* showEdgeDirection = new QCheckBox("Show edge direction");
    vLayout->addWidget(showEdgeDirection);
    auto* showOldVertices = new QCheckBox("Show old vertices");
    vLayout->addWidget(showOldVertices);
    auto* showNewVertices = new QCheckBox("Show new vertices");
    vLayout->addWidget(showNewVertices);
    auto* showNewControlPoints = new QCheckBox("Show new control points");
    vLayout->addWidget(showNewControlPoints);
    auto* showDebugInfo = new QCheckBox("Show debug info");
    vLayout->addWidget(showDebugInfo);

    /*auto* queueInfo = new QLabel("<h3>Queue</h3>");
    vLayout->addWidget(queueInfo);
    auto* scrollArea = new QScrollArea();
    auto* queueLabel = new QLabel();
    scrollArea->setWidget(queueLabel);
    scrollArea->setWidgetResizable(true);
    scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    scrollArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    scrollArea->setAlignment(Qt::AlignTop);
    queueLabel->setAlignment(Qt::AlignTop);
    vLayout->addWidget(scrollArea);*/

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

    connect(m_renderer, &GeometryWidget::clicked, [this](Point<Inexact> pt) {
        double minDist2Edge = std::numeric_limits<double>::infinity();
        std::optional<Graph::Edge_handle> closest;
        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            double minDist2 = CGAL::squared_distance(pt, eit->curve().nearest(pt).point);
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

    auto updateQueueInfo = [this/*, queueLabel*/]() {
        /*std::stringstream ss;
        if (!m_collapse.m_q.empty()) {
            auto top = m_collapse.m_q.pop();
            ss << "> " << top->data().collapse->cost << "\n";
            m_collapse.m_q.push(top);
            std::vector<double> costs;
            for (auto eit = m_graph.edges_begin(); eit != m_graph.edges_end(); ++eit) {
                auto clps = eit->data().collapse;
                if (clps.has_value()) {
                    costs.push_back(clps->cost);
                }
            }
            std::sort(costs.begin(), costs.end());
            for (const auto& cost : costs) {
                ss << cost << "\n";
            }
        }
        queueLabel->setText(QString::fromStdString(ss.str()));*/
    };

    auto updateComplexityInfo = [complexityLabel, complexity, complexityLog, desiredComplexity, this]() {
        m_backup = std::nullopt;
        complexity->setMinimum(std::min((int) m_graph.number_of_edges(), complexity->minimum()));
        complexity->setValue(m_graph.number_of_edges());
        complexityLabel->setText(QString::fromStdString("#Edges: " + std::to_string(m_graph.number_of_edges())));
        complexityLog->setMinimum(std::min(log(m_graph.number_of_edges() + 0.5), complexityLog->minimum()));
        complexityLog->setValue(log(m_graph.number_of_edges() + 0.5));
    };

    auto resetEdits = [editControlPoints, this]() {
        editControlPoints->setChecked(false);

        // todo fix
//        if (m_backup.has_value()) {
//            m_baseGraph = *m_backup;
//        }
    };

    connect(loadFileButton, &QPushButton::clicked, [this, complexity, updateComplexityInfo, complexityLog, desiredComplexity]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getOpenFileName(this, tr("Select input file (.ipe or .shp)"), startDir).toStdString();
        if (filePath == "") return;
        loadInput(filePath);

        complexity->setMaximum(m_graph.number_of_edges());
        complexity->setMinimum(m_graph.number_of_edges());
        complexity->setValue(m_graph.number_of_edges());
        complexityLog->setMaximum(log(m_graph.number_of_edges() + 0.5));
        complexityLog->setMinimum(log(m_graph.number_of_edges() + 0.5));
        complexityLog->setValue(log(m_graph.number_of_edges() + 0.5));
        desiredComplexity->setMaximum(m_graph.number_of_edges());

        updateComplexityInfo();
    });

    connect(loadReferenceDataButton, &QPushButton::clicked, [this]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getOpenFileName(this, tr("Select reference file (.tif or .shp)"), startDir).toStdString();
        if (filePath == "") return;
        if (filePath.extension() == ".tif") {
            m_referenceData.push_back(QImage(filePath.c_str()));
        } else {
            auto [regionSet, proj] = readRegionSetUsingGDAL(filePath);
            m_referenceData.push_back(regionSet);
        }
    });

    connect(loadReferencePolygonButton, &QPushButton::clicked, [this]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getOpenFileName(this, tr("Select reference polygon for image (.shp)"), startDir).toStdString();
        if (filePath == "") return;
        auto regionSet = readRegionSetUsingGDAL(filePath);
        std::vector<PolygonWithHoles<Inexact>> pgns;
        regionSet.first[0].geometry.polygons_with_holes(std::back_inserter(pgns));
        Rectangle<Inexact> rect = pgns[0].bbox();
        auto rectT = rect.transform(m_transform);
        m_referencePolygon = Box(rectT.xmin(), rectT.ymin(), rectT.xmax(), rectT.ymax());
    });

    connect(exportButton, &QPushButton::clicked, [this, stackPolygons]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getSaveFileName(this, tr("Select isolines"), startDir).toStdString();
        if (filePath.empty()) return;
        saveGraphIntoTopoSet(m_baseGraph, m_toposet);
        exportTopoSetUsingGDAL(filePath, m_toposet, m_transform.inverse(), m_spatialRef, stackPolygons->isChecked());
    });

    connect(complexity, &QSlider::valueChanged, [this, updateComplexityInfo, resetEdits](int value) {
        resetEdits();
        m_graph.recallComplexity(value);
        m_renderer->repaint();
        updateComplexityInfo();
    });

    connect(complexityLog, &DoubleSlider::valueChanged, [this, updateComplexityInfo, resetEdits, complexityLog](double value) {
        resetEdits();
        m_graph.recallComplexity(std::exp(value));
        m_renderer->repaint();
        updateComplexityInfo();
    });

	connect(stepButton, &QPushButton::clicked, [this, updateQueueInfo, updateComplexityInfo, resetEdits]() {
        resetEdits();
		m_collapse.step();
		m_renderer->repaint();
        updateQueueInfo();
        updateComplexityInfo();
	});

    connect(step10Button, &QPushButton::clicked, [this, updateComplexityInfo, resetEdits, updateQueueInfo]() {
        resetEdits();
        for (int i = 0; i < 10; ++i) m_collapse.step();
		m_renderer->repaint();
        updateQueueInfo();
        updateComplexityInfo();
	});

    connect(step100Button, &QPushButton::clicked, [this, updateComplexityInfo, resetEdits, updateQueueInfo]() {
        resetEdits();
        for (int i = 0; i < 100; ++i) m_collapse.step();
		m_renderer->repaint();
        updateQueueInfo();
        updateComplexityInfo();
	});

    connect(undoButton, &QPushButton::clicked, [this, complexity, resetEdits]() {
        resetEdits();
        m_graph.backInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(undo10Button, &QPushButton::clicked, [this, complexity, resetEdits]() {
        resetEdits();
        for (int i = 0; i < 10; ++i)
            m_graph.backInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(undo100Button, &QPushButton::clicked, [this, complexity, resetEdits]() {
        resetEdits();
        for (int i = 0; i < 100; ++i)
            m_graph.backInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(redoButton, &QPushButton::clicked, [this, complexity, resetEdits]() {
        resetEdits();
        m_graph.forwardInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(redo10Button, &QPushButton::clicked, [this, complexity, resetEdits]() {
        resetEdits();
        for (int i = 0; i < 10; ++i)
            m_graph.forwardInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(redo100Button, &QPushButton::clicked, [this, complexity, resetEdits]() {
        resetEdits();
        for (int i = 0; i < 100; ++i)
            m_graph.forwardInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(m_minDist, &DoubleSlider::valueChanged, [this]() {
        m_renderer->repaint();
        m_forcer.m_requiredMinDist = m_minDist->value();
        m_forcer.recomputeAuxiliary();
    });

    connect(mdInitializeButton, &QPushButton::clicked, [this]() {
        m_forcer.m_g = approximateBezierGraph(m_baseGraph, std::min(20.0, 1000000.0 / m_baseGraph.number_of_edges()));
        m_forcer.initialize();
        m_renderer->repaint();
    });

    connect(mdStepButton, &QPushButton::clicked, [this]() {
        m_forcer.step();
        m_renderer->repaint();
    });

    connect(mdStep100Button, &QPushButton::clicked, [this]() {
        for (int i = 0; i < 100; ++i) {
            m_forcer.step();
        }
        m_renderer->repaint();
    });

    connect(mdReconstructButton, &QPushButton::clicked, [this]() {
        m_baseGraph = reconstructBezierGraph(m_forcer.m_g);
        m_forcer.m_g.clear();
        m_forcer.m_withinDistEdges.clear();
        m_forcer.m_delaunay.clear();
        m_renderer->repaint();
    });

    connect(editControlPoints, &QCheckBox::stateChanged, [this]() {
        // todo
//        if (!m_backup.has_value())
//            m_backup = m_baseGraph;
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

    connect(runToComplexityButton, &QPushButton::clicked, [this, desiredComplexity, updateComplexityInfo]() {
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

    connect(m_renderer, &GeometryWidget::dragStarted, [this](const Point<Inexact>& p) {
        m_dragging = nullptr;

        std::optional<std::shared_ptr<ControlPoint>> closest;
        double minDist = std::numeric_limits<double>::infinity();

        for (auto editable : m_editables) {
            auto diff = m_renderer->convertPoint(editable->point) - m_renderer->convertPoint(p);
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

    connect(m_renderer, &GeometryWidget::dragMoved, [this, editAlignTangents](const Point<Inexact>& p) {
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

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
		renderer.setStroke(Color(200, 200, 200), 2.0);
        for (const auto& refData : m_referenceData) {
            if (auto qImageP = std::get_if<QImage>(&refData)) {
                auto& qImage = *qImageP;
                const auto formats = QImageReader::supportedImageFormats();
                if (!formats.contains("tiff")) {
                    std::cerr << "No support for tiff" << std::endl;
                }
                if (auto gw = dynamic_cast<GeometryWidget*>(&renderer)) {
                    std::cout << m_referencePolygon << std::endl;
//                    gw->drawImage(qImage.size())
                    gw->drawImage(m_referencePolygon, qImage);
                }
            } else if (auto regionSetP = std::get_if<RegionSet<Inexact>>(&refData)) {
                auto& regionSet = *regionSetP;
                for (const auto& region : regionSet) {
                    renderer.draw(region.geometry);
                }
            }
        }
	}, "Reference data");

	m_renderer->addPainting([showOldVertices, this](GeometryRenderer& renderer) {
		renderer.setStroke(Color(200, 200, 200), 2.0);
        for (auto eit = m_original.edges_begin(); eit != m_original.edges_end(); ++eit) {
            renderer.draw(eit->curve());
        }
        if (showOldVertices->isChecked()) {
            for (auto vit = m_original.vertices_begin(); vit != m_original.vertices_end(); ++vit) {
                renderer.draw(vit->point());
            }
        }
	}, "Original");

//    m_renderer->addPainting([this](GeometryRenderer& renderer) {
//        renderer.setMode(GeometryRenderer::stroke);
//        renderer.setStroke(Color(200, 0, 200), 2.0);
//        renderer.setStrokeOpacity(10);
//        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
//            int tSteps = 1000;
//            auto& curve = eit->curve();
//            for (int tStep = 0; tStep <= tSteps; ++tStep) {
//                double t = static_cast<double>(tStep) / tSteps;
//                auto c = curve.curvature(t);
//                auto n = curve.normal(t);
//                n /= sqrt(n.squared_length());
//                auto p = curve.position(t);
//                renderer.draw(Segment<Inexact>(p, p + 10 * c * n));
//            }
//        }
//    }, "Curvature");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::fill);
        renderer.setFill(Color(230, 230, 230));
        renderer.draw(Circle<Inexact>(m_renderer->mousePosition(), m_minDist->value() * m_minDist->value()));
    }, "Min. dist. disk");

	m_renderer->addPainting([this, showNewVertices, showNewControlPoints, showEdgeDirection, showDebugInfo](GeometryRenderer& renderer) {
	  	renderer.setMode(GeometryRenderer::stroke);
		renderer.setStroke(Color(0, 0, 0), 2.0);
		for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            if (showDebugInfo->isChecked() && eit->data().collapse.has_value()) {
                renderer.setStroke(Color(50, 200, 50), 2.0);
            } else {
                renderer.setStroke(Color(0, 0, 0), 2.0);
            }
			renderer.draw(eit->curve());
            if (showEdgeDirection->isChecked()) {
                renderer.setStroke(Color(0, 0, 0), 6.0);
                renderer.draw(eit->curve().split(0.25).first);
            }
		}
        // Control polylines
        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            if (showNewControlPoints->isChecked()) {
                renderer.setStroke(Color(0, 255, 0), 1.0);
                Polyline<Inexact> pl;
                for (int c = 0; c < 4; ++c) pl.push_back(eit->curve().control(c));
                renderer.draw(pl);
            }
        }
        // Control points
        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            if (showNewControlPoints->isChecked()) {
                renderer.setStroke(Color(0, 0, 255), 2.0);
                renderer.draw(eit->curve().sourceControl());
                renderer.draw(eit->curve().targetControl());
                renderer.setStroke(Color(255, 0, 255), 2.0);
                renderer.draw(eit->curve().source());
                renderer.draw(eit->curve().target());
            }
        }
        if (showNewVertices->isChecked()) {
            for (auto vit = m_baseGraph.vertices_begin(); vit != m_baseGraph.vertices_end(); ++vit) {
                renderer.draw(vit->point());
            }
        }
	}, "Simplification");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::stroke);
        auto voronoiDrawer = VoronoiDrawer<MinimumDistanceForcer<std::monostate, std::monostate>::Gt>(&renderer);
        for (auto eit = m_forcer.m_delaunay.finite_edges_begin(); eit != m_forcer.m_delaunay.finite_edges_end(); ++eit) {
            auto [ps, qs] = m_forcer.defining_storage_sites(*eit);
            if (m_forcer.withinDistanceAlongIsoline(ps, qs, 2 * m_forcer.m_requiredMinDist)) {
                continue;
            }
//
//            CGAL::Object o = m_forcer.m_delaunay.primal(eit);
            using Exact_SDG_traits = CGAL::Segment_Delaunay_graph_traits_2<Exact>;
            CGAL::Object o = exact_primal(*eit, m_forcer.m_delaunay);

            Segment<Exact> s;
            Line<Exact> l;
            Ray<Exact> r;
            CGAL::Parabola_segment_2<Exact_SDG_traits> parab;
            if (CGAL::assign(s, o)) {
//                if (!isfinite(s.source().x()) || !isfinite(s.target().x())) continue;
//                std::cout << "Segment: " << s.source() << " -> " << s.target() << std::endl;
//                if (s.source().)
                if (!do_overlap(m_forcer.m_bbox, s.bbox())) continue;
            }
            if (CGAL::assign(l, o)) {
//                continue;
//                if (!isfinite(l.a()) || !isfinite(l.b()) || !isfinite(l.c())) continue;
//                std::cout << "Line: " << l.a() << " " << l.b() << " " << l.c() << std::endl;
            }
            if (CGAL::assign(r, o)) {
//                continue;
//                if (!isfinite(r.source().x())) continue;
//                std::cout << "Ray: " << r.source() << " " << r.direction() << std::endl;

            }
            if (CGAL::assign(parab, o)) {
//                if (!isfinite(parab.p1.x()) || !isfinite(parab.p2.x()) || abs(parab.p1.x()) > 1E9 || abs(parab.p1.y()) > 1E9 || abs(parab.p2.x()) > 1E9 || abs(parab.p2.y()) > 1E9) continue;
//                std::cout << "Parabola segment: " << parab.p1 << " " << parab.p2 << " " << parab.center() << " " << parab.line() << std::endl;
//                if (!do_intersect(m_forcer.m_bbox, parab.p1) || !do_intersect(m_forcer.m_bbox, parab.p2)) continue;
            }

            renderer.setStroke(Color{210, 210, 210}, 1.0);
            draw_dual_edge_exact(m_forcer.m_delaunay, *eit, voronoiDrawer);
        }

        for (const auto& [vEdge, dEdge] : m_forcer.m_withinDistEdges) {
            std::visit([&](const auto &geom) {
                renderer.setStroke(Color{255, 50, 50}, 2.0);
                voronoiDrawer << geom;
            }, vEdge);
        }
    }, "Segment Voronoi diagram");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setStroke(Color{0, 0, 0}, 2.0);
        for (auto eit = m_forcer.m_g.edges_begin(); eit != m_forcer.m_g.edges_end(); ++eit) {
            renderer.draw(eit->curve());
        }
        for (auto vit = m_forcer.m_g.vertices_begin(); vit != m_forcer.m_g.vertices_end(); ++vit) {
            renderer.draw(vit->point());
        }
    }, "Linearized");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        if (!m_debugEdge.has_value()) return;
        renderer.setStroke(Color(255, 0, 0), 2.0);
        auto eh = *m_debugEdge;
        renderer.draw(eh->curve());
        const auto& d = eh->data();
        renderer.setStroke(Color(50, 200, 50), 2.0);
        if (d.collapse.has_value()) {
            const auto& col = *d.collapse;
            renderer.draw(col.before);
            renderer.draw(col.after);
            auto bb = col.before.bbox();
            auto bbC = Point<Inexact>((bb.xmin() + bb.xmax()) / 2, (bb.ymin() + bb.ymax())  / 2);
            renderer.drawText(bbC, std::to_string(col.cost));
        }
        if (d.hist == nullptr) return;
        auto op = std::dynamic_pointer_cast<curved_simplification::detail::CollapseEdgeOperation<BaseGraph>>(d.hist);
        renderer.setStroke(Color(200, 50, 200), 2.0);
        renderer.draw(op->m_c0);
        renderer.draw(op->m_c1);
        renderer.draw(op->m_c2);
    }, "Debug edge");

    m_renderer->addPainting([this, editControlPoints](GeometryRenderer& renderer) {
        if (!editControlPoints->isChecked()) return;

        // Control polylines
        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            renderer.setStroke(Color(0, 255, 0), 1.0);
            Polyline<Inexact> pl;
            for (int c = 0; c < 4; ++c) pl.push_back(eit->curve().control(c));
            renderer.draw(pl);
        }

        for (auto editable : m_editables) {
            if (auto vhP = std::get_if<Vertex_handle>(&editable->type)) {
                renderer.setStroke(Color(255, 0, 255), 2.0);
                renderer.draw(editable->point);
            } else if (auto eiP = std::get_if<std::pair<Edge_handle, bool>>(&editable->type)) {
                renderer.setStroke(Color(0, 0, 255), 2.0);
                auto [eh, b] = *eiP;
                if (b) {
                    renderer.draw(editable->point);
                } else {
                    renderer.draw(editable->point);
                }
            }
        }
    }, "Control points");
}

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
	BezierSimplificationDemo demo;
	demo.show();
	app.exec();
}
