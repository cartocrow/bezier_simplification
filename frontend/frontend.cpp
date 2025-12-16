#include "frontend.h"

#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QPushButton>
#include <QCheckBox>
#include <QSpinBox>
#include <QLabel>
#include <QScrollArea>
#include <QFileDialog>

#include "read_ipe_bezier_spline.h"

void BezierSimplificationDemo::loadInput(const std::filesystem::path& path) {
    m_splines = ipeSplinesToIsolines(path);

    m_baseGraph.clear();
    m_debugEdge = std::nullopt;

    std::unordered_map<Point<Inexact>, BaseGraph::Vertex_handle> pToV;

    auto getVertex = [&pToV, this](const Point<Inexact>& p) {
        if (pToV.contains(p)) {
            return pToV.at(p);
        } else {
            return m_baseGraph.insert_vertex(p);
        }
    };

    for (const auto& spline : m_splines) {
        if (spline.empty()) continue;
        auto lastV = getVertex(spline.curves()[0].source());
        pToV[spline.curves()[0].source()] = lastV;
        for (const auto &curve: spline.curves()) {
            auto targetV = getVertex(curve.target());
            pToV[curve.target()] = targetV;
            m_baseGraph.add_edge(lastV, targetV, curve);
            lastV = targetV;
        }
    }

    m_baseGraph.orient();
    m_collapse.initialize();

    std::vector<Point<Inexact>> points;
    for (auto vit = m_baseGraph.vertices_begin(); vit != m_baseGraph.vertices_end(); ++vit) {
        points.push_back(vit->point());
    }
    m_renderer->fitInView(CGAL::bbox_2(points.begin(), points.end()));
}

BezierSimplificationDemo::BezierSimplificationDemo() : m_graph(m_baseGraph), m_collapse(m_graph, Traits()) {
	setWindowTitle("Bézier simplification");
	m_renderer = new GeometryWidget();
	m_renderer->setDrawAxes(false);
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

    auto* redoButton = new QPushButton("Redo");
    vLayout->addWidget(redoButton);

	loadInput("splines.ipe");

    auto* complexityLabel = new QLabel("#Edges: ");
    auto* complexity = new QSlider();
    complexity->setMaximum(m_graph.number_of_edges());
    complexity->setMinimum(m_graph.number_of_edges());
    complexity->setValue(m_graph.number_of_edges());
    desiredComplexity->setMinimum(1);
    desiredComplexity->setMaximum(m_graph.number_of_edges());
    desiredComplexity->setValue(1);

    complexity->setOrientation(Qt::Horizontal);
    vLayout->addWidget(complexityLabel);
    vLayout->addWidget(complexity);

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

    auto* queueInfo = new QLabel("<h3>Queue</h3>");
    vLayout->addWidget(queueInfo);
    auto* scrollArea = new QScrollArea();
    auto* queueLabel = new QLabel();
    scrollArea->setWidget(queueLabel);
    scrollArea->setWidgetResizable(true);
    scrollArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    scrollArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    scrollArea->setAlignment(Qt::AlignTop);
    queueLabel->setAlignment(Qt::AlignTop);
    vLayout->addWidget(scrollArea);

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

    auto updateQueueInfo = [this, queueLabel]() {
        std::stringstream ss;
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
        queueLabel->setText(QString::fromStdString(ss.str()));
    };

    auto updateComplexityInfo = [complexityLabel, complexity, desiredComplexity, this]() {
        complexity->setMinimum(std::min((int) m_graph.number_of_edges(), complexity->minimum()));
        complexity->setValue(m_graph.number_of_edges());
        complexityLabel->setText(QString::fromStdString("#Edges: " + std::to_string(m_graph.number_of_edges())));
    };

    connect(loadFileButton, &QPushButton::clicked, [this, complexity, updateComplexityInfo]() {
        QString startDir = ".";
        std::filesystem::path filePath = QFileDialog::getOpenFileName(this, tr("Select isolines"), startDir).toStdString();
        if (filePath == "") return;
        loadInput(filePath);
        complexity->setMaximum(m_graph.number_of_edges());
        updateComplexityInfo();
    });

    connect(complexity, &QSlider::valueChanged, [this, updateComplexityInfo](int value) {
        m_graph.recallComplexity(value);
        m_renderer->repaint();
        updateComplexityInfo();
    });

	connect(stepButton, &QPushButton::clicked, [this, updateQueueInfo, updateComplexityInfo]() {
		m_collapse.step();
		m_renderer->repaint();
        updateQueueInfo();
        updateComplexityInfo();
	});

    connect(step10Button, &QPushButton::clicked, [this, updateComplexityInfo, updateQueueInfo]() {
        for (int i = 0; i < 10; ++i) m_collapse.step();
		m_renderer->repaint();
        updateQueueInfo();
        updateComplexityInfo();
	});

    connect(step100Button, &QPushButton::clicked, [this, updateComplexityInfo, updateQueueInfo]() {
        for (int i = 0; i < 100; ++i) m_collapse.step();
		m_renderer->repaint();
        updateQueueInfo();
        updateComplexityInfo();
	});

    connect(undoButton, &QPushButton::clicked, [this, complexity]() {
        m_graph.backInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(redoButton, &QPushButton::clicked, [this, complexity]() {
        m_graph.forwardInTime();
        complexity->setValue(m_graph.number_of_edges());
        m_renderer->repaint();
    });

    connect(runToComplexityButton, &QPushButton::clicked, [this, desiredComplexity, updateComplexityInfo]() {
        m_collapse.runToComplexity(desiredComplexity->value());
        updateComplexityInfo();
        m_renderer->repaint();
    });

	m_renderer->addPainting([showOldVertices, this](GeometryRenderer& renderer) {
        if (showOldVertices->isChecked()) {
            renderer.setMode(GeometryRenderer::stroke | GeometryRenderer::vertices);
        } else {
            renderer.setMode(GeometryRenderer::stroke);
        }
		renderer.setStroke(Color(200, 200, 200), 3.0);
        for (const auto& spline : m_splines) {
            renderer.draw(spline);
        }
	}, "Original");

//    m_renderer->addPainting([this](GeometryRenderer& renderer) {
//        renderer.setMode(GeometryRenderer::stroke);
//        renderer.setStroke(Color(200, 0, 200), 3.0);
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

	m_renderer->addPainting([this, showNewVertices, showNewControlPoints, showEdgeDirection, showDebugInfo](GeometryRenderer& renderer) {
	  	renderer.setMode(GeometryRenderer::stroke);
		renderer.setStroke(Color(0, 0, 0), 3.0);
		for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            if (showDebugInfo->isChecked() && eit->data().collapse.has_value()) {
                renderer.setStroke(Color(50, 200, 50), 3.0);
            } else {
                renderer.setStroke(Color(0, 0, 0), 3.0);
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
                renderer.setStroke(Color(0, 0, 255), 3.0);
                renderer.draw(eit->curve().sourceControl());
                renderer.draw(eit->curve().targetControl());
                renderer.setStroke(Color(255, 0, 255), 3.0);
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
        if (!m_debugEdge.has_value()) return;
        renderer.setStroke(Color(255, 0, 0), 3.0);
        auto eh = *m_debugEdge;
        renderer.draw(eh->curve());
        const auto& d = eh->data();
        renderer.setStroke(Color(50, 200, 50), 3.0);
        if (d.collapse.has_value()) {
            const auto& col = *d.collapse;
            renderer.draw(col.before);
            renderer.draw(col.after);
            auto bb = col.before.bbox();
            auto bbC = Point<Inexact>((bb.xmin() + bb.xmax()) / 2, (bb.ymin() + bb.ymax())  / 2);
            renderer.drawText(bbC, std::to_string(col.cost));
        }
        if (d.hist == nullptr) return;
        auto op = std::dynamic_pointer_cast<detail::CollapseEdgeOperation<BaseGraph>>(d.hist);
        renderer.setStroke(Color(200, 50, 200), 3.0);
        renderer.draw(op->m_c0);
        renderer.draw(op->m_c1);
        renderer.draw(op->m_c2);
    }, "Debug edge");
}

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
	BezierSimplificationDemo demo;
	demo.show();
	app.exec();
}
