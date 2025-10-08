#include "bezier_simplification_demo.h"

#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QPushButton>
#include <QCheckBox>
#include <QLabel>

#include "read_ipe_bezier_spline.h"

Point<Inexact> nearest(const Segment<Inexact> seg, const Point<Inexact>& q)
{
    const auto& a = seg.source();
    const auto& b = seg.target();
    auto ab = b - a;
    auto aq = q - a;
    auto ba = a - b;
    auto bq = q - b;

    if (ab * aq <= 0) {
        return a;
    }

    if (ba * bq <= 0) {
        return b;
    }

    return seg.supporting_line().projection(q);
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

    auto* simplificationSettings = new QLabel("<h3>Simplification</h3>");
    vLayout->addWidget(simplificationSettings);

	auto* stepButton = new QPushButton("Step");
	vLayout->addWidget(stepButton);

    auto* undoButton = new QPushButton("Undo");
    vLayout->addWidget(undoButton);

    auto* redoButton = new QPushButton("Redo");
    vLayout->addWidget(redoButton);

	auto splines = ipeSplinesToIsolines("splines.ipe");

    for (const auto& spline : splines) {
        auto lastV = m_baseGraph.insert_vertex(spline.curves()[0].source());
        for (const auto &curve: spline.curves()) {
            auto targetV = m_baseGraph.insert_vertex(curve.target());
            m_baseGraph.add_edge(lastV, targetV, curve);
            lastV = targetV;
        }
    }

	m_baseGraph.orient();

    auto* complexityLabel = new QLabel("#Edges");
    auto* complexity = new QSlider();
    complexity->setMaximum(m_graph.num_edges());
	m_collapse.initialize();
	m_collapse.runToComplexity(10);
    complexity->setMinimum(m_graph.num_edges());
    complexity->setValue(m_graph.num_edges());

    complexity->setOrientation(Qt::Horizontal);
    vLayout->addWidget(complexityLabel);
    vLayout->addWidget(complexity);

    connect(complexity, &QSlider::valueChanged, [this](int value) {
        m_graph.recallComplexity(value);
        m_renderer->repaint();
    });

    auto* drawSettings = new QLabel("<h3>Drawing</h3>");
    vLayout->addWidget(drawSettings);

    auto* showEdgeDirection = new QCheckBox("Show edge direction");
    vLayout->addWidget(showEdgeDirection);
    auto* showVertices = new QCheckBox("Show vertices");
    vLayout->addWidget(showVertices);
    auto* showDebugInfo = new QCheckBox("Show debug info");
    vLayout->addWidget(showDebugInfo);

    connect(showEdgeDirection, &QCheckBox::stateChanged, [this]() {
        m_renderer->repaint();
    });
    connect(showVertices, &QCheckBox::stateChanged, [this]() {
        m_renderer->repaint();
    });
    connect(showDebugInfo, &QCheckBox::stateChanged, [this]() {
        m_renderer->repaint();
    });

    connect(m_renderer, &GeometryWidget::clicked, [this](Point<Inexact> pt) {
        double minDist2Edge = std::numeric_limits<double>::infinity();
        std::optional<Graph::Edge_handle> closest;
        for (auto eit = m_baseGraph.edges_begin(); eit != m_baseGraph.edges_end(); ++eit) {
            // find approximate distance to eit->curve()
            double minDist2 = std::numeric_limits<double>::infinity();
            auto pl = eit->curve().polyline(10);
            for (auto sit = pl.edges_begin(); sit != pl.edges_end(); ++sit) {
                auto q = nearest(*sit, pt);
                auto dist2 = CGAL::squared_distance(pt, q);
                if (dist2 < minDist2) {
                    minDist2 = dist2;
                }
            }
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
        Traits debugTraits(10, 20, true);
        debugTraits.determineCollapse(*m_debugEdge);
        repaint();
    });

	connect(stepButton, &QPushButton::clicked, [this, complexity]() {
		m_collapse.step();
        complexity->setValue(m_graph.num_edges());
		m_renderer->repaint();
	});

    connect(undoButton, &QPushButton::clicked, [this, complexity]() {
        m_graph.backInTime();
        complexity->setValue(m_graph.num_edges());
        m_renderer->repaint();
    });

    connect(redoButton, &QPushButton::clicked, [this, complexity]() {
        m_graph.forwardInTime();
        complexity->setValue(m_graph.num_edges());
        m_renderer->repaint();
    });

	m_renderer->addPainting([splines](GeometryRenderer& renderer) {
		renderer.setMode(GeometryRenderer::stroke | GeometryRenderer::vertices);
		renderer.setStroke(Color(200, 200, 200), 3.0);
        for (const auto& spline : splines) {
            renderer.draw(spline);
        }
	}, "Original");

	m_renderer->addPainting([this, showVertices, showEdgeDirection, showDebugInfo](GeometryRenderer& renderer) {
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
        if (showVertices->isChecked()) {
            for (auto vit = m_baseGraph.vertices_begin(); vit != m_baseGraph.vertices_end(); ++vit) {
                renderer.draw(vit->point());
            }
        }
	}, "Graph");

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
            auto [_1, _2, _3, top] = col.before.extrema();
            renderer.drawText(top.point, std::to_string(col.cost));
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
