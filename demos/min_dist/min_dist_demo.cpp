#include "min_dist_demo.h"

#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QPushButton>

#include "library/read_graph_gdal.h"

#include <cartocrow/renderer/voronoi_drawer.h>
#include <cartocrow/core/transform_helpers.h>
#include <cartocrow/circle_segment_helpers/cs_types.h>

#include <cartocrow/renderer/ipe_renderer.h>

MinDistDemo::MinDistDemo() : m_forcer(m_g, 1.0) {
    setWindowTitle("Minimum distance");
    m_renderer = new GeometryWidget();
    m_renderer->setDrawAxes(false);
    setCentralWidget(m_renderer);

    auto *dockWidget = new QDockWidget();
    addDockWidget(Qt::RightDockWidgetArea, dockWidget);
    auto *vWidget = new QWidget();
    auto *vLayout = new QVBoxLayout(vWidget);
    vLayout->setAlignment(Qt::AlignTop);
    dockWidget->setWidget(vWidget);

    auto* stepButton = new QPushButton("Step");
    vLayout->addWidget(stepButton);

    auto* step100Button = new QPushButton("Step (x100)");
    vLayout->addWidget(step100Button);

    auto* resetButton = new QPushButton("Reset");
    vLayout->addWidget(resetButton);

    m_minDist = new DoubleSlider(Qt::Horizontal);
    m_minDist->setMaximum(30);
    m_minDist->setMinimum(0);
    m_minDist->setValue(6);

    m_forcer.m_requiredMinDist = m_minDist->value();
    vLayout->addWidget(m_minDist);

    connect(stepButton, &QPushButton::clicked, [this]() {
        m_forcer.step();
        m_renderer->repaint();
    });

    connect(step100Button, &QPushButton::clicked, [this]() {
        for (int i = 0; i < 100; ++i) {
            m_forcer.step();
        }
        m_renderer->repaint();
    });

    connect(resetButton, &QPushButton::clicked, [this]() {
        m_g = m_ogg;
        m_forcer.recomputeDelaunay();
        m_forcer.recomputeAuxiliary();
        m_renderer->repaint();
    });

    connect(m_minDist, &DoubleSlider::valueChanged, [this]() {
        m_renderer->repaint();
        m_forcer.m_requiredMinDist = m_minDist->value();
        m_forcer.recomputeAuxiliary();
    });

    m_g = readGraphUsingGDAL("data/small_min_dist_test/small_min_dist_test.shp");

    std::vector<Point<Inexact>> points;
    for (auto vit = m_g.vertices_begin(); vit != m_g.vertices_end(); ++vit) {
        points.push_back(vit->point());
    }
    Box bbox = CGAL::bbox_2(points.begin(), points.end());
    auto trans = fitInto(bbox, Box(0, 0, 1000, 1000));

    auto t = m_g.transform(trans);
    m_g = t;
    m_g.orient();
    m_ogg = m_g;

    m_renderer->setMaxZoom(1000);
    m_renderer->setMinZoom(0.01);

    m_forcer.m_g = m_g;
    m_forcer.recomputeAuxiliary();
    m_forcer.initialize();

    m_renderer->fitInView(Box(0, 0, 1000, 1000));

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::fill);
        renderer.setFill(Color(230, 230, 230));
        renderer.draw(Circle<Inexact>(m_renderer->mousePosition(), m_minDist->value() * m_minDist->value()));
    }, "Disk");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::stroke);
        auto voronoiDrawer = VoronoiDrawer<MinimumDistanceForcer<std::monostate, std::monostate>::Gt>(&renderer);
        for (auto eit = m_forcer.m_delaunay.finite_edges_begin(); eit != m_forcer.m_delaunay.finite_edges_end(); ++eit) {
            auto [ps, qs] = m_forcer.defining_storage_sites(*eit);
            if (m_forcer.withinDistanceAlongIsoline(ps, qs, 2 * m_forcer.m_requiredMinDist)) {
                continue;
            }

            renderer.setStroke(Color{210, 210, 210}, 1.0);
            draw_dual_edge(m_forcer.m_delaunay, *eit, voronoiDrawer);
        }

        //for (const auto& [vEdge, dEdge] : m_forcer.m_withinDistEdges) {
        //    std::visit([&](const auto &geom) {
        //        renderer.setStroke(Color{255, 50, 50}, 2.0);
        //        voronoiDrawer << geom;
        //    }, vEdge);
        //}
    }, "Segment Voronoi diagram");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::stroke);
        renderer.setStroke(Color(0, 0, 0), 3.0);
        for (auto vit = m_g.vertices_begin(); vit != m_g.vertices_end(); ++vit) {
            renderer.draw(vit->point());
        }
        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            renderer.draw(Segment<Inexact>(eit->source()->point(), eit->target()->point()));
        }
    }, "Graph");

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::stroke);
        renderer.setStroke(Color(100, 255, 100), 3.0);
        for (const auto& part : m_forcer.m_withinDistIsolineParts) {
            std::visit([&](const auto& geom) {
                renderer.draw(geom);
            }, part);
        }
    }, "Isoline parts");
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    MinDistDemo demo;
    demo.show();
    app.exec();
}