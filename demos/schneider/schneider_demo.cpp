#include "schneider_demo.h"

#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>

#include "frontend/read_ipe_bezier_spline.h"

SchneiderDemo::SchneiderDemo() {
    setWindowTitle("Demo for Schneider's algorithm to fit Béziers");
    m_renderer = new GeometryWidget();
    m_renderer->setDrawAxes(false);
    setCentralWidget(m_renderer);

    auto *dockWidget = new QDockWidget();
    addDockWidget(Qt::RightDockWidgetArea, dockWidget);
    auto *vWidget = new QWidget();
    auto *vLayout = new QVBoxLayout(vWidget);
    vLayout->setAlignment(Qt::AlignTop);
    dockWidget->setWidget(vWidget);

    std::vector<Point<Inexact>> points = {
        {0.0, 1.0},
        {0.5, 2.0},
        {1.0, 3.0},
        {2.0, 3.5},
        {3.0, 3.0},
        {4.0, 2.0},
        {5.0, 1.0},
        {6.0, 0.5},
        {7.0, 1.0},
    };

    auto cubic = fitCurve(points, 1000);
    auto spline = fitSpline(points, 0.001, 1000);

    m_renderer->fitInView(cubic.bbox());

    m_renderer->addPainting([points](GeometryRenderer& renderer) {
        for (const auto& pt : points) {
            renderer.draw(pt);
        }
    }, "Points");

    m_renderer->addPainting([cubic](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::stroke);
        renderer.setStroke(Color(0, 0, 0), 3.0);
        renderer.draw(cubic);
    }, "Cubic curve");

    m_renderer->addPainting([spline](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::stroke);
        renderer.setStroke(Color(0, 0, 0), 3.0);
        renderer.draw(spline);
    }, "Cubic spline");
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    SchneiderDemo demo;
    demo.show();
    app.exec();
}