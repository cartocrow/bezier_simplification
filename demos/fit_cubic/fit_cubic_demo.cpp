#include "fit_cubic_demo.h"

#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>

#include "../read_ipe_bezier_spline.h"

FitCubicDemo::FitCubicDemo() {
    setWindowTitle("Fit cubic");
    m_renderer = new GeometryWidget();
    m_renderer->setDrawAxes(false);
    setCentralWidget(m_renderer);

    auto *dockWidget = new QDockWidget();
    addDockWidget(Qt::RightDockWidgetArea, dockWidget);
    auto *vWidget = new QWidget();
    auto *vLayout = new QVBoxLayout(vWidget);
    vLayout->setAlignment(Qt::AlignTop);
    dockWidget->setWidget(vWidget);

//    std::vector<Point<Inexact>> points = {
//        {0.0, 1.0},
//        {0.5, 2.0},
//        {1.0, 3.0},
//        {2.0, 3.5},
//        {3.0, 3.0},
//        {4.0, 2.0},
//        {5.0, 1.0},
//        {6.0, 0.5},
//        {7.0, 1.0},
//    };

    std::vector<Point<Inexact>> points = {
            {143.272, 394.30},
            {143.231, 394.45},
            {143.189, 394.60},
            {143.148, 394.74},
            {143.106, 394.89},
            {143.065, 395.03},
            {143.024, 395.17},
            {142.982, 395.31},
            {142.941, 395.46},
            {142.9,   395.608},
            {142.859, 395.75},
            {142.819, 395.91},
            {142.778, 396.07},
            {142.738, 396.23},
            {142.698, 396.40},
            {142.658, 396.58},
            {142.619, 396.77},
            {142.58,  396.97},
            {142.542, 397.17},
            {142.504, 397.39},
            {142.466, 397.62},
            {142.41,  397.993},
            {142.368, 398.32},
            {142.338, 398.60},
            {142.319, 398.85},
            {142.308, 399.07},
            {142.305, 399.26},
            {142.307, 399.43},
            {142.313, 399.57},
            {142.321, 399.71},
            {142.331, 399.83},
            {142.339, 399.94},
            {142.344, 400.06},
            {142.345, 400.17},
            {142.34,  400.3},
            {142.328, 400.43},
            {142.307, 400.58},
            {142.274, 400.75},
            {142.23,  400.946},
            {142.171, 401.16},
            {142.096, 401.42}
    };

    auto spline = fitCurve(points);

    m_renderer->addPainting([points, spline](GeometryRenderer& renderer) {
        for (const auto& pt : points) {
            renderer.draw(pt);
        }
        renderer.setMode(GeometryRenderer::stroke);
        renderer.setStroke(Color(0, 0, 0), 3.0);
        renderer.draw(spline);
    }, "Fit cubic");
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    FitCubicDemo demo;
    demo.show();
    app.exec();
}