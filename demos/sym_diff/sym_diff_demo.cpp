#include "sym_diff_demo.h"

#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QCheckBox>

#include <cartocrow/core/cubic_bezier.h>

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    SymDiffDemo demo;
    demo.show();
    app.exec();
}

SymDiffDemo::SymDiffDemo() {
    setWindowTitle("Symmetric diff demo");
    m_renderer = new GeometryWidget();
    m_renderer->setDrawAxes(false);
    setCentralWidget(m_renderer);

    auto p1 = std::make_shared<Point<Inexact>>(0, 0);
    auto p2 = std::make_shared<Point<Inexact>>(1, 3);
    auto p3 = std::make_shared<Point<Inexact>>(2, 0);
    auto p4 = std::make_shared<Point<Inexact>>(3, 1);

    m_renderer->registerEditable(p1);
    m_renderer->registerEditable(p2);
    m_renderer->registerEditable(p3);
    m_renderer->registerEditable(p4);

    auto q1 = std::make_shared<Point<Inexact>>(0, 0);
    auto q2 = std::make_shared<Point<Inexact>>(0, 2);
    auto q3 = std::make_shared<Point<Inexact>>(1, 5);
    auto q4 = std::make_shared<Point<Inexact>>(3, 1);

    m_renderer->registerEditable(q1);
    m_renderer->registerEditable(q2);
    m_renderer->registerEditable(q3);
    m_renderer->registerEditable(q4);

    m_renderer->fitInView(Box(-1, -1, 4, 2));

    m_renderer->addPainting([p1, p2, p3, p4, q1, q2, q3, q4](GeometryRenderer& renderer) {
        renderer.setStroke(Color(0, 0, 0), 2.0);
        CubicBezierSpline spline1;
        spline1.appendCurve(*p1, *p2, *p3, *p4);
        CubicBezierSpline spline2;
        spline2.appendCurve(*q1, *q2, *q3, *q4);

        std::vector<CubicBezierIntersectionResult<CubicBezierSpline, CubicBezierSpline>> inters;
        spline1.intersections(spline2, std::back_inserter(inters), M_EPSILON);
        
        renderer.setStroke(Color(255, 0, 0), 1.0);
        for (const auto& inter : inters) {
            renderer.draw(inter.point);
        }

        double sym_diff = 0;
        for (int i = 0; i < inters.size() - 1; ++i) {
            auto start = inters[i];
            auto end = inters[i + 1];
            
            auto part1 = spline1.sub(start.parameter1, end.parameter1);
            auto part2 = spline2.sub(start.parameter2, end.parameter2);
            part1.appendSpline(part2.reversed());

            renderer.setFill(Color(255, 0, 255));
            renderer.setMode(GeometryRenderer::fill);
            renderer.draw(part1);
            sym_diff += std::abs(part1.signedArea());
        }

        std::stringstream ss;
        ss << "Symmetric difference: " << sym_diff;
        renderer.setStroke(Color(0, 0, 0), 1.0);
        renderer.setMode(GeometryRenderer::stroke);
        renderer.drawText({ 5, 0 }, ss.str());
    }, "Intersections");

    m_renderer->addPainting([p1, p2, p3, p4](GeometryRenderer& renderer) {
        renderer.setStroke(Color(0, 0, 0), 2.0);
        CubicBezierSpline spline;
        spline.appendCurve(*p1, *p2, *p3, *p4);
        renderer.draw(spline);
        Polyline<Inexact> pl;
        pl.push_back(*p1);
        pl.push_back(*p2);
        pl.push_back(*p3);
        pl.push_back(*p4);
        renderer.setStroke(Color(100, 100, 100), 1.0);
        renderer.setMode(GeometryRenderer::stroke | GeometryRenderer::vertices);
        renderer.draw(pl);
    }, "Spline 1");

    m_renderer->addPainting([q1, q2, q3, q4](GeometryRenderer& renderer) {
        renderer.setStroke(Color(0, 0, 0), 2.0);
        CubicBezierSpline spline;
        spline.appendCurve(*q1, *q2, *q3, *q4);
        renderer.draw(spline);
        Polyline<Inexact> pl;
        pl.push_back(*q1);
        pl.push_back(*q2);
        pl.push_back(*q3);
        pl.push_back(*q4);
        renderer.setStroke(Color(100, 100, 100), 1.0);
        renderer.setMode(GeometryRenderer::stroke | GeometryRenderer::vertices);
        renderer.draw(pl);
    }, "Spline 2");
}