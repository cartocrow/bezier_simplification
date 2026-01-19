#include "intersection_helpers_demo.h"

#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>

#include "library/intersection_helpers.h"

#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>

#include <cartocrow/core/segment_delaunay_graph_helpers.h>

typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<Inexact, CGAL::Field_with_sqrt_tag> Gt;

IntersectionHelpersDemo::IntersectionHelpersDemo() {
    setWindowTitle("Intersection helpers");
    m_renderer = new GeometryWidget();
    m_renderer->setDrawAxes(false);
    setCentralWidget(m_renderer);

    auto *dockWidget = new QDockWidget();
    addDockWidget(Qt::RightDockWidgetArea, dockWidget);
    auto *vWidget = new QWidget();
    auto *vLayout = new QVBoxLayout(vWidget);
    vLayout->setAlignment(Qt::AlignTop);
    dockWidget->setWidget(vWidget);

    CGAL::Aff_transformation_2<Inexact> trans(CGAL::ROTATION, sin(std::numbers::pi/4), cos(std::numbers::pi/4));
    CGAL::Parabola_segment_2<Gt> paraSeg(Point<Inexact>(0, 0.25).transform(trans), Line<Inexact>(0, 1, 0.25).transform(trans),
                                         Point<Inexact>(-2, 4).transform(trans), Point<Inexact>(-1, 1).transform(trans));

    Segment<Inexact> seg({-3, -1}, {4, 2});

    m_circle = std::make_shared<Circle<Inexact>>(Point<Inexact>(1.5, 0.5), 1.0);

    m_renderer->registerEditable(m_circle);

    m_renderer->addPainting([paraSeg, seg, this](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::stroke);
        renderer.setStroke(Color(0, 0, 0), 1.0);
        renderer.draw(parabolaSegmentToBezier(paraSeg));
        renderer.draw(paraSeg.center());
        renderer.draw(paraSeg.line());
        renderer.draw(seg);
        renderer.setMode(GeometryRenderer::stroke | GeometryRenderer::fill);
        renderer.setFill(Color(200, 200, 200));
        renderer.setMode(GeometryRenderer::stroke);
        auto piece = intersection(paraSeg, *m_circle);
        auto segPiece = intersection(seg, *m_circle);
        renderer.draw(*m_circle);
        if (piece.has_value()) {
            renderer.setStroke(Color(255, 0, 0), 2.0);
            renderer.draw(parabolaSegmentToBezier(*piece));
        }
        if (segPiece.has_value()) {
            renderer.setStroke(Color(255, 0, 0), 2.0);
            renderer.draw(*segPiece);
        }

        std::vector<Point<Inexact>> inters;
        intersections(paraSeg, *m_circle, std::back_inserter(inters));
        intersections(seg, *m_circle, std::back_inserter(inters));
        for (const auto& inter : inters) {
            renderer.draw(inter);
        }
    }, "Parabola");
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    IntersectionHelpersDemo demo;
    demo.show();
    app.exec();
}