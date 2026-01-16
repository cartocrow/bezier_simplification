#include "min_dist_demo.h"

#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QPushButton>

#include "library/read_graph_gdal.h"

#include <cartocrow/renderer/voronoi_drawer.h>
#include <cartocrow/core/segment_delaunay_graph_helpers.h>
#include <cartocrow/core/delaunay_voronoi_helpers.h>
#include <cartocrow/circle_segment_helpers/cs_types.h>

#include "library/conic_types.h"
#include "library/intersection_helpers.h"

VoronoiEdge
objectToVoronoiEdge(CGAL::Object o) {
    Segment<Inexact> s;
    Line<Inexact> l;
    Ray<Inexact> r;
    CGAL::Parabola_segment_2<Gt> ps;

    if (CGAL::assign(s, o)) {
        return s;
    }
    if (CGAL::assign(l, o)) {
        return l;
    }
    if (CGAL::assign(r, o)) {
        return r;
    }
    if (CGAL::assign(ps, o)) {
        return ps;
    }

    throw std::runtime_error("Object is not a segment, line, ray, or parabola segment");
}

/// Returns the smallest distance to a site attained by a point on this edge (multiplied by 2 and squared
double minDist(const SDG& delaunay, const SDG::Edge& edge) {
    VoronoiEdge vEdge = objectToVoronoiEdge(delaunay.primal(edge));
    if (auto sP = std::get_if<Segment<Inexact>>(&vEdge)) {
        auto& s = *sP;
        auto [p, q] = defining_sites<SDG>(edge);

        if (p.is_point() && q.is_point()) {
            // Check if midpoint of p and q is part of the segment.
            auto m = CGAL::midpoint(p.point(), q.point());
            auto v =  s.target() - s.source();
            auto dot = (m - s.source()) * v;
            if (dot >= 0 || dot <= v * v) {
                // midpoint lies on segment
                return CGAL::squared_distance(p.point(), q.point());
            }
            return CGAL::sqrt(std::min(CGAL::squared_distance(s.source(), p.point()), CGAL::squared_distance(s.target(), p.point())));
        }

        auto sProj = projection(p.segment(), s.source());
        auto tProj = projection(p.segment(), s.target());
        return CGAL::sqrt(std::min(CGAL::squared_distance(s.source(), sProj), CGAL::squared_distance(s.target(), tProj)));
    }
    if (auto psP = std::get_if<CGAL::Parabola_segment_2<Gt>>(&vEdge)) {
        auto& ps = *psP;
        const auto& directrix = ps.line();
        const auto& focus = ps.center();
        auto focusProj = directrix.projection(focus);
        auto dirV = directrix.to_vector();
        auto psStartProj = directrix.projection(ps.p1);
        auto psEndProj   = directrix.projection(ps.p2);
        auto psStartDot  = (psStartProj - directrix.point()) * dirV;
        auto psEndDot    = (psEndProj - directrix.point()) * dirV;
        auto psFocusDot  = (focusProj - directrix.point()) * dirV;
        if ((psStartDot > psFocusDot || psEndDot > psFocusDot) && (psStartDot < psFocusDot || psEndDot < psFocusDot)) {
            // focus is part of parabolic segment.
            return CGAL::sqrt(CGAL::squared_distance(focus, focusProj)) / 2;
        } else {
            // focus is not part of parabolic segment.
            return CGAL::sqrt(std::min(CGAL::squared_distance(ps.p1, psStartProj), CGAL::squared_distance(ps.p2, psEndProj)));
        }
    }
    if (auto lP = std::get_if<Line<Inexact>>(&vEdge)) {
        auto [p, q] = defining_sites<SDG>(edge);
        return CGAL::sqrt(CGAL::squared_distance(p.point(), q.point())) / 2;
    }
    if (auto rP = std::get_if<Ray<Inexact>>(&vEdge)) {
        auto& r = *rP;
        // The point with min. dist. is the midpoint of p and q if it lies on the ray, otherwise the source of ray
        auto [p, q] = defining_sites<SDG>(edge);
        auto m = CGAL::midpoint(p.point(), q.point());
        if ((m - r.source()) * r.to_vector() >= 0) {
            return CGAL::sqrt(CGAL::squared_distance(m, p.point()));
        } else {
            return CGAL::sqrt(CGAL::squared_distance(r.source(), p.point()));
        }
    }

    throw std::runtime_error("Unexpected Segment Delaunay Graph edge type!");
}

std::optional<VoronoiEdge>
withinDist(const SDG& delaunay, const SDG::Edge& edge, double dist) {
    auto vEdge = objectToVoronoiEdge(delaunay.primal(edge));

    if (minDist(delaunay, edge) > dist) {
        return std::nullopt;
    }

    Conic_Arrangement arr;
    auto [p, q] = defining_sites<SDG>(edge);

    if (auto sP = std::get_if<Segment<Inexact>>(&vEdge)) {
        auto& s = *sP;
        if (p.is_point() && q.is_point()) {
            CSArrangement arr;
            CSCurve circ(pretendExact(p.point()), dist * dist, CGAL::COUNTERCLOCKWISE);
            CGAL::insert(arr, circ);
            CGAL::insert(arr, pretendExact(s));
            std::vector<CSArrangement::Vertex_handle> inters;
            for (auto vh : arr.vertex_handles()) {
                if (vh->degree() > 2) {
                    inters.push_back(vh);
                }
            }

            auto point = p.point();
            auto sKept = CGAL::squared_distance(s.source(), point) < dist * dist;
            auto tKept = CGAL::squared_distance(s.target(), point) < dist * dist;

            if (inters.empty()) {
                assert(sKept && tKept);
                return s;
            }

            if (inters.size() == 1) {
                if (sKept) {
                    assert(!tKept);
                    return Segment<Inexact>(s.source(), approximateOneRootPoint(inters[0]->point()));
                }
                if (tKept) {
                    assert(!sKept);
                    return Segment<Inexact>(approximateOneRootPoint(inters[0]->point()), s.target());
                }
            }

            if (inters.size() == 2) {
                auto i0 = approximateOneRootPoint(inters[0]->point());
                auto i1 = approximateOneRootPoint(inters[1]->point());

                auto newSource = CGAL::squared_distance(i0, s.source()) < CGAL::squared_distance(i0, s.target()) ? i0 : i1;
                auto newTarget = newSource == i0 ? i1 : i0;

                return Segment<Inexact>(newSource, newTarget);
            }

            throw std::runtime_error("Unexpected number of intersections!");
        }

        // p and q are segments.
        // distance changes linearly
        auto sProj = projection(p.segment(), s.source());
        auto tProj = projection(p.segment(), s.target());
        auto sDist2 = CGAL::squared_distance(s.source(), sProj);
        auto tDist2 = CGAL::squared_distance(s.target(), tProj);
        auto largerDist = sDist2 > tDist2 ? sDist2 : tDist2;
        auto smallerDist = sDist2 > tDist2 ? tDist2 : sDist2;

        if (sDist2 < dist * dist && tDist2 < dist * dist) {
            return s;
        }

        auto sourceClosest = sDist2 < tDist2;

        auto narrow = sourceClosest ? s.source() : s.target();
        auto wide = sourceClosest ? s.target() : s.source();
        auto vec = wide - narrow;
        vec *= (dist - CGAL::sqrt(smallerDist)) / (CGAL::sqrt(largerDist) - CGAL::sqrt(smallerDist));
        return Segment<Inexact>(narrow, narrow + vec);
    }
    if (auto psP = std::get_if<CGAL::Parabola_segment_2<Gt>>(&vEdge)) {
        auto& ps = *psP;

        auto pointSite = p.is_point() ? p : q;
        auto point = pointSite.point();

        Circle<Inexact> circle({point.x(), point.y()}, dist * dist);

        return *intersection(ps, circle);
    }
    if (auto lP = std::get_if<Line<Inexact>>(&vEdge)) {
        auto& l = *lP;
        // todo
        return l;
    }
    if (auto rP = std::get_if<Ray<Inexact>>(&vEdge)) {
        auto& r = *rP;
        // todo
        return r;
    }
}

void MinDistDemo::recompute() {
    m_withinDistEdges.clear();
    for (auto eit = m_delaunay.finite_edges_begin(); eit != m_delaunay.finite_edges_end(); ++eit) {
        auto [p, q] = defining_sites<SDG>(*eit);
        // Skip edges that are defined by consecutive isoline edges
        if (p.is_segment() && q.is_segment() && (p.source() == q.source() || p.target() == q.target()  || p.source() == q.target() || p.target() == q.source()) ||
            p.is_point()   && q.is_segment() && (p.point()  == q.source() || p.point()  == q.target()) ||
            p.is_segment() && q.is_point()   && (p.source() == q.point()  || p.target() == q.point())) continue;

        auto vorEdge = withinDist(m_delaunay, *eit, m_minDist->value());
        if (vorEdge.has_value()) {
            m_withinDistEdges.push_back(*vorEdge);
        }
    }
}

MinDistDemo::MinDistDemo() {
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

    auto* resetButton = new QPushButton("Reset");
    vLayout->addWidget(resetButton);

    m_minDist = new DoubleSlider(Qt::Horizontal);
    m_minDist->setMaximum(100000);
    m_minDist->setMinimum(0);
    m_minDist->setValue(60000);
    vLayout->addWidget(m_minDist);

    connect(stepButton, &QPushButton::clicked, [this]() {
        forceDirectedMinimumDistance(m_g, m_minDist->value());
        m_renderer->repaint();
    });

    connect(resetButton, &QPushButton::clicked, [this]() {
        m_g = m_ogg;
        m_renderer->repaint();
    });

    connect(m_minDist, &DoubleSlider::valueChanged, [this]() {
        m_renderer->repaint();
        recompute();
    });

    m_g = readGraphUsingGDAL("data/small_min_dist_test.shp");
    m_g.orient();
    m_ogg = m_g;

    std::vector<Segment<Inexact>> segs;
    for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
        segs.push_back(Segment<Inexact>(eit->source()->point(), eit->target()->point()));
    }

    m_delaunay.insert_segments(segs.begin(), segs.end());

    m_renderer->setMaxZoom(1000000000);
    m_renderer->setMinZoom(0.00000001);

    recompute();

    std::vector<Point<Inexact>> points;
    for (auto vit = m_g.vertices_begin(); vit != m_g.vertices_end(); ++vit) {
        points.push_back(vit->point());
    }
    m_renderer->fitInView(CGAL::bbox_2(points.begin(), points.end()));

    m_renderer->addPainting([this](GeometryRenderer& renderer) {
        renderer.setMode(GeometryRenderer::fill);
        renderer.setFill(Color(240, 240, 240));
//        renderer.setFillOpacity(10);
//        for (auto vit = m_g.vertices_begin(); vit != m_g.vertices_end(); ++vit) {
//            renderer.draw(Circle<Inexact>(vit->point(), m_minDist->value() * m_minDist->value()));
//        }

        renderer.draw(Circle<Inexact>(m_renderer->mousePosition(), m_minDist->value() * m_minDist->value()));
    }, "Disk");

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
        auto voronoiDrawer = VoronoiDrawer<Gt>(&renderer);
        for (auto eit = m_delaunay.finite_edges_begin(); eit != m_delaunay.finite_edges_end(); ++eit) {
            renderer.setStroke(Color{210, 210, 210}, 1.0);
//            if (minDist(m_delaunay, *eit) <= m_minDist->value()) {
//                renderer.setStroke(Color{255, 50, 50}, 2.0);
//            }
            auto [p, q] = defining_sites<SDG>(*eit);
            // Skip edges that are defined by consecutive isoline edges
            if (p.is_segment() && q.is_segment() && (p.source() == q.source() || p.target() == q.target()  || p.source() == q.target() || p.target() == q.source()) ||
                p.is_point()   && q.is_segment() && (p.point()  == q.source() || p.point()  == q.target()) ||
                p.is_segment() && q.is_point()   && (p.source() == q.point()  || p.target() == q.point())) continue;

            draw_dual_edge(m_delaunay, *eit, voronoiDrawer);



//            auto projP = site_projection<SDG>(m_delaunay, *eit, p);
//            auto projQ = site_projection<SDG>(m_delaunay, *eit, q);

//            renderer.setStroke(Color{255, 50, 50}, 2.0);
//            std::visit([&renderer](const auto& segOrPoint) {
//                renderer.draw(segOrPoint);
//            }, projP);
//            std::visit([&renderer](const auto& segOrPoint) {
//                renderer.draw(segOrPoint);
//            }, projQ);
        }

        for (const auto& vorEdge : m_withinDistEdges) {
            std::visit([&](const auto &geom) {
                renderer.setStroke(Color{255, 50, 50}, 2.0);
                voronoiDrawer << geom;
            }, vorEdge);
        }

    }, "Segment Voronoi diagram");
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    MinDistDemo demo;
    demo.show();
    app.exec();
}