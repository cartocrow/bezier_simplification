#pragma once

#include <CGAL/Parabola_segment_2.h>
#include <cartocrow/core/core.h>

#include "conic_types.h"

/// Streams Point<Inexact> to out.
/// Returns the (approximate) points at which the parabola segment and circle intersect.
template <class Gt, class OutputIterator>
void intersections(const CGAL::Parabola_segment_2<Gt>& ps, const Circle<Inexact>& circleInexact, OutputIterator out) {
    std::cout << "!" << std::endl;
    Conic_Traits traits;
    auto ctr_cv = traits.construct_curve_2_object();

    auto circle = Rat_circle({circleInexact.center().x(), circleInexact.center().y()}, circleInexact.squared_radius());

    auto& directrix = ps.line();
    auto& focus = ps.center();
    auto div = CGAL::sqrt(directrix.a() * directrix.a() + directrix.b() * directrix.b());
    auto a = directrix.a() / div;
    auto b = directrix.b() / div;
    auto c = directrix.c() / div;
    auto r = 1 - a * a;
    auto s = 1 - b * b;
    auto t = -2 * a * b;
    auto u = -2 * focus.x() - 2 * a * c;
    auto v = -2 * focus.y() - 2 * b * c;
    auto w = focus.x() * focus.x() + focus.y() * focus.y() - c * c;

    auto ls = Line<Inexact>(ps.p1, focus);
    auto lt = Line<Inexact>(ps.p2, focus);

    auto parabolicSegment =
            ctr_cv(r, s, t, u, v, w, directrix.oriented_side(focus),
                   Conic_Point(ps.p1.x(), ps.p1.y()),
                   0, 0, 0, ls.a(), ls.b(), ls.c(),
                   Conic_Point(ps.p2.x(), ps.p2.y()),
                   0, 0, 0, lt.a(), lt.b(), lt.c()
            );

    Conic_Arrangement arrangement;
    CGAL::insert(arrangement, parabolicSegment);
    CGAL::insert(arrangement, circle);

    for (auto vh : arrangement.vertex_handles()) {
        if (vh->degree() > 2) {
            *out++ = approximate(vh->point());
        }
    }
    std::cout << ":)" << std::endl;
}

/// Gt must use Point<Inexact>
/// Returns the part of the parabola segment that lies inside the circle.
template <class Gt>
std::optional<CGAL::Parabola_segment_2<Gt>> intersection(const CGAL::Parabola_segment_2<Gt>& ps, const Circle<Inexact>& circle) {
    std::vector<Point<Inexact>> inters;
    intersections(ps, circle, std::back_inserter(inters));

    auto p1Kept = circle.has_on_bounded_side(ps.p1);
    auto p2Kept = circle.has_on_bounded_side(ps.p2);

    if (inters.empty()) {
        if (!p1Kept) return std::nullopt;
        assert(p1Kept && p2Kept);
        return ps;
    }

    if (inters.size() == 1) {
        if (p1Kept) {
            assert(!p2Kept);
            return CGAL::Parabola_segment_2<Gt>(ps.center(), ps.line(), ps.p1, inters[0]);
        }
        if (p2Kept) {
            assert(!p1Kept);
            return CGAL::Parabola_segment_2<Gt>(ps.center(), ps.line(), inters[0], ps.p2);
        }
    }

    if (inters.size() == 2) {
        auto& directrix = ps.line();
        auto proj1I = directrix.projection(inters[0]);
        auto proj2I = directrix.projection(inters[1]);
        auto proj1E = directrix.projection(approximate(ps.p1));
        auto proj2E = directrix.projection(approximate(ps.p2));
        auto v = directrix.to_vector();
        bool needToSwap = (proj2I - proj1I) * (proj2E - proj1E) < 0;
        auto newP1 = approximate(inters[0]);
        auto newP2 = approximate(inters[1]);
        if (needToSwap) {
            std::swap(newP1, newP2);
        }

        return CGAL::Parabola_segment_2<Gt>(ps.center(), ps.line(), newP1, newP2);
    }

    throw std::runtime_error("Unexpected number of intersections!");
}