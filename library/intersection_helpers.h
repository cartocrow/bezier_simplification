#pragma once

#include <CGAL/Parabola_segment_2.h>
#include <cartocrow/core/core.h>
#include <cartocrow/circle_segment_helpers/cs_types.h>

#include "conic_types.h"

namespace cartocrow {
template <class Gt>
Conic_arc
parabolaSegmentToConicArc(const CGAL::Parabola_segment_2<Gt>& ps) {
    Conic_Traits traits;
    auto ctr_cv = traits.construct_curve_2_object();

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

    auto ls = directrix.perpendicular(ps.p1);
    auto lt = directrix.perpendicular(ps.p2);
//    CGAL::Orientation orient = directrix.oriented_side(focus);

    return ctr_cv(r, s, t, u, v, w, CGAL::orientation(focus, ps.p1, ps.p2),
                  Conic_Point(ps.p1.x(), ps.p1.y()),
                  0, 0, 0, ls.a(), ls.b(), ls.c(),
                  Conic_Point(ps.p2.x(), ps.p2.y()),
                  0, 0, 0, lt.a(), lt.b(), lt.c()
           );
}

template <class Arr, class OutputIterator>
void intersectionsArrCurves(const typename Arr::Geometry_traits_2::Curve_2& curve1, const typename Arr::Geometry_traits_2::Curve_2& curve2, OutputIterator out) {
    Arr arrangement;
    CGAL::insert(arrangement, curve1);
    CGAL::insert(arrangement, curve2);

    for (auto vh : arrangement.vertex_handles()) {
        if (vh->degree() > 2) {
            *out++ = vh->point();
        }
    }
}

/// Streams Rat_point to out.
/// Returns the points at which the parabola segment and circle intersect.
template <class Gt, class OutputIterator>
void intersections(const CGAL::Parabola_segment_2<Gt>& ps, const Rat_circle& circle, OutputIterator out) {
    auto conicArc = parabolaSegmentToConicArc(ps);
    intersectionsArrCurves<Conic_Arrangement>(conicArc, circle, out);
}

/// Streams Point<Inexact> to out.
/// Returns the (approximate) points at which the parabola segment and circle intersect.
template <class Gt, class OutputIterator>
void intersections(const CGAL::Parabola_segment_2<Gt>& ps, const Circle<Inexact>& circleInexact, OutputIterator out) {
    auto circle = Rat_circle({circleInexact.center().x(), circleInexact.center().y()}, circleInexact.squared_radius());
    std::vector<Conic_Point> inters;
    intersections(ps, circle, std::back_inserter(inters));
    for (const auto& pt : inters) {
        *out++ = approximate(pt);
    }
}

template <class OutputIterator>
void intersections(const Segment<Inexact>& s, const Circle<Inexact>& circle, OutputIterator out) {
    std::vector<OneRootPoint> inters;
    CSCurve seg(pretendExact(s));
    intersectionsArrCurves<CSArrangement>(seg, pretendExact(circle), std::back_inserter(inters));
    for (const auto &pt: inters) {
        *out++ = approximateOneRootPoint(pt);
    }
}

template <class Gt, class OutputIterator>
void approximateIntersections(const CGAL::Parabola_segment_2<Gt>& ps, const Circle<Inexact>& circleInexact, OutputIterator out, double step) {
    std::vector<Point<Inexact>> points;
    ps.generate_points(points, step);
    for (int i = 0; i < points.size() - 1; ++i) {
        Segment<Inexact> seg(points[i], points[i+1]);
        intersections(seg, circleInexact, out);
    }
}

std::optional<Segment<Inexact>> intersection(const Segment<Inexact>& s, const Circle<Inexact>& circle);

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
        bool needToSwap = (proj2I - proj1I) * (proj2E - proj1E) < 0;
        auto newP1 = approximate(inters[0]);
        auto newP2 = approximate(inters[1]);
        if (needToSwap) {
            std::swap(newP1, newP2);
        }

        return CGAL::Parabola_segment_2<Gt>(ps.center(), ps.line(), newP1, newP2);
    }

    std::cerr << "Unexpected number of intersections!" << std::endl;
    return ps;
//    throw std::runtime_error("Unexpected number of intersections!");
}

/// Gt must use Point<Inexact>
/// Returns the part of the parabola segment that lies inside the circle.
template <class Gt>
std::optional<CGAL::Parabola_segment_2<Gt>> approximateIntersection(const CGAL::Parabola_segment_2<Gt>& ps, const Circle<Inexact>& circle, double step) {
    std::vector<Point<Inexact>> inters;
    approximateIntersections(ps, circle, std::back_inserter(inters), step);

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
        bool needToSwap = (proj2I - proj1I) * (proj2E - proj1E) < 0;
        auto newP1 = approximate(inters[0]);
        auto newP2 = approximate(inters[1]);
        if (needToSwap) {
            std::swap(newP1, newP2);
        }

        return CGAL::Parabola_segment_2<Gt>(ps.center(), ps.line(), newP1, newP2);
    }

    std::cerr << "Unexpected number of intersections!" << std::endl;
    return ps;
}
}