#include "intersection_helpers.h"

namespace cartocrow {
std::optional<Segment<Inexact>> intersection(const Segment<Inexact>& s, const Circle<Inexact>& circle) {
    std::vector<Point<Inexact>> inters;
    intersections(s, circle, std::back_inserter(inters));

    auto sKept = circle.has_on_bounded_side(s.source());
    auto tKept = circle.has_on_bounded_side(s.target());

    if (inters.empty()) {
        if (sKept && tKept) {
            return s;
        } else {
            assert(!sKept && !tKept);
            return std::nullopt;
        }
    }

    if (inters.size() == 1) {
        if (sKept) {
            assert(!tKept);
            return Segment<Inexact>(s.source(), inters[0]);
        }
        if (tKept) {
            assert(!sKept);
            return Segment<Inexact>(inters[0], s.target());
        }
    }

    if (inters.size() == 2) {
        auto i0 = inters[0];
        auto i1 = inters[1];

        auto newSource = CGAL::squared_distance(i0, s.source()) < CGAL::squared_distance(i0, s.target()) ? i0 : i1;
        auto newTarget = newSource == i0 ? i1 : i0;

        return Segment<Inexact>(newSource, newTarget);
    }

    throw std::runtime_error("Unexpected number of intersections!");
}
}