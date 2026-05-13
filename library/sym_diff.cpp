#include "sym_diff.h"

#include <cartocrow/core/arrangement_helpers.h>

namespace cartocrow::curved_simplification {
double symmetricDifference(const CubicBezierSpline& spline1, const CubicBezierSpline& spline2, double intersectionThreshold) {
    std::vector<CubicBezierIntersectionResult<CubicBezierSpline, CubicBezierSpline>> inters;
    spline1.intersections(spline2, std::back_inserter(inters), intersectionThreshold);

    double sym_diff = 0;
    for (int i = 0; i < inters.size() - 1; ++i) {
        auto start = inters[i];
        auto end = inters[i + 1];

        auto part1 = spline1.sub(start.parameter1, end.parameter1);
        auto part2 = spline2.sub(start.parameter2, end.parameter2);
        part1.appendSpline(part2.reversed());

        sym_diff += std::abs(part1.signedArea());
    }

    return sym_diff;
}

double symmetricDifferenceArrangement(const CubicBezierSpline& spline1, const CubicBezierSpline& spline2, int symDiffSegs) {
    auto beforePlE = pretendExact(spline1.polyline(symDiffSegs));
    auto afterPlE = pretendExact(spline2.polyline(symDiffSegs));
    Arrangement<Exact> arr;
    std::vector<Arrangement<Exact>::X_monotone_curve_2> beforePlXMCurves;
    for (auto eit = beforePlE.edges_begin(); eit != beforePlE.edges_end(); ++eit) {
        beforePlXMCurves.emplace_back(*eit);
    }
    CGAL::insert_non_intersecting_curves(arr, beforePlXMCurves.begin(), beforePlXMCurves.end());
    CGAL::insert(arr, afterPlE.edges_begin(), afterPlE.edges_end());

    double symDiffErr = 0;
    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (fit->is_unbounded()) continue;
        auto pwh = approximate(face_to_polygon_with_holes<Exact>(fit));
        symDiffErr += abs(pwh.outer_boundary().area());
        for (const auto& h : pwh.holes()) {
            symDiffErr -= abs(h.area());
        }
    }
    return symDiffErr;
}
}