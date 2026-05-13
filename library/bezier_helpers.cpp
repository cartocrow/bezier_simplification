#include "bezier_helpers.h"

#include <cartocrow/core/vector_helpers.h>

namespace cartocrow {
bool isStraight(const CubicBezierCurve& curve) {
    auto a1 = abs(CGAL::area(curve.control(0), curve.control(1), curve.control(2)));
    auto a2 = abs(CGAL::area(curve.control(1), curve.control(2), curve.control(3)));
    return a1 + a2 < M_EPSILON;
}

bool connectsSmoothlyTo(const CubicBezierCurve& c1, const CubicBezierCurve& c2) {
    auto control1 = CGAL::squared_distance(c1.targetControl(), c1.target()) > M_EPSILON ? c1.targetControl() :
        (CGAL::squared_distance(c1.sourceControl(), c1.target()) > M_EPSILON ? c1.sourceControl() : c1.source());
    auto control2 = CGAL::squared_distance(c2.source(), c2.sourceControl()) > M_EPSILON ? c2.sourceControl() :
        (CGAL::squared_distance(c2.source(), c2.targetControl()) > M_EPSILON ? c2.targetControl() : c2.target());
    return smallestAngleBetween(c1.target() - control1, control2 - c1.target()) < M_EPSILON;
}

bool isSmooth(const CubicBezierSpline& spline) {
    for (int i = 0; i < spline.numCurves() - 1; ++i) {
        if (!connectsSmoothlyTo(spline.curve(i), spline.curve(i + 1))) return false;
    }
    return true;
}
}