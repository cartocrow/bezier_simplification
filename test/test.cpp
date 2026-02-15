#include <catch2/catch_test_macros.hpp>
#include "library/steven_bezier_collapse.h"

using namespace cartocrow;
using namespace cartocrow::curved_simplification;

TEST_CASE("getAreaD") {
    CubicBezierCurve curve({0, 0}, {3, 1}, {2, 6}, {5, 6});
    auto a = curve.signedArea();
    auto p0 = curve.source();
    auto p3 = curve.target();
    auto t0 = curve.sourceControl() - curve.source();
    auto t1 = curve.targetControl() - curve.target();
    auto d0 = sqrt(t0.squared_length());
    auto d1 = sqrt(t1.squared_length());

    CHECK(abs(getAreaD0(a, p0, t0.direction(), t1.direction(), d1, p3) - d0) < M_EPSILON);
    CHECK(abs(getAreaD1(a, p0, t0.direction(), d0, t1.direction(), p3) - d1) < M_EPSILON);
}