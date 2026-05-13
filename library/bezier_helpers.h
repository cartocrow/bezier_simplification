#pragma once

#include "bezier_helpers.h"
#include <cartocrow/core/cubic_bezier.h>

namespace cartocrow {
bool isStraight(const CubicBezierCurve& curve);
// Whether curve c1 connects smoothly to curve c2: whether the end tangent of c1 aligns approximately with the
// start tangents of c2
bool connectsSmoothlyTo(const CubicBezierCurve& c1, const CubicBezierCurve& c2);
bool isSmooth(const CubicBezierSpline& spline);
}