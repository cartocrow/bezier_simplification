#pragma once

#include <cartocrow/core/cubic_bezier.h>

namespace cartocrow::curved_simplification {
double symmetricDifference(const CubicBezierSpline& spline1, const CubicBezierSpline& spline2, double intersectionThreshold=0.001);
double symmetricDifferenceArrangement(const CubicBezierSpline& spline1, const CubicBezierSpline& spline2, int symDiffSegs);
}