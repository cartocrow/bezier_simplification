#pragma once

#include <vector>
#include <cartocrow/core/cubic_bezier.h>

namespace cartocrow::curved_simplification {
CubicBezierSpline fitSpline(const std::vector<Point<Inexact>>& points, double maxSquaredError);
CubicBezierCurve fitCurve(const std::vector<Point<Inexact>>& points);
CubicBezierCurve fitCurve(const std::vector<Point<Inexact>>& points, const Vector<Inexact>& startTangent, const Vector<Inexact>& endTangent);
CubicBezierSpline fitTwoCurves(const std::vector<Point<Inexact>>& points);
}
