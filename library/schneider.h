#pragma once

#include <vector>
#include <cartocrow/core/cubic_bezier.h>

namespace cartocrow::curved_simplification {
CubicBezierSpline fitSpline(const std::vector<Point<Inexact>>& points, double maxSquaredError, int reparameterizationIterations=10);
CubicBezierSpline fitSpline(const std::vector<Point<Inexact>>& points, double maxSquaredError, const Vector<Inexact>& startTangent, const Vector<Inexact>& endTangent, int reparameterizationIterations=10);
CubicBezierCurve fitCurve(const std::vector<Point<Inexact>>& points, int reparameterizationIterations=10);
CubicBezierCurve fitCurve(const std::vector<Point<Inexact>>& points, const Vector<Inexact>& startTangent, const Vector<Inexact>& endTangent, int reparameterizationIterations=10);
CubicBezierSpline fitTwoCurves(const std::vector<Point<Inexact>>& points, int reparameterizationIterations=10);
}
