#pragma once

#include <cartocrow/core/cubic_bezier.h>

namespace cartocrow::curved_simplification {
struct BezierSimplificationTraits {
	using Curve = CubicBezierCurve;
	using Pt = Point<Inexact>;
};
}

