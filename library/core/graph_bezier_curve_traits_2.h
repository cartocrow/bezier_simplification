#pragma once

#include <cartocrow/core/cubic_bezier.h>
#include "graph_curve_traits_2.h"

namespace cartocrow {
struct Graph_Bezier_curve_traits_2 {
	using Curve_2 = CubicBezierCurve;
	using Point_2 = Point<Inexact>;
    using Kernel = Inexact;

    static Curve_2 reversed(const Curve_2& curve) {
        return curve.reversed();
    }

    static Curve_2 transform(const Curve_2& curve, const CGAL::Aff_transformation_2<Inexact>& trans) {
        return curve.transform(trans);
    }

    static_assert(GraphCurveTraits_2<Graph_Bezier_curve_traits_2>);
};
}
