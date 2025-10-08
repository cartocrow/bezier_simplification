#pragma once

namespace cartocrow::curve_simplification {
template <class CST>
concept CurveSimplificationTraits =
	requires(CST cst, typename CST::Curve curve) {
        typename CST::Curve;
	    typename CST::Pt;

	    {
		    curve.source()
	    } -> std::same_as<CST::Pt>;

		{
			curve.target()
		} -> std::same_as<CST::Pt>;

	    { // todo: do not require the curve class to have this function but follow cgal's approach
		    curve.reverse()
	    }
    };
}