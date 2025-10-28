#pragma once

#include <concepts>

namespace cartocrow {
template <class CST>
concept GraphCurveTraits_2 =
	requires(CST cst, typename CST::Curve_2 curve) {
        typename CST::Curve_2;
	    typename CST::Point_2;
        typename CST::Kernel;

	    {
		    curve.source()
	    } -> std::convertible_to<typename CST::Point_2>;

		{
			curve.target()
		} -> std::convertible_to<typename CST::Point_2>;

	    {
		    CST::reversed(curve)
	    } -> std::convertible_to<typename CST::Curve_2>;
    };
}