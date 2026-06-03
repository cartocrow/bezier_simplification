#pragma once

#include <cartocrow/core/cubic_bezier.h>
#include <cartocrow/data_structures/graph_curve_traits_2.h>

namespace cartocrow {
struct Graph_Bezier_curve_traits_2 {
	using Kernel = Inexact;
	using Point_2 = CGAL::Point_2<Kernel>;
	using Curve_2 = CubicBezierCurve;
	struct Curve_representation_2 {
		Point_2 source_control;
		Point_2 target_control;
	};

	static Curve_representation_2 representation(const Curve_2& curve) {
		return { curve.sourceControl(), curve.targetControl() };
	}
	static Curve_2 curve(const Point_2& start, const Point_2& end, const Curve_representation_2& rep) {
		return Curve_2(start, rep.source_control, rep.target_control, end);
	}
	static Curve_2 curve_reversed(const Point_2& start, const Point_2& end,
		const Curve_representation_2& rep) {
		return { end, rep.target_control, rep.source_control, start };
	}
	static CGAL::Bbox_2 bbox(const Point_2& start, const Point_2& end,
		const Curve_representation_2& rep) {
		return curve(start, end, rep).bbox();
	}
	static void reverse_representation(const Point_2& start, const Point_2& end,
		Curve_representation_2& rep) {
		std::swap(rep.source_control, rep.target_control);
	}
	static void move_start(const Point_2& start, const Point_2& end,
		Curve_representation_2& rep) {
		// skip
	}
	static void move_end(const Point_2& start, const Point_2& end, Curve_representation_2& rep) {
		// skip
	}
	static void transform(const CGAL::Aff_transformation_2<Kernel>& trans, const Point_2& start,
		const Point_2& end, Curve_representation_2& rep) {
		rep.source_control = trans.transform(rep.source_control);
		rep.target_control = trans.transform(rep.target_control);
	}
};
}
