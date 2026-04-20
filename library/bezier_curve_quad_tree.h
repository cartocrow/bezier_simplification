#pragma once

#include <cartocrow/data_structures/quad_tree.h>
#include "utils.h"

namespace cartocrow::curved_simplification {
template<class Graph>
struct BezierCurveQuadTreeTraits {
	using Element = Graph::Edge_handle;
	using Kernel = Graph::Kernel;

	static Rectangle<Kernel> get_bounding_box(Element& elt) {
		CubicBezierCurve curve = elt->curve();
		CGAL::Bbox_2 box = curve.bbox();
		return { box.xmin(), box.ymin(), box.xmax(), box.ymax() };
	}

	static bool element_overlaps_rectangle(Element& elt, Rectangle<Kernel>& rect) {
		CubicBezierCurve curve = elt->curve();
		return utils::overlaps(rect, curve);
	}
};

template<class Graph>
using BezierCurveQuadTree = cartocrow::data_structures::QuadTree<BezierCurveQuadTreeTraits<Graph>>;
}