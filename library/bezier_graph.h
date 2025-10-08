#pragma once

#include "curved_graph.h"
#include "bezier_simplification_traits.h"

namespace cartocrow::curved_simplification {
template <class VD, class ED>
using BezierGraph = CurvedGraph<VD, ED, BezierSimplificationTraits>;
}