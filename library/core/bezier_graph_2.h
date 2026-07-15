#pragma once

#include <cartocrow/data_structures/graph_2.h>
#include "graph_bezier_curve_traits_2.h"

namespace cartocrow {
template <class VD, class ED>
using Bezier_graph_2 = Graph_2<VD, ED, Graph_Bezier_curve_traits_2, DecomposedGraph<false, int>>;

template <class VD, class ED>
using Bezier_graph_with_history_2 = Graph_2<VD, ED, Graph_Bezier_curve_traits_2, DecomposedGraph<true, int>>;
}