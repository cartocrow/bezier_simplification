#pragma once

#include <cartocrow/core/cubic_bezier.h>
#include "core/topo_set.h"
#include "bezier_helpers.h"

namespace cartocrow {
using BezierTopoSet = TopoSet<CubicBezierSpline>;

BezierTopoSet bezierTopoSetFromStraightTopoSet(const StraightTopoSet<Inexact>& sTopoSet);
StraightTopoSet<Inexact> approximate(const BezierTopoSet& bTopoSet);
}