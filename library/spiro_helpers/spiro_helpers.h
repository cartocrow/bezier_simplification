#pragma once

#include "bezctx_cbs.h"
#include <cartocrow/core/cubic_bezier.h>

namespace cartocrow {
template <class InputIterator>
CubicBezierSpline spiro(InputIterator begin, InputIterator end, bool closed) {
    std::vector<spiro_cp> points;
    for (auto pit = begin; pit != end; ++pit) {
        points.emplace_back(pit->x(), pit->y(), 'c');
    }
//    points.front().ty = '{';
//    points.back().ty = '}';

    auto bc = new_bezctx_cbs();
    auto ncq = 0;
    SpiroCPsToBezier2(points.data(), points.size(), ncq, closed, bc);
    return bezctx_cbs_close(bc);
}
}