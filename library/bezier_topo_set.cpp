#include "bezier_topo_set.h"

namespace cartocrow {
BezierTopoSet bezierTopoSetFromStraightTopoSet(const StraightTopoSet<Inexact>& sTopoSet) {
    BezierTopoSet bTopoSet;
    bTopoSet.features = sTopoSet.features;
    for (const auto& arc : sTopoSet.arcs) {
        CubicBezierSpline spline;
        for (int i = 0; i < arc.size() - 1; ++i) {
            spline.appendCurve(arc.vertex(i), arc.vertex(i + 1));
        }
        bTopoSet.arcs.push_back(spline);
    }
    return bTopoSet;
}


StraightTopoSet<Inexact> approximate(const BezierTopoSet& bTopoSet) {
    StraightTopoSet<Inexact> sTopoSet;
    sTopoSet.features = bTopoSet.features;

    for (auto bArc : bTopoSet.arcs) {
        int nPoints = std::max(std::min(1000.0, 31998.0 / bArc.numCurves()), 5.0);

        Polyline<Inexact> sArc;

        if (!bArc.empty()) {
            bArc.curve(0).samplePoints(isStraight(bArc.curve(0)) ? 2 : nPoints, std::back_inserter(sArc));
            for (int i = 1; i < bArc.numCurves(); ++i) {
                sArc.pop_back();
                bArc.curve(i).samplePoints(isStraight(bArc.curve(i)) ? 2 : nPoints, std::back_inserter(sArc));
            }
        }
        sTopoSet.arcs.push_back(sArc);
    }

    return sTopoSet;
}
}