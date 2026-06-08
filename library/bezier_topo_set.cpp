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

std::vector<std::vector<CubicBezierSpline>> getGeometry(const TopoSet<CubicBezierSpline>::PolygonSetTopology& pgnSetTopology, const TopoSet<CubicBezierSpline>& topoSet) {
    std::vector<std::vector<CubicBezierSpline>> polygonSet;

    for (const auto& polygonWithHolesArcs : pgnSetTopology.arcs) {
        std::vector<CubicBezierSpline> polygonWithHoles;
        for (const auto& polygonArcs : polygonWithHolesArcs) {
            CubicBezierSpline polygon;
            Point<Inexact> last;
            for (int arcIx : polygonArcs) {
                const auto& arc = topoSet.arcs[arcIx];

                // We currently don't distinguish reverse arcs so just check which one we need
                bool reverse = false;
                if (polygon.empty()) {
                    if (polygonArcs.size() > 1) {
                        reverse = topoSet.arcs[polygonArcs[1]].source() != topoSet.arcs[polygonArcs[0]].target() && topoSet.arcs[polygonArcs[1]].target() != topoSet.arcs[polygonArcs[0]].target();

                    }
                }
                else {
                    reverse = arc.source() != last;
                }

                if (!reverse) {
                    if (!polygon.empty() && last != arc.source()) {
                        std::cerr << "Arcs do not connect!" << std::endl;
                    }
                    polygon.appendSpline(arc);
                    last = polygon.target();
                }
                else {
                    if (!polygon.empty() && last != arc.target()) {
                        std::cerr << "Arcs do not connect!" << std::endl;
                    }
                    polygon.appendSpline(arc.reversed());
                    last = polygon.target();
                }
            }
            polygonWithHoles.push_back(polygon);
        }
        polygonSet.push_back(polygonWithHoles);
    }

    return polygonSet;
}
}