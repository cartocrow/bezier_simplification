#include "topo_set.h"

namespace cartocrow {
TopoSet<Exact> pretendExact(const TopoSet<Inexact>& topoSet) {
    TopoSet<Exact> exact;
    for (const auto& arc : topoSet.arcs) {
        exact.arcs.push_back(pretendExact(arc));
    }
    for (const auto& feature : topoSet.features) {
        TopoSet<Exact>::Geometry geometry;
        if (auto gP = std::get_if<TopoSet<Inexact>::PolygonSetGeometry>(&feature.geometry)) {
            geometry = TopoSet<Exact>::PolygonSetGeometry(gP->arcs);
        } else {
            std::cerr << "Did not handle type of geometry!" << std::endl;
        }
        exact.features.push_back(TopoSet<Exact>::Feature{geometry, feature.attributes});
    }
    return exact;
}

TopoSet<Inexact> approximate(const TopoSet<Exact>& topoSet) {
    TopoSet<Inexact> inexact;
    for (const auto& arc : topoSet.arcs) {
        inexact.arcs.push_back(approximate(arc));
    }
    for (const auto& feature : topoSet.features) {
        TopoSet<Inexact>::Geometry geometry;
        if (auto gP = std::get_if<TopoSet<Exact>::PolygonSetGeometry>(&feature.geometry)) {
            geometry = TopoSet<Inexact>::PolygonSetGeometry(gP->arcs);
        } else {
            std::cerr << "Did not handle type of geometry!" << std::endl;
        }
        inexact.features.push_back(TopoSet<Inexact>::Feature{geometry, feature.attributes});
    }
    return inexact;
}
}