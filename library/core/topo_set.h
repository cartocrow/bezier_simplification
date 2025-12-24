#pragma once

#include <cartocrow/core/core.h>
#include <cartocrow/core/polyline.h>
#include "region_set.h"

#include <CGAL/General_polygon_set_2.h>

namespace cartocrow {
template <class K>
struct TopoSet {
    std::vector<Polyline<K>> arcs;

    using Arc_handle = int; // we ignore order for now

    struct PolylineGeometry {
        std::vector<Arc_handle> arcs;
    };

    struct PolygonGeometry {
        std::vector<std::vector<Arc_handle>> arcs;
    };

    struct PolygonSetGeometry {
        std::vector<std::vector<std::vector<Arc_handle>>> arcs;

        PolygonSet<K> getGeometry(const TopoSet<K>& topoSet) {
            PolygonSet<K> polygonSet;

            for (const auto& polygonWithHolesArcs : arcs) {
                PolygonWithHoles<K> polygonWithHoles;
                for (const auto& polygonArcs : polygonWithHolesArcs) {
                    Polygon<K> polygon;
                    Point<K> last;
                    for (int arcIx : polygonArcs) {
                        const auto& arc = topoSet.arcs[arcIx];
                        // We currently don't distinguish reverse arcs so just check which one we need
                        if (polygon.is_empty() || *arc.vertices_begin() == last) {
                            std::copy(arc.vertices_begin(), --arc.vertices_end(), std::back_inserter(polygon));
                            last = *(--arc.vertices_end());
                        } else {
                            std::copy(arc.vertices_rbegin(), --arc.vertices_rend(), std::back_inserter(polygon));
                            last = *(--arc.vertices_rend());
                        }
                    }
                    if (polygonWithHoles.is_empty()) {
                        if (polygon.is_clockwise_oriented()) {
                            polygon.reverse_orientation();
                        }
                        polygonWithHoles.outer_boundary() = polygon;
                    } else {
                        if (polygon.is_counterclockwise_oriented()) {
                            polygon.reverse_orientation();
                        }
                        polygonWithHoles.add_hole(polygon);
                    }
                }
                polygonSet.insert(polygonWithHoles);
            }

            return polygonSet;
        }
    };

    using Geometry = std::variant<PolylineGeometry, PolygonGeometry, PolygonSetGeometry>;

    struct Feature {
        Geometry geometry;
        RegionAttributes attributes;
    };

    std::vector<Feature> features;

    TopoSet() = default;
    TopoSet(const RegionSet<K>& regionSet) {
        // Structure as in https://bost.ocks.org/mike/topology/.
        // 1. Extract.

        for (const auto& region : regionSet) {
            std::vector<PolygonWithHoles<K>> pgns;
            region.geometry.polygons_with_holes(std::back_inserter(pgns));


            features.emplace_back(PolygonSetGeometry{}, region.attributes);
            // Copy attributes

            // Extract arcs and store references to them
            for (const PolygonWithHoles<K>& pgn : pgns) {
                PolygonSetGeometry& geom = std::get<PolygonSetGeometry>(features.back().geometry);
                auto& pgnArcs = geom.arcs.emplace_back();
                const Polygon<K>& outer = pgn.outer_boundary();
                arcs.emplace_back(outer.vertices_begin(), outer.vertices_end());
                arcs.back().push_back(*outer.vertices_begin());
                auto& outerVec = pgnArcs.emplace_back();
                outerVec.push_back(arcs.size()-1);
                for (const Polygon<K>& hole : pgn.holes()) {
                    arcs.emplace_back(hole.vertices_begin(), hole.vertices_end());
                    arcs.back().push_back(*hole.vertices_begin());
                    auto& holeVec = pgnArcs.emplace_back();
                    holeVec.push_back(arcs.size()-1);
                }
            }
        }

        // 2. Join (identify junctions).
        std::unordered_map<Point<K>, std::vector<Point<K>>> neighbours;
        std::unordered_set<Point<K>> junctions;
        for (const auto& arc : arcs) {
            bool closed = *arc.vertices_begin() == *(--arc.vertices_end());
            for (int i = 0; i < arc.size(); ++i) {
                std::vector<Point<K>> currentNeighbours;
                if (i > 0) {
                    currentNeighbours.push_back(arc.vertex(i-1));
                } else if (closed) {
                    currentNeighbours.push_back(arc.vertex(arc.size() - 2));
                }
                if (i < arc.size() - 1) {
                    currentNeighbours.push_back(arc.vertex(i+1));
                } else if (closed) {
                    currentNeighbours.push_back(arc.vertex(1));
                }

                auto v = arc.vertex(i);
                if (!neighbours.contains(v)) {
                    for (const auto& n : currentNeighbours) {
                        neighbours[v].push_back(n);
                    }
                } else {
                    const auto& ns = neighbours[v];
                    for (const auto& n : currentNeighbours) {
                        if (std::find(ns.begin(), ns.end(), n) == ns.end()) {
                            // n is a new neighbour => v is a junction
                            neighbours[v].push_back(n);
                            junctions.insert(v);
                            continue;
                        }
                    }
                }
            }
        }

        // 3. Cut arcs.
        std::vector<Polyline<K>> newArcs;
        std::vector<std::list<std::vector<Arc_handle>*>> newArcToUsages;
        std::vector<std::vector<int>> arcIxToNewArcIxs(arcs.size());
        arcIxToNewArcIxs.resize(arcs.size());
        for (int arcIx = 0; arcIx < arcs.size(); ++arcIx) {
            std::vector<typename Polyline<K>::iterator> cuts;
            const auto& arc = arcs[arcIx];

            bool closed = *arc.vertices_begin() == *(--arc.vertices_end());

            bool needsToBeCut = false;
            for (auto vit = arc.vertices_begin(); vit != arc.vertices_end(); ++vit) {
                if (junctions.contains(*vit)) {
                    needsToBeCut = true;
                    break;
                }
            }

            if (!closed || !needsToBeCut) {
                cuts.push_back(arc.vertices_begin());
            }
            for (auto vit = arc.vertices_begin(); vit != arc.vertices_end(); ++vit) {
                if (junctions.contains(*vit)) {
                    cuts.push_back(vit);
                }
            }
            if (!closed || !needsToBeCut) {
                cuts.push_back(arc.vertices_end());
            }

            int newArcStart = newArcs.size();
            for (int i = 0; i < cuts.size() - 1; ++i) {
                if (cuts[i] == cuts[i+1]) continue;

                auto end = cuts[i+1];
                if (end != arc.vertices_end()) {
                    ++end;
                } else {
                    auto start = cuts[i];
                    ++start;
                    if (start == end) continue;
                }

                newArcs.emplace_back(cuts[i], end);
            }

            if (closed && needsToBeCut) {
                Polyline<K> lastArc;
                if (cuts.back() != arc.vertices_end()) {
                    auto start = cuts.back();
                    ++start;
                    if (start != arc.vertices_end()) {
                        for (auto vit = cuts.back(); vit != arc.vertices_end(); ++vit) {
                            lastArc.push_back(*vit);
                        }
                    }
                }
                if (cuts.front() != arc.vertices_begin()) {
                    auto end = cuts.front();
                    ++end;
                    auto begin = arc.vertices_begin();
                    if (lastArc.size() > 0) {
                        ++begin;
                    }
                    for (auto vit = begin; vit != end; ++vit) {
                        lastArc.push_back(*vit);
                    }
                }
                if (lastArc.size() > 0) {
                    newArcs.push_back(lastArc);
                }
            }

            int newArcEnd = newArcs.size();

            if (newArcEnd > newArcToUsages.size()) {
                newArcToUsages.resize(newArcEnd);
            }
            for (int i = newArcStart; i < newArcEnd; ++i) {
                arcIxToNewArcIxs[arcIx].push_back(i);
            }

            // Can be optimized by changing the start of the closed arc to a cut point.
        }

        for (const auto& arc : arcs) {
            std::cout << "===" << std::endl;
            for (auto vit = arc.vertices_begin(); vit != arc.vertices_end(); ++vit) {
                std::cout << *vit << std::endl;
            }
        }

        for (auto& feature : features) {
            PolygonSetGeometry& polygonSetGeometry = get<PolygonSetGeometry>(feature.geometry);
            for (auto& polygonWithHoles : polygonSetGeometry.arcs) {
                for (auto& polygon : polygonWithHoles) {
                    if (polygon.empty()) continue;
                    assert(polygon.size() == 1);
                    auto& newIxs = arcIxToNewArcIxs[polygon[0]];
                    polygon.clear();
                    polygon = newIxs;
                }
            }
        }

        arcs = newArcs;

        // 4. Dedup

        // First modify all closed arcs such that they start at the left-most (then bottom-most) vertex.
        for (auto& arc : arcs) {
            if (*arc.vertices_begin() == *(--arc.vertices_end())) {
                // The arc is closed (it is a "ring")
                auto leftVertexIt = std::min_element(arc.vertices_begin(), arc.vertices_end(), [](const Point<K>& p1, const Point<K>& p2) {
                    if (p1.x() == p2.x()) return p1.y() < p2.y();
                    return p1.x() < p2.x();
                });
                if (leftVertexIt != arc.vertices_begin()) {
                    Polyline<K> newArc(leftVertexIt, arc.vertices_end());
                    for (auto vit = ++arc.vertices_begin(); vit != leftVertexIt; ++vit) {
                        newArc.push_back(*vit);
                    }
                    newArc.push_back(*leftVertexIt);
                    arc = newArc;
                }
            }
        }

        newArcs.clear();
//        newArcToUsages.clear();
//        newArcToUsages.resize(0);
        std::unordered_map<Point<K>, std::vector<int>> possibleDuplicates;
        for (int arcIx = 0; arcIx < arcs.size(); ++arcIx) {
            const auto& arc = arcs[arcIx];
            auto& start = *arc.vertices_begin();
            auto& end = *(--arc.vertices_end());
            if (start == end) {
                // Arc is closed
                possibleDuplicates[*arc.vertices_begin()].push_back(arcIx);
            } else {
                auto leftBottomMost = (start.x() < end.x() || (start.x() == end.x() && start.y() < end.y())) ? start : end;
                possibleDuplicates[leftBottomMost].push_back(arcIx);
            }
        }

        std::vector<int> newArcIndex(arcs.size());
        newArcIndex.resize(arcs.size());
        for (auto& [_, vec] : possibleDuplicates) {
            std::vector<std::optional<int>> duplicate(vec.size());
            for (int i = 0; i < vec.size(); ++i) {
                duplicate[i] = std::nullopt;
            }
            for (int i = 0; i < vec.size(); ++i) {
                int arcIx1 = vec[i];
                for (int j = i+1; j < vec.size(); ++j) {
                    if (duplicate[j].has_value()) continue;
                    int arcIx2 = vec[j];
                    auto& arc1 = arcs[arcIx1];
                    auto& arc2 = arcs[arcIx2];
                    if (arc1.size() != arc2.size()) continue;
                    bool same = true;
                    if (*arc1.vertices_begin() == *arc2.vertices_begin() && arc1.vertex(1) == arc2.vertex(1)) { // same orientation
                        for (int vIx = 0; vIx < arc1.size(); ++vIx) {
                            if (arc1.vertex(vIx) != arc2.vertex(vIx)) {
                                same = false;
                                break;
                            }
                        }
                    } else { // different orientation
                        for (int vIx = 0; vIx < arc1.size(); ++vIx) {
                            if (arc1.vertex(vIx) != arc2.vertex(arc2.size() - 1 - vIx)) {
                                same = false;
                                break;
                            }
                        }
                    }

                    if (!same) continue;

                    // Found duplicate.
                    duplicate[j] = i;
                }
            }

            for (int i = 0; i < vec.size(); ++i) {
                int arcIx = vec[i];
                const auto& arc = arcs[arcIx];
                if (!duplicate[i]) {
                    newArcs.push_back(arc);
                    newArcIndex[arcIx] = newArcs.size() - 1;
                } else {
                    // Index of the replacement arc in the newArcs vector.
                    int replacementArcIx = newArcIndex[vec[*duplicate[i]]];
                    newArcIndex[arcIx] = replacementArcIx;
                }
            }
        }

        for (auto& feature : features) {
            PolygonSetGeometry& polygonSetGeometry = get<PolygonSetGeometry>(feature.geometry);
            for (auto& polygonWithHoles : polygonSetGeometry.arcs) {
                for (auto& polygon : polygonWithHoles) {
                    for (auto& v : polygon) {
                        v = newArcIndex[v];
                    }
                }
            }
        }

        arcs = newArcs;
    }
};

TopoSet<Exact> pretendExact(const TopoSet<Inexact>& topoSet);
TopoSet<Inexact> approximate(const TopoSet<Exact>& topoSet);
}