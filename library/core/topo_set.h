#pragma once

#include <cartocrow/core/core.h>
#include <cartocrow/core/polyline.h>
#include <cartocrow/renderer/ipe_renderer.h>

#include "region_set.h"

#include <CGAL/General_polygon_set_2.h>

namespace cartocrow {
    template <class Arc>
    struct TopoSet;

    inline std::pair<int, bool> decodeArc(int v) {
        if (v < 0) {
            return { -v - 1, true };
        }
        return { v, false };
    }

    namespace detail {
        struct PolylineTopology {
            using Arc_handle = int;
            std::vector<Arc_handle> arcs;

            template <class K>
            Polyline<K> getGeometry(const TopoSet<Polyline<K>>& topoSet) const {
                Polyline<K> polyline;

                for (int h : arcs) {
                    auto [arcIx, reversed] = decodeArc(h);
                    const auto& pl = topoSet.arcs[arcIx];

                    if (!reversed) {
                        std::copy(pl.vertices_begin(), pl.vertices_end(),
                            std::back_inserter(polyline));
                    }
                    else {
                        std::copy(pl.vertices_rbegin(), pl.vertices_rend(),
                            std::back_inserter(polyline));
                    }
                }

                return polyline;
            }
        };

        struct PolygonTopology {
            using Arc_handle = int;
            std::vector<std::vector<Arc_handle>> arcs;


            template <class K>
            PolygonWithHoles<K> getGeometry(const TopoSet<Polyline<K>>& topoSet) const {
                PolygonWithHoles<K> polygonWithHoles;

                for (const auto& polygonArcs : arcs) {
                    Polygon<K> polygon;
                    Point<K> last;

                    for (int h : polygonArcs) {
                        auto [arcIx, reversed] = decodeArc(h);
                        const auto& arc = topoSet.arcs[arcIx];

                        if (!reversed) {
                            std::copy(arc.vertices_begin(), arc.vertices_end() - 1,
                                std::back_inserter(polygon));

                            last = *(arc.vertices_end() - 1);
                        }
                        else {
                            std::copy(arc.vertices_rbegin(), arc.vertices_rend() - 1,
                                std::back_inserter(polygon));

                            last = *(arc.vertices_rend() - 1);
                        }
                    }

                    if (polygonWithHoles.is_empty()) {
                        if (polygon.is_clockwise_oriented()) {
                            polygon.reverse_orientation();
                        }
                        polygonWithHoles.outer_boundary() = polygon;
                    }
                    else {
                        if (polygon.is_counterclockwise_oriented()) {
                            polygon.reverse_orientation();
                        }
                        polygonWithHoles.add_hole(polygon);
                    }
                }

                return polygonWithHoles;
            }
        };

        struct PolygonSetTopology {
            using Arc_handle = int;
            std::vector<std::vector<std::vector<Arc_handle>>> arcs;


            template <class K>
            PolygonSetRaw<K> getGeometry(const TopoSet<Polyline<K>>& topoSet) const {
                PolygonSetRaw<K> polygonSet;

                for (const auto& polygonWithHolesArcs : arcs) {
                    PolygonWithHoles<K> polygonWithHoles;

                    for (const auto& polygonArcs : polygonWithHolesArcs) {
                        Polygon<K> polygon;
                        Point<K> last;

                        for (int h : polygonArcs) {
                            auto [arcIx, _] = decodeArc(h);
                            //auto [arcIx, reversed] = decodeArc(h);
                            const auto& arc = topoSet.arcs[arcIx];

                            //if (!polygon.is_empty()) {
                            //    if (*arc.vertices_begin() != last) {
                            //        reversed = false;
                            //        //std::cout << "Arcs do not connect!" << std::endl;
                            //        //std::cout << "Last " << last << std::endl;
                            //        //std::cout << "== arc ==" << std::endl;
                            //        //for (auto vit = arc.vertices_begin(); vit != arc.vertices_end(); ++vit) {
                            //        //    std::cout << *vit << std::endl;
                            //        //}

                            //    }
                            //}

                            //if (!reversed) {
                            //    std::copy(arc.vertices_begin(), arc.vertices_end() - 1,
                            //        std::back_inserter(polygon));

                            //    last = *(arc.vertices_end() - 1);
                            //}
                            //else {
                            //    if (!polygon.is_empty()) {
                            //        if (*arc.vertices_rbegin() != last) {
                            //            reversed = !reversed;
                            //            std::cout << "Arcs do not connect!" << std::endl;
                            //            std::cout << "Last " << last << std::endl;
                            //            std::cout << "== arc ==" << std::endl;
                            //            for (auto vit = arc.vertices_begin(); vit != arc.vertices_end(); ++vit) {
                            //                std::cout << *vit << std::endl;
                            //            }
                            //        }
                            //    }

                            //    std::copy(arc.vertices_rbegin(), arc.vertices_rend() - 1,
                            //        std::back_inserter(polygon));

                            //    last = *(arc.vertices_rend() - 1);
                            //}

                            bool reverse = false;
                            if (polygon.is_empty()) {
                                if (polygonArcs.size() > 1) {
                                    auto i1 = decodeArc(polygonArcs[1]).first;
                                    auto i0 = decodeArc(polygonArcs[0]).first;
                                    reverse = topoSet.arcs[i1].source() != topoSet.arcs[i0].target() && topoSet.arcs[i1].target() != topoSet.arcs[i0].target();

                                }
                            }
                            else {
                                reverse = *arc.vertices_begin() != last;
                            }

                            if (!reverse) {
                                if (!polygon.is_empty() && last != *arc.vertices_begin()) {
                                    std::cerr << "Arcs do not connect!" << std::endl;
                                }
                                std::copy(arc.vertices_begin(), arc.vertices_end() - 1, std::back_inserter(polygon));
                                last = *(arc.vertices_end() - 1);
                            }
                            else {
                                if (!polygon.is_empty() && last != *arc.vertices_rbegin()) {
                                    std::cerr << "Arcs do not connect!" << std::endl;
                                }
                                std::copy(arc.vertices_rbegin(), arc.vertices_rend() - 1, std::back_inserter(polygon));
                                last = *(arc.vertices_rend() - 1);
                            }
                        }

                        if (polygonWithHoles.is_empty()) {
                            if (!polygon.is_empty() && polygon.is_simple() && polygon.is_clockwise_oriented()) {
                                polygon.reverse_orientation();
                            }
                            polygonWithHoles.outer_boundary() = polygon;
                        }
                        else {
                            if (!polygon.is_empty() && polygon.is_simple() && polygon.is_counterclockwise_oriented()) {
                                polygon.reverse_orientation();
                            }
                            polygonWithHoles.add_hole(polygon);
                        }
                    }

                    polygonSet.polygons_with_holes.push_back(polygonWithHoles);
                }

                return polygonSet;
            }
        };

        using Topology = std::variant<PolylineTopology, PolygonTopology, PolygonSetTopology>;

        struct Feature {
            Topology topology;
            RegionAttributes attributes;
        };
    }

    struct TopoSetTopology {
        using Feature = detail::Feature;
        std::vector<Feature> features;
    };

    template <class Arc>
    struct TopoSet {
        std::vector<Arc> arcs;
        using Arc_handle = int;
        using PolylineTopology = detail::PolylineTopology;
        using PolygonTopology = detail::PolygonTopology;
        using PolygonSetTopology = detail::PolygonSetTopology;
        using Topology = detail::Topology;
        using Feature = detail::Feature;

        std::vector<Feature> features;

        TopoSet transform(const CGAL::Aff_transformation_2<Inexact>& trans) const {
            TopoSet transformed;
            transformed.features = features;
            for (const auto& arc : arcs) {
                transformed.arcs.push_back(arc.transform(trans));
            }
            return transformed;
        }

        Box bbox() const {
            return CGAL::bbox_2(arcs.begin(), arcs.end());
        }

        TopoSetTopology topology() const {
            return { features };
        }

        TopoSet() = default;
    };

    template <class K>
    using StraightTopoSet = TopoSet<Polyline<K>>;

    template <class K>
    StraightTopoSet<K> regionSetToTopoSet(const RegionSet<K>& regionSet) {
        StraightTopoSet<K> ts;
        std::vector<Polyline<K>>& arcs = ts.arcs;
        std::vector<typename StraightTopoSet<K>::Feature>& features = ts.features;
        using PolygonSetTopology = typename StraightTopoSet<K>::PolygonSetTopology;
        using Arc_handle = typename StraightTopoSet<K>::Arc_handle;

        auto encodeArc = [](int idx, bool reversed) {
            return reversed ? -(idx + 1) : idx;
            };

        // 1. Extract
        for (const auto& region : regionSet.regions) {
            const std::vector<PolygonWithHoles<K>>& pgns = region.geometry.polygons_with_holes;

            features.emplace_back(PolygonSetTopology{}, region.attributes);

            for (const PolygonWithHoles<K>& pgn : pgns) {
                PolygonSetTopology& geom = std::get<PolygonSetTopology>(features.back().topology);
                auto& pgnArcs = geom.arcs.emplace_back();

                const Polygon<K>& outer = pgn.outer_boundary();
                arcs.emplace_back(outer.vertices_begin(), outer.vertices_end());
                arcs.back().push_back(*outer.vertices_begin());

                auto& outerVec = pgnArcs.emplace_back();
                outerVec.push_back(arcs.size() - 1);

                for (const Polygon<K>& hole : pgn.holes()) {
                    arcs.emplace_back(hole.vertices_begin(), hole.vertices_end());
                    arcs.back().push_back(*hole.vertices_begin());

                    auto& holeVec = pgnArcs.emplace_back();
                    holeVec.push_back(arcs.size() - 1);
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

                if (i > 0) currentNeighbours.push_back(arc.vertex(i - 1));
                else if (closed) currentNeighbours.push_back(arc.vertex(arc.size() - 2));

                if (i < arc.size() - 1) currentNeighbours.push_back(arc.vertex(i + 1));
                else if (closed) currentNeighbours.push_back(arc.vertex(1));

                auto v = arc.vertex(i);

                if (!neighbours.contains(v)) {
                    neighbours[v] = currentNeighbours;
                }
                else {
                    for (const auto& n : currentNeighbours) {
                        if (std::find(neighbours[v].begin(), neighbours[v].end(), n) == neighbours[v].end()) {
                            neighbours[v].push_back(n);
                            junctions.insert(v);
                        }
                    }
                }
            }
        }

        // 3. Cut arcs
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
                if (cuts[i] == cuts[i + 1]) continue;

                auto end = cuts[i + 1];
                if (end != arc.vertices_end()) {
                    ++end;
                }
                else {
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

        for (auto& feature : features) {
            PolygonSetTopology& polygonSetTopology = get<PolygonSetTopology>(feature.topology);
            for (auto& polygonWithHoles : polygonSetTopology.arcs) {
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

        std::vector<Polyline<K>> dedupArcs;
        std::unordered_map<Point<K>, std::vector<int>> possibleDuplicates;

        for (int i = 0; i < arcs.size(); ++i) {
            auto& arc = arcs[i];
            auto& start = *arc.vertices_begin();
            auto& end = *(--arc.vertices_end());

            Point<K> key = (start == end)
                ? start
                : (start.x() < end.x() || (start.x() == end.x() && start.y() < end.y()) ? start : end);

            possibleDuplicates[key].push_back(i);
        }

        std::vector<int> newArcIndex(arcs.size());

        for (auto& [_, vec] : possibleDuplicates) {

            std::vector<std::optional<int>> duplicate(vec.size(), std::nullopt);

            for (int i = 0; i < vec.size(); ++i) {
                for (int j = i + 1; j < vec.size(); ++j) {
                    if (duplicate[j]) continue;

                    auto& a1 = arcs[vec[i]];
                    auto& a2 = arcs[vec[j]];
                    if (a1.size() != a2.size()) continue;

                    bool same = true;
                    bool reversed = false;

                    if (*a1.vertices_begin() == *a2.vertices_begin() &&
                        a1.vertex(1) == a2.vertex(1)) {
                        for (int k = 0; k < a1.size(); ++k) {
                            if (a1.vertex(k) != a2.vertex(k)) { same = false; break; }
                        }
                    }
                    else {
                        reversed = true;
                        for (int k = 0; k < a1.size(); ++k) {
                            if (a1.vertex(k) != a2.vertex(a2.size() - 1 - k)) {
                                same = false;
                                break;
                            }
                        }
                    }

                    if (same) {
                        duplicate[j] = i;
                    }
                }
            }

            for (int i = 0; i < vec.size(); ++i) {
                int arcIx = vec[i];

                if (!duplicate[i]) {
                    dedupArcs.push_back(arcs[arcIx]);
                    newArcIndex[arcIx] = encodeArc(dedupArcs.size() - 1, false);
                }
                else {
                    int refIx = vec[*duplicate[i]];
                    int baseIdx = std::abs(newArcIndex[refIx]);

                    auto& a1 = arcs[arcIx];
                    auto& a2 = arcs[refIx];

                    bool reversed = !(
                        *a1.vertices_begin() == *a2.vertices_begin() &&
                        a1.vertex(1) == a2.vertex(1)
                        );

                    newArcIndex[arcIx] = encodeArc(baseIdx, reversed);
                }
            }
        }

        for (auto& feature : features) {
            PolygonSetTopology& topology = std::get<PolygonSetTopology>(feature.topology);

            for (auto& pwh : topology.arcs) {
                for (auto& polygon : pwh) {
                    for (auto& v : polygon) {
                        v = newArcIndex[v];
                    }
                }
            }
        }

        arcs = dedupArcs;

        return ts;
    }

    StraightTopoSet<Exact> pretendExact(const StraightTopoSet<Inexact>& topoSet);
    StraightTopoSet<Inexact> approximate(const StraightTopoSet<Exact>& topoSet);
}