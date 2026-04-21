#include "core/bezier_graph_2.h"
#include "core/straight_graph_2.h"

#include "conic_types.h"
#include "intersection_helpers.h"
#include "schneider.h"
#include "steven_bezier_collapse.h"

#include <cartocrow/core/core.h>
#include <cartocrow/core/vector_helpers.h>
#include <cartocrow/core/segment_delaunay_graph_helpers.h>
#include <cartocrow/core/delaunay_voronoi_helpers.h>

#include <cartocrow/renderer/ipe_renderer.h>

#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>

#include <chrono>
#include <cmath>

#define DEBUG_MIN_DIST 0

namespace cartocrow::curved_simplification {
template <class BezierGraph>
struct ApproximatedBezierGraphVertexData {
    std::optional<typename BezierGraph::Vertex_const_handle> originalVertex;
};

template <class BezierGraph>
struct ApproximatedBezierGraphEdgeData {
    typename BezierGraph::Edge_const_handle originalEdge;
    bool changed;
};

template<typename T>
concept HasChanged = requires(T t) {
    { t.changed } -> std::convertible_to<bool>;
};

static_assert(HasChanged<ApproximatedBezierGraphEdgeData<Bezier_graph_2<std::monostate, std::monostate>>>);

template <class BezierGraph>
using ApproximatedBezierGraph = Straight_graph_2<ApproximatedBezierGraphVertexData<BezierGraph>, ApproximatedBezierGraphEdgeData<BezierGraph>, Inexact>;

/// Approximates the given bezier graph by a straight-line graph.
/// nPoints is the number of points each curve is approximated by and should be at least 2.
template <class BezierGraph> // Expectation is that BG is a Bezier_graph_2
ApproximatedBezierGraph<BezierGraph>
approximateBezierGraph(const BezierGraph& bg, int nPoints) {
    using BG = BezierGraph;
    using SG = ApproximatedBezierGraph<BezierGraph>;

    // First copy over all vertices from bg.
    std::unordered_map<const typename BG::Vertex*, typename SG::Vertex_handle> vmap;

    SG g;

    for (auto vit = bg.vertices_begin(); vit != bg.vertices_end(); ++vit) {
        auto vh = g.add_vertex(vit->point());
        vh->data().originalVertex = vit;
        vmap[&*vit] = vh;
    }

    for (auto eit = bg.edges_begin(); eit != bg.edges_end(); ++eit) {
        auto curveStartVertex = vmap[&*eit->source()];
        auto curveEndVertex = vmap[&*eit->target()];
        
        auto lastVertex = curveStartVertex;

        CubicBezierCurve curve = eit->curve();

        std::vector<Point<Inexact>> pts;
        if (isStraight(curve)) {
            pts.push_back(curve.source());
            pts.push_back(curve.target());
        } else {
            curve.samplePoints(nPoints - 1, std::back_inserter(pts));
        }

        for (int i = 1; i < pts.size(); ++i) {
            auto endVertex = i == pts.size() - 1 ? curveEndVertex : g.add_vertex(pts[i]);
            auto eh = g.add_edge(lastVertex, endVertex, {lastVertex->point(), endVertex->point()});
            eh->data().originalEdge = eit;
            eh->data().changed = false;
            lastVertex = endVertex;
        }
    }

    return g;
}

template <class BezierGraph>
BezierGraph reconstructBezierGraph(const ApproximatedBezierGraph<BezierGraph>& sg, double maxSquaredError) {
    using BG = BezierGraph;
    using SG = ApproximatedBezierGraph<BezierGraph>;

    // First copy over all endpoint vertices from sg.
    std::unordered_map<const typename SG::Vertex*, typename BG::Vertex_handle> vmap;
    std::unordered_map<const typename BG::Vertex*, typename BG::Vertex_handle> vmapbb;

    BG bg;

    for (auto vit = sg.vertices_begin(); vit != sg.vertices_end(); ++vit) {
        if (vit->data().originalVertex.has_value()) {
            auto vh = bg.add_vertex(vit->point());
            auto& og = *vit->data().originalVertex;
            vh->data() = og->data();
            vmap[&*vit] = vh;
            vmapbb[&*og] = vh;
        }
    }

    for (auto vit = sg.vertices_begin(); vit != sg.vertices_end(); ++vit) {
        auto& og = vit->data().originalVertex;
        if (og.has_value()) {
            for (auto eit = vit->incident_edges_begin(); eit != vit->incident_edges_end(); ++eit) {
                auto eh = *eit;
                if (eh->source() != vit) continue;
                auto& ogEdge = eh->data().originalEdge;
                auto otherVertex = ogEdge->source() == *og ? vmapbb[&*ogEdge->target()] : vmapbb[&*ogEdge->source()];

                std::vector<Point<Inexact>> points;
                points.push_back(eh->source()->point());

                const CubicBezierCurve& ogCurve = ogEdge->curve();

                auto current = eh;
                Vector<Inexact> t1;
                if (eh->source()->degree() != 2) {
                    t1 = ogCurve.tangent(0);
                } else {
                    auto v1 = eh->target()->point() - eh->source()->point();
                    v1 /= sqrt(v1.squared_length());
                    auto v2 = eh->source()->point() - eh->source()->prev()->point();
                    v2 /= sqrt(v2.squared_length());
                    t1 = v1 + v2;
                }
                bool changed = current->data().changed;
                while (!current->target()->data().originalVertex.has_value()) {
                    current = current->next();
                    points.push_back(current->target()->point());
                    if (current->data().changed) {
                        changed = true;
                    }
                }
                Vector<Inexact> t2;
                if (current->target()->degree() != 2) {
                    t2 = -ogCurve.tangent(1);
                } else {
                    auto v1 = current->source()->point() - current->target()->point();
                    v1 /= sqrt(v1.squared_length());
                    auto v2 = current->target()->point() - current->target()->next()->point();
                    v2 /= sqrt(v2.squared_length());
                    t2 = v1 + v2;
                }
                points.push_back(current->target()->point());

                if (!changed) {
                    auto eh = bg.add_edge(vmap[&*vit], otherVertex, ogCurve);
                    eh->data() = ogEdge->data();
                } else {
                    if (points.size() == 2) {
                        auto eh = bg.add_edge(vmap[&*vit], otherVertex, CubicBezierCurve(points[0], points[1]));
                        eh->data() = ogEdge->data();
                        continue;
                    }

                    std::vector<Point<Inexact>> subsampled;
                    for (int i = 0; i < points.size()-1; ++i) {
                        auto v = points[i+1] - points[i];

                        subsampled.push_back(points[i]);
                        subsampled.push_back(points[i] + 0.25 * v);
                        subsampled.push_back(points[i] + 0.5 * v);
                        subsampled.push_back(points[i] + 0.75 * v);
                    }
                    subsampled.push_back(points.back());

                    auto spline = fitSpline(subsampled, maxSquaredError, t1, t2, 500);

                    auto lastVertex = vmap[&*vit];
                    for (int i = 0; i < spline.numCurves(); ++i) {
                        auto c = spline.curve(i);
                        auto nextVertex = i == spline.numCurves() - 1 ? otherVertex : bg.add_vertex(c.target());
                        typename BezierGraph::Edge_handle eh;
                        if (!std::isfinite(c.source().x()) || !std::isfinite(c.source().y()) ||
                            !std::isfinite(c.sourceControl().x()) || !std::isfinite(c.sourceControl().y()) ||
                            !std::isfinite(c.targetControl().x()) || !std::isfinite(c.targetControl().y()) ||
                            !std::isfinite(c.target().x()) || !std::isfinite(c.target().y())) {
                            eh = bg.add_edge(lastVertex, nextVertex, CubicBezierCurve(lastVertex->point(), nextVertex->point()));
                        } else {
                            eh = bg.add_edge(lastVertex, nextVertex, c);
                        }
                        eh->data() = ogEdge->data();
                        lastVertex = nextVertex;
                    }
                }
            }
        }
    }

    return bg;
}

template <class K>
Point<K> projection(const Segment<K>& seg, const Point<K>& p) {
    Line<K> l = seg.supporting_line();
    auto q = l.projection(p);
    auto s = (q - seg.source()) * (seg.target() - seg.source()) / seg.squared_length();
    if (s < 0) return seg.source();
    if (s > 1) return seg.target();
    return q;
}

template <class VD, class ED>
class MinimumDistanceForcer {
    public:
    using StraightGraph = Straight_graph_2<VD, ED, Inexact>;

    typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<Inexact, CGAL::Field_with_sqrt_tag> Gt;

    using SiteInfo = std::optional<typename StraightGraph::Edge_handle>;

    struct Converter {
        SiteInfo operator()(const SiteInfo& info, bool is_src) {
            return info;
        }
        SiteInfo operator()(const SiteInfo& info1, const SiteInfo& info2, bool first) {
            if (info1.has_value()) {
                return *info1;
            }
            return info2;
        }
    };

    struct Merger {
        SiteInfo operator()(const SiteInfo& info1, const SiteInfo& info2) {
            if (info1.has_value()) {
                return *info1;
            }
            return info2;
        }
    };

    typedef CGAL::Segment_Delaunay_graph_2<Gt, CGAL::Segment_Delaunay_graph_storage_traits_with_info_2<Gt, SiteInfo, Converter, Merger>> SDG;

    using VoronoiEdge = std::variant<Segment<Inexact>, Line<Inexact>, Ray<Inexact>, CGAL::Parabola_segment_2<Gt>>;

    double length(const VoronoiEdge& vEdge) {
        if (auto sP = std::get_if<Segment<Inexact>>(&vEdge)) {
            auto& s = *sP;
            return CGAL::sqrt(s.squared_length());
        }
        if (auto psP = std::get_if<CGAL::Parabola_segment_2<Gt>>(&vEdge)) {
            auto& ps = *psP;
            std::vector<Point<Inexact>> points;
            ps.generate_points(points, 0.01);

            double length = 0;
            for (int i = 0; i < points.size() - 1; ++i) {
                length += CGAL::sqrt(CGAL::squared_distance(points[i], points[i+1]));
            }
            return length;
        }
        if (auto lP = std::get_if<Line<Inexact>>(&vEdge)) {
            auto& l = *lP;
            return std::numeric_limits<double>::infinity();
        }
        if (auto rP = std::get_if<Ray<Inexact>>(&vEdge)) {
            auto& r = *rP;
            return std::numeric_limits<double>::infinity();
        }
    }

    VoronoiEdge
    objectToVoronoiEdge(CGAL::Object o) {
        Segment<Inexact> s;
        Line<Inexact> l;
        Ray<Inexact> r;
        CGAL::Parabola_segment_2<Gt> ps;

        if (CGAL::assign(s, o)) {
            return s;
        }
        if (CGAL::assign(l, o)) {
            return l;
        }
        if (CGAL::assign(r, o)) {
            return r;
        }
        if (CGAL::assign(ps, o)) {
            return ps;
        }

        throw std::runtime_error("Object is not a segment, line, ray, or parabola segment");
    }

    /// Computes the smallest distance to a site attained by a point on this edge.
    /// Returns both the distance and the point.
    std::pair<double, Point<Inexact>> minDist(const SDG& delaunay, const SDG::Edge& edge) {
        VoronoiEdge vEdge = objectToVoronoiEdge(delaunay.primal(edge));
        if (auto sP = std::get_if<Segment<Inexact>>(&vEdge)) {
            auto& s = *sP;
            auto [p, q] = defining_sites<SDG>(edge);

            if (p.is_point() && q.is_point()) {
                // Check if midpoint of p and q is part of the segment.
                auto m = CGAL::midpoint(p.point(), q.point());
                auto v =  s.target() - s.source();
                auto dot = (m - s.source()) * v;
                if (dot >= 0 && dot <= v * v) {
                    // midpoint lies on segment
                    return {CGAL::sqrt(CGAL::squared_distance(p.point(), q.point())) / 2, CGAL::midpoint(p.point(), q.point())};
                }
                auto sDist2 = CGAL::squared_distance(s.source(), p.point());
                auto tDist2 = CGAL::squared_distance(s.target(), p.point());

                if (sDist2 < tDist2) {
                    return {CGAL::sqrt(sDist2), s.source()};
                } else {
                    return {CGAL::sqrt(tDist2), s.target()};
                }
            }

            auto sProj = projection(p.segment(), s.source());
            auto tProj = projection(p.segment(), s.target());
            auto sDist2 = CGAL::squared_distance(s.source(), sProj);
            auto tDist2 = CGAL::squared_distance(s.target(), tProj);
            if (sDist2 < tDist2) {
                return {CGAL::sqrt(sDist2), s.source()};
            } else {
                return {CGAL::sqrt(tDist2), s.target()};
            }
        }
        if (auto psP = std::get_if<CGAL::Parabola_segment_2<Gt>>(&vEdge)) {
            auto& ps = *psP;
            const auto& directrix = ps.line();
            const auto& focus = ps.center();
            auto focusProj = directrix.projection(focus);
            auto dirV = directrix.to_vector();
            auto psStartProj = directrix.projection(ps.p1);
            auto psEndProj   = directrix.projection(ps.p2);
            auto psStartDot  = (psStartProj - directrix.point()) * dirV;
            auto psEndDot    = (psEndProj - directrix.point()) * dirV;
            auto psFocusDot  = (focusProj - directrix.point()) * dirV;
            if ((psStartDot > psFocusDot || psEndDot > psFocusDot) && (psStartDot < psFocusDot || psEndDot < psFocusDot)) {
                // vertex is part of parabolic segment.
                return {CGAL::sqrt(CGAL::squared_distance(focus, focusProj)) / 2, CGAL::midpoint(focus, focusProj)};
            } else {
                // vertex is not part of parabolic segment.
                auto distS2 = CGAL::squared_distance(ps.p1, psStartProj);
                auto distE2 = CGAL::squared_distance(ps.p2, psEndProj);
                if (distS2 < distE2) {
                    return {CGAL::sqrt(distS2), ps.p1};
                } else {
                    return {CGAL::sqrt(distE2), ps.p2};
                }
            }
        }
        if (auto lP = std::get_if<Line<Inexact>>(&vEdge)) {
            auto [p, q] = defining_sites<SDG>(edge);
            return {CGAL::sqrt(CGAL::squared_distance(p.point(), q.point())) / 2, CGAL::midpoint(p.point(), q.point())};
        }
        if (auto rP = std::get_if<Ray<Inexact>>(&vEdge)) {
            auto& r = *rP;
            // The point with min. dist. is the midpoint of p and q if it lies on the ray, otherwise the source of ray
            auto [p, q] = defining_sites<SDG>(edge);
            auto m = CGAL::midpoint(p.point(), q.point());
            if ((m - r.source()) * r.to_vector() >= 0) {
                return {CGAL::sqrt(CGAL::squared_distance(m, p.point())), m};
            } else {
                return {CGAL::sqrt(CGAL::squared_distance(r.source(), p.point())), r.source()};
            }
        }

        throw std::runtime_error("Unexpected Segment Delaunay Graph edge type!");
    }

    std::optional<VoronoiEdge>
    withinDist(const SDG& delaunay, const SDG::Edge& edge, double dist) {
        auto vEdge = objectToVoronoiEdge(delaunay.primal(edge));

        if (minDist(delaunay, edge).first > dist) {
            return std::nullopt;
        }

        Conic_Arrangement arr;
        auto [p, q] = defining_sites<SDG>(edge);

        if (auto sP = std::get_if<Segment<Inexact>>(&vEdge)) {
            auto& s = *sP;
            if (!isfinite(s.source().x()) || !isfinite(s.target().x())) return std::nullopt;
            if (p.is_point() && q.is_point()) {
                Circle<Inexact> circ(p.point(), dist * dist, CGAL::COUNTERCLOCKWISE);
                return intersection(s, circ);
            }

            // p and q are segments.
            // distance changes linearly
            auto sProj = projection(p.segment(), s.source());
            auto tProj = projection(p.segment(), s.target());
            auto sDist2 = CGAL::squared_distance(s.source(), sProj);
            auto tDist2 = CGAL::squared_distance(s.target(), tProj);
            auto largerDist = std::max(sDist2, tDist2);
            auto smallerDist = std::min(sDist2, tDist2);

            if (largerDist < dist * dist) {
                return s;
            }

            if (smallerDist > dist * dist) {
                return std::nullopt;
            }

            auto sourceClosest = sDist2 < tDist2;

            auto narrow = sourceClosest ? s.source() : s.target();
            auto wide = sourceClosest ? s.target() : s.source();
            auto vec = wide - narrow;
            vec *= (dist - CGAL::sqrt(smallerDist)) / (CGAL::sqrt(largerDist) - CGAL::sqrt(smallerDist));
            return Segment<Inexact>(narrow, narrow + vec);
        }
        if (auto psP = std::get_if<CGAL::Parabola_segment_2<Gt>>(&vEdge)) {
            auto& ps = *psP;

            auto pointSite = p.is_point() ? p : q;
            auto point = pointSite.point();

            Circle<Inexact> circle({point.x(), point.y()}, dist * dist);
            return approximateIntersection(ps, circle, 0.01);
        }
        if (auto lP = std::get_if<Line<Inexact>>(&vEdge)) {
            auto& l = *lP;
            // todo
            return l;
        }
        if (auto rP = std::get_if<Ray<Inexact>>(&vEdge)) {
            auto& r = *rP;
            // todo
            return r;
        }
        throw std::runtime_error("Impossible: Segment Voronoi edge is not a line, ray, line segment, or parabola segment.");
    }

    std::variant<Point<Inexact>, Segment<Inexact>>
    site_projection(const VoronoiEdge& vEdge, const typename SDG::Site_2& site) {
        if (site.is_point()) {
            return { site.point() };
        } else {
            // Ray and line cases cannot occur because they require both sites to be a point
            if (auto sP = std::get_if<Segment<Inexact>>(&vEdge)) {
                auto s = *sP;
                auto start = site.segment().supporting_line().projection(s.source());
                auto end = site.segment().supporting_line().projection(s.end());
                return { Segment<Inexact>(start, end) };
            }
            else if (auto psP = std::get_if<CGAL::Parabola_segment_2<Gt>>(&vEdge)) {
                auto ps = *psP;
                // Roundabout way to obtain start and end of parabolic segment because they are protected -_-
                auto p1 = ps.p1;
                auto p2 = ps.p2;

                if (std::isnan(p1.x()) || std::isnan(p2.x())) {
                    return {typename Gt::Segment_2(typename Gt::Point_2(0.0, 0.0), typename Gt::Point_2(1.0, 0.0))};
                }

                auto start = site.segment().supporting_line().projection(p1);
                auto end = site.segment().supporting_line().projection(p2);

                return {Segment<Inexact>(start, end)};
            }
            else if (auto lP = std::get_if<Line<Inexact>>(&vEdge)) {
                return site.segment();
            }
            else if (auto rP = std::get_if<Ray<Inexact>>(&vEdge)) {
                // todo
                return site.segment();
            }
            else {
                throw std::runtime_error("Impossible: a segment Voronoi edge is neither a line segment nor a parabolic "
                                         "segment, but at least one of its sites is a line segment.");
            }
        }
    }

    std::pair<typename SDG::Vertex::Storage_site_2, typename SDG::Vertex::Storage_site_2>
    defining_storage_sites(const typename SDG::Edge& edge) {
        return {edge.first->vertex(SDG::cw(edge.second))->storage_site(),
                edge.first->vertex(SDG::ccw(edge.second))->storage_site()};
    }

    using GraphFeature = std::variant<typename StraightGraph::Vertex_handle, typename StraightGraph::Edge_handle>;

    GraphFeature
    static defining_graph_feature(const typename SDG::Vertex::Storage_site_2& s) {
        auto edge = *s.info();

        if (s.is_segment()) {
            return edge;
        }

        Point<Inexact> pt(s.point()->x(), s.point()->y());
        if (pt == edge->source()->point()) {
            return edge->source();
        } else {
            assert(pt == edge->target()->point());
            return edge->target();
        }
    }

    std::pair<GraphFeature, GraphFeature>
    defining_graph_features(const typename SDG::Edge& edge) {
        auto [ps, qs] = defining_storage_sites(edge);
        return {defining_graph_feature(ps), defining_graph_feature(qs)};
    }

    bool
    withinDistanceAlongIsoline(typename SDG::Vertex::Storage_site_2 p, typename SDG::Vertex::Storage_site_2 q, double dist) {
        auto nextFeature = [](const GraphFeature& f) -> std::optional<GraphFeature> {
            if (auto vhP = std::get_if<typename StraightGraph::Vertex_handle>(&f)) {
                auto& vh = *vhP;
                if (vh->degree() == 2) {
                    return vh->outgoing();
                } else if (vh->degree() == 1) {
                    return *vh->incident_edges_begin();
                }
            } else if (auto ehP = std::get_if<typename StraightGraph::Edge_handle>(&f)) {
                auto& eh = *ehP;
                return eh->target();
            } else {
                throw std::runtime_error("Unknown graph feature!");
            }
            return std::nullopt;
        };

        auto prevFeature = [](const GraphFeature& f) -> std::optional<GraphFeature> {
            if (auto vhP = std::get_if<typename StraightGraph::Vertex_handle>(&f)) {
                auto& vh = *vhP;
                if (vh->degree() == 2) {
                    return vh->incoming();
                }
            } else if (auto ehP = std::get_if<typename StraightGraph::Edge_handle>(&f)) {
                auto& eh = *ehP;
                return eh->source();
            } else {
                throw std::runtime_error("Unknown graph feature!");
            }
            return std::nullopt;
        };

        auto featureLength = [](const GraphFeature& f) -> double {
            if (std::get_if<typename StraightGraph::Vertex_handle>(&f)) {
                return 0;
            } else if (auto ehP = std::get_if<typename StraightGraph::Edge_handle>(&f)) {
                auto& eh = *ehP;
                return CGAL::sqrt(CGAL::squared_distance(eh->source()->point(), eh->target()->point()));
            }
            throw std::runtime_error("Unknown graph feature!");
        };

        auto show = [](const GraphFeature& f) -> std::string {
            std::stringstream ss;
            if (auto vhP = std::get_if<typename StraightGraph::Vertex_handle>(&f)) {
                auto& vh = *vhP;
                ss << "Vertex: " << vh->point();
            } else if (auto ehP = std::get_if<typename StraightGraph::Edge_handle>(&f)) {
                auto& eh = *ehP;
                ss << "Edge: " << eh->source()->point() << " -> " << eh->target()->point();
            } else {
                throw std::runtime_error("Unknown graph feature!");
            }
            return ss.str();
        };

        auto f1 = defining_graph_feature(p);
        auto f2 = defining_graph_feature(q);

        if (f1 == f2) return true;

        double traversedLength = 0;
        std::optional<GraphFeature> current = f1;
        do {
            current = nextFeature(*current);
            if (!current.has_value()) break;
            if (*current == f2) return true;
            traversedLength += featureLength(*current);
        } while (traversedLength < dist);

        traversedLength = 0;
        current = f1;

        do {
            current = prevFeature(*current);
            if (!current.has_value()) break;
            if (*current == f2) return true;
            traversedLength += featureLength(*current);
        } while (traversedLength < dist);

        return false;
    }

    using Component = std::vector<std::pair<VoronoiEdge, typename SDG::Edge>>;
    std::vector<Component> m_withinDistEdgeComponents;
    std::vector<std::variant<Point<Inexact>, Segment<Inexact>>> m_withinDistIsolineParts;
    SDG m_delaunay;
    double m_requiredMinDist = 0;
    double m_requiredLength = 0;
    double m_minLoopArea = 0;
    double m_minAdjDist;
    double m_minAngle;
    double m_averageEdgeLength = 0;
    bool m_ignoreBbox = false;
    StraightGraph& m_g;
    using VertexP = const typename StraightGraph::Vertex*;
    using EdgeP = const typename StraightGraph::Edge*;
    using Loop = std::vector<EdgeP>;
    std::vector<Loop> m_loops;
    std::vector<Segment<Inexact>> m_originalSegments;
    std::vector<Loop*> m_edgeToLoop;
    std::vector<Loop*> m_vertexToLoop;

    Box m_bbox;

    bool liesOnBbox(const Point<Inexact>& p) {
        return abs(p.x() - m_bbox.xmin()) < M_EPSILON ||
               abs(p.x() - m_bbox.xmax()) < M_EPSILON ||
               abs(p.y() - m_bbox.ymin()) < M_EPSILON ||
               abs(p.y() - m_bbox.ymax()) < M_EPSILON;
    }

    void recomputeDelaunay() {
        m_delaunay.clear();
        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            Segment<Inexact> seg(eit->source()->point(), eit->target()->point());
            if (m_ignoreBbox && (
                liesOnBbox(seg.source()) &&
                liesOnBbox(seg.target()))
//                abs(seg.source().x() - seg.target().x()) < M_EPSILON ||
//                abs(seg.source().y() - seg.target().y()) < M_EPSILON)
                ) {
                continue;
            }
            typename SDG::Site_2 site = Gt::Site_2::construct_site_2(seg.source(), seg.target());
            SiteInfo info = eit;
            m_delaunay.insert(site, info);
        }
    }

    // Compute a vector of minimum length from one feature to the other
    static Vector<Inexact> minimumDistanceVector(const typename SDG::Vertex::Site_2& p, const typename SDG::Vertex::Site_2& q) {
        if (p.is_point()) {
            if (q.is_point()) {
                return q.point() - p.point();
            } else {
                auto qS = q.segment();
                return projection(qS, p.point()) - p.point();
            }
        }
        else {
            auto pS = p.segment();
            if (q.is_point()) {
                return q.point() - projection(pS, q.point());
            }
            else {
                auto qS = q.segment();
                auto projS = projection(qS, pS.source());
                auto projT = projection(qS, pS.target());
                if (CGAL::squared_distance(projS, pS.source()) < CGAL::squared_distance(projT, pS.target())) {
                    return projS - pS.source();
                }
                else {
                    return projT - pS.target();
                }
            }
        }
    }

    bool angleSmallerThan(const typename SDG::Vertex::Storage_site_2& p, const typename SDG::Vertex::Storage_site_2& q, double maxAngle) {
        auto getVector = [this](const typename SDG::Vertex::Storage_site_2& site) {
            auto siteF = defining_graph_feature(site);
            if (site.is_segment()) {
                return std::get<typename StraightGraph::Edge_handle>(siteF)->curve().to_vector();
            } else {
                auto vh = std::get<typename StraightGraph::Vertex_handle>(siteF);
                if (vh->degree() != 2) {
                    return Vector<Inexact>(0, 0);
                }
                auto prev = vh->prev();
                auto next = vh->next();
                return next->point() - prev->point();;
            }
        };

        auto pv = getVector(p);
        auto qv = getVector(q);
        if (pv == Vector<Inexact>(0, 0) || qv == Vector<Inexact>(0, 0)) return true;

        auto v = minimumDistanceVector(p.site(), q.site());

        bool tooSmall = std::min(std::min(smallestAngleBetween(pv, v),
                                          smallestAngleBetween(qv, v)),
                                 std::min(smallestAngleBetween(pv, -v),
                                          smallestAngleBetween(qv, -v))) < maxAngle;

#if DEBUG_MIN_DIST
//        if (!tooSmall) {
//            using namespace renderer;
//            IpeRenderer ipeRenderer;
//            ipeRenderer.addPainting([p, q, pv, qv, tooSmall, v](GeometryRenderer &renderer) {
//                if (p.is_point()) {
//                    renderer.draw(Point<Inexact>(p.point()->x(), p.point()->y()));
//                    renderer.draw(Segment<Inexact>(Point<Inexact>(p.point()->x(), p.point()->y()),
//                                                   Point<Inexact>(p.point()->x(), p.point()->y()) + v));
//                }
//                if (p.is_segment()) {
//                    auto s = p.source_of_supporting_site();
//                    auto e = p.target_of_supporting_site();
//                    renderer.draw(Segment<Inexact>({s->x(), s->y()}, {e->x(), e->y()}));
//                    renderer.draw(Segment<Inexact>({e->x(), e->y()}, Point<Inexact>{e->x(), e->y()} + v));
//                }
//                if (q.is_point()) {
//                    renderer.draw(Point<Inexact>(q.point()->x(), q.point()->y()));
//                }
//                if (q.is_segment()) {
//                    auto s = q.source_of_supporting_site();
//                    auto e = q.target_of_supporting_site();
//                    renderer.draw(Segment<Inexact>({s->x(), s->y()}, {e->x(), e->y()}));
//                }
//
////            renderer.draw(Point<Inexact>{pv.x(), pv.y()});
////            renderer.draw(Point<Inexact>{qv.x(), qv.y()});
////            renderer.draw(Point<Inexact>{v.x(), v.y()});
////                std::stringstream text;
////                text << "Too small: " << tooSmall;
////                renderer.drawText({0, 0}, text.str());
//            });
//            std::stringstream fn;
//            using namespace std::chrono;
//            milliseconds ms = duration_cast<milliseconds>(
//                    system_clock::now().time_since_epoch()
//            );
//            fn << "debug_min_dist_" << ms << ".ipe";
//            std::cout << fn.str() << std::endl;
//            ipeRenderer.save(fn.str());
//        }
#endif

        return tooSmall;
    }

    Loop* getLoop(GraphFeature& f) {
        if (auto vhP = std::get_if<typename StraightGraph::Vertex_handle>(&f)) {
            if (m_vertexToLoop[(*vhP)->graph_index()] != nullptr) {
                return m_vertexToLoop.at((*vhP)->graph_index());
            }
        } else if (auto ehP = std::get_if<typename StraightGraph::Edge_handle>(&f)) {
            if (m_edgeToLoop[(*ehP)->graph_index()] != nullptr) {
                return m_edgeToLoop.at((*ehP)->graph_index());
            }
        }
        return nullptr;
    }

    bool partOfSameSmallLoop(GraphFeature& f, GraphFeature& g) {
        auto* fl = getLoop(f);
        auto* fg = getLoop(g);

        if (fl != nullptr && fg != nullptr) {
            if (fl == fg) return true;
        }
        return false;
    }

    bool filterVoronoiEdge(const typename SDG::Edge& e) {
        auto v1 = e.first->vertex(SDG::cw(e.second));
        auto v2 = e.first->vertex(SDG::ccw(e.second));
        if (!v1->storage_site().is_defined() || !v2->storage_site().is_defined()) return true;
        auto [p, q] = defining_sites<SDG>(e);
        // Skip edges that are defined by consecutive isoline edges
        if (p.is_segment() && q.is_segment() && (p.source() == q.source() || p.target() == q.target()  || p.source() == q.target() || p.target() == q.source()) ||
            p.is_point()   && q.is_segment() && (p.point()  == q.source() || p.point()  == q.target()) ||
            p.is_segment() && q.is_point()   && (p.source() == q.point()  || p.target() == q.point())) return true;
        // Skip edges that are defined by sites that lie on the bounding box.
        auto connectsToBbox = [this](const typename SDG::Site_2& site) {
            if (site.is_point()) {
                return liesOnBbox(site.point());
            } else {
                const auto& s = site.segment();
                return liesOnBbox(s.source()) || liesOnBbox(s.target());
            }
        };
        if (m_ignoreBbox && (connectsToBbox(p) || connectsToBbox(q))) return true;
        auto [ps, qs] = defining_storage_sites(e);

        auto pf = defining_graph_feature(ps);
        auto qf = defining_graph_feature(qs);

        if (partOfSameSmallLoop(pf, qf)) return true;
        if (withinDistanceAlongIsoline(ps, qs, m_minAdjDist) && angleSmallerThan(ps, qs, m_minAngle)) return true;
        return false;
    }

    bool filterComponent(const Component& comp) {
        double totalLength = 0;
        for (const auto& [vEdge, dEdge] : comp) {
            totalLength += length(vEdge);
        }
        return totalLength < m_requiredLength;
    }

    void recomputeAuxiliary() {
        m_withinDistIsolineParts.clear();
        m_withinDistEdgeComponents.clear();

        std::vector<std::vector<typename SDG::Edge>> components;
        bfsOnVoronoiEdges(m_delaunay, [&](const typename SDG::Edge& e) {
            if (filterVoronoiEdge(e)) return false;
            return minDist(m_delaunay, e).first <= m_requiredMinDist;
        }, std::back_inserter(components));

        for (const auto& comp : components) {
            std::vector<std::pair<VoronoiEdge, typename SDG::Edge>> component;
            for (const auto& e : comp) {
                auto vorEdge = *withinDist(m_delaunay, e, m_requiredMinDist);
                component.push_back({vorEdge, e});

                // project vorEdge on the sites.
                auto [p, q] = defining_sites<SDG>(e);
                m_withinDistIsolineParts.push_back(site_projection(vorEdge, p));
                m_withinDistIsolineParts.push_back(site_projection(vorEdge, q));
            }
            m_withinDistEdgeComponents.push_back(component);
        }
    }

    void computeAverageEdgeLength() {
        double totalLength = 0;
        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            auto d = CGAL::sqrt(CGAL::squared_distance(eit->source()->point(), eit->target()->point()));
            totalLength += d;
        }
        m_averageEdgeLength = totalLength / m_g.number_of_edges();
    }

    void identifyLoops() {
        m_loops.clear();
        m_edgeToLoop.clear();
        m_vertexToLoop.clear();

        std::unordered_set<const typename StraightGraph::Edge*> seen;

        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            if (seen.contains(&*eit)) continue;
            auto start = eit;
            auto current = start;

            std::vector<const typename StraightGraph::Edge*> component;
            do {
                seen.insert(&*current);
                component.push_back(&*current);
                if (current->target()->degree() != 2) break;
                current = current->next();
            } while (current != start);
            // now either current == start, or there is no next edge because the target vertex has degree != 2

            if (current->target()->degree() == 2) {
                Polygon<Inexact> polygon;
                for (auto edgeP : component) {
                    polygon.push_back(edgeP->source()->point());
                }
                if (polygon.area() < m_minLoopArea) {
                    m_loops.push_back(component);
                }
            }
        }

        for (auto& component : m_loops) {
            for (auto edgeP: component) {
                m_edgeToLoop[edgeP->graph_index()] = &m_loops.back();
                m_vertexToLoop[edgeP->source()->graph_index()] = &component;
            }
        }
    }

    MinimumDistanceForcer() = default;

    MinimumDistanceForcer(StraightGraph& graph, double minDist) : m_g(graph), m_requiredMinDist(minDist) {};

    void initialize() {
        m_originalSegments.resize(m_g.number_of_edges());
        m_edgeToLoop.resize(m_g.number_of_edges());
        m_vertexToLoop.resize(m_g.number_of_vertices());
        if (!m_g.is_indexed()) {
            m_g.index_vertices_and_edges();
        }
        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            m_originalSegments[eit->graph_index()] = eit->curve();
        }
        m_bbox = m_g.bbox();
        computeAverageEdgeLength();
        identifyLoops();
        recomputeDelaunay();
        recomputeAuxiliary();
    }

    bool step() {
        if (m_withinDistEdgeComponents.empty()) return false;

        // Iterate over all m_withinDistEdges.
        // Apply a force to the defining sites, its strength proportional to the length of the edge.

        // todo map vertices properly once they have indices.
        std::unordered_map<Point<Inexact>, Vector<Inexact>> forces;

        auto addForce = [&forces, this](typename StraightGraph::Vertex_handle vh, const Vector<Inexact>& force) {
            if (m_vertexToLoop[vh->graph_index()] != nullptr) {
                auto loopP = m_vertexToLoop.at(vh->graph_index());
                auto& loop = *loopP;
                for (auto eh : loop) {
                    if (!forces.contains(eh->source()->point())) {
                        forces[eh->source()->point()] = force / 8;
                    } else {
                        forces[eh->source()->point()] += force / 8;
                    }
                }
            } else {
                if (!forces.contains(vh->point())) {
                    forces[vh->point()] = force;
                } else {
                    forces[vh->point()] += force;
                }
            }
        };

        auto repel = [&](const GraphFeature& feature, const Point<Inexact>& repeller, double magnitude) {
            if (auto vhP = std::get_if<typename StraightGraph::Vertex_handle>(&feature)) {
                auto& vh = *vhP;
                Vector<Inexact> force = vh->point() - repeller;
                force /= CGAL::sqrt(force.squared_length());
                force *= magnitude;

                addForce(vh, force);
            } else if (auto ehP = std::get_if<typename StraightGraph::Edge_handle>(&feature)) {
                auto& eh = *ehP;

                const auto& seg = eh->curve();

                auto repellerProj = projection(seg, repeller);
                Vector<Inexact> force = repellerProj - repeller;

                auto targetPercentage = (repellerProj - seg.source()) * (seg.target() - seg.source()) / seg.squared_length();

                force /= CGAL::sqrt(force.squared_length());
                force *= magnitude;

                addForce(eh->source(), (1-targetPercentage) * force);
                addForce(eh->target(), targetPercentage * force);
            }
        };

        auto applyForce = [&](typename StraightGraph::Vertex_handle vh, const Vector<Inexact>& force) {
            if (liesOnBbox(vh->point())) return;
            vh->point() = vh->point() + force;
            for (auto eit = vh->incident_edges_begin(); eit != vh->incident_edges_end(); ++eit) {
                (*eit)->curve() = {(*eit)->source()->point(), (*eit)->target()->point()};
                if constexpr (HasChanged<ED>) {
                    auto& ed = (*eit)->data();
                    ed.changed = true;
                }
            }
        };

        for (const auto& comp : m_withinDistEdgeComponents) {
            if (filterComponent(comp)) continue;
            for (const auto& [vEdge, dEdge] : comp) {
                auto [pf, qf] = defining_graph_features(dEdge);

                auto repeller = minDist(m_delaunay, dEdge).second;
                auto magnitude = length(vEdge) * 0.025 + m_averageEdgeLength / 40;

                repel(pf, repeller, magnitude);
                repel(qf, repeller, magnitude);
            }
        }

        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
//             auto lengthRatio = eit->curve().squared_length() / m_originalSegments[&*eit].squared_length();
             auto diff = sqrt(eit->curve().squared_length()) - sqrt(m_originalSegments[eit->graph_index()].squared_length());
             auto v = eit->curve().to_vector();
             v /= sqrt(v.squared_length());
             addForce(eit->source(), diff / 2 * v / 2);
             addForce(eit->target(), -diff / 2 * v / 2);
        }

        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            auto s = eit->curve();
            auto v = s.to_vector();
            auto ogS = m_originalSegments[eit->graph_index()];
            auto ogV = ogS.to_vector();
            auto angleCw = orientedAngleBetween(v, ogV, CGAL::CLOCKWISE);
            auto angleCcw = orientedAngleBetween(v, ogV, CGAL::COUNTERCLOCKWISE);
            auto angle = abs(angleCw) < abs(angleCcw) ? -angleCw : angleCcw;
            CGAL::Aff_transformation_2<Inexact> translate(CGAL::TRANSLATION, midpoint(s) - CGAL::ORIGIN);
            CGAL::Aff_transformation_2<Inexact> rotate(CGAL::ROTATION, sin(angle / 4), cos(angle / 4));
            CGAL::Aff_transformation_2<Inexact> transform = translate * rotate * translate.inverse();

            if (abs(angle) < M_EPSILON) continue;
            auto rotatedS = eit->curve().transform(transform);
//#if DEBUG_MIN_DIST
//            using namespace renderer;
//            IpeRenderer ipeRenderer;
//            ipeRenderer.addPainting([&](GeometryRenderer& renderer) {
//                renderer.setStroke(Color(0, 0, 0), 1.0);
//                renderer.draw(ogS);
//                renderer.setStroke(Color(255, 0, 0), 1.0);
//                renderer.draw(Segment<Inexact>({0, 0}, Point<Inexact>{0, 0} + ogV));
//                renderer.setStroke(Color(0, 255, 0), 1.0);
//                renderer.draw(s);
//                renderer.setStroke(Color(0, 0, 255), 1.0);
//                renderer.draw(rotatedS);
//            });
//            std::stringstream fn;
//            using namespace std::chrono;
//            milliseconds ms = duration_cast<milliseconds>(
//                    system_clock::now().time_since_epoch()
//            );
//            fn << "debug_min_dist_rotate_" << ms << ".ipe";
//            std::cout << fn.str() << std::endl;
//            ipeRenderer.save(fn.str());
////
//#endif

            addForce(eit->source(), rotatedS.source() - s.source());
            addForce(eit->target(), rotatedS.target() - s.target());
        }

        for (auto vit = m_g.vertices_begin(); vit != m_g.vertices_end(); ++vit) {
            if (forces.contains(vit->point())) {
                Vector<Inexact> force = forces.at(vit->point());
//                force /= sqrt(force.squared_length());
//                force *= m_averageEdgeLength / 10;
                if (std::isfinite(force.x()) && !isnan(force.x()) && std::isfinite(force.y()) && !isnan(force.y())) {
                    applyForce(vit, force);
                }
            }
        }

        recomputeDelaunay();
        recomputeAuxiliary();

        return true;
    }
};
}