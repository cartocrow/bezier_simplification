#include "core/bezier_graph_2.h"
#include "core/straight_graph_2.h"

#include <cartocrow/core/core.h>
#include <cartocrow/core/vector_helpers.h>

namespace cartocrow::curved_simplification {
/// Approximates the given bezier graph by a straight-line graph.
/// nPoints is the number of points each curve is approximated by and should be at least 2.
template <class VD, class ED> // Expectation is that BG is a Bezier_graph_2
Straight_graph_2<VD, ED, Inexact> approximateBezierGraph(const Bezier_graph_2<VD, ED>& bg, int nPoints) {
    using BG = Bezier_graph_2<VD, ED>;
    using SG = Straight_graph_2<VD, ED, Inexact>;

    // First copy over all vertices from bg.
    std::unordered_map<const typename BG::Vertex*, typename BG::Vertex_handle> vmap;

    SG g;

    for (auto vit = bg.vertices_begin(); vit != bg.vertices_end(); ++vit) {
        auto vh = g.insert_vertex(vit->point());
        vmap[&*vit] = vh;
    }

    for (auto eit = bg.edges_begin(); eit != bg.edges_end(); ++eit) {
        auto curveStartVertex = vmap[&*eit->source()];
        auto curveEndVertex = vmap[&*eit->target()];
        
        auto lastVertex = curveStartVertex;

        CubicBezierCurve curve = eit->curve();

        std::vector<Point<Inexact>> pts;
        curve.samplePoints(nPoints-1, std::back_inserter(pts));

        for (int i = 0; i < pts.size(); ++i) {
            auto endVertex = i == nPoints - 1 ? curveEndVertex : g.insert_vertex(pts[i]);
            g.insert_edge(lastVertex, endVertex);
            lastVertex = endVertex;
        }
    }

    return g;
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

/// Modifies the given straight graph such that its edges have a minimum distance from each other
template <class VD, class ED>
void forceDirectedMinimumDistance(Straight_graph_2<VD, ED, Inexact>& g, double minDist) {
	// For each vertex v, determine the force that acts on it.
	// To do is, find all edges that lie within minDist.
	// Find the point p on this edge closest to the vertex.
	// Apply a force proportional to v - p to v.
	// 
	// Edge case: edges are almost consecutive.
	// We still want to allow their vertices to push each other away, but not "along" the polyline. 
	// So here, compute the angle the vector makes with the polyline.

    if (!g.is_oriented()) {
        g.orient();
    }

    double stepSize = minDist / 10;

    for (auto vit = g.vertices_begin(); vit != g.vertices_end(); ++vit) {
        auto& v = *vit;
        auto p = v.point();
        Vector<Inexact> force;
        for (auto eit = g.edges_begin(); eit != g.edges_end(); ++eit) {
            if (eit->source() == vit || eit->target() == vit) continue;
            // Find closest point on eit to vit
            Segment<Inexact> seg = eit->curve();
            auto nearest = projection(seg, p);
            if (CGAL::squared_distance(p, nearest) < minDist * minDist) {
                // We want to push the vertex away
                auto vector = p - nearest;
                auto angle = smallestAngleBetween(vector, seg.to_vector());

                auto vertex = eit->target();
                bool badForce = false;
                for (int i = 0; i < 15; ++i) {
                    vertex = vertex->next();
                    if (vertex == vit) {// && angle < std::numbers::pi / 6) {
                        badForce = true;
                        break;
                    }
                }
                if (badForce) continue;
                vertex = eit->source();
                for (int i = 0; i < 15; ++i) {
                    vertex = vertex->prev();
                    if (vertex == vit) {// && angle < std::numbers::pi / 6) {
                        badForce = true;
                        break;
                    }
                }
                if (badForce) continue;

                force += vector;
            }
        }

        // Let's normalize and use a fixed small step size.
        if (force.squared_length() == 0) {
            return;
        }

        force /= sqrt(force.squared_length());
        force *= stepSize;

        // Set new point
        auto newP = p + force;
        v.point() = newP;
    }
}
}