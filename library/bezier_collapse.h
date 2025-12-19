#pragma once

#include <cartocrow/core/core.h>
#include <cartocrow/core/arrangement_helpers.h>
#include "core/bezier_graph_2.h"
#include "common.h"
#include "indexed_priority_queue.h"
#include "curved_graph_with_history.h"
#include "bezier_curve_quad_tree.h"

#include "utils.h"

#include <future>

namespace cartocrow::curved_simplification {
namespace detail {
struct Collapse {
    /// Cost of collapse
    Number<Inexact> cost;
    /// The Bézier after collapse that has point as a target
    CubicBezierCurve before;
    /// The Bézier after collapse that has point as a source
    CubicBezierCurve after;
};

template <class BG, class BCT>
concept BCSetup = requires(typename BG::Edge_handle e) {
	{ e->data().collapse } -> std::same_as<std::optional<Collapse>>;
	{ e->data().blocked_by } -> std::same_as<std::vector<typename BG::Edge_handle>>;
	{ e->data().blocking } -> std::same_as<std::vector<typename BG::Edge_handle>>;
	{ BCT::determineCollapse(e) };
};
struct ECData;
struct HECData;
using HECGraph = Bezier_graph_2<std::monostate, HECData>;
using BezierCollapseGraph = Bezier_graph_2<std::monostate, ECData>;
using BezierCollapseGraphWithHistory = CollapseHistoryGraphAdaptor<HECGraph>;

template <class Edge_handle>
struct ECBase {
    std::optional<Collapse> collapse;
	int qid;

	// Blocked information
	bool blocked_by_degzero; // todo: implement
	std::vector<Edge_handle> blocked_by;
	std::vector<Edge_handle> blocking;
};
struct ECData : ECBase<typename BezierCollapseGraph::Edge_handle> {
};
struct HECData : ECBase<typename BezierCollapseGraphWithHistory::Edge_handle> {
    std::shared_ptr<Operation<HECGraph>> hist;
    std::shared_ptr<Operation<HECGraph>> futr;
};
}
using BezierCollapseGraph = detail::BezierCollapseGraph;
using BezierCollapseGraphWithHistory = detail::BezierCollapseGraphWithHistory;
template <class BG, class BCT>// requires detail::BCSetup<BG, BCT>
class BezierCollapse {
  private:
	using Vertex_handle = BG::Vertex_handle;
	using Edge_handle = BG::Edge_handle;
	BG& m_g;
	BCT m_traits;
    std::unique_ptr<BezierCurveQuadTree<Edge_handle, Exact>> m_bcqt; // optional as we create it only on initialize() call

  public:
    IndexedPriorityQueue<GraphQueueTraits<Edge_handle, Inexact>> m_q;

  private:
	void update(Edge_handle e) {
        auto& edata = e->data();

        // clear topology
        for (Edge_handle b : edata.blocked_by) {
            utils::listRemove(e, b->data().blocking);
        }
        edata.blocked_by.clear();

		m_traits.determineCollapse(e);
		if (m_q.contains(e)) {
			m_q.update(e);
		} else {
			m_q.push(e);
		}
	}

  public:
	/// Traits are passed as an object to allow the user to specify parameters.
	BezierCollapse(BG& graph, BCT traits) : m_g(graph), m_traits(std::move(traits)) {};
	void initialize() {
        std::vector<CubicBezierCurve> curves;
        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            curves.push_back(eit->curve());
        }
        Box bigBbox = CGAL::bbox_2(curves.begin(), curves.end());
        Rectangle<Exact> bigRect(bigBbox.xmin(), bigBbox.ymin(), bigBbox.xmax(), bigBbox.ymax());
        m_bcqt = std::make_unique<BezierCurveQuadTree<Edge_handle, Exact>>(bigRect, 10, 0.05, [](Edge_handle e) { return e->curve(); });

        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            m_bcqt->insert(eit);
        }

        m_q = {};
        std::vector<std::future<void>> futures;

        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            futures.emplace_back(std::async(std::launch::async, [eit, this]() {
                m_traits.determineCollapse(eit);
            }));
        }

        for (auto& f : futures) f.get();
        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            m_q.push(eit);
        }
	}

    bool blocks(Edge_handle edge, Edge_handle collapse) {
        Edge_handle prev = collapse->prev();
        Edge_handle next = collapse->next();
        const auto& clps = *collapse->data().collapse;

        if (edge == collapse || edge == prev || edge == next) {
            // involved in collapse
            return false;
        }

        int nSegs = 100; // todo make parameter
        CubicBezierSpline beforeSpline;
        beforeSpline.appendCurve(prev->curve());
        beforeSpline.appendCurve(collapse->curve());
        beforeSpline.appendCurve(next->curve());

        CubicBezierSpline afterSpline;
        afterSpline.appendCurve(clps.before);
        afterSpline.appendCurve(clps.after);

        // Because CGAL's Bézier functionality is a bit buggy (https://github.com/CGAL/cgal/issues/9176),
        // we use polyline approximations here for now.
        auto beforePl = beforeSpline.polyline(nSegs);
        auto afterPl = afterSpline.polyline(nSegs);
        std::vector<Arrangement<Exact>::X_monotone_curve_2> xmCurvesBefore;
        for (auto eit = beforePl.edges_begin(); eit != beforePl.edges_end(); ++eit) {
            xmCurvesBefore.emplace_back(pretendExact(*eit));
        }
        std::vector<Arrangement<Exact>::X_monotone_curve_2> xmCurvesAfter;
        for (auto eit = afterPl.edges_begin(); eit != afterPl.edges_end(); ++eit) {
            xmCurvesAfter.emplace_back(pretendExact(*eit));
        }
        Arrangement<Exact> arr;
        CGAL::insert_non_intersecting_curves(arr, xmCurvesBefore.begin(), xmCurvesBefore.end());
        CGAL::insert(arr, xmCurvesAfter.begin(), xmCurvesAfter.end());

        Polyline<Exact> testPl = pretendExact(edge->curve().polyline(nSegs));

        for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
            if (fit->is_unbounded()) continue;

            Arrangement<Exact> polygonPolylineArr;
            auto polygon = face_to_polygon_with_holes<Exact>(fit);
            auto& outer = polygon.outer_boundary();
            for (auto eit = outer.edges_begin(); eit != outer.edges_end(); ++eit) {
                CGAL::insert(polygonPolylineArr, *eit);
            }
            for (auto& h : polygon.holes()) {
                for (auto eit = h.edges_begin(); eit != h.edges_end(); ++eit) {
                    CGAL::insert(polygonPolylineArr, *eit);
                }
            }
            for (auto eit = testPl.edges_begin(); eit != testPl.edges_end(); ++eit) {
                CGAL::insert(polygonPolylineArr, *eit);
            }
            for (auto vit = polygonPolylineArr.vertices_begin(); vit != polygonPolylineArr.vertices_end(); ++vit) {
                if (vit->degree() > 2 && vit->point() != pretendExact(beforeSpline.source()) && vit->point() != pretendExact(beforeSpline.target())) {
                    return true;
                }
            }
            auto arbitraryPoint = CGAL::midpoint(*testPl.edges_begin());
            bool arbitraryPointInPolygon = outer.has_on_bounded_side(arbitraryPoint);
            for (auto& h : polygon.holes()) {
                if (h.has_on_bounded_side(arbitraryPoint)) {
                    arbitraryPointInPolygon = false;
                    break;
                }
            }
            if (arbitraryPointInPolygon) return true;
        }

        return false;
    }

	bool step() {
		while (!m_q.empty()) {
			Edge_handle e = m_q.pop();
			auto& edata = e->data();

			if (!edata.collapse.has_value()) continue;

			Edge_handle prev = e->prev();
			Edge_handle next = e->next();
            if (prev == next) continue;

            const auto& clps = *edata.collapse;
			const auto& c0 = clps.before;
			const auto& c1 = clps.after;

            // possibly blocked?
            Box box = c0.bbox() + c1.bbox();
            Rectangle<Exact> rect(box.xmin(), box.ymin(), box.xmax(), box.ymax());

            edata.blocked_by_degzero = false;

            // todo: degzero vertices
//            pqt.findContained(rect, [&edata](Vertex& b) {
//                if (!edata.T1.has_on_unbounded_side(b.getPoint()) ||
//                    !edata.T2.has_on_unbounded_side(b.getPoint())) {
//                    // blocked, by an unmovable vertex
//                    edata.blocked_by_degzero = true;
//                }
//            });

//            if (!edata.blocked_by_degzero) {

                m_bcqt->findOverlapped(rect, [this, e](Edge_handle b) {
                    if (blocks(b, e)) {
                        b->data().blocking.push_back(e);
                        e->data().blocked_by.push_back(b);
                        std::cout << "Edge: " << e->curve().source() << " -> " << e->curve().target()
                                  << " is blocked by: " << b->curve().source() << " -> " << b->curve().target();

                    }
                });
//            }

            if (//edata.blocked_by_degzero ||
                !edata.blocked_by.empty()) continue;

            m_q.remove(prev);
            m_q.remove(next);

            m_bcqt->remove(prev);
            m_bcqt->remove(e);
            m_bcqt->remove(next);

            auto v = m_g.collapse_edge(e, c0, c1);
            auto eh1 = v->incoming();
            auto eh2 = v->outgoing();

            m_bcqt->insert(eh1);
            m_bcqt->insert(eh2);

            if (eh1->source()->degree() == 2) {
				update(eh1->prev());
			}
			update(eh1);
			update(eh2);
			if (eh2->target()->degree() == 2) {
				update(eh2->next());
			}

			return true;
		}
		return false;
	}
	bool runToComplexity(int k, std::optional<std::function<void(int)>> progress = std::nullopt,
                         std::optional<std::function<bool()>> cancelled = std::nullopt) {
		while (m_g.number_of_edges() > k) {
            if (m_g.number_of_edges() % 100 == 0) {
                std::cout << m_g.number_of_edges() << "\n";
            }
            if (progress.has_value()) {
                (*progress)(m_g.number_of_edges());
            }
            if (cancelled.has_value() && (*cancelled)()) {
                break;
            }
			if (!step()) {
				return false;
			}
		}
		return true;
	}
};
}