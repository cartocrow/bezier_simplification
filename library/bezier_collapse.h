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
/// Collapse the edge, the one that comes before, and the one that comes after
struct Collapse3 {
    /// The Bézier after collapse that has point as a target
    CubicBezierCurve before;
    /// The Bézier after collapse that has point as a source
    CubicBezierCurve after;
};

/// Collapse the edge and the one that comes before
struct Collapse2 {
    /// The Bézier after collapse
    CubicBezierCurve replacement;
};

struct Collapse {
    /// Cost of collapse
    Number<Inexact> cost;

    std::variant<Collapse2, Collapse3> result;
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
struct HECIData;

using HECGraph = Bezier_graph_2<std::monostate, HECData>;
using HECIGraph = Bezier_graph_2<std::monostate, HECIData>;
using BezierCollapseGraph = Bezier_graph_2<std::monostate, ECData>;
using BezierCollapseGraphWithHistory = CollapseHistoryGraphAdaptor<HECGraph>;
using BezierCollapseGraphWithHistoryAndIndex = CollapseHistoryGraphAdaptor<HECIGraph>;

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
struct HECIData : ECBase<typename BezierCollapseGraphWithHistoryAndIndex::Edge_handle> {
    std::shared_ptr<Operation<HECIGraph>> hist;
    std::shared_ptr<Operation<HECIGraph>> futr;
    int index;
};
}
using BezierCollapseGraph = detail::BezierCollapseGraph;
using BezierCollapseGraphWithHistory = detail::BezierCollapseGraphWithHistory;
using BezierCollapseGraphWithHistoryAndIndex = detail::BezierCollapseGraphWithHistoryAndIndex;
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
#ifndef __EMSCRIPTEN__
        std::vector<std::future<void>> futures;
#endif

        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            eit->data().blocked_by.clear();
            eit->data().blocking.clear();
            eit->data().blocked_by_degzero = false;
#ifndef __EMSCRIPTEN__
            futures.emplace_back(std::async(std::launch::async, [eit, this]() {
#endif
                m_traits.determineCollapse(eit);
#ifndef __EMSCRIPTEN__
            }));
#endif
        }

#ifndef __EMSCRIPTEN__
        for (auto& f : futures) f.get();
#endif
        for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
            m_q.push(eit);
        }
	}

    std::pair<CubicBezierSpline, CubicBezierSpline>
    beforeAndAfterCollapse(Edge_handle collapse) {
        CubicBezierSpline beforeSpline;
        CubicBezierSpline afterSpline;
        auto& clps = *collapse->data().collapse;

        if (auto* clps3P = std::get_if<detail::Collapse3>(&clps.result)) {
            auto& clps3 = *clps3P;
            Edge_handle prev = collapse->prev();
            Edge_handle next = collapse->next();

            beforeSpline.appendCurve(prev->curve());
            beforeSpline.appendCurve(collapse->curve());
            beforeSpline.appendCurve(next->curve());

            afterSpline.appendCurve(clps3.before);
            afterSpline.appendCurve(clps3.after);
        } else if (auto* clps2P = std::get_if<detail::Collapse2>(&clps.result)) {
            auto& clps2 = *clps2P;

            Edge_handle prev = collapse->prev();

            beforeSpline.appendCurve(prev->curve());
            beforeSpline.appendCurve(collapse->curve());

            afterSpline.appendCurve(clps2.replacement);
        }

        return {beforeSpline, afterSpline};
    }

    bool blocks(Edge_handle edge, Edge_handle collapse) {
        auto& clps = *collapse->data().collapse;
        if (std::holds_alternative<detail::Collapse3>(clps.result)) {
            Edge_handle prev = collapse->prev();
            Edge_handle next = collapse->next();

            if (edge == collapse || edge == prev || edge == next) {
                // involved in collapse
                return false;
            }
        } else if (std::holds_alternative<detail::Collapse2>(clps.result)) {
            Edge_handle prev = collapse->prev();

            if (edge == collapse || edge == prev) {
                // involved in collapse
                return false;
            }
        }

        auto [beforeSpline, afterSpline] = beforeAndAfterCollapse(collapse);

        CubicBezierCurve curve = edge->curve();
        Box curveBox = curve.bbox();
        Box combBox = beforeSpline.bbox() + afterSpline.bbox();
        if (!do_overlap(combBox, curveBox)) return false;

        // Blocks if edge intersects afterSpline
        if (curve.sub(0.01, 0.99).intersects(afterSpline, 0.01)) {
            return true;
        }
        Rectangle<Inexact> combRect(combBox);
        Rectangle<Inexact> curveRect(curveBox);

        // Blocks if edge is in between the beforeSpline and afterSpline
        if (!utils::encloses(combRect, curveRect)) return false;

        // Because CGAL's Bézier functionality is a bit buggy (https://github.com/CGAL/cgal/issues/9176),
        // we use polyline approximations here for now.

        int nSegs = 20; // todo make parameter

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

            auto polygon = face_to_polygon_with_holes<Exact>(fit);
            auto& outer = polygon.outer_boundary();

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

    std::optional<Edge_handle> findNextStep() {
        // Note that the queue contains only collapses for which we have not checked yet if its collapse violated topology.
        while (!m_q.empty()) {
            Edge_handle e = m_q.pop();
            auto& edata = e->data();

            if (!edata.collapse.has_value()) continue;

            const auto& clps = edata.collapse->result;

            if (std::holds_alternative<detail::Collapse2>(clps)) {
                if (e->source()->degree() != 2) continue;
            } else if (std::holds_alternative<detail::Collapse3>(clps)) {
                if (e->source()->degree() != 2 || e->target()->degree() != 2) continue;
            }

            Edge_handle prev = e->prev();
            Edge_handle next = e->next();

            if (prev == next) continue;

            auto [beforeSpline, afterSpline] = beforeAndAfterCollapse(e);
            Box box = beforeSpline.bbox() + afterSpline.bbox();
            Rectangle<Exact> rect(box.xmin(), box.ymin(), box.xmax(), box.ymax());

            edata.blocked_by_degzero = false;

            // Check whether this collapse is blocked
            m_bcqt->findOverlapped(rect, [this, e](Edge_handle b) {
                if (blocks(b, e)) {
                    b->data().blocking.push_back(e);
                    e->data().blocked_by.push_back(b);
                    return true;
                }
                return false;
            });

            // If it is blocked, then we do not collapse (keep it removed from the queue).
            // The edges that block this collapse store a reference to this collapse, so that it can be
            // added back in the queue once the blocker edge has changed.
            if (!edata.blocked_by.empty()) continue;

            // The edge is not blocked
            return e;
        }

        return std::nullopt;
    }

    void clearBlockInfo(Edge_handle e) {
        // Add back collapses unblocked by removal of this edge.
        for (Edge_handle b: e->data().blocking) {
            if (utils::listRemove(e, b->data().blocked_by)) {
                if (b->data().blocked_by.empty() && !b->data().blocked_by_degzero) {
                    m_q.push(b);
                }
            }
        }
        e->data().blocking.clear();

        // Remove references to this edge in other edges block info.
        for (Edge_handle b: e->data().blocked_by) {
            utils::listRemove(e, b->data().blocking);
        }
        e->data().blocked_by.clear();
    }

    void performStep(Edge_handle e) {
        auto& edata = e->data();

        if (!edata.collapse.has_value()) return;

        const auto& clps = edata.collapse->result;

        if (auto* clps2P = std::get_if<detail::Collapse2>(&clps)) {
            auto& clps2 = *clps2P;
            Edge_handle prev = e->prev();

            const auto& newCurve = clps2.replacement;

            // Update queue
            m_q.remove(prev);

            // Update quadtree
            m_bcqt->remove(prev);
            m_bcqt->remove(e);

            // Add back collapses unblocked by removal of this and prev edge, and remove references to them.
            clearBlockInfo(e);
            clearBlockInfo(prev);

            // Perform the collapse
            Edge_handle eh;
            if constexpr (std::is_same_v<BG, BezierCollapseGraphWithHistoryAndIndex>) {
                auto index = e->data().index;
                eh = m_g.merge_edge_with_prev(e, newCurve);
                eh->data().index = index;
            } else {
                eh = m_g.merge_edge_with_prev(e, newCurve);
            }

            // Update quadtree
            m_bcqt->insert(eh);

            // Compute new collapses for new edge, and ones whose neighbor is the new edge (then their collapse changes)
            if (eh->source()->degree() == 2) {
                update(eh->prev());
            }
            update(eh);
            if (eh->target()->degree() == 2) {
                update(eh->next());
            }
        } else if (auto* clps3P = std::get_if<detail::Collapse3>(&clps)) {
            auto& clps3 = *clps3P;

            Edge_handle prev = e->prev();
            Edge_handle next = e->next();
            const auto& c0 = clps3.before;
            const auto& c1 = clps3.after;

            // Update queue
            m_q.remove(prev);
            m_q.remove(next);

            // Update quadtree
            m_bcqt->remove(prev);
            m_bcqt->remove(e);
            m_bcqt->remove(next);

            // Add back collapses unblocked by removal of this edge, and remove references to this edge.
            // Do the same for prev and next as they are also removed.
            clearBlockInfo(e);
            clearBlockInfo(prev);
            clearBlockInfo(next);

            // Perform the collapse
            Edge_handle eh1;
            Edge_handle eh2;
            if constexpr (std::is_same_v<BG, BezierCollapseGraphWithHistoryAndIndex>) {
                auto index = e->data().index;
                auto v = m_g.collapse_edge(e, c0, c1);
                eh1 = v->incoming();
                eh2 = v->outgoing();
                eh1->data().index = index;
                eh2->data().index = index;
            } else {
                auto v = m_g.collapse_edge(e, c0, c1);
                eh1 = v->incoming();
                eh2 = v->outgoing();
            }

            // Update quadtree
            m_bcqt->insert(eh1);
            m_bcqt->insert(eh2);

            // Compute new collapses for new edges
            if (eh1->source()->degree() == 2) {
                update(eh1->prev());
            }
            update(eh1);
            update(eh2);
            if (eh2->target()->degree() == 2) {
                update(eh2->next());
            }
        } else {
            throw std::runtime_error("Unknown collapse");
        }
    }

	bool step() {
        auto eh = findNextStep();
        if (!eh.has_value()) {
            return false;
        }

        performStep(*eh);
        return true;
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