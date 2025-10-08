#pragma once

#include <cartocrow/core/core.h>
#include "bezier_graph.h"
#include "common.h"
#include "indexed_priority_queue.h"
#include "curved_graph_with_history.h"

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
	{ e->data().collapse } -> std::same_as<Collapse>;
	{ e->data().blocked_by } -> std::same_as<std::vector<typename BG::Edge_handle>>;
	{ e->data().blocking } -> std::same_as<std::vector<typename BG::Edge_handle>>;
	{ BCT::determineCollapse(e) };
};
struct ECData;
struct HECData;
using HECGraph = BezierGraph<std::monostate, HECData>;
using BezierCollapseGraph = BezierGraph<std::monostate, ECData>;
using BezierCollapseGraphWithHistory = CurvedGraphWithHistoryAdaptor<HECGraph>;

struct ECBase {
    std::optional<Collapse> collapse;
	int qid;

	// Blocked information; todo: use this
	bool blocked_by_degzero;
	std::vector<typename BezierCollapseGraph::Edge_handle> blocked_by;
	std::vector<typename BezierCollapseGraph::Edge_handle> blocking;
};
struct ECData : ECBase {
};
struct HECData : ECBase {
    std::shared_ptr<Operation<HECGraph>> hist;
    std::shared_ptr<Operation<HECGraph>> futr;
};
}
using BezierCollapseGraph = detail::BezierCollapseGraph;
using BezierCollapseGraphWithHistory = detail::BezierCollapseGraphWithHistory;
template <class BG, class BCT>
class BezierCollapse {
  private:
	using Vertex_handle = BG::Vertex_handle;
	using Edge_handle = BG::Edge_handle;
	BG& m_g;
	BCT m_traits;
	IndexedPriorityQueue<GraphQueueTraits<Edge_handle, Inexact>> m_q;

	void update(Edge_handle e) {
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
		for (auto eit = m_g.edges_begin(); eit != m_g.edges_end(); ++eit) {
			update(eit);
		}
	}
	bool step() {

		while (!m_q.empty()) {
			Edge_handle e = m_q.pop();
			auto& edata = e->data();

			if (!edata.collapse.has_value()) continue;

			Edge_handle prev = e->prev();
			Edge_handle next = e->next();
			m_q.remove(prev);
			m_q.remove(next);

            const auto& clps = *edata.collapse;
			const auto& c0 = clps.before;
			const auto& c1 = clps.after;

            auto v = m_g.collapse_edge(e, c0, c1);
            auto eh1 = v->incoming();
            auto eh2 = v->outgoing();

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
	bool runToComplexity(int k) {
		while (m_g.num_edges() > k) {
            if (m_g.num_edges() % 1000 == 0) {
                std::cout << m_g.num_edges() << "\n";
            }
			if (!step()) {
				return false;
			}
		}
		return true;
	}
};
}