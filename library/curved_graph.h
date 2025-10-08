#pragma once

#include <cartocrow/core/core.h>

#include "utils.h"

namespace cartocrow::curved_simplification {
template <class VD, class ED, class CST> class CurvedGraphVertex;
template <class VD, class ED, class CST> class CurvedGraphEdge;

template <class VD, class ED, class CST>
class CurvedGraph {
	friend class CurvedGraphVertex<VD, ED, CST>;
	friend class CurvedGraphEdge<VD, ED, CST>;

  public:
	using Vertex = CurvedGraphVertex<VD, ED, CST>;
	using Edge = CurvedGraphEdge<VD, ED, CST>;
	using Pt = CST::Pt;
    using Curve = CST::Curve;
    using EdgeData = ED;

  private:
	using Vertex_container = std::list<Vertex>;
	using Edge_container = std::list<Edge>;
	Vertex_container m_vertices;
	Edge_container m_edges;
	bool m_oriented = false;
	size_t m_num_edges = 0;

  public:
	using Vertex_iterator = std::list<Vertex>::iterator;
	using Edge_iterator = std::list<Edge>::iterator;
	using Vertex_const_iterator = std::list<Vertex>::const_iterator;
	using Edge_const_iterator = std::list<Edge>::const_iterator;
	using Vertex_handle = std::list<Vertex>::iterator;
	using Edge_handle = std::list<Edge>::iterator;
	using Vertex_const_handle = Vertex_const_iterator;
	using Edge_const_handle = Edge_const_iterator;

	Vertex_iterator vertices_begin() { return m_vertices.begin(); }
	Vertex_iterator vertices_end() { return m_vertices.end(); }
	Edge_iterator edges_begin() { return m_edges.begin(); }
	Edge_iterator edges_end() { return m_edges.end(); }

	Vertex_const_iterator vertices_begin() const { return m_vertices.cbegin(); }
	Vertex_const_iterator vertices_end() const { return m_vertices.cend(); }
	Edge_const_iterator edges_begin() const { return m_edges.cbegin(); }
	Edge_const_iterator edges_end() const { return m_edges.cend(); }

	size_t num_edges() {
		return m_num_edges;
	}

	Vertex_handle insert_vertex(Pt p) {
		Vertex& v = m_vertices.emplace_back();
		v.m_point = std::move(p);
		return (--m_vertices.end());
	}
	void remove_vertex(Vertex_handle v) {
		m_vertices.erase(v);
	}
	Edge_handle add_edge(Vertex_handle source, Vertex_handle target, Curve curve) {
		Edge& e = m_edges.emplace_back(std::move(source), std::move(target), std::move(curve));
		Edge_handle eh = --m_edges.end();
		Vertex_handle es = e.m_source;
		Vertex_handle et = e.m_target;

		std::vector<Edge_handle>& esInc = es->m_incident;
		std::vector<Edge_handle>& etInc = et->m_incident;
		// Make sure orientation is preserved
		esInc.push_back(eh);
		if (etInc.empty()) {
			etInc.push_back(eh);
		} else {
			etInc.push_back(etInc.front());
			etInc[0] = eh;
		}
		++m_num_edges;
		return eh;
	}
	void remove_edge(Edge_handle edge) {
		auto& sInc = edge->m_source->m_incident;
		auto& tInc = edge->m_target->m_incident;
		utils::vectorRemove(edge, sInc);
		utils::vectorRemove(edge, tInc);
		m_edges.erase(edge);
		--m_num_edges;
	}
    /// Split a vertex into two. Or equivalently, replace two curves by three new ones.
    /// Connect the two new vertices with three new curves.
    Edge_handle split_vertex(Vertex_handle v, Curve c0, Curve c1, Curve c2) {
        assert(c0.target() == c1.source());
        assert(c1.target() == c2.source());

        auto e0 = v->incoming();
        auto v0 = e0->source();
        auto e1 = v->outgoing();
        auto v3 = e1->target();
        remove_edge(e0);
        remove_edge(e1);
        remove_vertex(v);
        auto v1 = insert_vertex(c0.target());
        auto v2 = insert_vertex(c1.target());
        add_edge(v0, v1, c0);
        auto eh = add_edge(v1, v2, c1);
        add_edge(v2, v3, c2);

        return eh;
    }
    /// Collapse an edge, removing it together with the edge that comes before and after it.
    /// Replace it with a vertex connected by two new curves.
    /// This is equivalent to replacing three curves by two new ones.
    Vertex_handle collapse_edge(Edge_handle e, Curve toNewPoint, Curve fromNewPoint) {
        assert(toNewPoint.target() == fromNewPoint.source());
        Edge_handle prev = e->prev();
        Edge_handle next = e->next();
        Vertex_handle eSource = e->source();
        Vertex_handle eTarget = e->target();
        Vertex_handle prevSource = prev->source();
        Vertex_handle nextTarget = next->target();

        // todo: reuse edge/vertex objects and move them

        // Remove old edges and vertices
        remove_edge(e);
        remove_edge(prev);
        remove_edge(next);
        remove_vertex(eSource);
        remove_vertex(eTarget);

        // Add new vertex
        auto v = insert_vertex(toNewPoint.target());
        // Add the two new edges
        auto eh1 = add_edge(prevSource, v, toNewPoint);
        auto eh2 = add_edge(v, nextTarget, fromNewPoint);

        return v;
    }
	void orient() {
		for (Vertex_handle v = m_vertices.begin(); v != m_vertices.end(); ++v) {
			if (v->degree() != 2) continue; // irrelevant for orientation

			Edge_handle bwd = v->m_incident[0];
			Edge_handle fwd = v->m_incident[1];

			if (bwd->m_target == v && fwd->m_source == v) continue; // already satisfies orientation

			if (bwd->m_target != v) {
				bwd->reverse();
			}

			if (fwd->m_source != v) {
				fwd->reverse();
			}

			Vertex_handle fv = fwd->m_target;
			while (fv->degree() == 2 && fv != v) {
				if (fv->m_incident[0] != fwd) {
					fv->m_incident[1] = fv->m_incident[0];
					fv->m_incident[0] = fwd;
				}

				fwd = fv->m_incident[1];
				if (fwd->m_source != fv) {
					fwd->reverse();
				}
				fv = fwd->m_target;
			}

			if (fv != v) {
				Vertex_handle bv = bwd->m_source;
				while (bv->degree() == 2) {
					if (bv->m_incident[1] != bwd) {
						bv->m_incident[0] = bv->m_incident[1];
						bv->m_incident[1] = bwd;
					}

					bwd = bv->m_incident[0];
					if (bwd->m_target != bv) {
						bwd->reverse();
					}
					bv = bwd->m_source;
				}
			}
		}

		m_oriented = true;
	}
	bool verifyOriented() {
		for (Vertex_handle v = m_vertices.begin(); v != m_vertices.end(); ++v) {
			if (v->degree() != 2) continue; // irrelevant for orientation

			auto bwd = v->m_incident[0];
			auto fwd = v->m_incident[1];

			if (bwd->m_target == v && fwd->m_source == v) continue;

			return false;
		}

		return true;
	}
};

template <class VD, class ED, class CST>
class CurvedGraphVertex {
	friend class CurvedGraph<VD, ED, CST>;
	friend class CurvedGraphEdge<VD, ED, CST>;

  public:
	using Vertex = CurvedGraphVertex<VD, ED, CST>;
	using Edge = CurvedGraphEdge<VD, ED, CST>;
	using Graph = CurvedGraph<VD, ED, CST>;
	using Vertex_handle = Graph::Vertex_handle;
	using Edge_handle = Graph::Edge_handle;
	using Pt = CST::Pt;

  private:
	using Edge_container = std::vector<Edge_handle>;
	Pt m_point;
	Edge_container m_incident = Edge_container();

  public:
	const Pt& point() const {
		return m_point;
	}
	int degree() {
		return m_incident.size();
	}
	Edge_handle incoming() {
		assert(degree() == 2);
		return m_incident[0];
	}
	Edge_handle outgoing() {
		assert(degree() == 2);
		return m_incident[1];
	}
	Vertex_handle prev() {
		return incoming()->source();
	}
	Vertex_handle next() {
		return outgoing()->target();
	}
};

template <class VD, class ED, class CST>
class CurvedGraphEdge {
	friend class CurvedGraph<VD, ED, CST>;
	friend class CurvedGraphVertex<VD, ED, CST>;

  public:
	using Vertex = CurvedGraphVertex<VD, ED, CST>;
	using Edge = CurvedGraphEdge<VD, ED, CST>;
	using Graph = CurvedGraph<VD, ED, CST>;
	using Vertex_handle = Graph::Vertex_handle;
	using Vertex_const_handle = Graph::Vertex_const_handle;
	using Edge_handle = Graph::Edge_handle;
	using Edge_const_handle = Graph::Edge_const_handle;
	using Pt = CST::Pt;
	using Curve = CST::Curve;

  private:
	Vertex_handle m_source;
	Vertex_handle m_target;
	Curve m_curve;
	ED m_d;

  public:
	CurvedGraphEdge(Vertex_handle source, Vertex_handle target, Curve curve)
		: m_source(std::move(source)), m_target(std::move(target)), m_curve(std::move(curve)) {
	    assert(m_curve.source() == m_source->m_point);
		assert(m_curve.target() == m_target->m_point);
	};

	CurvedGraphEdge(Vertex_handle source, Vertex_handle target, Curve curve, ED d)
	    : m_source(std::move(source)), m_target(std::move(target)), m_curve(std::move(curve)), m_d(std::move(d)) {};

	Vertex_handle source() { return m_source; }
	Vertex_handle target() { return m_target; }
	Vertex_const_handle source() const { return m_source; }
	Vertex_const_handle target() const { return m_target; }
	const Curve& curve() const { return m_curve; }
	void reverse() {
		std::swap(m_source, m_target);
		m_curve.reverse();
	}
	Edge_handle prev() {
		return m_source->incoming();
	}
	Edge_handle next() {
		return m_target->outgoing();
	}
	ED& data() {
		return m_d;
	}
};
}
