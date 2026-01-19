#pragma once

#include "graph_curve_traits_2.h"

namespace utils {
    template<typename T>
    bool vectorRemove(T elt, std::vector<T>& vec) {
        auto pos = std::find(vec.begin(), vec.end(), elt);
        if (pos != vec.end()) {
            vec.erase(pos);
            return true;
        }
        else {
            return false;
        }
    }
}

namespace cartocrow {
template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits> class Graph_2_vertex;
template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits> class Graph_2_edge;

template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits>
class Graph_2 {
	friend class Graph_2_vertex<VertexData, EdgeData, CurveTraits>;
	friend class Graph_2_edge<VertexData, EdgeData, CurveTraits>;

  public:
	using Vertex = Graph_2_vertex<VertexData, EdgeData, CurveTraits>;
	using Edge = Graph_2_edge<VertexData, EdgeData, CurveTraits>;
	using Point_2 = CurveTraits::Point_2;
    using Curve_2 = CurveTraits::Curve_2;
    using Vertex_data = VertexData;
    using Edge_data = EdgeData;
    using Curve_traits = CurveTraits;
    using Kernel = Curve_traits::Kernel;

  private:
	using Vertex_container = std::list<Vertex>;
	using Edge_container = std::list<Edge>;
	Vertex_container m_vertices;
	Edge_container m_edges;
	bool m_oriented = false;
    bool m_sorted = false;

  public:
	using Vertex_iterator = Vertex_container::iterator;
    using Vertex_const_iterator = Vertex_container::const_iterator;
    using Vertex_handle = Vertex_iterator;
    using Vertex_const_handle = Vertex_const_iterator;
	using Edge_iterator = Edge_container::iterator;
	using Edge_const_iterator = Edge_container::const_iterator;
	using Edge_handle = Edge_iterator;
	using Edge_const_handle = Edge_const_iterator;

    Graph_2& operator=(const Graph_2& other) {
        if (this == &other) return *this;

        m_vertices.clear();
        m_edges.clear();

        m_oriented = other.m_oriented;
        m_sorted   = other.m_sorted;

        std::unordered_map<const Vertex*, Vertex_handle> vmap;

        for (auto vit = other.m_vertices.begin(); vit != other.m_vertices.end(); ++vit) {
            m_vertices.emplace_back(*vit);
            auto new_vit = std::prev(m_vertices.end());
            new_vit->m_incident.clear();

            vmap[&*vit] = new_vit;
        }

        std::unordered_map<const Edge*, Edge_handle> emap;

        for (auto eit = other.m_edges.begin(); eit != other.m_edges.end(); ++eit) {
            Vertex_handle new_source = vmap.at(&*eit->m_source);
            Vertex_handle new_target = vmap.at(&*eit->m_target);

            m_edges.emplace_back(
                    new_source,
                    new_target,
                    eit->m_curve,
                    eit->m_d
            );

            auto new_eit = std::prev(m_edges.end());
            emap[&*eit] = new_eit;
        }

        for (auto vit = other.m_vertices.begin(); vit != other.m_vertices.end(); ++vit) {
            Vertex_handle new_v = vmap.at(&*vit);
            for (auto old_e : vit->m_incident) {
                new_v->m_incident.push_back(emap.at(&*old_e));
            }
        }

        return *this;
    }

    Graph_2 transform(CGAL::Aff_transformation_2<Kernel> trans) const {
        Graph_2 transformed;

        transformed.m_oriented = m_oriented;
        transformed.m_sorted   = m_sorted;

        std::unordered_map<const Vertex*, Vertex_handle> vmap;

        for (auto vit = m_vertices.begin(); vit != m_vertices.end(); ++vit) {
            transformed.m_vertices.emplace_back(*vit);
            auto new_vit = std::prev(transformed.m_vertices.end());
            new_vit->point() = new_vit->point().transform(trans);
            new_vit->m_incident.clear();

            vmap[&*vit] = new_vit;
        }

        std::unordered_map<const Edge*, Edge_handle> emap;

        for (auto eit = m_edges.begin(); eit != m_edges.end(); ++eit) {
            Vertex_handle new_source = vmap.at(&*eit->m_source);
            Vertex_handle new_target = vmap.at(&*eit->m_target);

            transformed.m_edges.emplace_back(
                    new_source,
                    new_target,
                    Curve_traits::transform(eit->m_curve, trans),
                    eit->m_d
            );

            auto new_eit = std::prev(transformed.m_edges.end());
            emap[&*eit] = new_eit;
        }

        for (auto vit = m_vertices.begin(); vit != m_vertices.end(); ++vit) {
            Vertex_handle new_v = vmap.at(&*vit);
            for (auto old_e : vit->m_incident) {
                new_v->m_incident.push_back(emap.at(&*old_e));
            }
        }

        return transformed;
    }

	Vertex_iterator vertices_begin() { return m_vertices.begin(); }
	Vertex_iterator vertices_end() { return m_vertices.end(); }
	Edge_iterator edges_begin() { return m_edges.begin(); }
	Edge_iterator edges_end() { return m_edges.end(); }

	Vertex_const_iterator vertices_begin() const { return m_vertices.cbegin(); }
	Vertex_const_iterator vertices_end() const { return m_vertices.cend(); }
	Edge_const_iterator edges_begin() const { return m_edges.cbegin(); }
	Edge_const_iterator edges_end() const { return m_edges.cend(); }

    size_t number_of_vertices() const {
        return m_vertices.size();
    }
	size_t number_of_edges() const {
        return m_edges.size();
	}
    void clear() {
        m_vertices.clear();
        m_edges.clear();
        m_oriented = false;
        m_sorted = false;
    }
	Vertex_handle insert_vertex(Point_2 p) {
		Vertex& v = m_vertices.emplace_back();
		v.m_point = std::move(p);
		return (--m_vertices.end());
	}
	void remove_vertex(Vertex_handle v) {
		m_vertices.erase(v);
	}
	Edge_handle add_edge(Vertex_handle source, Vertex_handle target, const Curve_2& curve) {
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
		return eh;
	}
	void remove_edge(Edge_handle edge) {
		auto& sInc = edge->m_source->m_incident;
		auto& tInc = edge->m_target->m_incident;
		utils::vectorRemove(edge, sInc);
		utils::vectorRemove(edge, tInc);
		m_edges.erase(edge);
	}
    /// Split a vertex into two. Or equivalently, replace two curves by three new ones.
    /// Connect the two new vertices with three new curves.
    Edge_handle split_vertex(Vertex_handle v, Curve_2 c0, Curve_2 c1, Curve_2 c2) {
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
    Vertex_handle collapse_edge(Edge_handle e, Curve_2 toNewPoint, Curve_2 fromNewPoint) {
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
    /// <summary>
    /// Ensures that all degree-2 vertices v have edge(0) = (u,v) and edge(1) = (v,w).
    /// </summary>
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
	bool verify_oriented() {
		for (Vertex_handle v = m_vertices.begin(); v != m_vertices.end(); ++v) {
			if (v->degree() != 2) continue; // irrelevant for orientation

			auto bwd = v->m_incident[0];
			auto fwd = v->m_incident[1];

			if (bwd->m_target == v && fwd->m_source == v) continue;

			return false;
		}

		return true;
	}
    bool is_oriented() {
        assert(!m_oriented || verify_oriented());
        return m_oriented;
    }
    void sort_incident_edges() {
        for (Vertex_handle& v : m_vertices) {
            if (v->degree() > 2) {
                std::ranges::sort(v->incident, [&v](Edge* e, Edge* f) {
                    CGAL::Direction_2<Kernel> dir_e = CGAL::Direction_2<Kernel>(e->other(v)->getPoint() - v->getPoint());
                    CGAL::Direction_2<Kernel> dir_f = CGAL::Direction_2<Kernel>(f->other(v)->getPoint() - v->getPoint());
                    return dir_e < dir_f;
                });
            }
        }
        m_sorted = true;
    }
    bool verify_sorted() {
        for (const Vertex_const_handle& v : m_vertices) {
            if (v->degree() > 2) {
                CGAL::Direction_2<Kernel> dir_prev = CGAL::Direction_2<Kernel>(v->neighbor(0)->getPoint() - v->getPoint());
                for (int i = 1; i < v->degree(); i++) {
                    CGAL::Direction_2<Kernel> dir = CGAL::Direction_2<Kernel>(v->neighbor(i)->getPoint() - v->getPoint());
                    if (dir < dir_prev) {
                        return false;
                    }
                    dir_prev = dir;
                }
            }
        }
        return true;
    }
    bool is_sorted() {
        assert(!m_sorted || verify_sorted());
        return m_sorted;
    }
};

template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits>
class Graph_2_vertex {
	friend class Graph_2<VertexData, EdgeData, CurveTraits>;
	friend class Graph_2_edge<VertexData, EdgeData, CurveTraits>;

  public:
	using Vertex = Graph_2_vertex<VertexData, EdgeData, CurveTraits>;
	using Edge = Graph_2_edge<VertexData, EdgeData, CurveTraits>;
	using Graph = Graph_2<VertexData, EdgeData, CurveTraits>;
	using Vertex_handle = Graph::Vertex_handle;
    using Vertex_const_handle = Graph::Vertex_const_handle;
	using Edge_handle = Graph::Edge_handle;
    using Edge_const_handle = Graph::Edge_const_handle;
	using Point_2 = CurveTraits::Point_2;
    using Vertex_data = VertexData;
    using Edge_data = EdgeData;
    using Curve_traits = CurveTraits;

  private:
	using Edge_container = std::vector<Edge_handle>;
	Point_2 m_point;
	Edge_container m_incident = Edge_container();

  public:
	const Point_2& point() const {
		return m_point;
	}
    Point_2& point() {
        return m_point;
    }
	int degree() const {
		return m_incident.size();
	}
	Edge_handle incoming() {
		assert(degree() == 2);
		return m_incident[0];
	}
    Edge_const_handle incoming() const {
        assert(degree() == 2);
        return m_incident[0];
    }
	Edge_handle outgoing() {
		assert(degree() == 2);
		return m_incident[1];
	}
    Edge_const_handle outgoing() const {
        assert(degree() == 2);
        return m_incident[1];
    }
	Vertex_handle prev() {
		return incoming()->source();
	}
    Vertex_const_handle prev() const {
        return incoming()->source();
    }
	Vertex_handle next() {
		return outgoing()->target();
	}
    Vertex_const_handle next() const {
        return outgoing()->target();
    }
    Edge_container::iterator incident_edges_begin() {
        return m_incident.begin();
    }
    Edge_container::iterator incident_edges_end() {
        return m_incident.end();
    }
    Edge_container::const_iterator incident_edges_begin() const {
        return m_incident.cbegin();
    }
    Edge_container::const_iterator incident_edges_end() const {
        return m_incident.cend();
    }
};

template <class VertexData, class EdgeData, GraphCurveTraits_2 CurveTraits>
class Graph_2_edge {
	friend class Graph_2<VertexData, EdgeData, CurveTraits>;
	friend class Graph_2_vertex<VertexData, EdgeData, CurveTraits>;

  public:
	using Vertex = Graph_2_vertex<VertexData, EdgeData, CurveTraits>;
	using Edge = Graph_2_edge<VertexData, EdgeData, CurveTraits>;
	using Graph = Graph_2<VertexData, EdgeData, CurveTraits>;
	using Vertex_handle = Graph::Vertex_handle;
	using Vertex_const_handle = Graph::Vertex_const_handle;
	using Edge_handle = Graph::Edge_handle;
	using Edge_const_handle = Graph::Edge_const_handle;
	using Point_2 = CurveTraits::Point_2;
	using Curve_2 = CurveTraits::Curve_2;
    using Vertex_data = VertexData;
    using Edge_data = EdgeData;
    using Curve_traits = CurveTraits;

  private:
	Vertex_handle m_source;
	Vertex_handle m_target;
	Curve_2 m_curve;
	EdgeData m_d;

  public:
	Graph_2_edge(Vertex_handle source, Vertex_handle target, Curve_2 curve)
		: m_source(std::move(source)), m_target(std::move(target)), m_curve(std::move(curve)) {
	    assert(m_curve.source() == m_source->m_point);
		assert(m_curve.target() == m_target->m_point);
	};

	Graph_2_edge(Vertex_handle source, Vertex_handle target, Curve_2 curve, Edge_data d)
	    : m_source(std::move(source)), m_target(std::move(target)), m_curve(std::move(curve)), m_d(std::move(d)) {};

	Vertex_handle source() { return m_source; }
	Vertex_handle target() { return m_target; }
	Vertex_const_handle source() const { return m_source; }
	Vertex_const_handle target() const { return m_target; }
    Curve_2& curve() { return m_curve; }
	const Curve_2& curve() const { return m_curve; }
	void reverse() {
		std::swap(m_source, m_target);
		m_curve = Curve_traits::reversed(m_curve);
	}
	Edge_handle prev() {
		return m_source->incoming();
	}
    Edge_const_handle prev() const {
        return m_source->incoming();
    }
	Edge_handle next() {
		return m_target->outgoing();
	}
    Edge_const_handle next() const {
        return m_target->outgoing();
    }
	Edge_data& data() {
		return m_d;
	}
    const Edge_data& data() const {
        return m_d;
    }
};
}
