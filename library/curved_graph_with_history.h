#pragma once

#include <cartocrow/core/core.h>
#include "core/graph_2_like.h"

namespace cartocrow::curved_simplification {
namespace detail {
template <Graph2Like Graph>
struct Operation {
    using Edge_handle = Graph::Edge_handle;
    Edge_handle m_edge;

    Operation(Edge_handle e) : m_edge(e) {};

    virtual void undo(Graph& g) = 0;
    virtual void redo(Graph& g) = 0;
};

/// Note: this operation has not been tested.
template <Graph2Like Graph>
struct SplitVertexOperation : public Operation<Graph> {
  public:
    using Edge_handle = Graph::Edge_handle;
    using Vertex_handle = Graph::Vertex_handle;
    using Op = Operation<Graph>;
    using Curve = Graph::Curve_2;
    using ED = Graph::ED;
    using Operation<Graph>::m_edge;

    Curve m_incomingCurve;
    Curve m_outgoingCurve;
    Curve m_c0;
    Curve m_c1;
    Curve m_c2;

    ED m_incomingData;
    ED m_outgoingData;
    ED m_e0Data;
    ED m_e1Data;
    ED m_e2Data;

  public:
    SplitVertexOperation(Vertex_handle v, Curve c0, Curve c1, Curve c2) :
        Operation<Graph>(v->incoming()), m_c0(std::move(c0)), m_c1(std::move(c1)), m_c2(std::move(c2)),
        m_incomingCurve(v->incoming()->curve()), m_outgoingCurve(v->outgoing()->curve()) {
        m_incomingData = v->incoming()->data();
        m_outgoingData = v->outgoing()->data();
    }
    void undo(Graph& g) {
        auto e0 = m_edge->prev();
        auto& e1 = m_edge;
        auto e2 = m_edge->next();
        m_e0Data = e0->data();
        m_e1Data = e1->data();
        m_e2Data = e2->data();

        auto thisOp = e1->data().hist;

        auto v = g.collapse_edge(m_edge, m_incomingCurve, m_outgoingCurve);
        auto inc = v->incoming();
        auto out = v->outgoing();
        m_edge = inc;

        inc->data() = m_incomingData;
        out->data() = m_outgoingData;

        m_edge->data().futr = thisOp;

        if (m_incomingData.hist != nullptr) {
            m_incomingData.hist->m_edge = inc;
        }
        if (m_incomingData.futr != nullptr) {
            m_incomingData.futr->m_edge = inc;
        }
        if (m_outgoingData.hist != nullptr) {
            m_outgoingData.hist->m_edge = out;
        }
        if (m_outgoingData.futr != nullptr) {
            m_outgoingData.futr->m_edge = out;
        }
        return v;
    }
    Edge_handle perform(Graph& g) {
        m_edge = g.split_vertex(m_edge->target(), m_c0, m_c1, m_c2);

        auto e0 = m_edge->prev();
        auto& e1 = m_edge;
        auto e2 = m_edge->next();
        e0->data() = m_e0Data;
        e1->data() = m_e1Data;
        e2->data() = m_e2Data;

        if (m_e0Data.hist != nullptr) {
            m_e0Data.hist->m_edge = e0;
        }
        if (m_e0Data.futr != nullptr) {
            m_e0Data.futr->m_edge = e0;
        }
        if (m_e1Data.hist != nullptr) {
            m_e1Data.hist->m_edge = e1;
        }
        if (m_e1Data.futr != nullptr) {
            m_e1Data.futr->m_edge = e1;
        }
        if (m_e2Data.hist != nullptr) {
            m_e2Data.hist->m_edge = e2;
        }
        if (m_e2Data.futr != nullptr) {
            m_e2Data.futr->m_edge = e2;
        }

        return m_edge;
    }
    void redo(Graph& g) {
        perform(g);
    }
};
template <Graph2Like Graph>
struct CollapseEdgeOperation : public Operation<Graph> {
public:
    using Edge_handle = Graph::Edge_handle;
    using Vertex_handle = Graph::Vertex_handle;
    using Op = Operation<Graph>;
    using Curve = Graph::Curve_2;
    using ED = Graph::Edge_data;
    using Operation<Graph>::m_edge;

    Curve m_incomingCurve;
    Curve m_outgoingCurve;
    Curve m_c0;
    Curve m_c1;
    Curve m_c2;

    ED m_incomingData;
    ED m_outgoingData;
    ED m_e0Data;
    ED m_e1Data;
    ED m_e2Data;

public:
    CollapseEdgeOperation(Edge_handle e, Curve incomingCurve, Curve outgoingCurve) :
            Operation<Graph>(e), m_incomingCurve(incomingCurve), m_outgoingCurve(outgoingCurve),
            m_c0(e->prev()->curve()), m_c1(e->curve()), m_c2(e->next()->curve()) {
        auto& e1 = e;
        auto e0 = e->prev();
        auto e2 = e->next();
        m_e0Data = e0->data();
        m_e1Data = e1->data();
        m_e2Data = e2->data();
    }
    void undo(Graph& g) {
        auto inc = m_edge;
        auto out = m_edge->next();
        m_incomingData = inc->data();
        m_outgoingData = out->data();

        auto thisOp = inc->data().hist;
        m_edge = g.split_vertex(m_edge->target(), m_c0, m_c1, m_c2);

        auto e0 = m_edge->prev();
        auto& e1 = m_edge;
        auto e2 = m_edge->next();
        e0->data() = m_e0Data;
        e1->data() = m_e1Data;
        e2->data() = m_e2Data;

        m_edge->data().futr = thisOp;

        if (m_e0Data.hist != nullptr) {
            m_e0Data.hist->m_edge = e0;
        }
        if (m_e0Data.futr != nullptr) {
            m_e0Data.futr->m_edge = e0;
        }
        if (m_e1Data.hist != nullptr) {
            m_e1Data.hist->m_edge = e1;
        }
        if (m_e1Data.futr != nullptr) {
            m_e1Data.futr->m_edge = e1;
        }
        if (m_e2Data.hist != nullptr) {
            m_e2Data.hist->m_edge = e2;
        }
        if (m_e2Data.futr != nullptr) {
            m_e2Data.futr->m_edge = e2;
        }
    }
    Vertex_handle perform(Graph& g) {
        auto v = g.collapse_edge(m_edge, m_incomingCurve, m_outgoingCurve);
        auto inc = v->incoming();
        auto out = v->outgoing();
        m_edge = inc;

        inc->data() = m_incomingData;
        out->data() = m_outgoingData;

        if (m_incomingData.hist != nullptr) {
            m_incomingData.hist->m_edge = inc;
        }
        if (m_incomingData.futr != nullptr) {
            m_incomingData.futr->m_edge = inc;
        }
        if (m_outgoingData.hist != nullptr) {
            m_outgoingData.hist->m_edge = out;
        }
        if (m_outgoingData.futr != nullptr) {
            m_outgoingData.futr->m_edge = out;
        }
        return v;
    }
    void redo(Graph& g) {
        perform(g);
    }
};
template <Graph2Like Graph> struct OperationBatch;
template <class Graph>
concept EdgeStoredOperations = requires(typename Graph::Edge_data d) {
    requires(Graph2Like<Graph>);
    {
    d.hist
    } -> std::same_as<Operation<Graph>*&>;
    {
    d.futr
    } -> std::same_as<Operation<Graph>*&>;
};
}

template <Graph2Like Graph>
class CollapseHistoryGraphAdaptor {
  public:
    using Curve_traits = Graph::Curve_traits;
    using Vertex = Graph::Vertex;
    using Vertex_data = Graph::Vertex_data;
    using Vertex_handle = Graph::Vertex_handle;
    using Vertex_const_handle = Graph::Vertex_const_handle;
    using Vertex_iterator = Graph::Vertex_iterator;
    using Vertex_const_iterator = Graph::Vertex_const_iterator;
    using Edge = Graph::Edge;
    using Edge_data = Graph::Edge_data;
    using Edge_handle = Graph::Edge_handle;
    using Edge_const_handle = Graph::Edge_const_handle;
    using Edge_iterator = Graph::Edge_iterator;
    using Edge_const_iterator = Graph::Edge_const_iterator;
    using Kernel = Graph::Kernel;
    using Point_2 = Graph::Point_2;
    using Curve_2 = Graph::Curve_2;

    using BaseGraph = Graph;

    Graph& m_graph;
    std::stack<std::shared_ptr<detail::Operation<Graph>>> m_history;
    std::stack<std::shared_ptr<detail::Operation<Graph>>> m_undone;

    CollapseHistoryGraphAdaptor(Graph& graph) : m_graph(graph) {};

    Vertex_iterator vertices_begin() { return m_graph.vertices_begin(); }
    Vertex_iterator vertices_end() { return m_graph.vertices_end(); }
    Edge_iterator edges_begin() { return m_graph.edges_begin(); }
    Edge_iterator edges_end() { return m_graph.edges_end(); }

    Vertex_const_iterator vertices_begin() const { return m_graph.vertices_begin(); }
    Vertex_const_iterator vertices_end() const { return m_graph.vertices_end(); }
    Edge_const_iterator edges_begin() const { return m_graph.edges_begin(); }
    Edge_const_iterator edges_end() const { return m_graph.edges_end(); }

    size_t number_of_edges() { return m_graph.number_of_edges(); }
    size_t number_of_vertices() { return m_graph.number_of_vertices(); }

    bool atPresent() { return m_undone.empty(); }
    void recallComplexity(int c) {
        while (!m_history.empty() && number_of_edges() < c) {
            backInTime();
        }
        while (!m_undone.empty() && number_of_edges() > c) {
            forwardInTime();
        }
    }
    void backInTime() {
        if (m_history.empty()) return;
        auto lastOperation = m_history.top();
        lastOperation->undo(m_graph);
        m_undone.push(std::move(lastOperation));
        m_history.pop();
    }
    void forwardInTime() {
        if (m_undone.empty()) return;
        auto lastOperation = m_undone.top();
        lastOperation->redo(m_graph);
        m_history.push(std::move(lastOperation));
        m_undone.pop();
    }
    void goToPresent() {
        while (!atPresent()) {
            forwardInTime();
        }
    }
    void reset() {
        m_history = {};
        m_undone = {};
    }

//    void startBatch(Number<K> c);
//    void endBatch();

    /// Split a vertex into two. Or equivalently, replace two curves by three new ones.
    /// Connect the two new vertices with three new curves.
    Edge_handle split_vertex(Vertex_handle v, Curve_2 c0, Curve_2 c1, Curve_2 c2) {
        auto operation = std::make_shared<detail::SplitVertexOperation<Graph>>(v, std::move(c0), std::move(c1), std::move(c2));
        auto e = operation->perform(m_graph);
        e->data().hist = operation;
        m_history.push(std::move(operation));
        m_undone = {};
        return e;
    }
    /// Collapse an edge, removing it together with the edge that comes before and after it.
    /// Replace it with a vertex connected by two new curves.
    /// This is equivalent to replacing three curves by two new ones.
    Vertex_handle collapse_edge(Edge_handle e, Curve_2 toNewPoint, Curve_2 fromNewPoint) {
        auto operation = std::make_shared<detail::CollapseEdgeOperation<Graph>>(e, std::move(toNewPoint), std::move(fromNewPoint));
        auto v = operation->perform(m_graph);
        v->incoming()->data().hist = operation;
        m_history.push(std::move(operation));
        m_undone = {};
        return v;
    }
};
static_assert(Graph2Like<CollapseHistoryGraphAdaptor<Bezier_graph_2<std::monostate, std::monostate>>>);
}
