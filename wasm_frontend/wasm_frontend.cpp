#include <iostream>

#include "library/bezier_collapse.h"
#include "library/steven_bezier_collapse.h"

#include "cartocrow/renderer/svg_renderer.h"

#include <emscripten/bind.h>
#include <emscripten/val.h>

using namespace emscripten;

using namespace cartocrow;
using namespace cartocrow::curved_simplification;

class BezierSimplification {
private:
    using Graph = BezierCollapseGraphWithHistory;
    using BaseGraph = Graph::BaseGraph;
    using Traits = StevenBCTraits<Graph>;
    using Collapse = BezierCollapse<Graph, Traits>;
    BaseGraph m_bg;
    Graph m_g;
    Collapse m_collapse;
    int m_maxVid = 0;
    int m_maxEid = 0;

    std::unordered_map<const BaseGraph::Vertex*, int> m_vToVid;
    std::unordered_map<int, BaseGraph::Vertex_handle> m_vidToV;
    std::unordered_map<const BaseGraph::Edge*, int> m_eToEid;
    std::unordered_map<int, BaseGraph::Edge_handle> m_eidToE;


public:
    BezierSimplification() : m_g(m_bg), m_collapse(m_g, Traits()) {};

    int insert_vertex(double x, double y) {
        auto vh = m_bg.insert_vertex({x, y});
        int id = m_maxVid++;
        m_vToVid[&*vh] = id;
        m_vidToV[id] = vh;
        return id;
    }

    int add_edge(int id1, int id2, double c0x, double c0y, double c1x, double c1y) {
        auto vh1 = m_vidToV[id1];
        auto vh2 = m_vidToV[id2];

        CubicBezierCurve curve(vh1->point(), {c0x, c0y}, {c1x, c1y}, vh2->point());
        auto eh = m_bg.add_edge(vh1, vh2, curve);

        int id = m_maxEid++;
        m_eToEid[&*eh] = id;
        m_eidToE[id] = eh;

        return id;
    }

    val vertices() const {
        val arr = val::array();

        int i = 0;
        for (auto vit = m_bg.vertices_begin(); vit != m_bg.vertices_end(); ++vit) {
            arr.set(i++, m_vToVid.at(&*vit));
        }
        return arr;
    }

    val edges() const {
        val arr = val::array();

        int i = 0;
        for (auto eit = m_bg.edges_begin(); eit != m_bg.edges_end(); ++eit) {
            if (!m_eToEid.contains(&*eit)) continue;
            arr.set(i++, m_eToEid.at(&*eit));
        }
        return arr;
    }

    int number_of_vertices() const {
        return m_bg.number_of_vertices();
    }

    int number_of_edges() const {
        return m_bg.number_of_edges();
    }
private:
    val makePoint(const Point<Inexact>& pt) const {
        val obj = val::object();
        obj.set("x", pt.x());
        obj.set("y", pt.y());
        return obj;
    }

    void recomputeMaps() {
        m_vidToV.clear();
        m_eidToE.clear();
        m_vToVid.clear();
        m_eToEid.clear();
        m_maxVid = 0;
        m_maxEid = 0;

        int i = 0;
        for (auto vit = m_bg.vertices_begin(); vit != m_bg.vertices_end(); ++vit) {
            m_vidToV[i] = vit;
            m_vToVid[&*vit] = i++;
        }
        i = 0;
        for (auto eit = m_bg.edges_begin(); eit != m_bg.edges_end(); ++eit) {
            m_eidToE[i] = eit;
            m_eToEid[&*eit] = i++;
        }
    }
public:
    val get_vertex_point(int id) const {
        auto vh = m_vidToV.at(id);
        return makePoint(vh->point());
    }

    val get_edge_curve(int id) const {
        auto eh = m_eidToE.at(id);
        auto& c = eh->curve();
        val obj = val::object();
        obj.set("c0", makePoint(c.source()));
        obj.set("c1", makePoint(c.sourceControl()));
        obj.set("c2", makePoint(c.targetControl()));
        obj.set("c3", makePoint(c.target()));

        return obj;
    }

    val get_bbox() const {
        auto bb = m_bg.bbox();

        val obj = val::object();
        obj.set("xmin", bb.xmin());
        obj.set("ymin", bb.ymin());
        obj.set("xmax", bb.xmax());
        obj.set("ymax", bb.ymax());

        return obj;
    }

    void initialize() {
        m_collapse.initialize();
        m_collapse.runToComplexity(1);
        recomputeMaps();
    }

    void run_to_complexity(int k) {
        m_g.recallComplexity(k);
        recomputeMaps();
    }
};

EMSCRIPTEN_BINDINGS(bezier_simplification) {
    class_<BezierSimplification>("BezierSimplification")
            .constructor<>()
            .function("insert_vertex", &BezierSimplification::insert_vertex)
            .function("add_edge", &BezierSimplification::add_edge)
            .function("get_vertex_point", &BezierSimplification::get_vertex_point)
            .function("get_edge_curve", &BezierSimplification::get_edge_curve)
            .function("vertices", &BezierSimplification::vertices)
            .function("edges", &BezierSimplification::edges)
            .function("get_bbox", &BezierSimplification::get_bbox)
            .function("number_of_vertices", &BezierSimplification::number_of_vertices)
            .function("number_of_edges", &BezierSimplification::number_of_edges)
            .function("initialize", &BezierSimplification::initialize)
            .function("run_to_complexity", &BezierSimplification::run_to_complexity);
}