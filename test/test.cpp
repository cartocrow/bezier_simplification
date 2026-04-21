#include <catch2/catch_test_macros.hpp>
#include "library/steven_bezier_collapse.h"
#include "library/core/straight_graph_2.h"
#include "library/read_ipe_bezier_spline.h"

using namespace cartocrow;
using namespace cartocrow::curved_simplification;

TEST_CASE("getAreaD") {
    CubicBezierCurve curve({0, 0}, {3, 1}, {2, 6}, {5, 6});
    auto a = curve.signedArea();
    auto p0 = curve.source();
    auto p3 = curve.target();
    auto t0 = curve.sourceControl() - curve.source();
    auto t1 = curve.targetControl() - curve.target();
    auto d0 = sqrt(t0.squared_length());
    auto d1 = sqrt(t1.squared_length());

    CHECK(abs(getAreaD0(a, p0, t0.direction(), t1.direction(), d1, p3) - d0) < M_EPSILON);
    CHECK(abs(getAreaD1(a, p0, t0.direction(), d0, t1.direction(), p3) - d1) < M_EPSILON);
}

using BezierGraph = Bezier_graph_2<std::monostate, std::monostate>;

BezierGraph* readGraphFromIpe(std::filesystem::path path) {
    BezierGraph* g = new BezierGraph();

    std::unordered_map<Point<Inexact>, BezierGraph::Vertex_handle> pToV;
    auto getVertex = [&pToV, &g](const Point<Inexact> &p) {
        if (pToV.contains(p)) {
            return pToV.at(p);
        } else {
            auto newV = g->add_vertex(p);
            pToV[p] = newV;
            return newV;
        }
    };

    if (path.extension() == ".ipe") {
        auto splines = ipeSplinesToIsolines(path);

        for (const auto &spline: splines) {
            if (spline.empty()) continue;
            auto lastV = getVertex(spline.curves()[0].source());
            for (const auto &curve: spline.curves()) {
                auto targetV = getVertex(curve.target());
                g->add_edge(lastV, targetV, curve);
                lastV = targetV;
            }
        }

    }

    return g;
}

void checkGraph(BezierGraph& g) {
    double xSum = 0;
    for (auto vit = g.vertices_begin(); vit != g.vertices_end(); ++vit) {
        xSum += vit->point().x();
    }
    for (auto eit = g.edges_begin(); eit != g.edges_end(); ++eit) {
        if (eit->source()->degree() == 2) {
            CHECK(eit->source()->outgoing() == eit);
        }
        if (eit->target()->degree() == 2) {
            CHECK(eit->target()->incoming() == eit);
        }
        xSum += eit->curve().sourceControl().x();
        if (eit->target()->degree() == 2) {
            xSum += eit->next()->curve().sourceControl().x();
        }
    }
}

TEST_CASE("Read a graph") {
    readGraphFromIpe("splines.ipe");
}

TEST_CASE("Assign a graph") {
    auto* g = readGraphFromIpe("splines.ipe");
    BezierGraph assignmentCopy;
    assignmentCopy = *g;

    // We aim to verify that the copy doesn't have pointers to g
    delete g;

    checkGraph(assignmentCopy);
}

TEST_CASE("Copy a graph") {
    auto* g = readGraphFromIpe("splines.ipe");
    auto copy = *g;

    // We aim to verify that the copy doesn't have pointers to g
    delete g;

    checkGraph(copy);
}