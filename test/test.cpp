#include <catch2/catch_test_macros.hpp>
#include "library/steven_bezier_collapse.h"
#include <cartocrow/data_structures/straight_graph_2.h>
#include <cartocrow/reader/ipe_reader.h>
#include "library/read_ipe_bezier_spline.h"
#include "library/read_graph_gdal.h"
#include "library/sym_diff.h"
#include "library/vertex_snap.h"
#include <random>
#include <ipeshape.h>
#include <ipepath.h>

#include <cartocrow/renderer/ipe_renderer.h>

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

std::vector<PolygonSetRaw<Inexact>> polygonsInPage(ipe::Page* page) {
    auto pgns = std::vector<PolygonSetRaw<Inexact>>();

    for (int i = 0; i < page->count(); i++) {
        auto object = page->object(i);
        if (object->type() != ipe::Object::Type::EPath) continue;
        auto path = object->asPath();
        auto matrix = object->matrix();
        auto shape = path->shape();
        pgns.push_back(IpeReader::convertShapeToPolygonSetRaw(shape, matrix));
    }

    return pgns;
}

RegionSet<Inexact> readRegionSetFromIpe(std::filesystem::path path) {
    RegionSet<Inexact> rs;

    std::shared_ptr<ipe::Document> document = IpeReader::loadIpeFile(path);

    if (document->countPages() == 0) {
        throw std::runtime_error("Cannot read isolines from an Ipe file with no pages");
    }
    else if (document->countPages() > 1) {
        throw std::runtime_error("Cannot read isolines from an Ipe file with more than one page");
    }

    ipe::Page* page = document->page(0);

    auto pgns = polygonsInPage(page);
    for (const auto& pgn : pgns) {
        rs.regions.emplace_back(RegionAttributes{}, pgn);
    }

    return rs;
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

TEST_CASE("Symmetric difference of two Bézier splines with common endpoints") {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(0, 1);

    for (int i = 0; i < 10; ++i) {
        CubicBezierSpline spline1;
        spline1.appendCurve({ 0, 0 }, { dist6(rng), dist6(rng) }, { dist6(rng), dist6(rng) }, { 1.5, 1.5 });
        CubicBezierSpline spline2;
        spline2.appendCurve({ 0, 0 }, { dist6(rng), dist6(rng) }, { dist6(rng), dist6(rng) }, { 1.5, 1.5 });

        CHECK(abs(symmetricDifference(spline1, spline2) - symmetricDifferenceArrangement(spline1, spline2, 1000)) < 0.01);
    }
}

TEST_CASE("Remove duplicates") {
    Polygon<Inexact> p;
    p.push_back({ 0, 0 });
    p.push_back({ 1, 0 });
    p.push_back({ 1, 0 });
    p.push_back({ 1, 1 });
    p.push_back({ 0, 1 });
    
    removeDuplicates(p);

    CHECK(p.size() == 4);
    CHECK(p.vertex(0) == Point<Inexact>{0, 0});
    CHECK(p.vertex(1) == Point<Inexact>{1, 0});
    CHECK(p.vertex(2) == Point<Inexact>{1, 1});
    CHECK(p.vertex(3) == Point<Inexact>{0, 1});
}

TEST_CASE("Flatten polygons in bbox") {
    auto rs = readRegionSetFromIpe("data/test/flat_stacking/complex.ipe");
    auto bbox = rs.bbox();
    Rectangle<Inexact> rect(bbox.xmin(), bbox.ymin(), bbox.xmax(), bbox.ymax());
    std::vector<Polygon<Inexact>> polygons;
    for (const auto& r : rs.regions) {
        polygons.push_back(r.geometry.polygons_with_holes.front().outer_boundary());
        if (polygons.back().is_clockwise_oriented()) {
            polygons.back().reverse_orientation();
        }
    }
    std::vector<Polygon<Inexact>> flattened;
    flattenPolygonsInBbox(polygons.begin(), polygons.end(), std::back_inserter(flattened), rect, M_EPSILON);

    using namespace renderer;
    IpeRenderer renderer;

    renderer.addPainting([&flattened](GeometryRenderer& r) {
        for (const auto& flat : flattened) {
            r.draw(flat);
        }
    });
    renderer.save("debugging_flattening.ipe");
}