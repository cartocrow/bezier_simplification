#include "topojson_export.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace CGAL {
void to_json(json& j, const CGAL::Point_2<cartocrow::Inexact>& p) {
    j = json::array({ p.x(), p.y() });
}

void from_json(const json& j, cartocrow::Point<cartocrow::Inexact>& p) {
    p = { j.at(0).get<double>(), j.at(1).get<double>() };
}
}

namespace cartocrow {
void to_json(json& j, const CubicBezierSpline& s) {
    j = s.controlPoints();
}

void from_json(const json& j, CubicBezierSpline& s) {
	auto points = j.get<std::vector<Point<Inexact>>>();
	s = CubicBezierSpline{ points.begin(), points.end() };
}


}
//void to_json(nlohmann::json& j, const RegionAttributes& attrs)
//{
//    j = nlohmann::json::object();
//
//
//}

//void from_json(const nlohmann::json& j, RegionAttributes& attrs)
//{
//    attrs.clear();
//
//    for (auto it = j.begin(); it != j.end(); ++it) {
//        attrs[it.key()] = it.value().get<RegionAttribute>();
//    }
//}

namespace nlohmann {
using namespace cartocrow;


template<>
struct adl_serializer<RegionAttribute> {
static void to_json(json& j, const RegionAttribute& attr)
{
    std::visit([&j](const auto& value) {
        j = value;
        }, attr);
}

static void from_json(const json& j, RegionAttribute& attr)
{
    if (j.is_number_integer()) {
        attr = j.get<int64_t>();
    }
    else if (j.is_number_float()) {
        attr = j.get<double>();
    }
    else if (j.is_string()) {
        attr = j.get<std::string>();
    }
    else if (j.is_array()) {
        if (j.empty()) {
            attr = std::vector<int>{};
            return;
        }
        if (j[0].is_number_integer()) {
            attr = j.get<std::vector<int>>();
        }
        else if (j[0].is_number_float()) {
            attr = j.get<std::vector<double>>();
        }
        else if (j[0].is_string()) {
            attr = j.get<std::vector<std::string>>();
        }
    }
}
};

template<>
struct adl_serializer<BezierTopoSet::Feature> {
    static void to_json(json& j, const BezierTopoSet::Feature& f) {
        std::string type;
        json arcs;
        if (auto* plTP = std::get_if<BezierTopoSet::PolylineTopology>(&f.topology)) {
            type = "LineString";
            arcs = plTP->arcs;
        }
        else if (auto* pTP = std::get_if<BezierTopoSet::PolygonTopology>(&f.topology)) {
            type = "Polygon";
            arcs = pTP->arcs;
        }
        else if (auto* psTP = std::get_if<BezierTopoSet::PolygonSetTopology>(&f.topology)) {
            type = "MultiPolygon";
            arcs = psTP->arcs;
        }
        json attributes = f.attributes;
        j = { {"type", type}, {"properties", attributes}, {"arcs", arcs} };
    }

    static void from_json(const json& j, BezierTopoSet::Feature& f)
    {
        const auto& type = j.at("type");
        const auto& arcs = j.at("arcs");

        f.attributes =
            j.at("properties").get<RegionAttributes>();

        std::string t = type.get<std::string>();

        if (t == "LineString") {

            f.topology = BezierTopoSet::PolylineTopology{
                arcs.get<std::vector<int>>()
            };
        }
        else if (t == "Polygon") {

            f.topology = BezierTopoSet::PolygonTopology{
                arcs.get<std::vector<std::vector<int>>>()
            };
        }
        else if (t == "MultiPolygon") {

            f.topology = BezierTopoSet::PolygonSetTopology{
                arcs.get<std::vector<std::vector<std::vector<int>>>>()
            };
        }
        else {
            throw std::runtime_error(
                "Unknown topology type: " + t
            );
        }
    }
};
}

namespace cartocrow {
void exportTopoSetToJson(const std::filesystem::path& path, const BezierTopoSet& topoSet, std::optional<OGRSpatialReference> spatialReference) {
	json j;
	j["type"] = "Topology with cubic Bézier splines as arcs";
    json geometryCollection = {
        {"type", "GeometryCollection"},
        {"geometries", topoSet.features}
    };
    j["objects"] = { { path.stem(), geometryCollection } };
	j["arcs"] = topoSet.arcs;
    if (spatialReference.has_value()) {
        j["spatialReference"] = spatialReference->exportToWkt();
    }

    // Create parent directories if they don't exist
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream o(path);
    if (!o) {
        std::cerr << "Failed to open file for writing: " << path << '\n';
        return;
    }
    o << j << std::endl;
    if (!o.good()) {
        std::cerr << "Failed to write JSON to file: " << path << '\n';
    }
}

std::pair<BezierTopoSet, OGRSpatialReference> parseBezierTopoJson(const std::filesystem::path& path) {
    std::ifstream f(path);
    json j = json::parse(f);

    OGRSpatialReference spatialReference;
    if (j.contains("spatialReference")) { 
        std::string wkt = j["spatialReference"];
        spatialReference.importFromWkt(wkt.c_str());
    }
    std::vector<CubicBezierSpline> arcs = j["arcs"];

    auto& objects = j["objects"];
    // assume there is exactly 1.

    auto it = objects.begin();

    if (it == objects.end()) {
        throw std::runtime_error("No objects");
    }

    auto& geometries = it.value()["geometries"];
    std::vector<BezierTopoSet::Feature> features = geometries;

    BezierTopoSet ts;
    ts.arcs = arcs;
    ts.features = features;

    return { ts, spatialReference };

    //auto p2 = j.get<ns::person>();
}
}