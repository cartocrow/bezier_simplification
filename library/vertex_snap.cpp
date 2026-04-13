#include "vertex_snap.h"

namespace cartocrow::curved_simplification {
struct RegionSetVertex {
	RegionSet<Inexact>* regionSet;
	size_t regionIndex;
	size_t polygonIndex;
	std::optional<size_t> holeIndex;
	size_t vertexIndex;

    Point<Inexact>& point() const {
        auto& pgnWH = regionSet->regions[regionIndex].geometry.polygons_with_holes[polygonIndex];
        auto& pgn = holeIndex.has_value() ? pgnWH.holes()[*holeIndex] : pgnWH.outer_boundary();
        return pgn.vertex(vertexIndex);
    }
};

struct RegionSetVertexPQTTraits {
	using Element = RegionSetVertex;
	using Kernel = Inexact;

	static Point<Kernel>& get_point(Element& elt) {
		return elt.point();
	}
};

static_assert(datastructures::PointQuadTreeTraits<RegionSetVertexPQTTraits>);

void snapVertices(RegionSet<Inexact>& regionSet, std::optional<double> epsilon) {
	using namespace cartocrow::datastructures;
	auto bbox = regionSet.bbox();
	Rectangle<Inexact> box(bbox.xmin(), bbox.ymin(), bbox.xmax(), bbox.ymax());
	PointQuadTree<RegionSetVertexPQTTraits> pqt(box, 10);

    double avgLength = 0;
    int polygonCount = 0;
	
	std::vector<RegionSetVertex> rsVertices;
	for (int regionIndex = 0; regionIndex < regionSet.regions.size(); ++regionIndex) {
		auto& pgnWHs = regionSet.regions[regionIndex].geometry.polygons_with_holes;
		for (int polygonIndex = 0; polygonIndex < pgnWHs.size(); ++polygonIndex) {
			auto& pgnWH = pgnWHs[polygonIndex];

            double avgLengthOuter = 0;
			for (int vertexIndex = 0; vertexIndex < pgnWH.outer_boundary().size(); ++vertexIndex) {
				rsVertices.emplace_back(&regionSet, regionIndex, polygonIndex, std::nullopt, vertexIndex);
                avgLengthOuter += sqrt(CGAL::squared_distance(pgnWH.outer_boundary().vertex(vertexIndex), pgnWH.outer_boundary().vertex((vertexIndex + 1) % pgnWH.outer_boundary().size())));
			}
            avgLengthOuter /= pgnWH.outer_boundary().size();
            avgLength += avgLengthOuter;
			for (int holeIndex = 0; holeIndex < pgnWH.holes().size(); ++holeIndex) {
                double avgLengthHole = 0;
				auto& hole = pgnWH.holes()[holeIndex];
				for (int vertexIndex = 0; vertexIndex < hole.size(); ++vertexIndex) {
					rsVertices.emplace_back(&regionSet, regionIndex, polygonIndex, holeIndex, vertexIndex);
                    avgLengthHole += sqrt(CGAL::squared_distance(hole.vertex(vertexIndex), hole.vertex((vertexIndex + 1) % hole.size())));
				}
                avgLengthHole /= hole.size();
                avgLength += avgLengthHole;
                ++polygonCount;
			}
		}
	}
    avgLength /= polygonCount;

	for (auto& rsVertex : rsVertices) {
		pqt.insert(rsVertex);
	}

    for (auto& rsVertex : rsVertices) {
        Circle<Inexact> circle(rsVertex.point(), epsilon.has_value() ? *epsilon : avgLength * 0.0025);
        Box bbox = circle.bbox();
        Rectangle<Inexact> query(bbox.xmin(), bbox.ymin(), bbox.xmax(), bbox.ymax());

        std::vector<RegionSetVertex> candidates;
        pqt.findContained(query, [&rsVertex, &candidates](RegionSetVertex& vertex) {
            if (vertex.regionIndex == rsVertex.regionIndex) return;
            if (vertex.regionIndex != rsVertex.regionIndex ||
                vertex.holeIndex != rsVertex.holeIndex ||
                vertex.vertexIndex != rsVertex.vertexIndex ||
                vertex.polygonIndex != rsVertex.polygonIndex) {
                candidates.push_back(vertex);
            }
        });
        if (candidates.size() > 0) {
            auto closest = std::min_element(candidates.begin(), candidates.end(), [&rsVertex](const RegionSetVertex& v1, const RegionSetVertex& v2) {
                return CGAL::squared_distance(rsVertex.point(), v1.point()) < CGAL::squared_distance(rsVertex.point(), v2.point());
            });
            rsVertex.point() = closest->point();
        }
    }
}
}