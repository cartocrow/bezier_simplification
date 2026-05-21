#pragma once

#include <cartocrow/core/centroid.h>
#include <cartocrow/core/cubic_bezier.h>
#include <cartocrow/core/rectangle_helpers.h>
#include <cartocrow/data_structures/straight_graph_2.h>
#include <filesystem>
#include "library/core/region_set.h"
#include "library/core/topo_set.h"
#include "library/core/geometry_set.h"

#include <cartocrow/renderer/ipe_renderer.h>

#include <ogrsf_frmts.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>

namespace cartocrow::curved_simplification {
namespace read_graph_gdal_detail {
template <class K>
using Vb = CGAL::Triangulation_vertex_base_2<K>;
template <class K>
using Fb = CGAL::Constrained_triangulation_face_base_2<K>;
template <class K>
using TDS = CGAL::Triangulation_data_structure_2<Vb<K>, Fb<K>>;
typedef CGAL::Exact_predicates_tag                               Itag;
template <class K>
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TDS<K>, Itag>;
}

Straight_graph_2<std::monostate, std::monostate, Inexact> readGraphUsingGDAL(const std::filesystem::path& path);
std::pair<RegionSet<Inexact>, OGRSpatialReference> readRegionSetUsingGDAL(const std::filesystem::path& path);
void exportTopoSetUsingGDAL(const std::filesystem::path& path, const StraightTopoSet<Inexact> topoSet, const CGAL::Aff_transformation_2<Inexact>& trans, std::optional<OGRSpatialReference> spatialReference, bool stackPolygons);

GeometrySet<Inexact> readGeometrySetUsingGDAL(const std::filesystem::path& path);

// Returns whether the point lies on the boundary of the rectangle (under a given error threshold).
template <class K>
bool liesOn(const Point<K>& p, const Rectangle<K>& rect, const Number<K>& threshold = 0) {
	if (p.x() < rect.xmin() - threshold) return false;
	if (p.x() > rect.xmax() + threshold) return false;
	if (p.y() < rect.ymin() - threshold) return false;
	if (p.y() > rect.ymax() + threshold) return false;
	return abs(p.x() - rect.xmin()) < threshold || abs(p.x() - rect.xmax()) < threshold || abs(p.y() - rect.ymin()) < threshold || abs(p.y() - rect.ymax()) < threshold;
}

template <class K>
Point<K> pointInsidePolygon(const Polygon<K>& p) {
	using namespace read_graph_gdal_detail;
	CDT<K> cdt;
	cdt.insert_constraint(p.vertices_begin(), p.vertices_end(), true);
	using Face_handle = typename CDT<K>::Face_handle;
	std::unordered_map<Face_handle, bool> in_domain_map;
	boost::associative_property_map< std::unordered_map<Face_handle, bool> >
		in_domain(in_domain_map);

	//Mark facets that are inside the domain bounded by the polygon
	CGAL::mark_domain_in_triangulation(cdt, in_domain);

	for (Face_handle f : cdt.finite_face_handles())
	{
		if (get(in_domain, f)) {
			return CGAL::centroid(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
		}
	}
}

template <class PolygonTraits, class InputIterator, class OutputIterator, class K>
void flattenPolygonsInBbox(InputIterator begin, InputIterator end, OutputIterator out, const Rectangle<K>& bbox, const Number<K>& threshold = 0) {
	using PolygonLike = PolygonTraits::PolygonLike;
	std::vector<PolygonLike> bordering;
	std::vector<PolygonLike> within;
	for (auto it = begin; it != end; ++it) {
		bool borderingBbox = false;
		const auto& p = PolygonTraits::polygon(*it);
		for (auto vit = p.vertices_begin(); vit != p.vertices_end(); ++vit) {
			if (liesOn(*vit, bbox, threshold)) {
				borderingBbox = true;
				break;
			}
		}
		if (borderingBbox) {
			bordering.push_back(*it);
		}
		else {
			within.push_back(*it);
		}
	}

	// Sort on decreasing area.
	std::sort(bordering.begin(), bordering.end(), [](const PolygonLike& p1, const PolygonLike& p2) {
		return PolygonTraits::polygon(p1).area() > PolygonTraits::polygon(p2).area();
	});

	std::vector<PolygonLike> newPolygons;

	for (int bpI = 0; bpI < bordering.size(); ++bpI) {
		const Polygon<K>& bp = PolygonTraits::polygon(bordering[bpI]);
		auto circ = bp.vertices_circulator();
		
		std::deque<typename Polygon<K>::Vertex_circulator> d3s;
		
		auto curr = circ;
		
		auto initialPrev = curr;
		--initialPrev;
		bool prevLiesOnBbox = liesOn(*initialPrev, bbox, threshold);

		do {
			bool currLiesOnBbox = liesOn(*curr, bbox, threshold);
			auto next = curr;
			++next;
			bool nextLiesOnBbox = liesOn(*next, bbox, threshold);
			if (currLiesOnBbox && (!prevLiesOnBbox || !nextLiesOnBbox)) {
				d3s.push_back(curr);
			}
			prevLiesOnBbox = currLiesOnBbox;
		} while (++curr != circ);
		
		assert(d3s.size() % 2 == 0);

		// Ensure that the d3 vertices are ordered such that arc from i to i+1 is the one that may need to change.
		auto d3sF = d3s.front();
		auto d3sFN = d3sF;
		++d3sFN;
		if (liesOn(*d3sFN, bbox, threshold)) {
			d3s.pop_front();
			d3s.push_back(d3sF);
		}

		Polygon<K> newBp;

		for (int d3Index = 0; d3Index < d3s.size(); ++(++d3Index)) {
			// Replace non-bbox path by a bbox path.
			auto from = d3s[d3Index];
			auto to = d3s[d3Index + 1];

			//std::vector<Point<K>> border(from, to);
			std::vector<Point<K>> border;
			for (auto it = from; it != to + 1; ++it) {
				border.push_back(*it);
			}
			Polyline<K> newBorder;

			//newBp.push_back(*from);
			newBorder.push_back(*from);

			auto fromSide = closest_side(*from, bbox);
			auto toSide = closest_side(*to, bbox);
			
			auto currSide = fromSide;
			while (currSide != toSide) {
				auto nextSide = next_side(currSide);
				auto corner = get_corner(bbox, currSide, nextSide);
				//newBp.push_back(corner);
				newBorder.push_back(corner);
				currSide = nextSide;
			}

			// Vertex to is already added by the next sequence of operations
			//newBp.push_back(*to);
			newBorder.push_back(*to);

			// Add new border to newBp if it wouldn't cover already handled polygons.
			Polygon<K> coverPolygon(border.rbegin() + 1, border.rend() - 1);
			std::copy(newBorder.vertices_begin(), newBorder.vertices_end(), std::back_inserter(coverPolygon));

			//using namespace renderer;
			//IpeRenderer renderer;
			//renderer.addPainting([coverPolygon, newBorder, border, bp](GeometryRenderer& r) {
			//	r.setStroke(Color(0, 0, 0), 2.0);
			//	r.draw(bp);
			//	r.setStroke(Color(255, 0, 0), 2.0);
			//	Polyline<Inexact> borderPl(border.begin(), border.end());
			//	r.draw(borderPl);
			//	r.setStroke(Color(0, 255, 0), 2.0);
			//	r.draw(newBorder);
			//	r.setStroke(Color(0, 0, 255), 2.0);
			//	r.draw(coverPolygon);
			//});
			//std::stringstream ss;
			//ss << "debugging_fun_" << bpI << "_" << d3Index << ".ipe";
			//renderer.save(ss.str());
			
			bool covers = false;
			for (int handledBefore = 0; handledBefore < bpI; ++handledBefore) {
				const Polygon<K>& before = PolygonTraits::polygon(bordering.at(handledBefore));
				if (coverPolygon.has_on_bounded_side(pointInsidePolygon(before))) {
					covers = true;
					break;
				}
			}
			if (covers) {
				// Keep old border
				std::copy(border.begin(), border.end(), std::back_inserter(newBp));
			}
			else {
				std::copy(newBorder.vertices_begin(), newBorder.vertices_end(), std::back_inserter(newBp));
			}

			// Keep bbox path
			from = to;
			to = d3s[(d3Index + 2) % d3s.size()];
			auto keepCirc = from;
			do {
				newBp.push_back(*keepCirc);
			} while (++keepCirc != to);
		}

		newPolygons.push_back(PolygonTraits::changePolygon(bordering[bpI], newBp));
	}

	std::copy(within.begin(), within.end(), std::back_inserter(newPolygons));

	// Sort final polygons on decreasing area.
	std::sort(newPolygons.begin(), newPolygons.end(), [](const auto& p1, const auto& p2) {
		return PolygonTraits::polygon(p1).area() > PolygonTraits::polygon(p2).area();
	});

	for (const auto& pgn : newPolygons) {
		*out++ = pgn;
	}
}
}
