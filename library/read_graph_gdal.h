#pragma once

#include <cartocrow/core/cubic_bezier.h>
#include "library/core/straight_graph_2.h"
#include <filesystem>
#include "library/core/region_set.h"
#include "library/core/topo_set.h"
#include "library/core/geometry_set.h"

#include <ogrsf_frmts.h>

namespace cartocrow::curved_simplification {
Straight_graph_2<std::monostate, std::monostate, Inexact> readGraphUsingGDAL(const std::filesystem::path& path);
std::pair<RegionSet<Inexact>, OGRSpatialReference> readRegionSetUsingGDAL(const std::filesystem::path& path);
void exportTopoSetUsingGDAL(const std::filesystem::path& path, const TopoSet<Inexact> topoSet, const CGAL::Aff_transformation_2<Inexact>& trans, std::optional<OGRSpatialReference> spatialReference, bool stackPolygons);

GeometrySet<Inexact> readGeometrySetUsingGDAL(const std::filesystem::path& path);
}
