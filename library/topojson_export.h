#pragma once

#include <filesystem>

#include "bezier_topo_set.h"
#include <ogrsf_frmts.h>

namespace cartocrow {
void exportTopoSetToJson(const std::filesystem::path& path, const BezierTopoSet& topoSet, std::optional<OGRSpatialReference> spatialReference);
std::pair<BezierTopoSet, OGRSpatialReference> parseBezierTopoJson(const std::filesystem::path& path);
}