#include "core/region_set.h"
#include <cartocrow/data_structures/point_quad_tree.h>

namespace cartocrow::curved_simplification {
void removeDuplicates(Polygon<Inexact>& p);
void removeDuplicates(RegionSet<Inexact>& regionSet);
void snapVertices(RegionSet<Inexact>& regionSet, std::optional<double> epsilon = std::nullopt);
}