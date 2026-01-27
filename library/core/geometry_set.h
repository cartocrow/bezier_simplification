#pragma once

#include <cartocrow/core/core.h>
#include <cartocrow/core/polyline.h>
#include <cartocrow/core/transform_helpers.h>
#include <cartocrow/core/polyline_set.h"

namespace cartocrow {
template<class... Ts>
struct overloaded : Ts... {
    using Ts::operator()...;
};
template<class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

template<class K>
struct GeometrySet {
    using Geometry = std::variant<PolygonSet<K>, PolygonWithHoles<K>, Polygon<K>, PolylineSet<K>, Polyline<K>>;
    std::vector<Geometry> geometries;

    GeometrySet<K> transform(const CGAL::Aff_transformation_2<K>& trans) {
        GeometrySet<K> transformed;
        for (const auto& g : geometries) {
            transformed.geometries.push_back(std::visit(
                overloaded{
                    [&](const PolygonWithHoles<K>& pwh) -> Geometry {
                        return cartocrow::transform(trans, pwh);
                    },
                    [&](const PolygonSet<K>& ps) -> Geometry {
                        return cartocrow::transform(trans, ps);
                    },
                    [&](const Polygon<K>& p) -> Geometry {
                        return CGAL::transform(trans, p);
                    },
                    [&](const auto& someG) -> Geometry {
                        return someG.transform(trans);
                    }
                }, g));
        }
        return transformed;
    }
};
}