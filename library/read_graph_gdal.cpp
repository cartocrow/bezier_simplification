#include "read_graph_gdal.h"
#include <cartocrow/reader/gdal_conversion.h>
#include <cartocrow/core/transform_helpers.h>

#include <cartocrow/renderer/ipe_renderer.h>

using namespace cartocrow::renderer;

namespace cartocrow::curved_simplification {
OGRLinearRing polygonToOGRLinearRingT(const Polygon<Inexact>& polygon, const CGAL::Aff_transformation_2<Inexact>& trans) {
    OGRLinearRing ring;
    for (const auto& v : polygon.vertices()) {
        auto v_ = trans.transform(v);
        ring.addPoint(v_.x(), v_.y());
    }
    auto v_ = trans.transform(polygon.vertices().front());
    ring.addPoint(v_.x(), v_.y());

    return ring;
}

OGRPolygon polygonWithHolesToOGRPolygonT(const PolygonWithHoles<Inexact>& polygon, const CGAL::Aff_transformation_2<Inexact>& trans) {
    assert(!polygon.is_unbounded());
    auto outerRing = polygonToOGRLinearRingT(polygon.outer_boundary(), trans);
    std::vector<OGRLinearRing> holeRings;
    for (const auto& h : polygon.holes()) {
        holeRings.push_back(polygonToOGRLinearRingT(h, trans));
    }
    OGRPolygon ogrPolygon;
    ogrPolygon.addRing(&outerRing);
    for (auto& hr : holeRings) {
        ogrPolygon.addRing(&hr);
    }
    return ogrPolygon;
}

OGRMultiPolygon polygonSetToOGRMultiPolygonT(const PolygonSet<Inexact>& polygonSet, const CGAL::Aff_transformation_2<Inexact>& trans) {
    std::vector<PolygonWithHoles<Inexact>> pgns;
    polygonSet.polygons_with_holes(std::back_inserter(pgns));

    OGRMultiPolygon ogrMultiPolygon;
    for (const auto& pgn : pgns) {
        auto ogrPolygon = polygonWithHolesToOGRPolygonT(pgn, trans);
        ogrMultiPolygon.addGeometry(&ogrPolygon);
    }

    return ogrMultiPolygon;
}

OGRMultiPolygon polygonSetRawToOGRMultiPolygonT(const PolygonSetRaw<Inexact>& polygonSetRaw, const CGAL::Aff_transformation_2<Inexact>& trans) {
    const std::vector<PolygonWithHoles<Inexact>>& pgns = polygonSetRaw.polygons_with_holes;

    OGRMultiPolygon ogrMultiPolygon;
    for (const auto& pgn : pgns) {
        auto ogrPolygon = polygonWithHolesToOGRPolygonT(pgn, trans);
        ogrMultiPolygon.addGeometry(&ogrPolygon);
    }

    return ogrMultiPolygon;
}

Straight_graph_2<std::monostate, std::monostate, Inexact> readGraphUsingGDAL(const std::filesystem::path& path) {
    GDALAllRegister();
    GDALDataset *poDS;

    poDS = (GDALDataset*) GDALOpenEx( path.string().c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr );
    if( poDS == nullptr ) {
        printf( "GDAL open failed.\n" );
        exit( 1 );
    }
    OGRLayer* poLayer;
    poLayer = poDS->GetLayer(0);

    poLayer->ResetReading();

    using Graph = Straight_graph_2<std::monostate, std::monostate, Inexact>;

    Graph g;

    std::unordered_map<Point<Inexact>, Graph::Vertex_handle> pToV;

    auto getVertex = [&pToV, &g](const Point<Inexact>& p) {
        if (pToV.contains(p)) {
            return pToV.at(p);
        } else {
            auto newV = g.add_vertex(p);
            pToV[p] = newV;
            return newV;
        }
    };

    auto insertPolygon = [&](const Polygon<Inexact> polygon) {
        if (polygon.is_empty()) return;
        auto lastV = getVertex(polygon.vertex(0));
        pToV[polygon.vertex(0)] = lastV;
        for (const auto& seg : polygon.edges()) {
            auto targetV = getVertex(seg.target());

            bool oppositeExists = false;
            bool alreadyExists = false;
            for (auto ieit = targetV->incident_edges_begin(); ieit != targetV->incident_edges_end(); ++ieit) {
                if ((*ieit)->target() == lastV && (*ieit)->curve().opposite() == seg) {
                    oppositeExists = true;
                    break;
                }
                if ((*ieit)->source() == lastV && (*ieit)->curve() == seg) {
                    alreadyExists = true;
                    break;
                }
            }
            if (!oppositeExists && !alreadyExists) {
                g.add_edge(lastV, targetV, seg);
            }
            lastV = targetV;
        }
    };

    auto insertPolygonWithHoles = [&](const PolygonWithHoles<Inexact> polygonWithHoles) {
        insertPolygon(polygonWithHoles.outer_boundary());
        for (const auto& hole : polygonWithHoles.holes()) {
            insertPolygon(hole);
        }
    };

    auto insertPolygonSet = [&](const PolygonSetRaw<Inexact>& polygonSet) {
        const std::vector<PolygonWithHoles<Inexact>>& pgnsWithHoles = polygonSet.polygons_with_holes;
        for (const auto& pgnWithHoles : pgnsWithHoles) {
            insertPolygonWithHoles(pgnWithHoles);
        }
    };

    for (auto& poFeature : *poLayer) {
        OGRGeometry *poGeometry;

        poGeometry = poFeature->GetGeometryRef();
        switch(wkbFlatten(poGeometry->getGeometryType())) {
            case wkbMultiPolygon: {
                OGRMultiPolygon *poMultiPolygon = poGeometry->toMultiPolygon();
                insertPolygonSet(ogrMultiPolygonToPolygonSetRaw(*poMultiPolygon));
                break;
            }
            case wkbPolygon: {
                OGRPolygon* poly = poGeometry->toPolygon();
                insertPolygonWithHoles(ogrPolygonToPolygonWithHoles(*poly));
                break;
            }
//            case wkbLineString: {
//            }
            default: std::cout << "Did not handle this type of geometry: " << poGeometry->getGeometryName() << std::endl;
        }
    }

    return g;
}

std::pair<RegionSet<Inexact>, OGRSpatialReference> readRegionSetUsingGDAL(const std::filesystem::path& path) {
    GDALAllRegister();
    GDALDataset *poDS;

    poDS = (GDALDataset*) GDALOpenEx( path.string().c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr );
    if( poDS == nullptr ) {
        printf( "GDAL open failed.\n" );
        exit( 1 );
    }
    OGRLayer* poLayer;
    poLayer = poDS->GetLayer(0);

    poLayer->ResetReading();

    RegionSet<Inexact> regionSet;

    for (auto& poFeature : *poLayer) {
        OGRGeometry *poGeometry;

        poGeometry = poFeature->GetGeometryRef();

        PolygonSetRaw<Inexact> polygonSet;
        switch(wkbFlatten(poGeometry->getGeometryType())) {
            case wkbMultiPolygon: {
                OGRMultiPolygon *poMultiPolygon = poGeometry->toMultiPolygon();
                polygonSet = ogrMultiPolygonToPolygonSetRaw(*poMultiPolygon);
                break;
            }
            case wkbPolygon: {
                OGRPolygon* poly = poGeometry->toPolygon();
//                std::cout << "=====" << std::endl;
//                for (const auto& ring : poly) {
//                    std::cout << "===" << std::endl;
//                    for (const auto& pt : ring) {
//                        std::cout << pt.getX() << ", " << pt.getY() << std::endl;
//                    }
//                }
                polygonSet = ogrPolygonToPolygonSetRaw(*poly);
                break;
            }
            default: std::cout << "Did not handle this type of geometry: " << poGeometry->getGeometryName() << std::endl;
        }

        Region<Inexact> region;
        region.geometry = polygonSet;

        int i = 0;
        for( auto&& oField: *poFeature ) {
            std::string name = poFeature->GetDefnRef()->GetFieldDefn(i)->GetNameRef();
            switch(oField.GetType()) {
                case OFTInteger:
                    region.attributes[name] = static_cast<int>(oField.GetInteger());
                    break;
                case OFTReal:
                    region.attributes[name] = oField.GetDouble();
                    break;
                case OFTInteger64:
                    region.attributes[name] = static_cast<int64_t>(oField.GetInteger64());
                    break;
                default:
                    std::cout << "Did not handle this type of attribute: " << oField.GetType() << std::endl;
                    break;
            }
            ++i;
        }

        regionSet.regions.push_back(region);
    }

    return {regionSet, *poLayer->GetSpatialRef()};
}

// mixed geometries is a problem, we assume here toposet contains only PolygonSet geometries.
void exportTopoSetUsingGDAL(const std::filesystem::path& path, const StraightTopoSet<Inexact> topoSet, const CGAL::Aff_transformation_2<Inexact>& trans, std::optional<OGRSpatialReference> spatialReference, bool stackPolygons) {
    using PolygonSetTopology = StraightTopoSet<Inexact>::PolygonSetTopology;
    
    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;

    GDALAllRegister();

    poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName);
    if( poDriver == nullptr )
    {
        printf("%s driver not available.\n", pszDriverName);
        exit( 1 );
    }

    GDALDataset *poDS;

    poDS = poDriver->Create(path.string().c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if( poDS == nullptr )
    {
        printf( "Creation of output file failed.\n" );
        exit( 1 );
    }

    OGRLayer *poLayer;

    auto layerName = path.stem().string();
    poLayer = poDS->CreateLayer(layerName.c_str(), spatialReference.has_value() ? &*spatialReference : nullptr, wkbMultiPolygon, NULL);
    if( poLayer == NULL )
    {
        printf( "Layer creation failed.\n" );
        exit( 1 );
    }

//    auto topoSet = pretendExact(topoSetInexact);

    for (const auto& [attribute, value] : topoSet.features[0].attributes) {
        OGRFieldDefn oField = [&]() {
            if (std::holds_alternative<double>(value)) {
                return OGRFieldDefn(attribute.c_str(), OFTReal);
            } else if (std::holds_alternative<int>(value)) {
                return OGRFieldDefn(attribute.c_str(), OFTInteger);
            } else if (std::holds_alternative<int64_t>(value)) {
                return OGRFieldDefn(attribute.c_str(), OFTInteger64);
            } else {
                std::cerr << "Did not handle attribute type." << std::endl;
            }
            // todo: other attributes
        }();

        if (poLayer->CreateField(&oField) != OGRERR_NONE) {
            printf("Creating field failed.\n");
            exit(1);
        }
    }

    if (!stackPolygons) {
        for (const auto &feature: topoSet.features) {
            auto psg = get<PolygonSetTopology>(feature.topology);
            auto ps = psg.getGeometry<Inexact>(topoSet);
            auto mPgn = polygonSetRawToOGRMultiPolygonT(ps, trans);
            OGRFeature *poFeature;

            poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
            for (const auto &[attribute, data]: feature.attributes) {
                if (auto vDouble = std::get_if<double>(&data)) {
                    poFeature->SetField(attribute.c_str(), *vDouble);
                } else if (auto vInt = std::get_if<int>(&data)) {
                    poFeature->SetField(attribute.c_str(), static_cast<int>(*vInt));
                } else if (auto vInt64 = std::get_if<int64_t>(&data)) {
                    poFeature->SetField(attribute.c_str(), static_cast<GIntBig>(*vInt64));
                } else {
                    std::cout << "Did not handle attribute value." << std::endl;
                }
            }

            poFeature->SetGeometry(&mPgn);

            if (poLayer->CreateFeature(poFeature) != OGRERR_NONE) {
                printf("Failed to create feature in shapefile.\n");
                exit(1);
            }

            OGRFeature::DestroyFeature(poFeature);
        }
    } else {
        std::vector<std::tuple<PolygonSetRaw<Inexact>, RegionAttributes, int>> entries;

        //auto bbox = topoSet.bbox();
        for (const auto &feature: topoSet.features) {
            auto psg = get<PolygonSetTopology>(feature.topology);
            auto ps = psg.getGeometry<Inexact>(topoSet);
            PolygonSetRaw <Inexact> transformed;
            const std::vector <PolygonWithHoles<Inexact>>& pgsWHs = ps.polygons_with_holes;
            for (const auto &pgnWH: pgsWHs) {
                //auto pgn = closeAlongBox(pgnWH.outer_boundary(), bbox);
                auto pgn = pgnWH.outer_boundary();
                auto pgnT = transform(trans, pgn);
                PolygonWithHoles<Inexact> outer(pgnT);
                transformed.polygons_with_holes.push_back(outer);
            }
            entries.emplace_back(transformed, feature.attributes, 0);
        }

        for (int i = 0; i < entries.size(); ++i) {
            for (int j = 0; j < entries.size(); ++j) {
                if (i == j) continue;
                auto& [pgs1, _1, d1] = entries[i];
                auto& [pgs2, _2, d2] = entries[j];
                const std::vector<PolygonWithHoles<Inexact>>& pgns1 = pgs1.polygons_with_holes;
                const std::vector<PolygonWithHoles<Inexact>>& pgns2 = pgs2.polygons_with_holes;

                bool contained = false;
                for (const auto& pgn1 : pgns1) {
                    for (const auto& pgn2 : pgns2) {
                        if (pgn2.outer_boundary().has_on_bounded_side(pgn1.outer_boundary().vertex(0))) {
                            contained = true;
                            ++d1;
                            break;
                        }
                    }
                    if (contained) break;
                }
            }
        }

        std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
            return std::get<2>(a) < std::get<2>(b);
        });

        int index = 0;
        for (const auto& ps : entries) {
            auto mPgn = polygonSetRawToOGRMultiPolygonT(std::get<0>(ps), CGAL::IDENTITY);
            OGRFeature *poFeature;

            poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
            for (const auto &[attribute, data]: std::get<1>(ps)) {
                if (auto vDouble = std::get_if<double>(&data)) {
                    poFeature->SetField(attribute.c_str(), *vDouble);
                } else if (auto vInt = std::get_if<int>(&data)) {
                    poFeature->SetField(attribute.c_str(), static_cast<int>(*vInt));
                } else if (auto vInt64 = std::get_if<int64_t>(&data)) {
                    poFeature->SetField(attribute.c_str(), static_cast<GIntBig>(*vInt64));
                } else {
                    std::cout << "Did not handle attribute value." << std::endl;
                }
            }

            //poFeature->SetField("stacking", index);

            poFeature->SetGeometry(&mPgn);

            if (poLayer->CreateFeature(poFeature) != OGRERR_NONE) {
                printf("Failed to create feature in shapefile.\n");
                exit(1);
            }

            OGRFeature::DestroyFeature(poFeature);

            ++index;
        }
    }
    GDALClose( poDS );
}

GeometrySet<Inexact> readGeometrySetUsingGDAL(const std::filesystem::path& path) {
    GDALAllRegister();
    GDALDataset *poDS;

    poDS = (GDALDataset*) GDALOpenEx( path.string().c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr );
    if( poDS == nullptr ) {
        printf( "GDAL open failed.\n" );
        exit( 1 );
    }
    OGRLayer* poLayer;
    poLayer = poDS->GetLayer(0);

    poLayer->ResetReading();

    GeometrySet<Inexact> geometrySet;

    for (auto& poFeature : *poLayer) {
        OGRGeometry *poGeometry;
        poGeometry = poFeature->GetGeometryRef();

        // todo: extract into subroutine
        switch (wkbFlatten(poGeometry->getGeometryType())) {
            case wkbMultiPolygon: {
                OGRMultiPolygon *poMultiPolygon = poGeometry->toMultiPolygon();
                geometrySet.geometries.push_back(ogrMultiPolygonToPolygonSetRaw(*poMultiPolygon));
                break;
            }
            case wkbPolygon: {
                OGRPolygon *poly = poGeometry->toPolygon();
                geometrySet.geometries.push_back(ogrPolygonToPolygonWithHoles(*poly));
                break;
            }
            case wkbLinearRing: {
                OGRLinearRing *poly = poGeometry->toLinearRing();
                geometrySet.geometries.push_back(ogrLinearRingToPolygon(*poly));
                break;
            }
            case wkbLineString: {
                OGRLineString *pl = poGeometry->toLineString();
                geometrySet.geometries.push_back(ogrLineStringToPolyline(*pl));
                break;
            }
            case wkbMultiLineString: {
                OGRMultiLineString *pl = poGeometry->toMultiLineString();
                geometrySet.geometries.push_back(ogrMultiLineStringToPolylineSet(*pl));
                break;
            }
            case wkbPoint: {
                OGRPoint* p = poGeometry->toPoint();
                geometrySet.geometries.push_back(ogrPointToPoint(*p));
                break;
            }
            case wkbMultiPoint: {
                OGRMultiPoint* mp = poGeometry->toMultiPoint();
                geometrySet.geometries.push_back(ogrMultiPointToPointSet(*mp));
                break;
            }
            default: std::cout << "Did not handle this type of geometry: " << poGeometry->getGeometryName() << std::endl;
        }
    }

    return geometrySet;
}
}