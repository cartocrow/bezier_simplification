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
            auto newV = g.insert_vertex(p);
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

    auto insertPolygonSet = [&](const PolygonSet<Inexact>& polygonSet) {
        std::vector<PolygonWithHoles<Inexact>> pgnsWithHoles;
        polygonSet.polygons_with_holes(std::back_inserter(pgnsWithHoles));
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
                insertPolygonSet(approximate(ogrMultiPolygonToPolygonSet(*poMultiPolygon)));
                break;
            }
            case wkbPolygon: {
                OGRPolygon* poly = poGeometry->toPolygon();
                insertPolygonWithHoles(approximate(ogrPolygonToPolygonWithHoles(*poly)));
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

        PolygonSet<Inexact> polygonSet;
        switch(wkbFlatten(poGeometry->getGeometryType())) {
            case wkbMultiPolygon: {
                OGRMultiPolygon *poMultiPolygon = poGeometry->toMultiPolygon();
                polygonSet = approximate(ogrMultiPolygonToPolygonSet(*poMultiPolygon));
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
                polygonSet = approximate(ogrPolygonToPolygonSet(*poly));
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

        regionSet.push_back(region);
    }

    return {regionSet, *poLayer->GetSpatialRef()};
}

// mixed geometries is a problem, we assume here toposet contains only PolygonSet geometries.
void exportTopoSetUsingGDAL(const std::filesystem::path& path, const TopoSet<Inexact> topoSet, const CGAL::Aff_transformation_2<Inexact>& trans, std::optional<OGRSpatialReference> spatialReference, bool stackPolygons) {
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

    poLayer = poDS->CreateLayer( "regions", spatialReference.has_value() ? &*spatialReference : nullptr, wkbMultiPolygon, NULL );
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
            auto psg = get<TopoSet<Inexact>::PolygonSetGeometry>(feature.geometry);
            auto ps = psg.getGeometry(topoSet);
            auto mPgn = polygonSetToOGRMultiPolygonT(ps, trans);
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
        /*OGRFieldDefn oField("stacking", OFTInteger);
        if (poLayer->CreateField(&oField) != OGRERR_NONE) {
            printf("Creating field failed.\n");
            exit(1);
        }*/

        std::vector<std::pair<PolygonSet<Inexact>, RegionAttributes>> polygonSets;

        for (const auto &feature: topoSet.features) {
            auto psg = get<TopoSet<Inexact>::PolygonSetGeometry>(feature.geometry);
            auto ps = psg.getGeometry(topoSet);
            PolygonSet <Inexact> transformed;
            std::vector <PolygonWithHoles<Inexact>> pgsWHs;
            ps.polygons_with_holes(std::back_inserter(pgsWHs));
            for (const auto &pgnWH: pgsWHs) {
                auto pgnWHT = transform(trans, pgnWH);
                transformed.insert(pgnWHT.outer_boundary());
            }
            polygonSets.emplace_back(transformed, feature.attributes);
        }

        /*IpeRenderer ipeRenderer;*/

        /*std::sort(polygonSets.begin(), polygonSets.end(), [&ipeRenderer](const auto& a, const auto& b) {
            std::cout << "===" << std::endl;
            const PolygonSet<Inexact>& pgs1 = a.first;
            const PolygonSet<Inexact>& pgs2 = b.first;
            std::vector<PolygonWithHoles<Inexact>> pgns1;
            pgs1.polygons_with_holes(std::back_inserter(pgns1));
            std::vector<PolygonWithHoles<Inexact>> pgns2;
            pgs2.polygons_with_holes(std::back_inserter(pgns2));

            for (const auto& pgn1 : pgns1) {
                for (const auto& pgn2 : pgns2) {
                    ipeRenderer.addPainting([pgn1, pgn2](GeometryRenderer& renderer) {
                        renderer.setStroke(Color(255, 0, 0), 1.0);
                        renderer.draw(pgn1);
                        renderer.setStroke(Color(0, 0, 255), 1.0);
                        renderer.draw(pgn2);
                        std::stringstream ss;
                        ss << "inside: " << (pgn1.outer_boundary().has_on_bounded_side(pgn2.outer_boundary().vertex(0)));
                        renderer.drawText({ 0, 0 }, ss.str());
                     });
                    ipeRenderer.nextPage();
                    if (pgn1.outer_boundary().has_on_bounded_side(pgn2.outer_boundary().vertex(0))) return true;
                }
            }

            return false;
           
        });*/

        //ipeRenderer.save("debugging_nesting.ipe");

        int index = 0;
        for (const auto& ps : polygonSets) {
            auto mPgn = polygonSetToOGRMultiPolygon(ps.first);
            OGRFeature *poFeature;

            poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
            for (const auto &[attribute, data]: ps.second) {
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
}