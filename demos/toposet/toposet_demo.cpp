#include "toposet_demo.h"

#include "library/core/region_set.h"
#include "library/core/topo_set.h"

#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QCheckBox>

using namespace cartocrow;

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    TopoSetDemo demo;
    demo.show();
    app.exec();
}

TopoSetDemo::TopoSetDemo() {
    RegionSet<Inexact> regionSet;
    PolygonSet<Inexact> polygonSet1;
//    PolygonSet<Inexact> polygonSet2;
    Polygon<Inexact> polygon1;
    polygon1.push_back({0, 0});
    polygon1.push_back({1, 0});
    polygon1.push_back({1, 1});
    polygon1.push_back({0, 1});
    polygonSet1.insert(polygon1);
//    Polygon<Inexact> polygon2;
//    polygon2.push_back({-1, 0});
//    polygon2.push_back({0, 0});
//    polygon2.push_back({0, 1});
//    polygon2.push_back({-1, 1});
//    polygonSet2.insert(polygon2);
    Polygon<Inexact> polygon3;
    polygon3.push_back({-5, 5});
    polygon3.push_back({-4, 5});
    polygon3.push_back({-4, 6});
    polygon3.push_back({-5, 6});
    PolygonSet<Inexact> polygonSet3;
    polygonSet3.insert(polygon3);
    Polygon<Inexact> polygon4;
    polygon4.push_back({-2, -2});
    polygon4.push_back({2, -2});
    polygon4.push_back({2, 2});
    polygon4.push_back({-2, 2});
    PolygonWithHoles<Inexact> withHoles;
    withHoles.outer_boundary() = polygon4;
    auto hole = polygon1;
    hole.reverse_orientation();
    withHoles.add_hole(hole);
    polygonSet3.insert(withHoles);

    regionSet.push_back(Region<Inexact>({}, polygonSet1));
//    regionSet.push_back(Region<Inexact>({}, polygonSet2));
    regionSet.push_back(Region<Inexact>({}, polygonSet3));
    TopoSet<Inexact> topoSet(regionSet);

    setWindowTitle("TopoSet demo");
    m_renderer = new GeometryWidget();
    m_renderer->setDrawAxes(false);
    setCentralWidget(m_renderer);

    m_renderer->addPainting([topoSet](GeometryRenderer& renderer) {
//        for (const auto& feature : topoSet.features) {
//            auto polygonSet = get<TopoSet<Inexact>::PolygonSetGeometry>(feature.geometry);
//            auto ps = polygonSet.getGeometry(topoSet);
//            renderer.draw(ps);
//        }
        for (const auto& arc : topoSet.arcs) {
            renderer.draw(arc);
        }
    }, "TopoSet");

    std::cout << "done" << std::endl;
}