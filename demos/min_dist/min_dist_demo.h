#pragma once

#include <cartocrow/renderer/geometry_widget.h>
#include "library/minimum_distance.h"
#include <QMainWindow>

#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>

#include <cartocrow/widgets/double_slider.h>

using namespace cartocrow;
using namespace cartocrow::renderer;
using namespace cartocrow::curved_simplification;

typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<Inexact, CGAL::Field_with_sqrt_tag> Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt> SDG;

using VoronoiEdge = std::variant<Segment<Inexact>, Line<Inexact>, Ray<Inexact>, CGAL::Parabola_segment_2<Gt>>;

class MinDistDemo : public QMainWindow {
Q_OBJECT

private:
    GeometryWidget* m_renderer;
    Straight_graph_2<std::monostate, std::monostate, Inexact> m_g;
    Straight_graph_2<std::monostate, std::monostate, Inexact> m_ogg;
    SDG m_delaunay;

    std::vector<VoronoiEdge> m_withinDistEdges;

    DoubleSlider* m_minDist;
public:
    MinDistDemo();
    void recompute();
};
