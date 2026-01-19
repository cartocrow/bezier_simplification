#pragma once

#include <cartocrow/renderer/geometry_widget.h>
#include "library/minimum_distance.h"
#include <QMainWindow>

#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>

#include <cartocrow/widgets/double_slider.h>

using namespace cartocrow;
using namespace cartocrow::renderer;
using namespace cartocrow::curved_simplification;

using StraightGraph = Straight_graph_2<std::monostate, std::monostate, Inexact>;

using SiteInfo = std::optional<typename StraightGraph::Edge_handle>;

struct Converter {
    SiteInfo operator()(const SiteInfo& info, bool is_src) {
        return info;
    }
    SiteInfo operator()(const SiteInfo& info1, const SiteInfo& info2, bool first) {
        if (info1.has_value()) {
            return *info1;
        }
        return info2;
    }
};

struct Merger {
    SiteInfo operator()(const SiteInfo& info1, const SiteInfo& info2) {
        if (info1.has_value()) {
            return *info1;
        }
        return info2;
    }
};

typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<Inexact, CGAL::Field_with_sqrt_tag> Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt, CGAL::Segment_Delaunay_graph_storage_traits_with_info_2<Gt, SiteInfo, Converter, Merger>> SDG;

using VoronoiEdge = std::variant<Segment<Inexact>, Line<Inexact>, Ray<Inexact>, CGAL::Parabola_segment_2<Gt>>;

class MinDistDemo : public QMainWindow {
Q_OBJECT

private:
    GeometryWidget* m_renderer;
    StraightGraph m_g;
    StraightGraph m_ogg;
    SDG m_delaunay;

    std::vector<VoronoiEdge> m_withinDistEdges;
    std::vector<std::variant<Point<Inexact>, Segment<Inexact>>> m_withinDistIsolineParts;

    DoubleSlider* m_minDist;
public:
    MinDistDemo();
    void recompute();
};
