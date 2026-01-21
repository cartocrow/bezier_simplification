#pragma once

#include <cartocrow/renderer/geometry_widget.h>
//#include <cartocrow/widgets/double_slider.h>
#include "double_slider.h"
#include "library/core/bezier_graph_2.h"
#include "library/core/topo_set.h"
#include "library/bezier_collapse.h"
#include "library/steven_bezier_collapse.h"
#include "library/minimum_distance.h"

#include <QMainWindow>

#include <filesystem>

#include <ogrsf_frmts.h>

using namespace cartocrow;
using namespace cartocrow::renderer;
using namespace cartocrow::curved_simplification;

using Graph = BezierCollapseGraphWithHistoryAndIndex;
using BaseGraph = Graph::BaseGraph;
using ApproximatedGraph = ApproximatedBezierGraph<BaseGraph>;
using Traits = StevenBCTraits<Graph>;
using Collapse = BezierCollapse<Graph, Traits>;
using Edge_handle = Graph::Edge_handle;
using Vertex_handle = Graph::Vertex_handle;

class BezierSimplificationDemo : public QMainWindow {
	Q_OBJECT

  private:
	GeometryWidget* m_renderer;

	BaseGraph m_baseGraph;
    Graph m_graph;
	Collapse m_collapse;
    std::optional<Edge_handle> m_debugEdge;
    BaseGraph m_original;
    TopoSet<Inexact> m_toposet;
    CGAL::Aff_transformation_2<Inexact> m_transform;
    OGRSpatialReference m_spatialRef;
    DoubleSlider* m_minDist;
    MinimumDistanceForcer<typename ApproximatedGraph::Vertex_data, typename ApproximatedGraph::Edge_data> m_forcer;
    ApproximatedGraph m_approxGraph;

    std::optional<BaseGraph> m_backup;

    struct ControlPoint {
        Point<Inexact> point;
        std::variant<std::pair<Edge_handle, bool>, Vertex_handle> type;
    };

    std::vector<std::shared_ptr<ControlPoint>> m_editables;
    std::shared_ptr<ControlPoint> m_dragging;

    void loadInput(const std::filesystem::path& path);

  public:
	BezierSimplificationDemo();
};