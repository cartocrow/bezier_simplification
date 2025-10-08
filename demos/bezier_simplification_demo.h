#pragma once

#include <cartocrow/renderer/geometry_widget.h>
#include "../library/bezier_graph.h"
#include "../library/bezier_collapse.h"
#include "../library/steven_bezier_collapse.h"
#include <QMainWindow>

using namespace cartocrow;
using namespace cartocrow::renderer;
using namespace cartocrow::curved_simplification;

class BezierSimplificationDemo : public QMainWindow {
	Q_OBJECT

  private:
	GeometryWidget* m_renderer;
    using Graph = BezierCollapseGraphWithHistory;
    using BaseGraph = Graph::BaseGraph;
	using Traits = StevenBCTraits<Graph>;
	using Collapse = BezierCollapse<Graph, Traits>;
    using Edge_handle = Graph::Edge_handle;
	BaseGraph m_baseGraph;
    Graph m_graph;
	Collapse m_collapse;
    std::optional<Edge_handle> m_debugEdge;

  public:
	BezierSimplificationDemo();
    virtual ~BezierSimplificationDemo() {};
};