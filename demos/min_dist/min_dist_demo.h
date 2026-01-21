#pragma once

#include <cartocrow/renderer/geometry_widget.h>
#include "library/minimum_distance.h"
#include <QMainWindow>

#include <cartocrow/widgets/double_slider.h>

using namespace cartocrow;
using namespace cartocrow::renderer;
using namespace cartocrow::curved_simplification;

using StraightGraph = Straight_graph_2<std::monostate, std::monostate, Inexact>;

class MinDistDemo : public QMainWindow {
Q_OBJECT

private:
    GeometryWidget* m_renderer;
    MinimumDistanceForcer<std::monostate, std::monostate> m_forcer;
    StraightGraph m_g;
    StraightGraph m_ogg;

    DoubleSlider* m_minDist;
public:
    MinDistDemo();
};
