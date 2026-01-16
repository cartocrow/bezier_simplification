#pragma once

#include <cartocrow/renderer/geometry_widget.h>
#include "library/minimum_distance.h"
#include <QMainWindow>

using namespace cartocrow;
using namespace cartocrow::renderer;

class IntersectionHelpersDemo : public QMainWindow {
    Q_OBJECT

private:
    GeometryWidget* m_renderer;
    std::shared_ptr<Circle<Inexact>> m_circle;
public:
    IntersectionHelpersDemo();
};
