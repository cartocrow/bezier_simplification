#pragma once

#include <cartocrow/renderer/geometry_widget.h>
#include "../../library/fit_cubic.h"
#include <QMainWindow>

using namespace cartocrow;
using namespace cartocrow::renderer;
using namespace cartocrow::curved_simplification;

class FitCubicDemo : public QMainWindow {
Q_OBJECT

private:
    GeometryWidget* m_renderer;

public:
    FitCubicDemo();
};
