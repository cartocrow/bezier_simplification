#pragma once

#include <cartocrow/renderer/geometry_widget.h>
#include "library/schneider.h"
#include <QMainWindow>

using namespace cartocrow;
using namespace cartocrow::renderer;
using namespace cartocrow::curved_simplification;

class SchneiderDemo : public QMainWindow {
Q_OBJECT

private:
    GeometryWidget* m_renderer;

public:
    SchneiderDemo();
};
