#pragma once

#include <cartocrow/renderer/geometry_widget.h>
#include <QMainWindow>

using namespace cartocrow;
using namespace cartocrow::renderer;

class TopoSetDemo : public QMainWindow {
Q_OBJECT

private:
    GeometryWidget* m_renderer;

public:
    TopoSetDemo();
};
