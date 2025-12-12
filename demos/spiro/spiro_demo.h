#pragma once

#include <cartocrow/renderer/geometry_widget.h>
#include <QMainWindow>

using namespace cartocrow;
using namespace cartocrow::renderer;

class SpiroDemo : public QMainWindow {
Q_OBJECT

private:
    GeometryWidget* m_renderer;

public:
    SpiroDemo();
};
