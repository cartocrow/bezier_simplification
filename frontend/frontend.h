#pragma once

#include <cartocrow/renderer/geometry_widget.h>
#include <cartocrow/widgets/double_slider.h>
#include "library/core/bezier_graph_2.h"
#include "library/core/topo_set.h"
#include "library/core/geometry_set.h"
#include "library/bezier_collapse.h"
#include "library/steven_bezier_collapse.h"
#include "library/minimum_distance.h"

#include <cartocrow/renderer/painting_renderer.h>

#include <QMainWindow>
#include <QTabWidget>
#include <QCheckBox>
#include <QSpinBox>
#include <QSlider>
#include <QObject>
#include <QLabel>

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
using Forcer = MinimumDistanceForcer<typename ApproximatedGraph::Vertex_data, typename ApproximatedGraph::Edge_data>;

using ReferenceData = std::variant<QImage, GeometrySet<Inexact>>;

class DoubleSliderSpinBox : public QObject {
    Q_OBJECT

private:
    double m_value;
    DoubleSlider* m_slider;
    QDoubleSpinBox* m_spinBox;

public:
    DoubleSliderSpinBox() = default;
    DoubleSliderSpinBox& operator=(const DoubleSliderSpinBox& other) = default;

    DoubleSliderSpinBox(DoubleSlider* slider, QDoubleSpinBox* spinBox) : m_slider(slider), m_spinBox(spinBox) {
        connect(m_slider, &DoubleSlider::valueChanged, [this](double v) {
            setValue(v);
        });
        connect(m_spinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), [this](double v) {
            setValue(v);
        });
    }

    void setValue(double v) {
        m_value = v;
        m_slider->blockSignals(true);
        m_slider->setValue(v);
        m_slider->blockSignals(false);
        m_spinBox->blockSignals(true);
        m_spinBox->setValue(v);
        m_spinBox->blockSignals(false);

        emit valueChanged(v);
    }

    void setMinimum(double v) {
        m_slider->setMinimum(v);
        m_spinBox->setMinimum(v);
    }

    void setMaximum(double v) {
        m_slider->setMaximum(v);
        m_spinBox->setMaximum(v);
    }

    double value() const {
        return m_value;
    }

signals:
    void valueChanged(double newValue);
};

class BezierSimplificationDemo : public QMainWindow {
	Q_OBJECT

  private:
	GeometryWidget* m_renderer;

	BaseGraph m_baseGraph;
    Graph m_graph;
    BaseGraph m_beforeReconstruct;
	Collapse m_collapse;
    std::optional<Edge_handle> m_debugEdge;
    BaseGraph m_original;
    TopoSet<Inexact> m_toposet;
    CGAL::Aff_transformation_2<Inexact> m_transform;
    OGRSpatialReference m_spatialRef;
    std::unique_ptr<DoubleSliderSpinBox> m_minDist;
    std::unique_ptr<DoubleSliderSpinBox> m_minAdjDist;
    std::unique_ptr<DoubleSliderSpinBox> m_minComponentLength;
    Forcer m_forcer;
    ApproximatedGraph m_approxGraph;
    QTabWidget* m_tabs;
    QCheckBox* m_editControlPoints;
    QLabel* complexityLabel;
    QSlider* complexity;
    DoubleSlider* complexityLog;
    QSpinBox* desiredComplexity;
    QCheckBox* showEdgeDirection;
    QCheckBox* showOldVertices;
    QCheckBox* showNewVertices;
    QCheckBox* showNewControlPoints;
    QCheckBox* showDebugInfo;

    BaseGraph m_backup;

    struct ControlPoint {
        Point<Inexact> point;
        std::variant<std::pair<Edge_handle, bool>, Vertex_handle> type;
    };

    std::vector<std::shared_ptr<ControlPoint>> m_editables;
    std::shared_ptr<ControlPoint> m_dragging;
    std::vector<ReferenceData> m_referenceData;
    Box m_referencePolygon;

    PaintingRenderer m_voronoiPainting;

    std::vector<Color> m_colors = {
        {0xB9E1EE},
        {0x9AC019},
        {0xCD6814},
        {0xE53389},
        {0xC1BC56},
        {0x923B8B},
        {0xFBD2AA},
        {0x999999},
        {0xFECD0F},
        {0xCB9A03},
        {0xF3983B},
        {0x4B8EC7},
        {0x2E9A67},
        {0xE95937},
        {0xF8EE82},
        {0xE74646},
        {0xCBBC9D},
        {0x6699CD},
        {0x6FC4C6},
        {0xF1979A},
        {0x8F5A9C},
        {0xBB3087},
    };

    void loadInput(const std::filesystem::path& path);
    void repaintVoronoi();
    void addIOTab();
    void addSimplificationTab();
    void addMinimumDistanceTab();
    void addDrawingTab();
    void addPaintings();
    void updateComplexityInfo();
    void resetEdits();
    double getScale() const;

  public:
	BezierSimplificationDemo();
};