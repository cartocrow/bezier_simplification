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

class ComplexitySliders : public QObject {
    Q_OBJECT
private:
    QLabel* m_complexityLabel;
    QSlider* m_complexity;
    DoubleSlider* m_complexityLog;
    
    int m_value;

public:
    ComplexitySliders(QLayout* layout) {
        m_complexityLabel = new QLabel("#Edges: ");
        m_complexity = new QSlider();
        m_complexityLog = new DoubleSlider();

        m_complexity->setOrientation(Qt::Horizontal);
        layout->addWidget(m_complexityLabel);
        layout->addWidget(m_complexity);
        m_complexityLog->setOrientation(Qt::Horizontal);
        m_complexityLog->setPrecision(1000);
        layout->addWidget(m_complexityLog);

        connect(m_complexity, &QSlider::valueChanged, [this](int value) {
            setValue(value);
        });

        connect(m_complexityLog, &DoubleSlider::valueChanged, [this](double value) {
            setValue(std::exp(value));
        });
    }

    void setValue(int v) {
        m_value = v;
        m_complexity->blockSignals(true);
        m_complexity->setValue(v);
        m_complexity->blockSignals(false);
        m_complexityLog->blockSignals(true);
        m_complexityLog->setValue(log(v));
        m_complexityLog->blockSignals(false);

        m_complexityLabel->setText(QString::fromStdString("#Edges: " + std::to_string(v)));

        emit valueChanged(v);
    }

    void setMinimum(double v) {
        m_complexity->setMinimum(v);
        m_complexityLog->setMinimum(log(v));
    }

    void setMaximum(double v) {
        m_complexity->setMaximum(v);
        m_complexityLog->setMaximum(log(v));
    }

    int minimum() {
        return m_complexity->minimum();
    }

    int maximum() {
        return m_complexity->maximum();
    }

    int value() const {
        return m_value;
    }

    void setEnabled(bool enabled) {
        m_complexity->setEnabled(enabled);
        m_complexityLog->setEnabled(enabled);
    }

signals:
    void valueChanged(double newValue);
};

using InnerControlPoint = std::pair<Edge_handle, bool>;

struct ControlPoint {
    std::variant<InnerControlPoint, Vertex_handle> type;

    Point<Inexact> point() const {
        if (const auto icpP = std::get_if<InnerControlPoint>(&type)) {
            const auto& curve = icpP->first->curve();
            return icpP->second ? curve.targetControl() : curve.sourceControl();
        }
        else if (const auto vhP = std::get_if<Vertex_handle>(&type)) {
            return (*vhP)->point();
        }
        else {
            throw std::runtime_error("Unknown control point type!");
        }
    }
};

struct Dragging {
    ControlPoint controlPoint;
    Point<Inexact> startPoint;
};

class GraphPainting : public GeometryPainting {
public:
    BaseGraph& m_graph;

    struct DrawSettings {
        CGAL::Aff_transformation_2<Inexact> m_trans;
        bool m_showVertices;
        bool m_showEdgeDirection;
        bool m_showControlPoints;

        DrawSettings() : m_trans(CGAL::IDENTITY), m_showVertices(false), m_showEdgeDirection(false), m_showControlPoints(false) {};
    };

    DrawSettings m_drawSettings;

    GraphPainting(BaseGraph& graph, DrawSettings drawSettings = DrawSettings{}) : m_graph(graph), m_drawSettings(drawSettings) {};

    void paint(GeometryRenderer& renderer) const override {
        const auto& m_trans = m_drawSettings.m_trans;

        const auto& m_showVertices = m_drawSettings.m_showVertices;
        const auto& m_showEdgeDirection = m_drawSettings.m_showEdgeDirection;
        const auto& m_showControlPoints = m_drawSettings.m_showControlPoints;

        renderer.setMode(GeometryRenderer::stroke);
        renderer.setStroke(Color(0, 0, 0), 2.0);
        for (auto eit = m_graph.edges_begin(); eit != m_graph.edges_end(); ++eit) {
//            if (m_showDebugInfo->isChecked() && eit->data().collapse.has_value()) {
//                renderer.setStroke(Color(50, 200, 50), 2.0);
//            } else {
            renderer.setStroke(Color(0, 0, 0), 2.0);
//            }
            renderer.draw(eit->curve().transform(m_trans));
            if (m_showEdgeDirection) {
                renderer.setStroke(Color(0, 0, 0), 6.0);
                renderer.draw(eit->curve().split(0.25).first.transform(m_trans));
            }
        }
        if (m_showControlPoints) {
            // Control polylines
            for (auto eit = m_graph.edges_begin(); eit != m_graph.edges_end(); ++eit) {
                renderer.setStroke(Color(0, 255, 0), 1.0);
                Polyline<Inexact> pl;
                for (int c = 0; c < 4; ++c) pl.push_back(eit->curve().control(c));
                renderer.draw(pl.transform(m_trans));
            }
            // Control points
            for (auto eit = m_graph.edges_begin(); eit != m_graph.edges_end(); ++eit) {
                renderer.setStroke(Color(0, 0, 255), 2.0);
                renderer.draw(eit->curve().sourceControl().transform(m_trans));
                renderer.draw(eit->curve().targetControl().transform(m_trans));
                renderer.setStroke(Color(255, 0, 255), 2.0);
                renderer.draw(eit->curve().source().transform(m_trans));
                renderer.draw(eit->curve().target().transform(m_trans));
            }
        }
        if (m_showVertices) {
            for (auto vit = m_graph.vertices_begin(); vit != m_graph.vertices_end(); ++vit) {
                renderer.draw(vit->point().transform(m_trans));
            }
        }
    }
};

class BezierSimplificationDemo : public QMainWindow {
	Q_OBJECT

  private:
	GeometryWidget* m_renderer;

	BaseGraph m_baseGraph;
    Graph m_graph;
    BaseGraph m_editBaseGraph;
    Graph m_editGraph;
    BaseGraph m_beforeReconstruct;
	Collapse m_collapse;
    std::optional<Edge_handle> m_debugEdge;
    BaseGraph m_original;
    TopoSet<Inexact> m_toposet;
    CGAL::Aff_transformation_2<Inexact> m_transform;
    CGAL::Aff_transformation_2<Inexact> m_backupTransform;
    OGRSpatialReference m_spatialRef;
    std::unique_ptr<DoubleSliderSpinBox> m_minDist;
    std::unique_ptr<DoubleSliderSpinBox> m_minAdjDist;
    std::unique_ptr<DoubleSliderSpinBox> m_minComponentLength;
    std::unique_ptr<ComplexitySliders> m_complexitySliders;
    Forcer m_forcer;
    DoubleSlider* m_minAngle;
    ApproximatedGraph m_approxGraph;
    QTabWidget* m_tabs;
    QCheckBox* m_editControlPoints;
    QSpinBox* m_desiredComplexity;
    QCheckBox* m_showOldVertices;
    QCheckBox* m_showDebugInfo;
    QCheckBox* m_ignoreBbox;
    QCheckBox* m_showCurvature;
    DoubleSlider* m_curvatureScale;

   void baseModified(bool modified);

    std::vector<std::shared_ptr<ControlPoint>> m_editables;
    std::shared_ptr<Dragging> m_dragging;
    std::vector<ReferenceData> m_referenceData;
    Box m_referencePolygon;

    PaintingRenderer m_voronoiPainting;
    PaintingRenderer m_intersectionsPainting;

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

    std::shared_ptr<GraphPainting> m_graphPainting;
    std::shared_ptr<GraphPainting> m_editGraphPainting = nullptr;

    bool m_shiftDown = false;
    bool m_altDown = false;
    std::optional<std::pair<BaseGraph::Edge_handle, CubicBezierCurve::CurvePoint>> m_shiftNearest;
    std::optional<BaseGraph::Vertex_handle> m_altNearest;

    void loadInput(const std::filesystem::path& path);
    void repaintVoronoi();
    void addIOTab();
    void addSimplificationTab();
    void addMinimumDistanceTab();
    void addDrawingTab();
    void addPaintings();
    void updateComplexityInfo();
    void resetEdits();
    void checkIntersections();
    void updateEditables();
    double getScale() const;

    void keyPressEvent(QKeyEvent *event) override;
    void keyReleaseEvent(QKeyEvent *event) override;

  public:
	BezierSimplificationDemo();
};