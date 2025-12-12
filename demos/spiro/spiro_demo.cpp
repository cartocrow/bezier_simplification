#include "spiro_demo.h"

#include <spiroentrypoints.h>
#include <QApplication>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QCheckBox>

#include "../../library/spiro_helpers/spiro_helpers.h"

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    SpiroDemo demo;
    demo.show();
    app.exec();
}

SpiroDemo::SpiroDemo() {
    setWindowTitle("Spiro demo");
    m_renderer = new GeometryWidget();
    m_renderer->setDrawAxes(false);
    setCentralWidget(m_renderer);

    auto *dockWidget = new QDockWidget();
    addDockWidget(Qt::RightDockWidgetArea, dockWidget);
    auto *vWidget = new QWidget();
    auto *vLayout = new QVBoxLayout(vWidget);
    vLayout->setAlignment(Qt::AlignTop);
    dockWidget->setWidget(vWidget);
    auto* closed = new QCheckBox("Closed");
    vLayout->addWidget(closed);

    auto p1 = std::make_shared<Point<Inexact>>(0, 0);
    auto p2 = std::make_shared<Point<Inexact>>(1, 0);
    auto p3 = std::make_shared<Point<Inexact>>(2, 1);
    auto p4 = std::make_shared<Point<Inexact>>(3, 1);

    m_renderer->registerEditable(p1);
    m_renderer->registerEditable(p2);
    m_renderer->registerEditable(p3);
    m_renderer->registerEditable(p4);

    m_renderer->fitInView(Box(-1, -1, 4, 2));

    connect(closed, &QCheckBox::stateChanged, [this](){
        m_renderer->repaint();
    });

    m_renderer->addPainting([p1, p2, p3, p4, closed](GeometryRenderer& renderer) {
        renderer.setStroke(Color(0, 0, 0), 2.0);
        std::vector<Point<Inexact>> points({*p1, *p2, *p3, *p4});
        auto spline = spiro(points.begin(), points.end(), closed->isChecked());
        renderer.draw(spline);
        for (const auto& pt : points) {
            renderer.draw(pt);
        }
    }, "Spline");
}