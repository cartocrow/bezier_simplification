set(SOURCES
    fit_cubic_demo.cpp
    ../read_ipe_bezier_spline.cpp
)

set(HEADERS
    fit_cubic_demo.h
    ../read_ipe_bezier_spline.h
)

add_executable(fit_cubic_demo ${SOURCES})

target_link_libraries(
    fit_cubic_demo
    PRIVATE
    Qt::Widgets

    cartocrow::core
    cartocrow::renderer
    cartocrow::reader

    curved_simplification
)

install(TARGETS fit_cubic_demo DESTINATION ${INSTALL_BINARY_DIR})
