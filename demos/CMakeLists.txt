set(SOURCES
    bezier_simplification_demo.cpp
    read_ipe_bezier_spline.cpp
)

set(HEADERS
    bezier_simplification_demo.h
    read_ipe_bezier_spline.h
)

add_executable(bezier_simplification_demo ${SOURCES})

target_link_libraries(
    bezier_simplification_demo
    PRIVATE
    Qt::Widgets

    cartocrow::core
    cartocrow::renderer
    cartocrow::reader

    curved_simplification
)

install(TARGETS bezier_simplification_demo DESTINATION ${INSTALL_BINARY_DIR})
