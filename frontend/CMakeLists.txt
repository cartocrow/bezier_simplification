set(SOURCES
    frontend.cpp
    read_ipe_bezier_spline.cpp
    double_slider.cpp
)

set(HEADERS
    frontend.h
    read_ipe_bezier_spline.h
    double_slider.h
)

set(CMAKE_AUTOMOC ON)

add_executable(frontend ${SOURCES})

target_link_libraries(
    frontend
    PRIVATE
    Qt::Widgets

    cartocrow::core
    cartocrow::renderer
    cartocrow::reader

    curved_simplification
)

install(TARGETS frontend DESTINATION ${INSTALL_BINARY_DIR})
