set(SOURCES
    min_dist_demo.cpp
)

set(HEADERS
    min_dist_demo.h
)

add_executable(min_dist_demo ${SOURCES})

target_link_libraries(
    min_dist_demo
    PRIVATE
    Qt::Widgets

    cartocrow::core
    cartocrow::renderer
    cartocrow::reader
    cartocrow::widgets
    cartocrow::circle_segment_helpers

    curved_simplification
)

install(TARGETS min_dist_demo DESTINATION ${INSTALL_BINARY_DIR})
