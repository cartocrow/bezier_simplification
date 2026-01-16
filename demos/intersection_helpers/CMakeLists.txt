set(SOURCES
    intersection_helpers_demo.cpp
)

set(HEADERS
    intersection_helpers_demo.h
)

add_executable(intersection_helpers_demo ${SOURCES})

target_link_libraries(
    intersection_helpers_demo
    PRIVATE
    Qt::Widgets

    cartocrow::core
    cartocrow::renderer
    cartocrow::reader
    cartocrow::widgets
    cartocrow::circle_segment_helpers

    curved_simplification
)

install(TARGETS intersection_helpers_demo DESTINATION ${INSTALL_BINARY_DIR})
