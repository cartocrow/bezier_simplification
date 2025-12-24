set(SOURCES
    toposet_demo.cpp
)

set(HEADERS
    toposet_demo.h
)

add_executable(toposet_demo ${SOURCES})

target_link_libraries(
    toposet_demo
    PRIVATE
    Qt::Widgets

    cartocrow::core
    cartocrow::renderer
    cartocrow::reader

    curved_simplification
)

install(TARGETS spiro_demo DESTINATION ${INSTALL_BINARY_DIR})
