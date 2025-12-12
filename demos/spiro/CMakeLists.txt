set(SOURCES
    spiro_demo.cpp
)

set(HEADERS
    spiro_demo.h
)

add_executable(spiro_demo ${SOURCES})

target_link_libraries(
    spiro_demo
    PRIVATE
    Qt::Widgets

    cartocrow::core
    cartocrow::renderer
    cartocrow::reader

    spiro
    curved_simplification
)

install(TARGETS spiro_demo DESTINATION ${INSTALL_BINARY_DIR})
