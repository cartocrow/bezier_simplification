set(SOURCES
    schneider_demo.cpp
)

set(HEADERS
    schneider_demo.h
)

add_executable(schneider_demo ${SOURCES})

target_link_libraries(
    schneider_demo
    PRIVATE
    Qt::Widgets

    cartocrow::core
    cartocrow::renderer
    cartocrow::reader

    curved_simplification
)

install(TARGETS schneider_demo DESTINATION ${INSTALL_BINARY_DIR})
