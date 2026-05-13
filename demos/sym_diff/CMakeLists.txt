set(SOURCES
    sym_diff_demo.cpp
)

set(HEADERS
    sym_diff_demo.h
)

add_executable(sym_diff_demo ${SOURCES})

target_link_libraries(
    sym_diff_demo
    PRIVATE
    Qt::Widgets

    cartocrow::core
    cartocrow::renderer

    curved_simplification
)

install(TARGETS sym_diff_demo DESTINATION ${INSTALL_BINARY_DIR})
