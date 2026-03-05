set(SOURCES
    frontend.cpp
)

set(HEADERS
    frontend.h
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
    cartocrow::widgets

    curved_simplification
)

install(TARGETS frontend DESTINATION ${INSTALL_BINARY_DIR})
