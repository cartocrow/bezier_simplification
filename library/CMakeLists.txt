set(SOURCES
        steven_bezier_collapse.cpp
)
set(HEADERS
        bezier_graph.h
        bezier_simplification_traits.h
        common.h
        curved_graph.h
        curve_simplification_traits.h
        indexed_priority_queue.h
        steven_bezier_collapse.h
        curved_graph_with_history.h
        utils.h
)

add_library(curved_simplification ${SOURCES})
target_link_libraries(curved_simplification
        PUBLIC
        cartocrow::core
        cartocrow::renderer
)

install(FILES ${HEADERS} DESTINATION ${CARTOCROW_INSTALL_DIR}/curved_simplification)
