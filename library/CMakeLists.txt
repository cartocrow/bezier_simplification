set(SOURCES
    steven_bezier_collapse.cpp
    fit_cubic.cpp
    spiro_helpers/bezctx_cbs.cpp
)
set(HEADERS
    common.h
    core/bezier_graph_2.h
    core/graph_2.h
    core/graph_bezier_curve_traits_2.h
    core/graph_segment_traits_2.h
    core/graph_curve_traits_2.h
    core/graph_2_like.h
    core/region_set.h
    core/topo_set.h
    indexed_priority_queue.h
    steven_bezier_collapse.h
    curved_graph_with_history.h
    fit_cubic.h
    spiro_helpers/bezctx_cbs.h
    spiro_helpers/spiro_helpers.h
    bezier_curve_quad_tree.h
    utils.h

    spiro_helpers/spiro/bezctx.h
    spiro_helpers/spiro/bezctx_intf.h
    spiro_helpers/spiro/spiro-config.h
    spiro_helpers/spiro/spiro.h
    spiro_helpers/spiro/spiroentrypoints.h
    spiro_helpers/spiro/spiro_intf.h
)

add_library(curved_simplification ${SOURCES})
target_link_libraries(curved_simplification
        PUBLIC
        cartocrow::core
        cartocrow::renderer
)

install(FILES ${HEADERS} DESTINATION ${CARTOCROW_INSTALL_DIR}/curved_simplification)
