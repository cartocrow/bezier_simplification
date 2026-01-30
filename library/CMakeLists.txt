if (NOT EMSCRIPTEN)
set(SOURCES
    steven_bezier_collapse.cpp
    schneider.cpp
    read_graph_gdal.cpp
    intersection_helpers.cpp

    spiro_helpers/spiro/bezctx.c
    spiro_helpers/spiro/spiro.c
    spiro_helpers/spiro/spiroentrypoints.c

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
    core/geometry_set.h
    indexed_priority_queue.h
    steven_bezier_collapse.h
    curved_graph_with_history.h
    schneider.h
    spiro_helpers/bezctx_cbs.h
    spiro_helpers/spiro_helpers.h
    bezier_curve_quad_tree.h
    utils.h
    minimum_distance.h
    read_graph_gdal.h
    conic_types.h
    intersection_helpers.h

    spiro_helpers/spiro/bezctx.h
    spiro_helpers/spiro/bezctx_intf.h
    spiro_helpers/spiro/spiro-config.h
    spiro_helpers/spiro/spiro.h
    spiro_helpers/spiro/spiroentrypoints.h
    spiro_helpers/spiro/spiro_intf.h
)
else()
set(SOURCES
    spiro_helpers/spiro/bezctx.c
    spiro_helpers/spiro/spiro.c
    spiro_helpers/spiro/spiroentrypoints.c

    steven_bezier_collapse.cpp
    schneider.cpp
    spiro_helpers/bezctx_cbs.cpp
    intersection_helpers.cpp
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
    ../../core/cartocrow/core/polyline_set.h
    core/geometry_set.h
    indexed_priority_queue.h
    steven_bezier_collapse.h
    curved_graph_with_history.h
        schneider.h
    spiro_helpers/bezctx_cbs.h
    spiro_helpers/spiro_helpers.h
    bezier_curve_quad_tree.h
    utils.h
    minimum_distance.h
    conic_types.h
    intersection_helpers.h

    spiro_helpers/spiro/bezctx.h
    spiro_helpers/spiro/bezctx_intf.h
    spiro_helpers/spiro/spiro-config.h
    spiro_helpers/spiro/spiro.h
    spiro_helpers/spiro/spiroentrypoints.h
    spiro_helpers/spiro/spiro_intf.h
)
endif()

add_library(curved_simplification ${SOURCES})

if (NOT EMSCRIPTEN)
target_link_libraries(curved_simplification
        PUBLIC
        cartocrow::core
        cartocrow::renderer
        cartocrow::reader
        cartocrow::circle_segment_helpers
)
else()
target_link_libraries(curved_simplification
        PUBLIC
        cartocrow::core
        cartocrow::renderer
        cartocrow::circle_segment_helpers
)
endif()

install(FILES ${HEADERS} DESTINATION ${CARTOCROW_INSTALL_DIR}/curved_simplification)
