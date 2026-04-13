set(SHARED_SOURCES
    spiro_helpers/spiro/bezctx.c
    spiro_helpers/spiro/spiro.c
    spiro_helpers/spiro/spiroentrypoints.c

    steven_bezier_collapse.cpp
    schneider.cpp
    spiro_helpers/bezctx_cbs.cpp
    intersection_helpers.cpp
    vertex_snap.cpp

)
set(SHARED_HEADERS
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
    conic_types.h
    intersection_helpers.h
    vertex_snap.h

    spiro_helpers/spiro/bezctx.h
    spiro_helpers/spiro/bezctx_intf.h
    spiro_helpers/spiro/spiro-config.h
    spiro_helpers/spiro/spiro.h
    spiro_helpers/spiro/spiroentrypoints.h
    spiro_helpers/spiro/spiro_intf.h
)

set(NATIVE_SOURCES
    read_graph_gdal.cpp
    read_ipe_bezier_spline.cpp
)
set(NATIVE_HEADERS
    read_graph_gdal.h
    read_ipe_bezier_spline.h
)

if (EMSCRIPTEN)
set(SOURCES
    ${SHARED_SOURCES}
)
set(HEADERS
    ${SHARED_HEADERS}
)
else()
set(SOURCES
    ${SHARED_SOURCES}
    ${NATIVE_SOURCES}
)
set(HEADERS
    ${SHARED_HEADERS}
    ${NATIVE_HEADERS}
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

        GDAL::GDAL
)
else()
target_link_libraries(curved_simplification
        PUBLIC
        cartocrow::core
        cartocrow::renderer
        cartocrow::circle_segment_helpers

        GDAL::GDAL
)
endif()

install(FILES ${HEADERS} DESTINATION ${CARTOCROW_INSTALL_DIR}/curved_simplification)
