set(SOURCES
    wasm_frontend.cpp
)

add_executable(wasm_frontend ${SOURCES})

target_link_libraries(
    wasm_frontend
    PRIVATE
    curved_simplification
)