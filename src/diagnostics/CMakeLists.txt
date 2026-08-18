# CMake configuration for src/diagnostics.
add_library(MOM_src_diagnostics
  MOM_sum_output.cpp
)

target_include_directories(MOM_src_diagnostics
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(MOM_src_diagnostics
  PUBLIC
    MOM_src_core
)
add_amrex_to_protomom_target(MOM_src_diagnostics)
