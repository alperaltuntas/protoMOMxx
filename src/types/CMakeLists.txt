# CMake configuration for src/types
add_library(MOM_src_types
  MOM_grid.cpp
)

target_include_directories(MOM_src_types
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(MOM_src_types
  PUBLIC
    MOM_src_framework
)
add_amrex_to_protomom_target(MOM_src_types)
