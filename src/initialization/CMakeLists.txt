# CMake configuration for src/initialization
add_library(MOM_src_initialization
  MOM_fixed_initialization.cpp
  MOM_grid_initialize.cpp
  MOM_shared_initialization.cpp
  MOM_state_initialization.cpp
)

target_include_directories(MOM_src_initialization
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(MOM_src_initialization
  PUBLIC
    MOM_src_types
)
add_amrex_to_protomom_target(MOM_src_initialization)
