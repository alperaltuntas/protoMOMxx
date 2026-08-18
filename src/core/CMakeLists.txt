# CMake configuration for src/core
add_library(MOM_src_core
  MOM.cpp
  MOM_forcing_type.cpp
  MOM_grid.cpp
  MOM_state.cpp
  MOM_vertical_grid.cpp
)

target_include_directories(MOM_src_core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(MOM_src_core
  PUBLIC
    MOM_src_initialization
)
add_amrex_to_protomom_target(MOM_src_core)
