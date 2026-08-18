# CMake configuration for src/types.
#
# The passive model types and their construction precursors. See README.md for
# what belongs here and why it sits below src/initialization.
add_library(MOM_src_types
  MOM_forcing_type.cpp
  MOM_grid.cpp
  MOM_state.cpp
  MOM_vertical_grid.cpp
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
