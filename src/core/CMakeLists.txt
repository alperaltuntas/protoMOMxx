# CMake configuration for src/core.
#
# The library also carries src/parameterizations. MOM6's directory graph has
# cycles -- src/parameterizations uses MOM_grid from src/core, and src/core's
# dynamics drivers call the parameterizations -- which MOM6's single link unit
# hides. Splitting them into separate CMake libraries would need the cycle
# broken first; until then they share one library and keep MOM6's directory
# layout.
add_library(MOM_src_core
  MOM.cpp
  MOM_continuity_PPM.cpp
  MOM_CoriolisAdv.cpp
  MOM_dynamics_unsplit_RK2.cpp
  MOM_forcing_type.cpp
  MOM_grid.cpp
  MOM_PressureForce.cpp
  MOM_state.cpp
  MOM_vertical_grid.cpp
  ../parameterizations/lateral/MOM_hor_visc.cpp
  ../parameterizations/vertical/MOM_set_viscosity.cpp
  ../parameterizations/vertical/MOM_vert_friction.cpp
)

target_include_directories(MOM_src_core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../parameterizations/lateral
    ${CMAKE_CURRENT_SOURCE_DIR}/../parameterizations/vertical
)

target_link_libraries(MOM_src_core
  PUBLIC
    MOM_src_initialization
)
add_amrex_to_protomom_target(MOM_src_core)
