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
# Floating-point contraction has to match the Fortran side expression by
# expression, not file by file (see the PROTOMOM_NO_FP_CONTRACT comment in the
# top-level CMakeLists). gfortran contracts almost nothing in these three
# modules, because MOM6 writes their intermediates -- rel_vort, KEx, the
# strain rates, the PPM edge values -- into arrays before using them, which
# blocks fusion. Turning contraction off here reproduces that; the handful of
# expressions gfortran does contract are written out as std::fma in the source.
set_source_files_properties(MOM_CoriolisAdv.cpp MOM_continuity_PPM.cpp
  ../parameterizations/lateral/MOM_hor_visc.cpp
  PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")
