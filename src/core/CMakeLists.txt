# CMake configuration for src/core: the time stepping schemes and the term
# modules they sequence.
add_library(MOM_src_core
  MOM.cpp
  MOM_continuity_PPM.cpp
  MOM_CoriolisAdv.cpp
  MOM_dynamics_unsplit_RK2.cpp
  MOM_PressureForce.cpp
)

target_include_directories(MOM_src_core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(MOM_src_core
  PUBLIC
    MOM_src_initialization
    MOM_src_parameterizations
)
add_amrex_to_protomom_target(MOM_src_core)

# Floating-point contraction has to match the Fortran side expression by
# expression, not file by file (see the PROTOMOM_NO_FP_CONTRACT comment in the
# top-level CMakeLists). gfortran contracts almost nothing in these two
# modules, because MOM6 writes their intermediates -- rel_vort, KEx, the PPM
# edge values -- into arrays before using them, which blocks fusion. Turning
# contraction off here reproduces that; the expressions gfortran does contract
# are written out as std::fma in the source.
set_source_files_properties(MOM_CoriolisAdv.cpp MOM_continuity_PPM.cpp
  PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")
