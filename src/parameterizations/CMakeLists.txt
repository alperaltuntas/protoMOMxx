# CMake configuration for src/parameterizations.
#
# The parameterizations read the model types and are called by src/core's time
# stepping schemes, so they sit between src/types and src/core. In MOM6 the
# same relationship is a directory cycle -- src/parameterizations uses
# MOM_grid from src/core, and src/core calls the parameterizations -- which
# MOM6's single link unit hides. Moving the types below both breaks it.
add_library(MOM_src_parameterizations
  lateral/MOM_hor_visc.cpp
  vertical/MOM_set_viscosity.cpp
  vertical/MOM_vert_friction.cpp
)

target_include_directories(MOM_src_parameterizations
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lateral
    ${CMAKE_CURRENT_SOURCE_DIR}/vertical
)

target_link_libraries(MOM_src_parameterizations
  PUBLIC
    MOM_src_types
)
add_amrex_to_protomom_target(MOM_src_parameterizations)

# Floating-point contraction has to match the Fortran side expression by
# expression, not file by file (see the PROTOMOM_NO_FP_CONTRACT comment in the
# top-level CMakeLists). gfortran contracts nothing in MOM6's horizontal
# viscosity, because it writes its intermediates -- the strain rates, the
# stresses -- into arrays before using them, which blocks fusion. The one
# expression it does contract is written out as std::fma in the source.
set_source_files_properties(lateral/MOM_hor_visc.cpp
  PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")
