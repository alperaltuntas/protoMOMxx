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
