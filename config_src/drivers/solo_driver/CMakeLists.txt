# CMake configuration for the solo driver
add_executable(protoMOMxx MOM_driver.cpp MOM_surface_forcing.cpp)
set_target_properties(protoMOMxx PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

target_link_libraries(protoMOMxx PRIVATE
    MOM_src_framework
    MOM_src_core
    MOM_src_diagnostics
    MOM_config_src_infra
)
add_amrex_to_protomom_target(protoMOMxx)
