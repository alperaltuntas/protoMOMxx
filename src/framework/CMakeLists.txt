# CMake configuration file for src/framework
add_library(MOM_src_framework
  MOM_logger.cpp
  MOM_parser_utils.cpp
  MOM_param_table.cpp
  MOM_directories.cpp
  MOM_nml_parser.cpp
  MOM_file_parser.cpp
  MOM_document.cpp
  MOM_domains.cpp
  MOM_clock.cpp
  MOM_coms.cpp
  MOM_debug_dump.cpp
)

target_include_directories(MOM_src_framework
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(MOM_src_framework
  PUBLIC
    MOM_config_src_infra
)
add_amrex_to_protomom_target(MOM_src_framework)
