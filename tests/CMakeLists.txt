find_package(GTest QUIET)
if(NOT GTest_FOUND)
  if(PROTOMOM_FETCH_DEPS)
    include(FetchContent)
    FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(googletest)
  else()
    message(FATAL_ERROR
      "GTest not found. Install GoogleTest or set CMAKE_PREFIX_PATH "
      "to include its install location, or configure with "
      "-DPROTOMOM_FETCH_DEPS=ON to download it automatically.")
  endif()
endif()
enable_testing()

add_executable(test_MOM_logger MOM_logger/test_MOM_logger.cpp)
add_executable(test_MOM_nml_parser MOM_nml_parser/test_MOM_nml_parser.cpp)
add_executable(test_MOM_io MOM_io/test_MOM_io.cpp)
add_executable(test_MOM_directories MOM_directories/test_MOM_directories.cpp)
add_executable(test_MOM_file_parser MOM_file_parser/test_MOM_file_parser.cpp)
add_executable(test_MOM_document MOM_document/test_MOM_document.cpp)
add_executable(test_MOM_model MOM_model/test_MOM_model.cpp)
add_executable(test_MOM_domains MOM_domains/test_MOM_domains.cpp)
add_executable(test_MOM_grid MOM_grid/test_MOM_grid.cpp)
add_executable(test_MOM_vertical_grid MOM_vertical_grid/test_MOM_vertical_grid.cpp)
add_executable(test_MOM_state MOM_state/test_MOM_state.cpp)
add_executable(test_MOM_fixed_initialization MOM_fixed_initialization/test_MOM_fixed_initialization.cpp)
target_link_libraries(test_MOM_logger GTest::gtest_main MOM_src_framework)
target_link_libraries(test_MOM_nml_parser GTest::gtest_main MOM_src_framework)
target_link_libraries(test_MOM_io GTest::gtest_main MOM_src_framework)
target_link_libraries(test_MOM_directories GTest::gtest_main MOM_src_framework)
target_link_libraries(test_MOM_file_parser GTest::gtest_main MOM_src_framework)
target_link_libraries(test_MOM_document GTest::gtest_main MOM_src_framework)
target_link_libraries(test_MOM_model GTest::gtest MOM_src_core MOM_config_src_infra)
target_link_libraries(test_MOM_domains GTest::gtest MOM_src_framework MOM_config_src_infra)
target_link_libraries(test_MOM_grid GTest::gtest MOM_src_core MOM_config_src_infra)
target_link_libraries(test_MOM_vertical_grid GTest::gtest_main MOM_src_core)
target_link_libraries(test_MOM_state GTest::gtest MOM_src_core MOM_config_src_infra)
target_link_libraries(test_MOM_fixed_initialization GTest::gtest MOM_src_core MOM_config_src_infra)

include(GoogleTest)
gtest_discover_tests(test_MOM_logger PROPERTIES LABELS "unit")
gtest_discover_tests(test_MOM_nml_parser PROPERTIES LABELS "unit")
gtest_discover_tests(test_MOM_io PROPERTIES LABELS "unit")
gtest_discover_tests(test_MOM_directories PROPERTIES LABELS "unit")
gtest_discover_tests(test_MOM_file_parser PROPERTIES LABELS "unit")
gtest_discover_tests(test_MOM_document PROPERTIES LABELS "unit")
gtest_discover_tests(test_MOM_model PROPERTIES LABELS "unit")
gtest_discover_tests(test_MOM_domains PROPERTIES LABELS "unit")
gtest_discover_tests(test_MOM_grid PROPERTIES LABELS "unit")
gtest_discover_tests(test_MOM_vertical_grid PROPERTIES LABELS "unit")
gtest_discover_tests(test_MOM_state PROPERTIES LABELS "unit")
gtest_discover_tests(test_MOM_fixed_initialization PROPERTIES LABELS "unit")

add_test(NAME test_MOM_domains_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_MOM_domains> ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_MOM_domains_mpi PROPERTIES LABELS "mpi" PROCESSORS 2)

add_test(NAME test_MOM_grid_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_MOM_grid> ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_MOM_grid_mpi PROPERTIES LABELS "mpi" PROCESSORS 2)

add_test(NAME test_MOM_fixed_initialization_mpi
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_MOM_fixed_initialization> ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_MOM_fixed_initialization_mpi PROPERTIES LABELS "mpi" PROCESSORS 2)

add_subdirectory(config_src)

# System/regression test: run the full protoMOMxx executable against the double_gyre example
# and confirm its MOM_parameters_doc.* output matches the golden copies committed under
# double_gyre/.
set(DOUBLE_GYRE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/double_gyre)
set(DOUBLE_GYRE_RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/double_gyre)
file(COPY ${DOUBLE_GYRE_SOURCE_DIR}/ DESTINATION ${DOUBLE_GYRE_RUN_DIR}
     PATTERN "MOM_parameters_doc.*" EXCLUDE)

add_test(NAME double_gyre_smoke
         COMMAND $<TARGET_FILE:protoMOMxx>
         WORKING_DIRECTORY ${DOUBLE_GYRE_RUN_DIR})
set_tests_properties(double_gyre_smoke PROPERTIES
  LABELS doubleGyreFileSetup
  FIXTURES_SETUP DoubleGyreDataGen)

add_test(NAME double_gyre_regression
         COMMAND ${CMAKE_COMMAND}
                 -DEXPECTED_DIR=${DOUBLE_GYRE_SOURCE_DIR}
                 -DACTUAL_DIR=${DOUBLE_GYRE_RUN_DIR}
                 -P ${DOUBLE_GYRE_SOURCE_DIR}/check_regression.cmake)
set_tests_properties(double_gyre_regression PROPERTIES
  LABELS system
  FIXTURES_REQUIRED DoubleGyreDataGen)
