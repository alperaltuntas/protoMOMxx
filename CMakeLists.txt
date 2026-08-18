cmake_minimum_required(VERSION 3.24)
project(protoMOMxx)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(PROTOMOM_CUDA "Build with Cuda backend" OFF)
option(PROTOMOM_ADD_TEST_TARGETS "Build test targets" OFF)
option(PROTOMOM_FETCH_DEPS "Fetch missing dependencies automatically" ON)
# Floating-point contraction has to MATCH the Fortran side of a parity run,
# not simply be disabled, and the match is per expression rather than per
# build: contraction only ever happens inside one expression, so what has to
# agree is where each model rounds. The reference MOM6 build here contracts, so
# the default is to leave contraction on and to force rounding at the few
# places where gfortran does not fuse and g++ does (see MOM_fp_contract.h and
# the note in src/core/CMakeLists.txt). Set this ON when the Fortran side is
# built with contraction disabled throughout.
option(PROTOMOM_NO_FP_CONTRACT "Disable FMA contraction. Set to match the Fortran side of a bit-for-bit parity run." OFF)

if(PROTOMOM_CUDA)
  enable_language(CUDA)
  find_package(AMReX REQUIRED CUDA)
  include(AMReXTargetHelpers)
else()
  find_package(AMReX REQUIRED)
endif()

include(cmake/protomomtarget.cmake)

# TIM (tim/cpp) is a self-contained CMake subproject yielding the TIM::tim
# target. Added after AMReX is found, so tim/cpp's guarded find_package(AMReX)
# is skipped and the AMReX already found by protoMOMxx is used. TIM_SOURCE_DIR
# defaults to the bundled extern/TIM submodule, but can be overridden to point
# at an existing TIM checkout.
set(TIM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/extern/TIM" CACHE PATH
  "Path to the TIM checkout (defaults to the extern/TIM submodule).")
if(NOT EXISTS "${TIM_SOURCE_DIR}/tim/cpp/CMakeLists.txt")
  message(FATAL_ERROR
    "TIM sources not found under '${TIM_SOURCE_DIR}/tim/cpp'.\n"
    "If you cloned protoMOMxx without submodules, initialize them:\n"
    "    git submodule update --init --recursive\n"
    "or point -DTIM_SOURCE_DIR=<path> at an existing TIM checkout.")
endif()
add_subdirectory("${TIM_SOURCE_DIR}/tim/cpp" "${CMAKE_BINARY_DIR}/extern/TIM/tim_cpp")

add_subdirectory(src/framework)
add_subdirectory(src/initialization)
add_subdirectory(src/core)
add_subdirectory(src/diagnostics)
add_subdirectory(config_src/drivers/solo_driver)
add_subdirectory(config_src)

if(PROTOMOM_ADD_TEST_TARGETS)
  enable_testing()
  add_subdirectory(tests)
endif()