function(add_amrex_to_protomom_target)
  foreach(PROTOMOM_TARGET IN LISTS ARGN)
    add_mom_test("${PROTOMOM_TARGET}")
  endforeach()
endfunction()

function(add_amrex_to_protomom_target PROTOMOM_TARGET)
  if(PROTOMOM_CUDA)
    set_cpp_sources_to_cuda_language(${PROTOMOM_TARGET})
  endif()
  target_link_libraries(${PROTOMOM_TARGET} PUBLIC AMReX::amrex)

  # The build half of the bit-for-bit parity contract (DESIGN.md section 5):
  # floating-point contraction must match the Fortran side of the comparison.
  # GCC and Clang contract a*b+c into an FMA at -O by default, and so does
  # gfortran, so the default here is to leave contraction on. Turning it off
  # in C++ alone changes answers away from a stock MOM6 build.
  if(PROTOMOM_NO_FP_CONTRACT)
    target_compile_options(${PROTOMOM_TARGET} PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang>:-ffp-contract=off>
      $<$<CXX_COMPILER_ID:IntelLLVM,Intel>:-fp-model=precise>)
  endif()
endfunction()
