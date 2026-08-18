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
  # floating-point contraction has to match the Fortran side of the
  # comparison, because contraction happens inside a single expression and so
  # decides where each model rounds. Both sides disable it. That is already
  # what MOM6's intel and nvhpc mkmf templates do (-no-fma, -Mnofma); the gnu
  # template needs -ffp-contract=off adding to match them.
  #
  # Disabling it is not a numerical preference and not optional: with
  # contraction on, the C++ compiler fuses expressions the Fortran compiler
  # leaves alone and the answers separate in the last bit. Matching it here
  # rather than expression by expression is what makes the parity hold under
  # all three compilers instead of only under gcc.
  target_compile_options(${PROTOMOM_TARGET} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off>
    $<$<CXX_COMPILER_ID:IntelLLVM,Intel>:-fp-model=precise;-ffp-contract=off>
    $<$<CXX_COMPILER_ID:NVHPC>:-Kieee;-Mnofma>)
endfunction()
