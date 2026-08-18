#pragma once
/// @file MOM_fp_contract.h
/// @brief A barrier that keeps the compiler from fusing a multiply into a
///        following add, for bit-for-bit parity with legacy MOM6.

#include <AMReX_REAL.H>
#include <AMReX_Extension.H>

namespace MOM {

/// @brief Return a value the compiler must treat as already rounded.
///
/// Both compilers contract a*b+c into a fused multiply-add at -O, but they do
/// not make the same choices, because contraction only happens inside one
/// expression. Where MOM6 writes an intermediate into an array before using
/// it -- dudx and dvdy in the horizontal viscosity, for instance -- gfortran
/// cannot fuse across the store, while the same quantity written as a local in
/// C++ is fused into whatever comes next. Wrapping such an intermediate here
/// reproduces the Fortran rounding.
///
/// This is parity scaffolding, not a numerical preference: it does not change
/// the result of any expression, only which rounding the result carries. The
/// asm is empty; it exists to make the value opaque to the optimizer. On
/// targets without the constraint (or on a GPU build) the barrier is absent
/// and bit-for-bit agreement with a Fortran build is not guaranteed.
/// @param x The value to round.
/// @return The same value.
AMREX_FORCE_INLINE amrex::Real fp_rounded(amrex::Real x) {
#if defined(__GNUC__) && defined(__x86_64__) && !defined(AMREX_USE_GPU)
  __asm__("" : "+x"(x));
#endif
  return x;
}

} // namespace MOM
