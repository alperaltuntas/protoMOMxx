#pragma once
/// @file MOM_kernel_inline.h
/// @brief An attribute that makes inlining a ParallelFor kernel mandatory.

#include <AMReX_Config.H>

/// @brief Put after a ParallelFor lambda's parameter list to force it inline.
///
/// AMReX writes each kernel as a lambda handed to ParallelFor, and the
/// compiler decides per call site whether to inline it into the loop. When it
/// declines, every grid point costs a call, with the captured Array4
/// descriptors reloaded across it, which is more work than the arithmetic in
/// most of these kernels.
///
/// The decision is a size heuristic, so it turns on how large a kernel happens
/// to have grown rather than on anything about the kernel, and it does not
/// hold still: gcc inlines all of them at -O3 and declines for six of them at
/// -O2, where the main loop then takes twice as long. A kernel body is a loop
/// body and is always worth inlining, so the choice is taken away from the
/// heuristic rather than tuned with --param.
///
/// This is not `AMREX_FORCE_INLINE` for a language reason rather than a
/// stylistic one: that macro expands to `inline __attribute__((always_inline))`
/// and the `inline` keyword is not part of a lambda-declarator, so writing it
/// here is a syntax error ("'inline' invalid in lambda"). AMReX has no
/// attribute-only spelling of it. Its one attribute-only inlining macro,
/// `AMREX_FLATTEN`, does compile here, and put on the enclosing routine it does
/// get the kernels inlined -- but it inlines everything else those routines
/// call as well, which leaves icpx executing three times the instructions and
/// running 2.4x slower. The narrow attribute on the lambda is the tool that
/// fits.
///
/// GCC, Clang and icpx honour it. nvc++ accepts it and does essentially
/// nothing with it -- the same kernels stay out of line and the run time does
/// not move -- which is presumably why AMReX defines `AMREX_FORCE_INLINE` as
/// plain `inline` there. A GPU build is left alone, so the heuristic decides
/// as before.
#if defined(__GNUC__) && !defined(AMREX_USE_GPU)
#define MOM_KERNEL_INLINE __attribute__((always_inline))
#else
#define MOM_KERNEL_INLINE
#endif
