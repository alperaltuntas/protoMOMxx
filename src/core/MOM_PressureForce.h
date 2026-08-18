#pragma once
/// @file MOM_PressureForce.h
/// @brief The horizontal pressure gradient accelerations. The analogue of
///        MOM6's MOM_PressureForce dispatch plus the Boussinesq branch of
///        MOM_PressureForce_Montgomery.

#include <AMReX_MultiFab.H>

#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_grid.h"
#include "MOM_vertical_grid.h"

namespace MOM {

/// @class PressureForce
/// @brief Computes the pressure gradient accelerations PFu and PFv from the
/// layer thicknesses.
///
/// Only the Montgomery-potential form in Boussinesq mode without an equation
/// of state is implemented, which is the layered, adiabatic case: the
/// Montgomery potential of each layer follows from the interface heights and
/// the reduced gravities, and the accelerations are its horizontal gradient.
/// MOM6's default is the finite-volume form (PRESSUREFORCE = "FV"), which
/// needs the density integrals; that is the next step here.
class PressureForce {
public:
  /// @brief Read the pressure force configuration. The analogue of MOM6's
  /// PressureForce_init.
  /// @param params Runtime parameters.
  /// @throws logger::FatalError on an unsupported configuration.
  explicit PressureForce(RuntimeParams &params);

  /// @brief Compute the pressure gradient accelerations. The analogue of
  /// MOM6's PressureForce_Mont_Bouss.
  /// @param PFu The zonal acceleration, written [L T-2 ~> m s-2].
  /// @param PFv The meridional acceleration, written [L T-2 ~> m s-2].
  /// @param h The layer thicknesses [H ~> m].
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid (the reduced gravities).
  void calculate(amrex::MultiFab &PFu, amrex::MultiFab &PFv,
                 const amrex::MultiFab &h, const Domain &domain,
                 const Grid &grid, const VerticalGrid &vgrid) const;

private:
  /// @brief Build the interface heights and the Montgomery potential, one
  /// thread per column with the vertical recursions inside it.
  ///
  /// This is the form calculate() uses. It is also the form a GPU wants: one
  /// launch with a thread per column, the sequential dependence confined to a
  /// thread, and the fast thread index running along i so a warp's reads
  /// coalesce. On a CPU it is the slower of the two, because the innermost
  /// loop then steps by a whole horizontal plane and cannot be vectorised.
  /// montgomery_by_plane() is the same computation with the loops inverted;
  /// docs/loop_order.md has the measurements and the reasoning.
  ///
  /// @param columns The (i,j) footprint to build over, one k index deep.
  /// @param nk The number of layers.
  /// @param g_prime The reduced gravities, nk+1 of them, device-visible.
  /// @param hh The layer thicknesses [H ~> m].
  /// @param D The bottom depth [Z ~> m].
  /// @param ee The interface heights, written [Z ~> m].
  /// @param MM The Montgomery potential, written [L2 T-2 ~> m2 s-2].
  static void montgomery_by_column(const amrex::Box &columns, int nk,
                                   const amrex::Real *g_prime,
                                   const amrex::Array4<const amrex::Real> &hh,
                                   const amrex::Array4<const amrex::Real> &D,
                                   const amrex::Array4<amrex::Real> &ee,
                                   const amrex::Array4<amrex::Real> &MM);

  /// @brief The same, as a sweep of horizontal planes: the recursion in k is
  /// the outer loop and each of its steps is a kernel over the whole (i,j)
  /// footprint.
  ///
  /// Nothing calls this. It is kept because it is the CPU-shaped alternative
  /// to montgomery_by_column(), it is measurably faster there, and the two
  /// are exactly interchangeable -- the arithmetic and its order are
  /// identical, so ocean.stats does not move. Swapping the call in
  /// calculate() is the whole experiment.
  ///
  /// Why it is faster on a CPU: the innermost loop now runs along i, which is
  /// contiguous, so eight iterations share a cache line and four fit in an
  /// AVX2 register. icpx vectorises it four wide; the by-column form it
  /// cannot vectorise at all, and emits scalar vaddsd. Measured at -O2, best
  /// of three, in nanoseconds per gridpoint per timestep: 10.3 against 3.2 at
  /// 44x40x2, and 14.1 against 8.2 at 44x40x64.
  ///
  /// Why it is not the default: it is the wrong shape for the GPU. It trades
  /// one kernel launch for nk of them with a barrier between each, and at
  /// these footprints that is mostly launch latency. Which form should win is
  /// a decision about the backend, not about the physics, so it belongs in a
  /// loop-shape policy rather than in this file -- see docs/loop_order.md.
  ///
  /// Parameters are those of montgomery_by_column().
  static void montgomery_by_plane(const amrex::Box &columns, int nk,
                                  const amrex::Real *g_prime,
                                  const amrex::Array4<const amrex::Real> &hh,
                                  const amrex::Array4<const amrex::Real> &D,
                                  const amrex::Array4<amrex::Real> &ee,
                                  const amrex::Array4<amrex::Real> &MM);

  // defer: the FV/AFV form, the atmospheric pressure term, self-attraction
  //        and loading, tides, and the pbce/eta outputs the split scheme
  //        needs.
};

} // namespace MOM
