#pragma once
/// @file MOM_continuity_PPM.h
/// @brief The continuity solver: the thickness fluxes from a piecewise
///        parabolic reconstruction of the layer thicknesses, and the
///        thickness update they imply. The analogue of MOM6's
///        MOM_continuity_PPM.

#include <AMReX_MultiFab.H>

#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_grid.h"
#include "MOM_vertical_grid.h"

namespace MOM {

/// @class ContinuityPPM
/// @brief Advances the layer thicknesses and produces the thickness fluxes.
///
/// The scheme is directionally split: it advects zonally over the range the
/// following meridional pass needs, then meridionally using the updated
/// thicknesses, as MOM6 does when FIRST_DIRECTION selects x first.
///
/// Only the unsplit path is implemented: no barotropic transport target
/// (uhbt/vhbt), no viscous remnant, no barotropic continuity structure. Those
/// are what the split scheme's iteration needs, and they are the bulk of
/// MOM6's module.
class ContinuityPPM {
public:
  /// @brief Read the continuity configuration. The analogue of MOM6's
  /// continuity_PPM_init.
  /// @param params Runtime parameters.
  /// @throws logger::FatalError on an unsupported option.
  explicit ContinuityPPM(RuntimeParams &params);

  /// @brief The halo width the scheme reads beyond the computational domain.
  /// @return The stencil width.
  int stencil() const { return 3; }

  /// @brief Advance the thicknesses and compute the thickness fluxes. The
  /// analogue of MOM6's continuity_PPM.
  /// @param h The updated layer thicknesses, written [H ~> m].
  /// @param uh The zonal thickness flux, written [H L2 T-1 ~> m3 s-1].
  /// @param vh The meridional thickness flux, written [H L2 T-1 ~> m3 s-1].
  /// @param u The zonal velocity [L T-1 ~> m s-1].
  /// @param v The meridional velocity [L T-1 ~> m s-1].
  /// @param hin The layer thicknesses at the start of the step [H ~> m].
  /// @param dt The time increment [T ~> s].
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid.
  void solve(amrex::MultiFab &h, amrex::MultiFab &uh, amrex::MultiFab &vh,
             const amrex::MultiFab &u, const amrex::MultiFab &v,
             const amrex::MultiFab &hin, amrex::Real dt, const Domain &domain,
             const Grid &grid, const VerticalGrid &vgrid) const;

private:
  // defer: the barotropic-correction machinery (uhbt/vhbt targets, the
  //        Newton iteration in zonal_flux_adjust, BT_cont), the viscous
  //        remnant, the monotonic and upwind reconstructions, the
  //        volume-based CFL, porous barriers, and open boundaries.
};

} // namespace MOM
