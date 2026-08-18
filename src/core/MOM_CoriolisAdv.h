#pragma once
/// @file MOM_CoriolisAdv.h
/// @brief The Coriolis and momentum-advection accelerations. The analogue of
///        MOM6's MOM_CoriolisAdv (CorAdCalc).

#include <AMReX_MultiFab.H>

#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_grid.h"
#include "MOM_vertical_grid.h"

namespace MOM {

/// @class CoriolisAdv
/// @brief Computes CAu and CAv, the sum of the Coriolis acceleration and the
/// gradient of kinetic energy, from the velocities and the thickness fluxes.
///
/// The Sadourny 1975 energy-conserving vorticity scheme and the Arakawa &
/// Lamb kinetic energy are implemented, with the BOUND_CORIOLIS limiter.
/// Those are the defaults for the double_gyre configuration; the other
/// vorticity, kinetic-energy and PV-advection schemes abort.
class CoriolisAdv {
public:
  /// @brief Read the Coriolis and momentum-advection configuration. The
  /// analogue of MOM6's CoriolisAdv_init.
  /// @param params Runtime parameters.
  /// @throws logger::FatalError on an unsupported scheme.
  explicit CoriolisAdv(RuntimeParams &params);

  /// @brief Compute the Coriolis and momentum-advection accelerations. The
  /// analogue of MOM6's CorAdCalc.
  /// @param CAu The zonal acceleration, written [L T-2 ~> m s-2].
  /// @param CAv The meridional acceleration, written [L T-2 ~> m s-2].
  /// @param u The zonal velocity [L T-1 ~> m s-1].
  /// @param v The meridional velocity [L T-1 ~> m s-1].
  /// @param h The layer thicknesses [H ~> m].
  /// @param uh The zonal thickness flux [H L2 T-1 ~> m3 s-1].
  /// @param vh The meridional thickness flux [H L2 T-1 ~> m3 s-1].
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid.
  void calculate(amrex::MultiFab &CAu, amrex::MultiFab &CAv,
                 const amrex::MultiFab &u, const amrex::MultiFab &v,
                 const amrex::MultiFab &h, const amrex::MultiFab &uh,
                 const amrex::MultiFab &vh, const Domain &domain,
                 const Grid &grid, const VerticalGrid &vgrid) const;

private:
  bool bound_Coriolis_ = false;  ///< Bound CAu/CAv by the neighboring (f+rv)v estimates.
  bool no_slip_ = false;         ///< Use a no-slip vorticity boundary condition.

  // defer: the other vorticity schemes (Sadourny enstrophy, Arakawa & Hsu,
  //        Arakawa & Lamb, the blend, robust enstrophy, the WENO family), the
  //        energy-dissipating CORIOLIS_EN_DIS variant, the other KE schemes,
  //        Stokes vorticity, porous barriers, and open boundaries.
};

} // namespace MOM
