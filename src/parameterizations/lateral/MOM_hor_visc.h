#pragma once
/// @file MOM_hor_visc.h
/// @brief The horizontal viscosity accelerations. The analogue of MOM6's
///        MOM_hor_visc (horizontal_viscosity).

#include <AMReX_MultiFab.H>

#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_grid.h"
#include "MOM_vertical_grid.h"

namespace MOM {

/// @class HorizontalViscosity
/// @brief Computes diffu and diffv, the accelerations from the divergence of
/// the horizontal stress tensor.
///
/// The Laplacian branch is implemented, with the Smagorinsky and
/// velocity-scale viscosities and the BOUND_KH stability limiter. The
/// biharmonic branch, the Leith family, MEKE backscatter, GME and the
/// anisotropic viscosity all abort.
///
/// The metric combinations MOM6 precomputes in hor_visc_init (DY_dxT, dx2q,
/// the background and maximum viscosities, ...) are precomputed here too, as
/// member fields, because they depend on the timestep as well as the grid and
/// so do not belong on the Grid.
class HorizontalViscosity {
public:
  /// @brief Read the viscosity configuration and precompute the metric
  /// combinations. The analogue of MOM6's hor_visc_init.
  /// @param params Runtime parameters.
  /// @param dt The dynamics timestep, which sets the stability bound [T ~> s].
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @throws logger::FatalError on an unsupported option.
  HorizontalViscosity(RuntimeParams &params, amrex::Real dt, const Domain &domain,
                      const Grid &grid);

  /// @brief Compute the horizontal viscous accelerations. The analogue of
  /// MOM6's horizontal_viscosity.
  /// @param diffu The zonal acceleration, written [L T-2 ~> m s-2].
  /// @param diffv The meridional acceleration, written [L T-2 ~> m s-2].
  /// @param u The zonal velocity [L T-1 ~> m s-1].
  /// @param v The meridional velocity [L T-1 ~> m s-1].
  /// @param h The layer thicknesses [H ~> m].
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid.
  void calculate(amrex::MultiFab &diffu, amrex::MultiFab &diffv,
                 const amrex::MultiFab &u, const amrex::MultiFab &v,
                 const amrex::MultiFab &h, const Domain &domain,
                 const Grid &grid, const VerticalGrid &vgrid) const;

private:
  bool Smagorinsky_Kh_ = false;  ///< Use the Smagorinsky nonlinear viscosity.
  bool bound_Kh_ = true;         ///< Bound the viscosity for stability.
  amrex::Real Kh_bg_min_ = 0.0;  ///< The minimum viscosity [L2 T-1 ~> m2 s-1].

  // Metric combinations, precomputed as MOM6 precomputes them.
  amrex::MultiFab DY_dxT_;    ///< dyT/dxT at h points [nondim].
  amrex::MultiFab DX_dyT_;    ///< dxT/dyT at h points [nondim].
  amrex::MultiFab dx2h_;      ///< dxT squared [L2 ~> m2].
  amrex::MultiFab dy2h_;      ///< dyT squared [L2 ~> m2].
  amrex::MultiFab DY_dxBu_;   ///< dyBu/dxBu at q points [nondim].
  amrex::MultiFab DX_dyBu_;   ///< dxBu/dyBu at q points [nondim].
  amrex::MultiFab dx2q_;      ///< dxBu squared [L2 ~> m2].
  amrex::MultiFab dy2q_;      ///< dyBu squared [L2 ~> m2].
  amrex::MultiFab reduction_xx_;  ///< Blocked-face reduction factor at h points [nondim].
  amrex::MultiFab reduction_xy_;  ///< Blocked-face reduction factor at q points [nondim].
  amrex::MultiFab Kh_bg_xx_;      ///< Background viscosity at h points [L2 T-1 ~> m2 s-1].
  amrex::MultiFab Kh_bg_xy_;      ///< Background viscosity at q points [L2 T-1 ~> m2 s-1].
  amrex::MultiFab Laplac2_const_xx_;  ///< Smagorinsky constant times the squared grid scale
                                      ///< at h points [L2 ~> m2].
  amrex::MultiFab Laplac2_const_xy_;  ///< The same at q points [L2 ~> m2].
  amrex::MultiFab Kh_Max_xx_;     ///< The stability bound on the viscosity at h points
                                  ///< [L2 T-1 ~> m2 s-1].
  amrex::MultiFab Kh_Max_xy_;     ///< The same at q points [L2 T-1 ~> m2 s-1].
};

} // namespace MOM
