#pragma once
/// @file MOM_dynamics_unsplit_RK2.h
/// @brief The unsplit quasi-second-order Runge-Kutta time stepping of the
///        ocean dynamics. The analogue of MOM6's
///        MOM_dynamics_unsplit_RK2 (step_MOM_dyn_unsplit_RK2).

#include <AMReX_MultiFab.H>

#include "MOM_continuity_PPM.h"
#include "MOM_CoriolisAdv.h"
#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_forcing_type.h"
#include "MOM_grid.h"
#include "MOM_hor_visc.h"
#include "MOM_PressureForce.h"
#include "MOM_set_viscosity.h"
#include "MOM_state.h"
#include "MOM_vert_friction.h"
#include "MOM_vertical_grid.h"

namespace MOM {

/// @class DynamicsUnsplitRK2
/// @brief Advances u, v and h one dynamics step with the unsplit
/// quasi-second-order Runge-Kutta scheme.
///
/// The scheme takes a predictor step of length BE*dt using the accelerations
/// evaluated at the start of the step, then re-evaluates the Coriolis and
/// advection terms at the predicted state and takes the full step. The
/// pressure gradient and the horizontal viscosity are evaluated once, at the
/// start, which is what makes the scheme quasi- rather than fully
/// second-order.
///
/// This class owns the term modules, which is the C++ replacement for MOM6's
/// nest of control structures inside MOM_dyn_unsplit_RK2_CS.
class DynamicsUnsplitRK2 {
public:
  /// @brief Construct the dynamics and its term modules. The analogue of
  /// MOM6's initialize_dyn_unsplit_RK2.
  /// @param params Runtime parameters.
  /// @param dt The dynamics timestep [T ~> s].
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid.
  DynamicsUnsplitRK2(RuntimeParams &params, amrex::Real dt, const Domain &domain,
                     const Grid &grid, const VerticalGrid &vgrid);

  /// @brief The number of velocity truncations since the last call, and reset
  /// the count. MOM6 shares the same counter between set_visc and
  /// MOM_sum_output through a pointer.
  /// @return The truncation count.
  int take_truncations() { return vert_friction_.take_truncations(); }

  /// @brief Advance the state one dynamics step. The analogue of MOM6's
  /// step_MOM_dyn_unsplit_RK2.
  /// @param state The prognostic state, updated in place.
  /// @param forces The mechanical forcing over this step.
  /// @param dt The dynamics timestep [T ~> s].
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid.
  void step(State &state, const MechForcing &forces, amrex::Real dt,
            const Domain &domain, const Grid &grid, const VerticalGrid &vgrid);

  /// @brief The zonal thickness flux of the last step [H L2 T-1 ~> m3 s-1].
  /// @return The zonal thickness flux.
  const amrex::MultiFab &uh() const { return uh_; }

  /// @brief The meridional thickness flux of the last step [H L2 T-1 ~> m3 s-1].
  /// @return The meridional thickness flux.
  const amrex::MultiFab &vh() const { return vh_; }

private:
  amrex::Real BE_ = 0.6;    ///< The weighting between RK2 and backward Euler [nondim].
  amrex::Real begw_ = 0.0;  ///< The gravity-wave damping weight [nondim].

  PressureForce pressure_force_;
  CoriolisAdv coriolis_adv_;
  ContinuityPPM continuity_;
  HorizontalViscosity hor_visc_;
  SetViscosity set_visc_;
  VertFriction vert_friction_;
  VertVisc visc_;

  // The accelerations and the working state, held as members so they are
  // allocated once rather than each step.
  amrex::MultiFab CAu_, CAv_;      ///< Coriolis and advection [L T-2 ~> m s-2].
  amrex::MultiFab PFu_, PFv_;      ///< Pressure gradient [L T-2 ~> m s-2].
  amrex::MultiFab diffu_, diffv_;  ///< Horizontal viscosity [L T-2 ~> m s-2].
  amrex::MultiFab up_, vp_;        ///< The predicted velocities [L T-1 ~> m s-1].
  amrex::MultiFab hp_, h_av_;      ///< The predicted and averaged thicknesses [H ~> m].
  amrex::MultiFab uh_, vh_;        ///< The thickness fluxes [H L2 T-1 ~> m3 s-1].
};

} // namespace MOM
