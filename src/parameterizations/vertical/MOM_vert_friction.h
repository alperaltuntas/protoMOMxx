#pragma once
/// @file MOM_vert_friction.h
/// @brief The implicit vertical viscosity and the surface wind stress. The
///        analogue of MOM6's MOM_vert_friction (vertvisc_coef + vertvisc).

#include <AMReX_MultiFab.H>

#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_forcing_type.h"
#include "MOM_grid.h"
#include "MOM_set_viscosity.h"
#include "MOM_vertical_grid.h"

namespace MOM {

/// @class VertFriction
/// @brief Applies the vertical viscosity implicitly and adds the wind stress.
///
/// The coupling coefficients between layers are computed first
/// (vertvisc_coef), then a tridiagonal solve advances the velocities
/// (vertvisc). The coefficients depend on the thicknesses, so the two steps
/// are separate calls in MOM6 as well, with the coefficients held between
/// them.
///
/// The harmonic-mean thickness at velocity points and the direct application
/// of the wind stress over HMIX_STRESS are implemented; the law-of-the-wall
/// mixed layer, the GL90 scheme and the ice-shelf boundary layer abort or are
/// absent.
class VertFriction {
public:
  /// @brief Read the vertical friction configuration. The analogue of MOM6's
  /// vertvisc_init.
  /// @param params Runtime parameters.
  /// @param domain The computational domain.
  /// @param vgrid The vertical grid, which supplies the layer count.
  /// @throws logger::FatalError on an unsupported option.
  VertFriction(RuntimeParams &params, const Domain &domain, const VerticalGrid &vgrid);

  /// @brief Compute the coupling coefficients between layers. The analogue of
  /// MOM6's vertvisc_coef.
  /// @param u The zonal velocity [L T-1 ~> m s-1].
  /// @param v The meridional velocity [L T-1 ~> m s-1].
  /// @param h The layer thicknesses [H ~> m].
  /// @param visc The bottom boundary layer fields.
  /// @param Kv_interior The background interior viscosity [Z2 T-1 ~> m2 s-1].
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid.
  void coefficients(const amrex::MultiFab &u, const amrex::MultiFab &v,
                    const amrex::MultiFab &h, const VertVisc &visc,
                    amrex::Real Kv_interior, const Domain &domain,
                    const Grid &grid, const VerticalGrid &vgrid);

  /// @brief Advance the velocities under the vertical viscosity and the wind
  /// stress. The analogue of MOM6's vertvisc.
  /// @param u The zonal velocity, updated in place [L T-1 ~> m s-1].
  /// @param v The meridional velocity, updated in place [L T-1 ~> m s-1].
  /// @param h The layer thicknesses [H ~> m].
  /// @param forces The mechanical forcing (the wind stress).
  /// @param dt The time increment [T ~> s].
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid.
  void apply(amrex::MultiFab &u, amrex::MultiFab &v, const amrex::MultiFab &h,
             const MechForcing &forces, amrex::Real dt, const Domain &domain,
             const Grid &grid, const VerticalGrid &vgrid);

  /// @brief The number of velocity truncations since the last call, and reset
  /// the count. The analogue of MOM6's CS%ntrunc, which set_visc shares with
  /// MOM_sum_output through a pointer.
  /// @return The truncation count.
  int take_truncations() {
    const int n = ntrunc_;
    ntrunc_ = 0;
    return n;
  }

private:
  bool direct_stress_ = false;   ///< Spread the wind stress over HMIX_STRESS.
  amrex::Real Hmix_stress_ = 0.0;  ///< The depth the stress is spread over [H ~> m].
  amrex::Real Hmix_ = 0.0;         ///< The fixed mixed-layer depth [H ~> m].
  amrex::Real Kvml_invZ2_ = 0.0;   ///< The near-surface viscosity scale [Z2 T-1 ~> m2 s-1].
  amrex::Real Hbbl_ = 0.0;         ///< The nominal bottom boundary layer thickness [Z ~> m].
  amrex::Real CFL_trunc_ = 0.5;    ///< The CFL number above which velocities are truncated [nondim].
  amrex::Real vel_underflow_ = 0.0;  ///< Velocities smaller than this are set to zero [L T-1 ~> m s-1].
  int ntrunc_ = 0;                 ///< Truncations since the count was last taken.

  amrex::MultiFab a_u_;  ///< Coupling coefficient at u-point interfaces [Z T-1 ~> m s-1].
  amrex::MultiFab a_v_;  ///< The same at v-point interfaces [Z T-1 ~> m s-1].
  amrex::MultiFab h_u_;  ///< The thickness at u points used by the solve [H ~> m].
  amrex::MultiFab h_v_;  ///< The same at v points [H ~> m].

  /// @brief Truncate velocities that would exceed CFL_TRUNCATE and count the
  /// truncations. The analogue of MOM6's vertvisc_limit_vel, which vertvisc
  /// calls on its way out.
  /// @param u The zonal velocity, updated in place [L T-1 ~> m s-1].
  /// @param v The meridional velocity, updated in place [L T-1 ~> m s-1].
  /// @param dt The time increment [T ~> s].
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid.
  void limit_velocity(amrex::MultiFab &u, amrex::MultiFab &v, amrex::Real dt,
                      const Grid &grid, const VerticalGrid &vgrid);
};

} // namespace MOM
