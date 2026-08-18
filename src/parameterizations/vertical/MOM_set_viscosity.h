#pragma once
/// @file MOM_set_viscosity.h
/// @brief The bottom boundary layer thickness and viscosity. The analogue of
///        MOM6's MOM_set_viscosity (set_viscous_BBL).

#include <AMReX_MultiFab.H>

#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_grid.h"
#include "MOM_vertical_grid.h"

namespace MOM {

/// @class VertVisc
/// @brief The vertical viscosity fields the dynamics passes between
/// set_viscous_BBL and the vertical friction. The analogue of MOM6's
/// vertvisc_type, reduced to what the driving testcase uses.
struct VertVisc {
  amrex::MultiFab Kv_bbl_u;      ///< Bottom boundary layer viscosity at u points [Z2 T-1 ~> m2 s-1].
  amrex::MultiFab Kv_bbl_v;      ///< The same at v points [Z2 T-1 ~> m2 s-1].
  amrex::MultiFab bbl_thick_u;   ///< Bottom boundary layer thickness at u points [Z ~> m].
  amrex::MultiFab bbl_thick_v;   ///< The same at v points [Z ~> m].

  /// @brief Allocate the fields on the domain's decomposition.
  /// @param domain The computational domain.
  explicit VertVisc(const Domain &domain);
};

/// @class SetViscosity
/// @brief Computes the bottom boundary layer thickness and viscosity from the
/// bottom drag law. The analogue of MOM6's set_visc_CS + set_viscous_BBL.
///
/// Only the LINEAR_DRAG path is implemented: the bottom friction velocity is
/// the constant CDRAG^(1/2) * DRAG_BG_VEL rather than being derived from the
/// near-bottom flow, which is what makes the drag linear in the velocity.
class SetViscosity {
public:
  /// @brief Read the viscosity configuration. The analogue of MOM6's
  /// set_visc_init.
  /// @param params Runtime parameters.
  /// @throws logger::FatalError on an unsupported option.
  explicit SetViscosity(RuntimeParams &params);

  /// @brief The background interior viscosity [Z2 T-1 ~> m2 s-1].
  /// @return KV.
  amrex::Real Kv() const { return Kv_; }

  /// @brief Compute the bottom boundary layer thickness and viscosity. The
  /// analogue of MOM6's set_viscous_BBL.
  /// @param visc The viscosity fields to fill.
  /// @param u The zonal velocity [L T-1 ~> m s-1].
  /// @param v The meridional velocity [L T-1 ~> m s-1].
  /// @param h The layer thicknesses [H ~> m].
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid (the layer target densities).
  void set_viscous_BBL(VertVisc &visc, const amrex::MultiFab &u,
                       const amrex::MultiFab &v, const amrex::MultiFab &h,
                       const Domain &domain, const Grid &grid,
                       const VerticalGrid &vgrid) const;

private:
  amrex::Real cdrag_ = 0.003;       ///< The drag coefficient [nondim].
  amrex::Real drag_bg_vel_ = 0.0;   ///< The assumed bottom velocity [L T-1 ~> m s-1].
  amrex::Real BBL_thick_min_ = 0.0; ///< The minimum bottom boundary layer thickness [Z ~> m].
  amrex::Real Kv_ = 0.0;            ///< The background interior viscosity [Z2 T-1 ~> m2 s-1].
  amrex::Real Kv_BBL_min_ = 0.0;    ///< The minimum bottom boundary layer viscosity [Z2 T-1 ~> m2 s-1].
  amrex::Real Hbbl_ = 0.0;          ///< The nominal bottom boundary layer thickness [Z ~> m].

  // defer: the non-linear drag path (the near-bottom velocity average), the
  //        channel drag, the dynamic viscous mixed layer, the ice-shelf
  //        boundary layer, the equation-of-state density variables, and the
  //        body-force drag.
};

} // namespace MOM
