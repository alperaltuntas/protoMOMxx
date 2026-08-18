#pragma once
/// @file MOM_forcing_type.h
/// @brief The mechanical forcing applied to the ocean surface. The analogue of
///        MOM6's mech_forcing type (MOM_forcing_type.F90); the thermodynamic
///        forcing type of the same module is not carried yet, because the
///        driving testcase runs with ENABLE_THERMODYNAMICS false.

#include <AMReX_MultiFab.H>

#include "MOM_domain_infra.h"

namespace MOM {

/// @class MechForcing
/// @brief The mechanical forcing fields handed to the model each forcing
/// interval: the surface wind stress components on the C-grid velocity points.
///
/// The object is created once and rewritten each interval, so its accessors
/// are non-const, like State's. The driver's SurfaceForcing fills it; the
/// dynamics reads it.
class MechForcing {
public:
  /// @brief Allocate the forcing fields on the domain's decomposition and set
  /// them to zero.
  /// @param domain The computational domain the fields are created on.
  /// @pre The infrastructure layer (MOM::Infra) is initialized.
  explicit MechForcing(const Domain &domain);

  /// @brief The zonal wind stress at u points [R L Z T-2 ~> Pa].
  /// @return The zonal wind stress field.
  amrex::MultiFab &taux() { return taux_; }
  /// @brief The zonal wind stress at u points [R L Z T-2 ~> Pa].
  /// @return The zonal wind stress field.
  const amrex::MultiFab &taux() const { return taux_; }

  /// @brief The meridional wind stress at v points [R L Z T-2 ~> Pa].
  /// @return The meridional wind stress field.
  amrex::MultiFab &tauy() { return tauy_; }
  /// @brief The meridional wind stress at v points [R L Z T-2 ~> Pa].
  /// @return The meridional wind stress field.
  const amrex::MultiFab &tauy() const { return tauy_; }

private:
  amrex::MultiFab taux_;  ///< Zonal wind stress at u points [R L Z T-2 ~> Pa].
  amrex::MultiFab tauy_;  ///< Meridional wind stress at v points [R L Z T-2 ~> Pa].

  // defer: the rest of mech_forcing -- ustar and its gustiness, the surface
  //        pressure and its derivatives, the ice-shelf and wave fields, and
  //        the accumulated-forcing bookkeeping.
};

} // namespace MOM
