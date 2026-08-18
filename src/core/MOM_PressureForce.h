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
  // defer: the FV/AFV form, the atmospheric pressure term, self-attraction
  //        and loading, tides, and the pbce/eta outputs the split scheme
  //        needs.
};

} // namespace MOM
