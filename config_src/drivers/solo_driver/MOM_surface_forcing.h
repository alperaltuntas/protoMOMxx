#pragma once
/// @file MOM_surface_forcing.h
/// @brief The solo driver's surface forcing: the analytic wind-stress
///        configurations that stand in for a coupler. The analogue of MOM6's
///        solo-driver MOM_surface_forcing
///        (config_src/drivers/solo_driver/MOM_surface_forcing.F90).
///
/// This lives with the driver, not in src/core, because in a coupled run the
/// coupler supplies the same fields; only the solo driver invents them.

#include <string>

#include <AMReX_REAL.H>

#include "MOM_file_parser.h"
#include "MOM_forcing_type.h"
#include "MOM_grid.h"

namespace MOM {

/// @class SurfaceForcing
/// @brief Reads the wind configuration once and writes the mechanical forcing
/// fields on demand.
class SurfaceForcing {
public:
  /// @brief Read the surface forcing configuration (WIND_CONFIG and its
  /// parameters). The analogue of MOM6's surface_forcing_init.
  /// @param grid The horizontal grid, whose latitudes shape the analytic winds.
  /// @param params Runtime parameters.
  /// @throws logger::FatalError on an unsupported configuration.
  SurfaceForcing(const Grid &grid, RuntimeParams &params);

  /// @brief Write the wind stresses for the given time into the forcing
  /// object. The analogue of MOM6's set_forcing.
  /// @param forces The mechanical forcing to fill.
  /// @param time The time of the fluxes [T ~> s]. Unused while the winds are
  ///        steady (VARIABLE_WINDS false).
  void set_forcing(MechForcing &forces, amrex::Real time) const;

private:
  const Grid &grid_;             ///< The horizontal grid (not owned).
  std::string wind_config_;      ///< The selected wind configuration.
  amrex::Real taux_mag_ = 0.1;   ///< Peak zonal wind stress [R L Z T-2 ~> Pa].
  bool variable_winds_ = true;   ///< Whether the winds vary in time.
};

} // namespace MOM
