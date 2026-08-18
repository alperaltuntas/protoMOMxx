#pragma once
/// @file MOM_sum_output.h
/// @brief The globally summed diagnostics written to ocean.stats. The
///        analogue of MOM6's MOM_sum_output (write_energy and the depth list
///        it needs).

#include <string>
#include <vector>

#include <AMReX_MultiFab.H>

#include "MOM_coms.h"
#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_grid.h"
#include "MOM_state.h"
#include "MOM_vertical_grid.h"

namespace MOM {

/// @class SumOutput
/// @brief Writes the run's total energy, mass and other globally summed
/// quantities to an ASCII file, one line per save interval.
///
/// The sums go through MOM::reproducing_sum, so the numbers do not depend on
/// the decomposition. The available potential energy is measured against the
/// depth at which each layer's volume would sit if the ocean were at rest,
/// which needs the domain's hypsometry: the sorted list of bottom depths with
/// the area and the volume below each of them, built once at construction.
///
/// The file is byte-compatible with MOM6's ocean.stats. The NetCDF companion
/// (ocean.stats.nc) is not written; it waits on the I/O layer.
class SumOutput {
public:
  /// @brief Read the diagnostic parameters and build the depth list. The
  /// analogue of MOM6's MOM_sum_output_init.
  /// @param params Runtime parameters.
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid.
  /// @param directory The directory the energy file is written to.
  /// @throws logger::FatalError on an unsupported configuration.
  SumOutput(RuntimeParams &params, const Domain &domain, const Grid &grid,
            const VerticalGrid &vgrid, const std::string &directory);

  /// @brief Write one line of globally summed diagnostics, if the save
  /// interval has come round. The analogue of MOM6's write_energy.
  /// @param state The prognostic state.
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param vgrid The vertical grid.
  /// @param time The time since the start of the run [T ~> s].
  /// @param n_steps The number of dynamics steps taken so far.
  /// @param dt_forcing The forcing interval [T ~> s], which sets how close to
  ///        the save time a call has to be to count as being at it.
  /// @throws logger::FatalError if the energy is not finite, if it exceeds
  ///         MAX_ENERGY, or if the velocity was truncated more than MAXTRUNC
  ///         times since the last write.
  void write_energy(const State &state, const Domain &domain, const Grid &grid,
                    const VerticalGrid &vgrid, amrex::Real time, int n_steps,
                    amrex::Real dt_forcing);

  /// @brief Record velocity truncations, to be reported at the next write.
  /// @param n The number of truncations to add.
  void add_truncations(const int n) { ntrunc_ += n; }

  /// @brief The depths in the hypsometry list, deepest first [Z ~> m].
  /// @return The depth list, indexed from 1 as in MOM6.
  const std::vector<amrex::Real> &depth_list() const { return depth_; }

private:
  bool do_APE_calc_ = true;        ///< Calculate the available potential energy.
  amrex::Real dt_ = 0.0;           ///< The dynamics timestep the CFL is measured with [T ~> s].
  int maxtrunc_ = 0;               ///< Truncations tolerated between writes.
  amrex::Real max_energy_ = 0.0;   ///< The energy per unit mass that stops the run [L2 T-2 ~> m2 s-2].
  amrex::Real timeunit_ = 86400.0; ///< The time unit the Day column is in [s].
  amrex::Real energysavedays_ = 86400.0;  ///< The interval between writes [T ~> s].
  std::string energyfile_;         ///< The path of the ASCII output file.

  // The hypsometry, indexed from 1 so the search reads like MOM6's.
  std::vector<amrex::Real> depth_;      ///< Bottom depth of each entry [Z ~> m].
  std::vector<amrex::Real> area_;       ///< Open area at that depth [L2 ~> m2].
  std::vector<amrex::Real> vol_below_;  ///< Open volume below it [Z L2 ~> m3].
  int listsize_ = 0;                    ///< The number of entries.
  std::vector<int> lH_;                 ///< Per layer, the entry the last search landed on.

  int ntrunc_ = 0;                 ///< Velocity truncations since the last write.
  int previous_calls_ = 0;         ///< Writes made so far.
  amrex::Real next_write_time_ = 0.0;  ///< The time of the next write [T ~> s].
  EFP mass_prev_;                  ///< The total mass at the last write [R Z L2 ~> kg].

  /// @brief Build the sorted list of bottom depths with the area and the
  /// volume below each. The analogue of MOM6's create_depth_list.
  /// @param domain The computational domain.
  /// @param grid The horizontal grid.
  /// @param min_depth_inc The smallest depth difference that earns an entry.
  void create_depth_list(const Domain &domain, const Grid &grid, amrex::Real min_depth_inc);
};

} // namespace MOM
