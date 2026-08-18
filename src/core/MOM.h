#pragma once
/// @file MOM.h
/// @brief Main header for the Modular Ocean Model (MOM) core.

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_grid.h"
#include "MOM_vertical_grid.h"

namespace MOM {

/// @brief Scales input to be within necessary bounds
/// @param x Global index
/// @param x_min x min
/// @param x_max x max
/// @param xi_min grid x min
/// @param xi_max grid x max
/// @return Index within scaled space for local grid
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real LinearMapCoordinates(const amrex::Real x,
                                 const amrex::Real x_min, const amrex::Real x_max,
                                 const amrex::Real xi_min, const amrex::Real xi_max);

/// @brief The Model class is the main interface of the MOM core: it owns the
/// model state and subsystems and provides the entry points the driver calls.
///
/// The constructor performs the full initialization. It is the analogue of
/// MOM6's initialize_MOM (src/core/MOM.F90). The constructor is decomposed
/// into phases that mirror the topology of MOM6's initialize_MOM where each
/// remaining phase is currently a stub that will be filled in by upcoming PRs
/// (VerticalGrid, State, Dynamics). The initialize_state stub runs
/// the original psi (stream function) demo, which exercises the AMReX machinery
/// its real replacement will use.
class Model {
public:
  /// @brief Scalar configuration switches of the model (the analogue of the
  /// scalar members of MOM6's MOM_control_struct).
  struct Config {
    bool split = true;            ///< Use split time stepping.
    bool split_rk4 = false;       ///< Use the RK4 variant of the split scheme.
    bool use_RK2 = false;         ///< Use RK2 (not RK3) in unsplit stepping.
    bool fpmix = false;           ///< Use the FPMIX algorithm.
    bool debug = false;           ///< Write verbose debugging data.
  };

  /// @brief Initialize the model (the analogue of MOM6's initialize_MOM).
  /// @param params Runtime parameters. Injected by the driver.
  /// @pre The infrastructure layer (AMReX) is initialized.
  explicit Model(RuntimeParams &params);

  /// @brief Read-only access to the model configuration switches.
  /// @return Const reference to the model configuration switches.
  const Config &config() const { return config_; }

  /// @brief Read-only access to the computational domain.
  /// @return Const reference to the model's domain.
  const Domain &domain() const { return domain_; }

  /// @brief Read-only access to the horizontal grid.
  /// @return Const reference to the model's horizontal grid.
  const Grid &grid() const { return grid_; }

  /// @brief Read-only access to the vertical grid.
  /// @return Const reference to the model's vertical grid.
  const VerticalGrid &vertical_grid() const { return vgrid_; }

private:
  // config_ initialization must precede domain_: its initializer sets the log 
  // verbosity in effect for the later initializers' messages.
  Config config_;

  /// @brief The computational domain: global extents, connectivity, halo
  /// metadata, and the horizontal decomposition.
  Domain domain_;

  /// @brief The horizontal grid: metric fields at the h/q/u/v points and the
  /// Coriolis parameter, on domain_'s decomposition.
  Grid grid_;

  /// @brief The vertical grid: the layer count, the interface reduced
  /// gravities, and the layer target densities.
  VerticalGrid vgrid_;

  /// @brief Read the scalar configuration switches into a Config object.
  static Config read_config_switches(RuntimeParams &params);

  /// @brief Initialize the prognostic state (u, v, h, ...). Analogue of
  /// MOM6's MOM_initialize_state. (stub -- currently runs the original psi
  /// demo to exercise the AMReX machinery.)
  void initialize_state(RuntimeParams &params);

  /// @brief Initialize the dynamics subsystem for the configured time
  /// stepping scheme. Analogue of MOM6's register_restarts_dyn_* +
  /// initialize_dyn_* four-way dispatch. (stub)
  void initialize_dynamics(RuntimeParams &params);

  /// @brief tmp: fill the psi (stream function) demo field on the domain's
  /// decomposition. Retires with the demo when the real State arrives.
  void fill_psi_demo(amrex::MultiFab &psi) const;
};

} // namespace MOM
