#pragma once
/// @file MOM_vertical_grid.h
/// @brief The vertical ocean grid and the coordinate it carries: the number
///        of layers, the reduced gravities across interfaces, and the target
///        (coordinate) densities of layers. The analogue of MOM6's
///        MOM_verticalGrid (verticalGrid_type) as well as the coordinate values
///        set by MOM_coord_initialization (MOM_initialize_coord).

#include <vector>

#include <AMReX_REAL.H>

#include "MOM_file_parser.h"

namespace MOM {

// A note on unit descriptions in comments: MOM6 rescales units at runtime for
// dimensional consistency testing and annotates them like
// "[L2 Z-1 T-2 ~> m s-2]". protoMOMxx doesn't implement a unit scaling yet,
// so all values are in the MKS units on the right-hand side of "~>". The comments
// keep MOM6's notation so the annotations map one-to-one onto MOM6's.

/// @class VerticalGrid
/// @brief The vertical ocean grid of a model instance: the number of layers
/// (nk), the reduced gravity across each interface (g_prime), and the target
/// coordinate value (potential density) of each layer (Rlay), along with
/// the scalar constants they were built from.
class VerticalGrid {
public:
  /// @brief Construct the vertical grid from runtime parameters: read the
  /// vertical grid parameters (NK, G_EARTH, RHO_0, BOUSSINESQ, MAXIMUM_DEPTH) and the
  /// coordinate configuration (COORD_CONFIG, GFS, GINT, ...), and set up the
  /// coordinate values. The analogue of MOM6's verticalGridInit
  /// (MOM_verticalGrid.F90) + MOM_initialize_coord (MOM_coord_initialization.F90).
  /// @param params Runtime parameters.
  /// @throws logger::FatalError on an unsupported or an invalid configuration.
  explicit VerticalGrid(RuntimeParams &params);

  /// @brief Number of layers in the vertical (MOM6's GV%ke; parameter NK).
  /// @return The layer count.
  int nk() const { return nk_; }

  /// @brief Maximum depth of the ocean [Z ~> m].
  /// @return The maximum depth.
  amrex::Real max_depth() const { return max_depth_; }

  /// @brief Gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  /// @return The gravitational acceleration.
  amrex::Real g_Earth() const { return g_Earth_; }

  /// @brief Mean ocean density used in the Boussinesq approximation
  /// [R ~> kg m-3].
  /// @return The Boussinesq reference density.
  amrex::Real Rho0() const { return Rho0_; }

  /// @brief Whether the Boussinesq approximation is made.
  /// @return True if the Boussinesq approximation is made.
  bool Boussinesq() const { return Boussinesq_; }

  /// @brief Reduced gravity across each interface [L2 Z-1 T-2 ~> m s-2].
  /// The bottom value (index nk) does not matter physically and is set only
  /// to avoid an uninitialized value in output, as in MOM6.
  /// @return The nk+1 interface reduced gravities, surface first.
  const std::vector<amrex::Real> &g_prime() const { return g_prime_; }

  /// @brief Target coordinate value (potential density) of each layer
  /// [R ~> kg m-3].
  /// @return The nk layer target densities, surface first.
  const std::vector<amrex::Real> &Rlay() const { return Rlay_; }

private:
  int nk_ = 0;                        ///< Number of layers in the vertical.
  amrex::Real max_depth_ = 0.0;       ///< Maximum ocean depth [Z ~> m].
  amrex::Real g_Earth_ = 0.0;         ///< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  amrex::Real Rho0_ = 0.0;            ///< Boussinesq reference density [R ~> kg m-3].
  bool Boussinesq_ = true;            ///< Whether the Boussinesq approximation is made.
  std::vector<amrex::Real> g_prime_;  ///< Interface reduced gravities (nk+1) [L2 Z-1 T-2 ~> m s-2].
  std::vector<amrex::Real> Rlay_;     ///< Layer target densities (nk) [R ~> kg m-3].
};

} // namespace MOM
