#pragma once
/// @file MOM_state_fields.h
/// @brief The construction-phase counterpart of the prognostic state: the
///        values initialize_state needs, and the struct of fields it produces
///        for the State constructor.

#include <AMReX_MultiFab.H>

namespace MOM {

/// @brief What the state initialization needs from the rest of the model.
///
/// The values are passed rather than the Grid and VerticalGrid objects they
/// come from: those live in src/core, and the per-directory libraries run one
/// way (framework <- initialization <- core), so src/initialization cannot see
/// them. The bathymetry is passed alongside as a field.
struct StateSpec {
  int nk = 0;                    ///< The number of layers.
  amrex::Real max_depth = 0.0;   ///< The maximum depth of the ocean [Z ~> m].
  amrex::Real angstrom = 0.0;    ///< The minimum layer thickness [H ~> m].
};

/// @brief The prognostic fields as produced by the state initialization,
/// before the State takes ownership of them. The analogue of GridFields for
/// the grid.
struct StateFields {
  amrex::MultiFab h;  ///< Layer thickness at h points [H ~> m].
  amrex::MultiFab u;  ///< Zonal velocity at u points [L T-1 ~> m s-1].
  amrex::MultiFab v;  ///< Meridional velocity at v points [L T-1 ~> m s-1].
};

} // namespace MOM
