#pragma once
/// @file MOM_state_fields.h
/// @brief The construction-phase counterpart of the prognostic state: the
///        struct of fields initialize_state produces for the State
///        constructor.

#include <AMReX_MultiFab.H>

namespace MOM {

/// @brief The prognostic fields as produced by the state initialization,
/// before the State takes ownership of them. The analogue of GridFields for
/// the grid.
struct StateFields {
  amrex::MultiFab h;  ///< Layer thickness at h points [H ~> m].
  amrex::MultiFab u;  ///< Zonal velocity at u points [L T-1 ~> m s-1].
  amrex::MultiFab v;  ///< Meridional velocity at v points [L T-1 ~> m s-1].
};

} // namespace MOM
