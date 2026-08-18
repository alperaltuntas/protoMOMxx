#pragma once
/// @file MOM_state.h
/// @brief The prognostic state of a model instance: the layer thicknesses and
///        the horizontal velocity components. The analogue of the u/v/h
///        members of MOM6's MOM_control_struct.

#include <AMReX_MultiFab.H>

#include "MOM_state_fields.h"

namespace MOM {

/// @class State
/// @brief The prognostic ocean state: h at cell centers, u on x faces, v on y
/// faces, all with nk layers on the domain's decomposition.
///
/// Unlike Grid, State is mutable by design: the time stepping writes it every
/// step, so the accessors hand out non-const references, the way MOM6's
/// step_MOM writes CS%u, CS%v and CS%h. The initial values are computed in
/// src/initialization and moved in, as they are for the grid.
class State {
public:
  /// @brief Construct the state: check that every field is created and take
  /// ownership.
  /// @param fields The initialized prognostic fields; moved from.
  /// @pre The infrastructure layer (MOM::Infra) is initialized.
  explicit State(StateFields &&fields);

  /// @brief The layer thickness at h points [H ~> m].
  /// @return The thickness field.
  amrex::MultiFab &h() { return fields_.h; }
  /// @brief The layer thickness at h points [H ~> m].
  /// @return The thickness field.
  const amrex::MultiFab &h() const { return fields_.h; }

  /// @brief The zonal velocity at u points [L T-1 ~> m s-1].
  /// @return The zonal velocity field.
  amrex::MultiFab &u() { return fields_.u; }
  /// @brief The zonal velocity at u points [L T-1 ~> m s-1].
  /// @return The zonal velocity field.
  const amrex::MultiFab &u() const { return fields_.u; }

  /// @brief The meridional velocity at v points [L T-1 ~> m s-1].
  /// @return The meridional velocity field.
  amrex::MultiFab &v() { return fields_.v; }
  /// @brief The meridional velocity at v points [L T-1 ~> m s-1].
  /// @return The meridional velocity field.
  const amrex::MultiFab &v() const { return fields_.v; }

private:
  StateFields fields_;  ///< The prognostic fields (owned).
};

} // namespace MOM
