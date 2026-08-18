#pragma once
/// @file MOM_state_initialization.h
/// @brief Runtime-parameter-driven initialization of the prognostic state.
///        The analogue of MOM6's MOM_initialize_state
///        (MOM_state_initialization.F90).

#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_state_fields.h"

namespace MOM {

/// @brief Read the state configuration (THICKNESS_CONFIG, VELOCITY_CONFIG),
/// create the prognostic fields, and fill them. The analogue of MOM6's
/// MOM_initialize_state.
/// @param domain The computational domain the fields are created on.
/// @param spec What the initialization needs from the grids; see StateSpec.
/// @param bathyT The ocean bottom depth at h points [Z ~> m].
/// @param params Runtime parameters.
/// @return The initialized prognostic fields, ready for the State constructor.
/// @throws logger::FatalError on an unsupported configuration.
StateFields initialize_state(const Domain &domain, const StateSpec &spec,
                             const amrex::MultiFab &bathyT, RuntimeParams &params);

} // namespace MOM
