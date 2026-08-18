#pragma once
/// @file MOM_state_initialization.h
/// @brief Runtime-parameter-driven initialization of the prognostic state.
///        The analogue of MOM6's MOM_initialize_state
///        (MOM_state_initialization.F90).

#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_grid.h"
#include "MOM_state_fields.h"
#include "MOM_vertical_grid.h"

namespace MOM {

/// @brief Read the state configuration (THICKNESS_CONFIG, VELOCITY_CONFIG),
/// create the prognostic fields, and fill them. The analogue of MOM6's
/// MOM_initialize_state.
/// @param domain The computational domain the fields are created on.
/// @param grid The horizontal grid, which supplies the bathymetry and the
///        maximum depth.
/// @param vgrid The vertical grid, which supplies the layer count and the
///        minimum layer thickness.
/// @param params Runtime parameters.
/// @return The initialized prognostic fields, ready for the State constructor.
/// @throws logger::FatalError on an unsupported configuration.
StateFields initialize_state(const Domain &domain, const Grid &grid,
                             const VerticalGrid &vgrid, RuntimeParams &params);

} // namespace MOM
