#pragma once
/// @file MOM_fixed_initialization.h
/// @brief Runtime-parameter-driven construction of the fixed (time-invariant)
///        aspects of the model: the horizontal grid metrics and the planetary
///        rotation, and later the topography and land/sea masks. The analogue
///        of MOM6's MOM_initialize_fixed (MOM_fixed_initialization.F90) and
///        the parameter-reading halves of MOM_grid_initialize.F90 /
///        MOM_shared_initialization.F90.

#include "MOM_domain_infra.h"
#include "MOM_file_parser.h"
#include "MOM_grid.h"

namespace MOM {

/// @brief Read the grid parameters (GRID_CONFIG and its configuration-specific
/// extents, ROTATION and its rate), compute the grid fields on the domain's
/// decomposition, and construct the Grid from them. The analogue of MOM6's
/// MOM_initialize_fixed.
/// @param domain The computational domain the grid fields are created on.
/// @param params Runtime parameters.
/// @return The constructed Grid.
/// @pre The infrastructure layer (MOM::Infra) is initialized.
/// @throws logger::FatalError on an unsupported or an invalid configuration.
Grid make_grid(const Domain &domain, RuntimeParams &params);

} // namespace MOM
