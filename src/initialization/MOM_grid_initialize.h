#pragma once
/// @file MOM_grid_initialize.h
/// @brief Computation of the horizontal grid metrics for the supported
///        GRID_CONFIG options. The analogue of MOM6's MOM_grid_initialize.F90
///        (set_grid_metrics_*); the land/sea mask initialization
///        (initialize_masks) lands here with the topography PR.

#include "MOM_domain_infra.h"
#include "MOM_grid_fields.h"

namespace MOM {

/// @brief Create the metric fields of a simple spherical grid
/// (GRID_CONFIG = "spherical") on the domain's decomposition and compute the
/// geographic locations, grid spacings, and cell areas, over the full grown
/// boxes (halos included). The analogue of MOM6's set_grid_metrics_spherical.
/// @param domain The computational domain the fields are created on.
/// @param spec The grid specification.
/// @return The computed grid fields (the Coriolis parameter is left to
///         planetary_rotation).
GridFields spherical_grid_fields(const Domain &domain, const GridSpec &spec);

/// @brief Create the land/sea masks at the h/u/v/q points from the bottom
/// depth: an h point is ocean where its depth exceeds the masking depth, and
/// the u/v/q masks are products of the surrounding h masks. The analogue of
/// MOM6's initialize_masks.
///
/// MOM6's MASKING_DEPTH branch is not carried: only its default path, where
/// MINIMUM_DEPTH is the land/sea threshold, is implemented.
/// @param domain The computational domain the fields are created on.
/// @param topo_spec The topography specification (the masking depth).
/// @param bathyT The bottom depth at h points, halos filled [Z ~> m].
/// @param fields The grid fields; the four mask fields are created and filled.
void initialize_masks(const Domain &domain, const TopoSpec &topo_spec,
                      const amrex::MultiFab &bathyT, GridFields &fields);

/// @brief Compute the reciprocal spacings, the q-cell area, and the masked
/// face lengths and u/v-cell areas. The analogue of MOM6's
/// set_derived_dyn_horgrid plus the block at the end of initialize_masks that
/// sets dy_Cu / dx_Cv / areaCu / areaCv.
///
/// Every reciprocal is an Adcroft reciprocal: 1/x where x is nonzero and 0
/// where it is zero, so that a zero length or area propagates as a zero
/// rather than an infinity.
/// @param domain The computational domain the fields are created on.
/// @param fields The grid fields; the derived fields are created and filled.
void set_derived_metrics(const Domain &domain, GridFields &fields);

} // namespace MOM
