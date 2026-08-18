#pragma once
/// @file MOM_shared_initialization.h
/// @brief Initialization code shared between configurations: currently the
///        planetary rotation, with the topography helpers to follow. The
///        analogue of MOM6's MOM_shared_initialization.F90.

#include <string>

#include <AMReX_MultiFab.H>

#include "MOM_domain_infra.h"
#include "MOM_grid_fields.h"

namespace MOM {

/// @brief Create the Coriolis parameter field at q points on the domain's
/// decomposition and compute f = 2 OMEGA sin(latitude) from the given q-point
/// latitudes, over the full grown boxes (halos included). The analogue of
/// MOM6's set_rotation_planetary.
/// @param domain The computational domain the field is created on.
/// @param spec The grid specification (the rotation rate).
/// @param geoLatBu The geographic latitude at q points [degrees_N].
/// @return The computed Coriolis parameter field [T-1 ~> s-1].
amrex::MultiFab planetary_rotation(const Domain &domain, const GridSpec &spec,
                                   const amrex::MultiFab &geoLatBu);

/// @brief Create the bottom depth field at h points and fill it with one of
/// the named analytic topographies. The analogue of MOM6's
/// initialize_topography_named followed by limit_topography.
///
/// The shape is evaluated on the valid boxes only, as in MOM6, and the halos
/// are then filled: points inside the global domain by a halo exchange, and
/// points outside it by the clamped zero MOM6 leaves there (bathyT is
/// allocated with source=0.0 and limit_topography runs over the data domain
/// before the exchange). Those outside points become land when the masks are
/// set, which is what closes a non-reentrant basin.
/// @param domain The computational domain the field is created on.
/// @param config The named topography: "flat" or "spoon".
/// @param grid_spec The grid specification (the geographic extents).
/// @param topo_spec The topography specification (the depths).
/// @param geoLonT The geographic longitude at h points [degrees_E].
/// @param geoLatT The geographic latitude at h points [degrees_N].
/// @return The computed bottom depth field [Z ~> m].
/// @throws logger::FatalError on an unimplemented named topography.
amrex::MultiFab named_topography(const Domain &domain, const std::string &config,
                                 const GridSpec &grid_spec, const TopoSpec &topo_spec,
                                 const amrex::MultiFab &geoLonT,
                                 const amrex::MultiFab &geoLatT);

} // namespace MOM
