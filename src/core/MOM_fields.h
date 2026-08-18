#pragma once
/// @file MOM_fields.h
/// @brief Factories for the model's distributed fields, naming the vertical
///        residency (2-D, layer, interface) instead of passing a level count.
///        In MOM6 the same information is the SZK_(GV) / SZK_(GV)+1 choice in
///        an array declaration.

#include <AMReX_MultiFab.H>

#include "MOM_domain_infra.h"
#include "MOM_vertical_grid.h"

namespace MOM {

/// @brief Create a 2-D field on the domain's decomposition.
/// @param domain The computational domain.
/// @param stagger Where the field's values sit within a grid cell.
/// @param ncomp Number of field components.
/// @return The newly created field, with the domain's halo widths.
inline amrex::MultiFab make_2d_field(const Domain &domain, const Stagger stagger,
                                     const int ncomp = 1) {
  return domain.make_field(stagger, 1, ncomp);
}

/// @brief Create a 3-D field with one value per layer (nk levels).
/// @param domain The computational domain.
/// @param vgrid The vertical grid, which supplies the layer count.
/// @param stagger Where the field's values sit within a grid cell.
/// @param ncomp Number of field components.
/// @return The newly created field, with the domain's halo widths.
inline amrex::MultiFab make_layer_field(const Domain &domain, const VerticalGrid &vgrid,
                                        const Stagger stagger, const int ncomp = 1) {
  return domain.make_field(stagger, vgrid.nk(), ncomp);
}

/// @brief Create a 3-D field with one value per interface (nk+1 levels).
/// @param domain The computational domain.
/// @param vgrid The vertical grid, which supplies the layer count.
/// @param stagger Where the field's values sit within a grid cell.
/// @param ncomp Number of field components.
/// @return The newly created field, with the domain's halo widths.
inline amrex::MultiFab make_interface_field(const Domain &domain, const VerticalGrid &vgrid,
                                            const Stagger stagger, const int ncomp = 1) {
  return domain.make_field(stagger, vgrid.nk() + 1, ncomp);
}

} // namespace MOM
