#pragma once
/// @file MOM_grid.h
/// @brief The horizontal ocean grid of a model instance: the geographic
///        locations, grid spacings, and cell areas at the h/q/u/v points of
///        the Arakawa C-grid, and the Coriolis parameter. The analogue of
///        MOM6's MOM_grid (ocean_grid_type). The field values are computed
///        in src/initialization and handed to the constructor as a
///        GridFields struct (MOM_grid_fields.h).

#include <AMReX_MultiFab.H>

#include "MOM_grid_fields.h"

namespace MOM {

/// @class Grid
/// @brief The horizontal grid: metric fields at the four C-grid point types
/// (h: cell centers, q: cell corners, u: east faces, v: north faces) and the
/// Coriolis parameter at q points.
///
/// GridFields is the construction-phase counterpart of this class. The
/// field values are computed in src/initialization (where the GRID_CONFIG
/// options are handled), but the grid itself should be immutable (nearly
/// every part of the model reads the grid throughout the run, and none of
/// them should be able to alter it), so the computed fields are handed to
/// the constructor as a plain struct, and the constructor validates them
/// and takes ownership. MOM6 does the same job
/// with two full grid types: it fills dyn_horgrid_type, copies it into
/// ocean_grid_type, and both remain mutable. Here the precursor (GridFields)
/// is a bare struct with no behavior of its own, the handoff is a move
/// rather than a copy, and the resulting grid is read-only. So a Grid
/// is complete by construction and cannot be altered afterwards.
class Grid {
public:
  /// @brief Construct the grid: check that every field is created and take
  /// ownership.
  /// @param fields The computed grid fields; moved from.
  /// @pre The infrastructure layer (MOM::Infra) is initialized.
  explicit Grid(GridFields &&fields);

  /// @brief The southern latitude of the domain [degrees_N].
  /// @return The southern latitude.
  amrex::Real south_lat() const { return fields_.south_lat; }

  /// @brief The latitudinal length of the domain [degrees_N].
  /// @return The latitudinal length.
  amrex::Real len_lat() const { return fields_.len_lat; }

  /// @brief The western longitude of the domain [degrees_E].
  /// @return The western longitude.
  amrex::Real west_lon() const { return fields_.west_lon; }

  /// @brief The longitudinal length of the domain [degrees_E].
  /// @return The longitudinal length.
  amrex::Real len_lon() const { return fields_.len_lon; }

  /// @brief The radius of the Earth [L ~> m].
  /// @return The Earth radius.
  amrex::Real rad_earth() const { return fields_.rad_earth; }

  /// @brief The geographic latitude at h (tracer) points [degrees_N].
  /// @return The h-point latitude field.
  const amrex::MultiFab &geoLatT() const { return fields_.geoLatT; }

  /// @brief The geographic longitude at h (tracer) points [degrees_E].
  /// @return The h-point longitude field.
  const amrex::MultiFab &geoLonT() const { return fields_.geoLonT; }

  /// @brief Delta x at h points [L ~> m].
  /// @return The h-point zonal grid spacing.
  const amrex::MultiFab &dxT() const { return fields_.dxT; }

  /// @brief Delta y at h points [L ~> m].
  /// @return The h-point meridional grid spacing.
  const amrex::MultiFab &dyT() const { return fields_.dyT; }

  /// @brief The area of an h-cell [L2 ~> m2].
  /// @return The h-cell area field.
  const amrex::MultiFab &areaT() const { return fields_.areaT; }

  /// @brief The geographic latitude at u points [degrees_N].
  /// @return The u-point latitude field.
  const amrex::MultiFab &geoLatCu() const { return fields_.geoLatCu; }

  /// @brief The geographic longitude at u points [degrees_E].
  /// @return The u-point longitude field.
  const amrex::MultiFab &geoLonCu() const { return fields_.geoLonCu; }

  /// @brief Delta x at u points [L ~> m].
  /// @return The u-point zonal grid spacing.
  const amrex::MultiFab &dxCu() const { return fields_.dxCu; }

  /// @brief Delta y at u points [L ~> m].
  /// @return The u-point meridional grid spacing.
  const amrex::MultiFab &dyCu() const { return fields_.dyCu; }

  /// @brief The geographic latitude at v points [degrees_N].
  /// @return The v-point latitude field.
  const amrex::MultiFab &geoLatCv() const { return fields_.geoLatCv; }

  /// @brief The geographic longitude at v points [degrees_E].
  /// @return The v-point longitude field.
  const amrex::MultiFab &geoLonCv() const { return fields_.geoLonCv; }

  /// @brief Delta x at v points [L ~> m].
  /// @return The v-point zonal grid spacing.
  const amrex::MultiFab &dxCv() const { return fields_.dxCv; }

  /// @brief Delta y at v points [L ~> m].
  /// @return The v-point meridional grid spacing.
  const amrex::MultiFab &dyCv() const { return fields_.dyCv; }

  /// @brief The geographic latitude at q (corner) points [degrees_N].
  /// @return The q-point latitude field.
  const amrex::MultiFab &geoLatBu() const { return fields_.geoLatBu; }

  /// @brief The geographic longitude at q (corner) points [degrees_E].
  /// @return The q-point longitude field.
  const amrex::MultiFab &geoLonBu() const { return fields_.geoLonBu; }

  /// @brief Delta x at q points [L ~> m].
  /// @return The q-point zonal grid spacing.
  const amrex::MultiFab &dxBu() const { return fields_.dxBu; }

  /// @brief Delta y at q points [L ~> m].
  /// @return The q-point meridional grid spacing.
  const amrex::MultiFab &dyBu() const { return fields_.dyBu; }

  /// @brief The Coriolis parameter at q points [T-1 ~> s-1].
  /// @return The q-point Coriolis parameter field.
  const amrex::MultiFab &CoriolisBu() const { return fields_.CoriolisBu; }

  /// @brief The maximum depth of the ocean [Z ~> m].
  /// @return The maximum depth.
  amrex::Real max_depth() const { return fields_.max_depth; }

  /// @brief The ocean bottom depth at h points, positive down [Z ~> m].
  /// @return The bathymetry field.
  const amrex::MultiFab &bathyT() const { return fields_.bathyT; }

  /// @brief The land/sea mask at h points: 1 for ocean, 0 for land [nondim].
  /// @return The h-point mask field.
  const amrex::MultiFab &mask2dT() const { return fields_.mask2dT; }

  /// @brief The land/sea mask at u points: 1 for ocean, 0 for land [nondim].
  /// @return The u-point mask field.
  const amrex::MultiFab &mask2dCu() const { return fields_.mask2dCu; }

  /// @brief The land/sea mask at v points: 1 for ocean, 0 for land [nondim].
  /// @return The v-point mask field.
  const amrex::MultiFab &mask2dCv() const { return fields_.mask2dCv; }

  /// @brief The land/sea mask at q points: 1 for ocean, 0 for land [nondim].
  /// @return The q-point mask field.
  const amrex::MultiFab &mask2dBu() const { return fields_.mask2dBu; }

  /// @brief 1/dxT at h points [L-1 ~> m-1].
  /// @return The IdxT field.
  const amrex::MultiFab &IdxT() const { return fields_.IdxT; }

  /// @brief 1/dyT at h points [L-1 ~> m-1].
  /// @return The IdyT field.
  const amrex::MultiFab &IdyT() const { return fields_.IdyT; }

  /// @brief 1/areaT at h points [L-2 ~> m-2].
  /// @return The IareaT field.
  const amrex::MultiFab &IareaT() const { return fields_.IareaT; }

  /// @brief 1/dxCu at u points [L-1 ~> m-1].
  /// @return The IdxCu field.
  const amrex::MultiFab &IdxCu() const { return fields_.IdxCu; }

  /// @brief 1/dyCu at u points [L-1 ~> m-1].
  /// @return The IdyCu field.
  const amrex::MultiFab &IdyCu() const { return fields_.IdyCu; }

  /// @brief 1/dxCv at v points [L-1 ~> m-1].
  /// @return The IdxCv field.
  const amrex::MultiFab &IdxCv() const { return fields_.IdxCv; }

  /// @brief 1/dyCv at v points [L-1 ~> m-1].
  /// @return The IdyCv field.
  const amrex::MultiFab &IdyCv() const { return fields_.IdyCv; }

  /// @brief 1/dxBu at q points [L-1 ~> m-1].
  /// @return The IdxBu field.
  const amrex::MultiFab &IdxBu() const { return fields_.IdxBu; }

  /// @brief 1/dyBu at q points [L-1 ~> m-1].
  /// @return The IdyBu field.
  const amrex::MultiFab &IdyBu() const { return fields_.IdyBu; }

  /// @brief The area of a q-cell [L2 ~> m2].
  /// @return The areaBu field.
  const amrex::MultiFab &areaBu() const { return fields_.areaBu; }

  /// @brief 1/areaBu at q points [L-2 ~> m-2].
  /// @return The IareaBu field.
  const amrex::MultiFab &IareaBu() const { return fields_.IareaBu; }

  /// @brief The unblocked length of a u face [L ~> m].
  /// @return The dy_Cu field.
  const amrex::MultiFab &dy_Cu() const { return fields_.dy_Cu; }

  /// @brief The unblocked length of a v face [L ~> m].
  /// @return The dx_Cv field.
  const amrex::MultiFab &dx_Cv() const { return fields_.dx_Cv; }

  /// @brief The area of a u-cell [L2 ~> m2].
  /// @return The areaCu field.
  const amrex::MultiFab &areaCu() const { return fields_.areaCu; }

  /// @brief The area of a v-cell [L2 ~> m2].
  /// @return The areaCv field.
  const amrex::MultiFab &areaCv() const { return fields_.areaCv; }

  /// @brief The masked 1/areaCu [L-2 ~> m-2].
  /// @return The IareaCu field.
  const amrex::MultiFab &IareaCu() const { return fields_.IareaCu; }

  /// @brief The masked 1/areaCv [L-2 ~> m-2].
  /// @return The IareaCv field.
  const amrex::MultiFab &IareaCv() const { return fields_.IareaCv; }

  /// @brief IdxCu masked by the u-point open-boundary mask [L-1 ~> m-1].
  /// @return The IdxCu_OBCmask field.
  const amrex::MultiFab &IdxCu_OBCmask() const { return fields_.IdxCu_OBCmask; }

  /// @brief IdyCv masked by the v-point open-boundary mask [L-1 ~> m-1].
  /// @return The IdyCv_OBCmask field.
  const amrex::MultiFab &IdyCv_OBCmask() const { return fields_.IdyCv_OBCmask; }

private:
  GridFields fields_;  ///< The grid fields (owned).
};

} // namespace MOM
