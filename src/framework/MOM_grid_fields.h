#pragma once
/// @file MOM_grid_fields.h
/// @brief The construction-phase counterpart of the horizontal grid: the
///        grid specification read from the runtime parameters, and the
///        struct of grid fields that src/initialization computes and the
///        Grid constructor (src/core) takes over. The analogue of MOM6's
///        MOM_dyn_horgrid (dyn_horgrid_type).

#include <AMReX_MultiFab.H>

namespace MOM {

// A note on unit descriptions in comments: MOM6 rescales units at runtime for
// dimensional consistency testing and annotates them like "[L ~> m]".
// protoMOMxx doesn't implement a unit scaling yet, so all values are in the
// MKS units on the right-hand side of "~>".

/// @brief The construction specification of a horizontal grid: the geographic
/// extents of a simple spherical grid (GRID_CONFIG = "spherical") and the
/// planetary rotation rate. Other grid configurations (mosaic, cartesian,
/// mercator) will extend this specification when they are implemented.
struct GridSpec {
  amrex::Real south_lat = 0.0;      ///< The southern latitude of the domain [degrees_N].
  amrex::Real len_lat = 0.0;        ///< The latitudinal length of the domain [degrees_N].
  amrex::Real west_lon = 0.0;       ///< The western longitude of the domain [degrees_E].
  amrex::Real len_lon = 0.0;        ///< The longitudinal length of the domain [degrees_E].
  amrex::Real rad_earth = 6.378e6;  ///< The radius of the Earth [L ~> m].
  amrex::Real omega = 7.2921e-5;    ///< The rotation rate of the Earth [T-1 ~> s-1].
};

/// @brief The construction specification of the bottom topography: the named
/// analytic shape (TOPO_CONFIG) and the depths that bound and shape it. The
/// masking depth of MOM6's MASKING_DEPTH is not carried: protoMOMxx only
/// implements MOM6's default path, where MINIMUM_DEPTH is both the clamp and
/// the land/sea threshold.
struct TopoSpec {
  amrex::Real max_depth = 0.0;        ///< The maximum depth of the ocean [Z ~> m].
  amrex::Real min_depth = 0.0;        ///< The minimum depth of the ocean [Z ~> m].
  amrex::Real edge_depth = 100.0;     ///< The depth at the edge of a named topography [Z ~> m].
  amrex::Real expdecay = 400000.0;    ///< The decay scale of the sloping boundaries [L ~> m].
};

/// @brief The grid fields (the metrics at the four C-grid point types and the
/// Coriolis parameter) and the geographic extents they were computed from.
/// The analogue of MOM6's dyn_horgrid_type.
///
/// This is the construction-phase counterpart of Grid (src/core): a plain
/// struct with no behavior of its own, so the setup functions in
/// src/initialization can fill it freely. Each field is empty until its
/// setup function creates and computes it. Once complete, the struct is
/// moved into the Grid constructor, which checks that every field is
/// created and becomes the read-only owner.
struct GridFields {
  amrex::Real south_lat = 0.0;  ///< The southern latitude of the domain [degrees_N].
  amrex::Real len_lat = 0.0;    ///< The latitudinal length of the domain [degrees_N].
  amrex::Real west_lon = 0.0;   ///< The western longitude of the domain [degrees_E].
  amrex::Real len_lon = 0.0;    ///< The longitudinal length of the domain [degrees_E].
  amrex::Real rad_earth = 0.0;  ///< The radius of the Earth [L ~> m].

  amrex::MultiFab geoLatT;   ///< The geographic latitude at h points [degrees_N].
  amrex::MultiFab geoLonT;   ///< The geographic longitude at h points [degrees_E].
  amrex::MultiFab dxT;       ///< Delta x at h points [L ~> m].
  amrex::MultiFab dyT;       ///< Delta y at h points [L ~> m].
  amrex::MultiFab areaT;     ///< The area of an h-cell [L2 ~> m2].

  amrex::MultiFab geoLatCu;  ///< The geographic latitude at u points [degrees_N].
  amrex::MultiFab geoLonCu;  ///< The geographic longitude at u points [degrees_E].
  amrex::MultiFab dxCu;      ///< Delta x at u points [L ~> m].
  amrex::MultiFab dyCu;      ///< Delta y at u points [L ~> m].

  amrex::MultiFab geoLatCv;  ///< The geographic latitude at v points [degrees_N].
  amrex::MultiFab geoLonCv;  ///< The geographic longitude at v points [degrees_E].
  amrex::MultiFab dxCv;      ///< Delta x at v points [L ~> m].
  amrex::MultiFab dyCv;      ///< Delta y at v points [L ~> m].

  amrex::MultiFab geoLatBu;  ///< The geographic latitude at q points [degrees_N].
  amrex::MultiFab geoLonBu;  ///< The geographic longitude at q points [degrees_E].
  amrex::MultiFab dxBu;      ///< Delta x at q points [L ~> m].
  amrex::MultiFab dyBu;      ///< Delta y at q points [L ~> m].

  amrex::MultiFab CoriolisBu;  ///< The Coriolis parameter at q points [T-1 ~> s-1].

  amrex::Real max_depth = 0.0;  ///< The maximum depth of the ocean [Z ~> m].
  amrex::MultiFab bathyT;    ///< The ocean bottom depth at h points, positive down [Z ~> m].

  amrex::MultiFab mask2dT;   ///< 1 for ocean, 0 for land at h points [nondim].
  amrex::MultiFab mask2dCu;  ///< 1 for ocean, 0 for land at u points [nondim].
  amrex::MultiFab mask2dCv;  ///< 1 for ocean, 0 for land at v points [nondim].
  amrex::MultiFab mask2dBu;  ///< 1 for ocean, 0 for land at q points [nondim].

  // The derived metrics of MOM6's set_derived_dyn_horgrid plus the masked
  // face lengths and areas set at the end of initialize_masks. Every
  // reciprocal is an Adcroft reciprocal: 1/x, or 0 where x is 0.
  amrex::MultiFab IdxT;      ///< 1/dxT at h points [L-1 ~> m-1].
  amrex::MultiFab IdyT;      ///< 1/dyT at h points [L-1 ~> m-1].
  amrex::MultiFab IareaT;    ///< 1/areaT at h points [L-2 ~> m-2].
  amrex::MultiFab IdxCu;     ///< 1/dxCu at u points [L-1 ~> m-1].
  amrex::MultiFab IdyCu;     ///< 1/dyCu at u points [L-1 ~> m-1].
  amrex::MultiFab IdxCv;     ///< 1/dxCv at v points [L-1 ~> m-1].
  amrex::MultiFab IdyCv;     ///< 1/dyCv at v points [L-1 ~> m-1].
  amrex::MultiFab IdxBu;     ///< 1/dxBu at q points [L-1 ~> m-1].
  amrex::MultiFab IdyBu;     ///< 1/dyBu at q points [L-1 ~> m-1].
  amrex::MultiFab areaBu;    ///< The area of a q-cell [L2 ~> m2].
  amrex::MultiFab IareaBu;   ///< 1/areaBu at q points [L-2 ~> m-2].
  amrex::MultiFab dy_Cu;     ///< The unblocked length of a u face [L ~> m].
  amrex::MultiFab dx_Cv;     ///< The unblocked length of a v face [L ~> m].
  amrex::MultiFab areaCu;    ///< The area of a u-cell [L2 ~> m2].
  amrex::MultiFab areaCv;    ///< The area of a v-cell [L2 ~> m2].
  amrex::MultiFab IareaCu;   ///< The masked 1/areaCu at u points [L-2 ~> m-2].
  amrex::MultiFab IareaCv;   ///< The masked 1/areaCv at v points [L-2 ~> m-2].
};

} // namespace MOM
