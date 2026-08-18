#include <cmath>

#include "MOM_grid_initialize.h"

namespace MOM {

GridFields spherical_grid_fields(const Domain &domain, const GridSpec &spec) {

  // The metric fields are 2-D (single-level) fields on the domain's
  // horizontal decomposition, created through the domain's field factory.
  const int n_levels = 1;
  const int ncomp = 1;

  GridFields fields;

  fields.south_lat = spec.south_lat;
  fields.len_lat = spec.len_lat;
  fields.west_lon = spec.west_lon;
  fields.len_lon = spec.len_lon;
  fields.rad_earth = spec.rad_earth;

  fields.geoLatT = domain.make_field(Stagger::Cell, n_levels, ncomp);
  fields.geoLonT = domain.make_field(Stagger::Cell, n_levels, ncomp);
  fields.dxT = domain.make_field(Stagger::Cell, n_levels, ncomp);
  fields.dyT = domain.make_field(Stagger::Cell, n_levels, ncomp);
  fields.areaT = domain.make_field(Stagger::Cell, n_levels, ncomp);

  fields.geoLatCu = domain.make_field(Stagger::XFace, n_levels, ncomp);
  fields.geoLonCu = domain.make_field(Stagger::XFace, n_levels, ncomp);
  fields.dxCu = domain.make_field(Stagger::XFace, n_levels, ncomp);
  fields.dyCu = domain.make_field(Stagger::XFace, n_levels, ncomp);

  fields.geoLatCv = domain.make_field(Stagger::YFace, n_levels, ncomp);
  fields.geoLonCv = domain.make_field(Stagger::YFace, n_levels, ncomp);
  fields.dxCv = domain.make_field(Stagger::YFace, n_levels, ncomp);
  fields.dyCv = domain.make_field(Stagger::YFace, n_levels, ncomp);

  fields.geoLatBu = domain.make_field(Stagger::Node, n_levels, ncomp);
  fields.geoLonBu = domain.make_field(Stagger::Node, n_levels, ncomp);
  fields.dxBu = domain.make_field(Stagger::Node, n_levels, ncomp);
  fields.dyBu = domain.make_field(Stagger::Node, n_levels, ncomp);

  const amrex::Real PI = 4.0 * std::atan(1.0);
  const amrex::Real PI_180 = PI / 180.0;

  // The change in longitude/latitude between successive grid points [degrees].
  const amrex::Real dLon = spec.len_lon / domain.ni_global();
  const amrex::Real dLat = spec.len_lat / domain.nj_global();
  // dLon rescaled from degrees to radians [radians]. MOM6 computes the zonal
  // spacings from this expression (rather than from dLon*PI_180 directly) to
  // reproduce the set_grid_metrics_mercator solution on a simple spherical
  // grid; kept identical here for future parity.
  const amrex::Real dL_di = (spec.len_lon * PI) / (180.0 * domain.ni_global());
  // The meridional spacing is uniform on a spherical grid [L ~> m].
  const amrex::Real dy = spec.rad_earth * dLat * PI_180;

  const amrex::Real south_lat = spec.south_lat;
  const amrex::Real west_lon = spec.west_lon;
  const amrex::Real rad_earth = spec.rad_earth;

  // All values are computed analytically from the global index (AMReX indices
  // are global), over the grown boxes: halo values (including those beyond
  // the global domain edges) extrapolate the same formulas, as in MOM6.
  // Cell centers sit at half-integer indices, corners at integer indices;
  // latitudes are clamped to [-90, 90] as in MOM6.

  // h points (cell centers)
  for (amrex::MFIter mfi(fields.geoLatT); mfi.isValid(); ++mfi) {
    const amrex::Box bx = mfi.growntilebox();
    const amrex::Array4<amrex::Real> geoLatT = fields.geoLatT.array(mfi);
    const amrex::Array4<amrex::Real> geoLonT = fields.geoLonT.array(mfi);
    const amrex::Array4<amrex::Real> dxT = fields.dxT.array(mfi);
    const amrex::Array4<amrex::Real> dyT = fields.dyT.array(mfi);
    const amrex::Array4<amrex::Real> areaT = fields.areaT.array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      geoLonT(i, j, k) = west_lon + dLon * (i + 0.5);
      geoLatT(i, j, k) = amrex::min(amrex::max(south_lat + dLat * (j + 0.5),
                                               amrex::Real(-90.0)), amrex::Real(90.0));
      dxT(i, j, k) = rad_earth * std::cos(geoLatT(i, j, k) * PI_180) * dL_di;
      dyT(i, j, k) = dy;
      areaT(i, j, k) = dxT(i, j, k) * dyT(i, j, k);
    });
  }

  // u points (east faces: corner longitudes, cell-center latitudes)
  for (amrex::MFIter mfi(fields.geoLatCu); mfi.isValid(); ++mfi) {
    const amrex::Box bx = mfi.growntilebox();
    const amrex::Array4<amrex::Real> geoLatCu = fields.geoLatCu.array(mfi);
    const amrex::Array4<amrex::Real> geoLonCu = fields.geoLonCu.array(mfi);
    const amrex::Array4<amrex::Real> dxCu = fields.dxCu.array(mfi);
    const amrex::Array4<amrex::Real> dyCu = fields.dyCu.array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      geoLonCu(i, j, k) = west_lon + dLon * i;
      geoLatCu(i, j, k) = amrex::min(amrex::max(south_lat + dLat * (j + 0.5),
                                                amrex::Real(-90.0)), amrex::Real(90.0));
      dxCu(i, j, k) = rad_earth * std::cos(geoLatCu(i, j, k) * PI_180) * dL_di;
      dyCu(i, j, k) = dy;
    });
  }

  // v points (north faces: cell-center longitudes, corner latitudes)
  for (amrex::MFIter mfi(fields.geoLatCv); mfi.isValid(); ++mfi) {
    const amrex::Box bx = mfi.growntilebox();
    const amrex::Array4<amrex::Real> geoLatCv = fields.geoLatCv.array(mfi);
    const amrex::Array4<amrex::Real> geoLonCv = fields.geoLonCv.array(mfi);
    const amrex::Array4<amrex::Real> dxCv = fields.dxCv.array(mfi);
    const amrex::Array4<amrex::Real> dyCv = fields.dyCv.array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      geoLonCv(i, j, k) = west_lon + dLon * (i + 0.5);
      geoLatCv(i, j, k) = amrex::min(amrex::max(south_lat + dLat * j,
                                                amrex::Real(-90.0)), amrex::Real(90.0));
      dxCv(i, j, k) = rad_earth * std::cos(geoLatCv(i, j, k) * PI_180) * dL_di;
      dyCv(i, j, k) = dy;
    });
  }

  // q points (cell corners)
  for (amrex::MFIter mfi(fields.geoLatBu); mfi.isValid(); ++mfi) {
    const amrex::Box bx = mfi.growntilebox();
    const amrex::Array4<amrex::Real> geoLatBu = fields.geoLatBu.array(mfi);
    const amrex::Array4<amrex::Real> geoLonBu = fields.geoLonBu.array(mfi);
    const amrex::Array4<amrex::Real> dxBu = fields.dxBu.array(mfi);
    const amrex::Array4<amrex::Real> dyBu = fields.dyBu.array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      geoLonBu(i, j, k) = west_lon + dLon * i;
      geoLatBu(i, j, k) = amrex::min(amrex::max(south_lat + dLat * j,
                                                amrex::Real(-90.0)), amrex::Real(90.0));
      dxBu(i, j, k) = rad_earth * std::cos(geoLatBu(i, j, k) * PI_180) * dL_di;
      dyBu(i, j, k) = dy;
    });
  }

  return fields;
}

namespace {

// 1/x, or 0 where x is 0. MOM6's Adcroft_reciprocal.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real Adcroft_reciprocal(const amrex::Real val) {
  return (val != 0.0) ? (1.0 / val) : 0.0;
}

} // namespace

void initialize_masks(const Domain &domain, const TopoSpec &topo_spec,
                      const amrex::MultiFab &bathyT, GridFields &fields) {

  const int n_levels = 1;
  const int ncomp = 1;

  fields.mask2dT = domain.make_field(Stagger::Cell, n_levels, ncomp);
  fields.mask2dCu = domain.make_field(Stagger::XFace, n_levels, ncomp);
  fields.mask2dCv = domain.make_field(Stagger::YFace, n_levels, ncomp);
  fields.mask2dBu = domain.make_field(Stagger::Node, n_levels, ncomp);

  // MOM6 zeroes all four masks before filling them; the points its loops skip
  // (the outermost face and corner rows of the data domain) stay zero here too.
  fields.mask2dT.setVal(0.0);
  fields.mask2dCu.setVal(0.0);
  fields.mask2dCv.setVal(0.0);
  fields.mask2dBu.setVal(0.0);

  // MOM6 masks with MASKING_DEPTH when it is set and with MINIMUM_DEPTH
  // otherwise; only the latter path is implemented.
  const amrex::Real Dmask = topo_spec.min_depth;

  for (amrex::MFIter mfi(fields.mask2dT); mfi.isValid(); ++mfi) {
    // The cell-centered box including halos: MOM6's data domain, isd:ied.
    const amrex::Box cells = amrex::grow(mfi.validbox(), domain.nghost());

    const amrex::Array4<amrex::Real> maskT = fields.mask2dT.array(mfi);
    const amrex::Array4<amrex::Real> maskCu = fields.mask2dCu.array(mfi);
    const amrex::Array4<amrex::Real> maskCv = fields.mask2dCv.array(mfi);
    const amrex::Array4<amrex::Real> maskBu = fields.mask2dBu.array(mfi);
    const amrex::Array4<const amrex::Real> D = bathyT.const_array(mfi);

    amrex::ParallelFor(cells, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      maskT(i, j, k) = (D(i, j, k) <= Dmask) ? 0.0 : 1.0;
    });

    // AMReX anchors face i on the low side of cell i, so the u point between
    // cells i-1 and i carries index i; MOM6 labels the same point I = i-1.
    // Its loop runs I = isd .. ied-1, which is i = isd+1 .. ied here.
    const amrex::Box u_faces =
        amrex::surroundingNodes(cells, 0).growLo(0, -1).growHi(0, -1);
    amrex::ParallelFor(u_faces, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      maskCu(i, j, k) = maskT(i - 1, j, k) * maskT(i, j, k);
    });

    const amrex::Box v_faces =
        amrex::surroundingNodes(cells, 1).growLo(1, -1).growHi(1, -1);
    amrex::ParallelFor(v_faces, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      maskCv(i, j, k) = maskT(i, j - 1, k) * maskT(i, j, k);
    });

    // The corner mask follows from the face masks, as in MOM6.
    amrex::Box corners = amrex::surroundingNodes(cells, 0);
    corners = amrex::surroundingNodes(corners, 1);
    corners.growLo(0, -1).growHi(0, -1).growLo(1, -1).growHi(1, -1);
    amrex::ParallelFor(corners, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      maskBu(i, j, k) = (maskCu(i, j - 1, k) * maskCu(i, j, k)) *
                        (maskCv(i - 1, j, k) * maskCv(i, j, k));
    });
  }

  // defer: the OBC-aware mask revisions of initialize_masks (OBC_dir_u/v and
  //        open_corner_OBCs); the double_gyre configuration has no open
  //        boundaries.
}

void set_derived_metrics(const Domain &domain, GridFields &fields) {

  const int n_levels = 1;
  const int ncomp = 1;

  fields.IdxT = domain.make_field(Stagger::Cell, n_levels, ncomp);
  fields.IdyT = domain.make_field(Stagger::Cell, n_levels, ncomp);
  fields.IareaT = domain.make_field(Stagger::Cell, n_levels, ncomp);
  fields.IdxCu = domain.make_field(Stagger::XFace, n_levels, ncomp);
  fields.IdyCu = domain.make_field(Stagger::XFace, n_levels, ncomp);
  fields.IdxCv = domain.make_field(Stagger::YFace, n_levels, ncomp);
  fields.IdyCv = domain.make_field(Stagger::YFace, n_levels, ncomp);
  fields.IdxBu = domain.make_field(Stagger::Node, n_levels, ncomp);
  fields.IdyBu = domain.make_field(Stagger::Node, n_levels, ncomp);
  fields.areaBu = domain.make_field(Stagger::Node, n_levels, ncomp);
  fields.IareaBu = domain.make_field(Stagger::Node, n_levels, ncomp);
  fields.dy_Cu = domain.make_field(Stagger::XFace, n_levels, ncomp);
  fields.dx_Cv = domain.make_field(Stagger::YFace, n_levels, ncomp);
  fields.areaCu = domain.make_field(Stagger::XFace, n_levels, ncomp);
  fields.areaCv = domain.make_field(Stagger::YFace, n_levels, ncomp);
  fields.IareaCu = domain.make_field(Stagger::XFace, n_levels, ncomp);
  fields.IareaCv = domain.make_field(Stagger::YFace, n_levels, ncomp);
  fields.IdxCu_OBCmask = domain.make_field(Stagger::XFace, n_levels, ncomp);
  fields.IdyCv_OBCmask = domain.make_field(Stagger::YFace, n_levels, ncomp);

  for (amrex::MFIter mfi(fields.IdxT); mfi.isValid(); ++mfi) {
    const amrex::Box cells = amrex::grow(mfi.validbox(), domain.nghost());
    const amrex::Box u_faces = amrex::surroundingNodes(cells, 0);
    const amrex::Box v_faces = amrex::surroundingNodes(cells, 1);
    const amrex::Box corners = amrex::surroundingNodes(amrex::surroundingNodes(cells, 0), 1);

    const amrex::Array4<const amrex::Real> dxT = fields.dxT.const_array(mfi);
    const amrex::Array4<const amrex::Real> dyT = fields.dyT.const_array(mfi);
    const amrex::Array4<const amrex::Real> areaT = fields.areaT.const_array(mfi);
    const amrex::Array4<amrex::Real> IdxT = fields.IdxT.array(mfi);
    const amrex::Array4<amrex::Real> IdyT = fields.IdyT.array(mfi);
    const amrex::Array4<amrex::Real> IareaT = fields.IareaT.array(mfi);
    amrex::ParallelFor(cells, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      IdxT(i, j, k) = Adcroft_reciprocal(dxT(i, j, k));
      IdyT(i, j, k) = Adcroft_reciprocal(dyT(i, j, k));
      IareaT(i, j, k) = Adcroft_reciprocal(areaT(i, j, k));
    });

    const amrex::Array4<const amrex::Real> dxCu = fields.dxCu.const_array(mfi);
    const amrex::Array4<const amrex::Real> dyCu = fields.dyCu.const_array(mfi);
    const amrex::Array4<const amrex::Real> maskCu = fields.mask2dCu.const_array(mfi);
    const amrex::Array4<amrex::Real> IdxCu = fields.IdxCu.array(mfi);
    const amrex::Array4<amrex::Real> IdyCu = fields.IdyCu.array(mfi);
    const amrex::Array4<amrex::Real> dy_Cu = fields.dy_Cu.array(mfi);
    const amrex::Array4<amrex::Real> areaCu = fields.areaCu.array(mfi);
    const amrex::Array4<amrex::Real> IareaCu = fields.IareaCu.array(mfi);
    const amrex::Array4<amrex::Real> IdxCu_OBCmask = fields.IdxCu_OBCmask.array(mfi);
    amrex::ParallelFor(u_faces, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      IdxCu(i, j, k) = Adcroft_reciprocal(dxCu(i, j, k));
      // With no open boundaries MOM6's OBCmaskCu is the land/sea mask.
      IdxCu_OBCmask(i, j, k) = maskCu(i, j, k) * IdxCu(i, j, k);
      IdyCu(i, j, k) = Adcroft_reciprocal(dyCu(i, j, k));
      dy_Cu(i, j, k) = maskCu(i, j, k) * dyCu(i, j, k);
      areaCu(i, j, k) = dxCu(i, j, k) * dy_Cu(i, j, k);
      IareaCu(i, j, k) = maskCu(i, j, k) * Adcroft_reciprocal(areaCu(i, j, k));
    });

    const amrex::Array4<const amrex::Real> dxCv = fields.dxCv.const_array(mfi);
    const amrex::Array4<const amrex::Real> dyCv = fields.dyCv.const_array(mfi);
    const amrex::Array4<const amrex::Real> maskCv = fields.mask2dCv.const_array(mfi);
    const amrex::Array4<amrex::Real> IdxCv = fields.IdxCv.array(mfi);
    const amrex::Array4<amrex::Real> IdyCv = fields.IdyCv.array(mfi);
    const amrex::Array4<amrex::Real> dx_Cv = fields.dx_Cv.array(mfi);
    const amrex::Array4<amrex::Real> areaCv = fields.areaCv.array(mfi);
    const amrex::Array4<amrex::Real> IareaCv = fields.IareaCv.array(mfi);
    const amrex::Array4<amrex::Real> IdyCv_OBCmask = fields.IdyCv_OBCmask.array(mfi);
    amrex::ParallelFor(v_faces, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      IdxCv(i, j, k) = Adcroft_reciprocal(dxCv(i, j, k));
      IdyCv_OBCmask(i, j, k) = maskCv(i, j, k) * Adcroft_reciprocal(dyCv(i, j, k));
      IdyCv(i, j, k) = Adcroft_reciprocal(dyCv(i, j, k));
      dx_Cv(i, j, k) = maskCv(i, j, k) * dxCv(i, j, k);
      areaCv(i, j, k) = dyCv(i, j, k) * dx_Cv(i, j, k);
      IareaCv(i, j, k) = maskCv(i, j, k) * Adcroft_reciprocal(areaCv(i, j, k));
    });

    const amrex::Array4<const amrex::Real> dxBu = fields.dxBu.const_array(mfi);
    const amrex::Array4<const amrex::Real> dyBu = fields.dyBu.const_array(mfi);
    const amrex::Array4<amrex::Real> IdxBu = fields.IdxBu.array(mfi);
    const amrex::Array4<amrex::Real> IdyBu = fields.IdyBu.array(mfi);
    const amrex::Array4<amrex::Real> areaBu = fields.areaBu.array(mfi);
    const amrex::Array4<amrex::Real> IareaBu = fields.IareaBu.array(mfi);
    amrex::ParallelFor(corners, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      IdxBu(i, j, k) = Adcroft_reciprocal(dxBu(i, j, k));
      IdyBu(i, j, k) = Adcroft_reciprocal(dyBu(i, j, k));
      // On a spherical grid set_grid_metrics_spherical already sets areaBu
      // to dxBu*dyBu; MOM6's fallback here computes the same value.
      areaBu(i, j, k) = dxBu(i, j, k) * dyBu(i, j, k);
      IareaBu(i, j, k) = Adcroft_reciprocal(areaBu(i, j, k));
    });
  }
}

} // namespace MOM
