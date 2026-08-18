#include <cmath>

#include "MOM_hor_visc.h"

#include "MOM_logger.h"
#include "MOM_loop_boxes.h"

namespace MOM {

namespace {

// A thickness so small it is lost in roundoff [H ~> m]: MOM6's
// GV%H_subroundoff at the default ANGSTROM.
constexpr amrex::Real H_NEGLECT = 1.0e-30;

// Abort if a parameter selects a branch that is not implemented.
void reject_if_set(RuntimeParams &params, const std::string &key, const bool default_value,
                   const std::string &desc) {
  bool value = default_value;
  params.get(key, value, {.default_value = default_value, .desc = desc});
  if (value != default_value) {
    logger::fatal("HorizontalViscosity: ", key, " is not implemented yet.");
  }
}

} // namespace

HorizontalViscosity::HorizontalViscosity(RuntimeParams &params, const amrex::Real dt,
                                         const Domain &domain, const Grid &grid) {

  params.doc_module("MOM_hor_visc", "");

  bool laplacian = false;
  params.get("LAPLACIAN", laplacian,
             {.default_value = false,
              .desc = "If true, use a Laplacian horizontal viscosity."});
  if (!laplacian) {
    // defer: running with no Laplacian viscosity at all, which is only
    //        sensible together with the biharmonic branch.
    logger::fatal("HorizontalViscosity: LAPLACIAN = False is not implemented yet.");
  }

  amrex::Real Kh = 0.0;
  params.get("KH", Kh,
             {.default_value = 0.0,
              .desc = "The background Laplacian horizontal viscosity.",
              .units = "m2 s-1"});

  amrex::Real Kh_bg_min = 0.0;
  params.get("KH_BG_MIN", Kh_bg_min,
             {.default_value = 0.0,
              .desc = "The minimum value allowed for Laplacian horizontal viscosity.",
              .units = "m2 s-1"});

  amrex::Real Kh_vel_scale = 0.0;
  params.get("KH_VEL_SCALE", Kh_vel_scale,
             {.default_value = 0.0,
              .desc = "The velocity scale which is multiplied by the grid spacing to calculate "
                      "the Laplacian viscosity. The final viscosity is the largest of this "
                      "scaled viscosity, the Smagorinsky and Leith viscosities, and KH.",
              .units = "m s-1"});

  params.get("SMAGORINSKY_KH", Smagorinsky_Kh_,
             {.default_value = false,
              .desc = "If true, use a Smagorinsky nonlinear eddy viscosity."});

  amrex::Real Smag_Lap_const = 0.0;
  if (Smagorinsky_Kh_) {
    params.get("SMAG_LAP_CONST", Smag_Lap_const,
               {.default_value = 0.0,
                .desc = "The nondimensional Laplacian Smagorinsky constant, often 0.15.",
                .units = "nondim"});
  }

  reject_if_set(params, "LEITH_KH", false,
                "If true, use a Leith nonlinear eddy viscosity.");

  params.get("BOUND_KH", bound_Kh_,
             {.default_value = true,
              .desc = "If true, the Laplacian coefficient is locally limited to be stable."});

  reject_if_set(params, "ANISOTROPIC_VISCOSITY", false,
                "If true, allow anisotropic viscosity in the Laplacian horizontal viscosity.");
  reject_if_set(params, "ADD_LES_VISCOSITY", false,
                "If true, adds the viscosity from Smagorinsky and Leith to the background "
                "viscosity instead of taking the maximum.");

  bool biharmonic = true;
  params.get("BIHARMONIC", biharmonic,
             {.default_value = true,
              .desc = "If true, use a biharmonic horizontal viscosity. BIHARMONIC may be used "
                      "with LAPLACIAN."});
  if (biharmonic) {
    // defer: the biharmonic branch, its Smagorinsky and Leith variants, and
    //        the Del2u/Del2v machinery they need.
    logger::fatal("HorizontalViscosity: BIHARMONIC is not implemented yet.");
  }

  bool use_land_mask = true;
  params.get("USE_LAND_MASK_FOR_HVISC", use_land_mask,
             {.default_value = true,
              .desc = "If true, use the land mask for the computation of thicknesses at velocity "
                      "locations. This eliminates the dependence on arbitrary values over land or "
                      "outside of the domain."});
  if (!use_land_mask) {
    // defer: the unmasked thickness interpolation.
    logger::fatal("HorizontalViscosity: USE_LAND_MASK_FOR_HVISC = False is not implemented yet.");
  }

  amrex::Real bound_coef = 0.8;
  params.get("HORVISC_BOUND_COEF", bound_coef,
             {.default_value = 0.8,
              .desc = "The nondimensional coefficient of the ratio of the viscosity bounds to "
                      "the theoretical maximum for stability without considering other terms.",
              .units = "nondim"});

  // --- The metric combinations of hor_visc_init. ---
  const int n = 1;
  DY_dxT_ = domain.make_field(Stagger::Cell, n, n);
  DX_dyT_ = domain.make_field(Stagger::Cell, n, n);
  dx2h_ = domain.make_field(Stagger::Cell, n, n);
  dy2h_ = domain.make_field(Stagger::Cell, n, n);
  DY_dxBu_ = domain.make_field(Stagger::Node, n, n);
  DX_dyBu_ = domain.make_field(Stagger::Node, n, n);
  dx2q_ = domain.make_field(Stagger::Node, n, n);
  dy2q_ = domain.make_field(Stagger::Node, n, n);
  reduction_xx_ = domain.make_field(Stagger::Cell, n, n);
  reduction_xy_ = domain.make_field(Stagger::Node, n, n);
  Kh_bg_xx_ = domain.make_field(Stagger::Cell, n, n);
  Kh_bg_xy_ = domain.make_field(Stagger::Node, n, n);
  Laplac2_const_xx_ = domain.make_field(Stagger::Cell, n, n);
  Laplac2_const_xy_ = domain.make_field(Stagger::Node, n, n);
  Kh_Max_xx_ = domain.make_field(Stagger::Cell, n, n);
  Kh_Max_xy_ = domain.make_field(Stagger::Node, n, n);
  Laplac2_const_xx_.setVal(0.0);
  Laplac2_const_xy_.setVal(0.0);
  Kh_Max_xx_.setVal(0.0);
  Kh_Max_xy_.setVal(0.0);

  const amrex::Real Idt = 1.0 / dt;
  const bool smag = Smagorinsky_Kh_;
  const bool bound = bound_Kh_;

  for (amrex::MFIter mfi(DY_dxT_); mfi.isValid(); ++mfi) {
    const amrex::Box valid = mfi.validbox();
    const amrex::Box cells = amrex::grow(valid, domain.nghost());
    const amrex::Box corners = amrex::surroundingNodes(amrex::surroundingNodes(cells, 0), 1);

    const amrex::Array4<const amrex::Real> dxT = grid.dxT().const_array(mfi);
    const amrex::Array4<const amrex::Real> dyT = grid.dyT().const_array(mfi);
    const amrex::Array4<const amrex::Real> IdxT = grid.IdxT().const_array(mfi);
    const amrex::Array4<const amrex::Real> IdyT = grid.IdyT().const_array(mfi);
    const amrex::Array4<amrex::Real> DY_dxT = DY_dxT_.array(mfi);
    const amrex::Array4<amrex::Real> DX_dyT = DX_dyT_.array(mfi);
    const amrex::Array4<amrex::Real> dx2h = dx2h_.array(mfi);
    const amrex::Array4<amrex::Real> dy2h = dy2h_.array(mfi);
    amrex::ParallelFor(cells, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      dx2h(i, j, k) = dxT(i, j, k) * dxT(i, j, k);
      dy2h(i, j, k) = dyT(i, j, k) * dyT(i, j, k);
      DX_dyT(i, j, k) = dxT(i, j, k) * IdyT(i, j, k);
      DY_dxT(i, j, k) = dyT(i, j, k) * IdxT(i, j, k);
    });

    const amrex::Array4<const amrex::Real> dxBu = grid.dxBu().const_array(mfi);
    const amrex::Array4<const amrex::Real> dyBu = grid.dyBu().const_array(mfi);
    const amrex::Array4<const amrex::Real> IdxBu = grid.IdxBu().const_array(mfi);
    const amrex::Array4<const amrex::Real> IdyBu = grid.IdyBu().const_array(mfi);
    const amrex::Array4<amrex::Real> DY_dxBu = DY_dxBu_.array(mfi);
    const amrex::Array4<amrex::Real> DX_dyBu = DX_dyBu_.array(mfi);
    const amrex::Array4<amrex::Real> dx2q = dx2q_.array(mfi);
    const amrex::Array4<amrex::Real> dy2q = dy2q_.array(mfi);
    amrex::ParallelFor(corners, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      dx2q(i, j, k) = dxBu(i, j, k) * dxBu(i, j, k);
      dy2q(i, j, k) = dyBu(i, j, k) * dyBu(i, j, k);
      DX_dyBu(i, j, k) = dxBu(i, j, k) * IdyBu(i, j, k);
      DY_dxBu(i, j, k) = dyBu(i, j, k) * IdxBu(i, j, k);
    });

    // The blocked-face reduction factors. They stay at 1 unless a face has
    // been narrowed by reset_face_lengths, which is not implemented.
    const amrex::Array4<const amrex::Real> dy_Cu = grid.dy_Cu().const_array(mfi);
    const amrex::Array4<const amrex::Real> dyCu = grid.dyCu().const_array(mfi);
    const amrex::Array4<const amrex::Real> dx_Cv = grid.dx_Cv().const_array(mfi);
    const amrex::Array4<const amrex::Real> dxCv = grid.dxCv().const_array(mfi);
    const amrex::Array4<amrex::Real> red_xx = reduction_xx_.array(mfi);
    amrex::ParallelFor(amrex::grow(valid, amrex::IntVect(1, 1, 0)), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      amrex::Real r = 1.0;
      if ((dy_Cu(i + 1, j, k) > 0.0) && (dy_Cu(i + 1, j, k) < dyCu(i + 1, j, k)) &&
          (dy_Cu(i + 1, j, k) < dyCu(i + 1, j, k) * r)) r = dy_Cu(i + 1, j, k) / dyCu(i + 1, j, k);
      if ((dy_Cu(i, j, k) > 0.0) && (dy_Cu(i, j, k) < dyCu(i, j, k)) &&
          (dy_Cu(i, j, k) < dyCu(i, j, k) * r)) r = dy_Cu(i, j, k) / dyCu(i, j, k);
      if ((dx_Cv(i, j + 1, k) > 0.0) && (dx_Cv(i, j + 1, k) < dxCv(i, j + 1, k)) &&
          (dx_Cv(i, j + 1, k) < dxCv(i, j + 1, k) * r)) r = dx_Cv(i, j + 1, k) / dxCv(i, j + 1, k);
      if ((dx_Cv(i, j, k) > 0.0) && (dx_Cv(i, j, k) < dxCv(i, j, k)) &&
          (dx_Cv(i, j, k) < dxCv(i, j, k) * r)) r = dx_Cv(i, j, k) / dxCv(i, j, k);
      red_xx(i, j, k) = r;
    });
    const amrex::Array4<amrex::Real> red_xy = reduction_xy_.array(mfi);
    amrex::ParallelFor(loops::q_points(valid), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      amrex::Real r = 1.0;
      if ((dy_Cu(i, j - 1, k) > 0.0) && (dy_Cu(i, j - 1, k) < dyCu(i, j - 1, k)) &&
          (dy_Cu(i, j - 1, k) < dyCu(i, j - 1, k) * r)) r = dy_Cu(i, j - 1, k) / dyCu(i, j - 1, k);
      if ((dy_Cu(i, j, k) > 0.0) && (dy_Cu(i, j, k) < dyCu(i, j, k)) &&
          (dy_Cu(i, j, k) < dyCu(i, j, k) * r)) r = dy_Cu(i, j, k) / dyCu(i, j, k);
      if ((dx_Cv(i - 1, j, k) > 0.0) && (dx_Cv(i - 1, j, k) < dxCv(i - 1, j, k)) &&
          (dx_Cv(i - 1, j, k) < dxCv(i - 1, j, k) * r)) r = dx_Cv(i - 1, j, k) / dxCv(i - 1, j, k);
      if ((dx_Cv(i, j, k) > 0.0) && (dx_Cv(i, j, k) < dxCv(i, j, k)) &&
          (dx_Cv(i, j, k) < dxCv(i, j, k) * r)) r = dx_Cv(i, j, k) / dxCv(i, j, k);
      red_xy(i, j, k) = r;
    });

    // The background and Smagorinsky viscosity constants, from the harmonic
    // mean of the squared grid spacings.
    const amrex::Array4<amrex::Real> Kh_bg_xx = Kh_bg_xx_.array(mfi);
    const amrex::Array4<amrex::Real> Lap2_xx = Laplac2_const_xx_.array(mfi);
    amrex::ParallelFor(amrex::grow(valid, amrex::IntVect(1, 1, 0)), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const amrex::Real grid_sp_h2 = (2.0 * dx2h(i, j, k) * dy2h(i, j, k)) /
                                     (dx2h(i, j, k) + dy2h(i, j, k));
      if (smag) Lap2_xx(i, j, k) = Smag_Lap_const * grid_sp_h2;
      Kh_bg_xx(i, j, k) = amrex::max(Kh, Kh_vel_scale * std::sqrt(grid_sp_h2));
    });
    const amrex::Array4<amrex::Real> Kh_bg_xy = Kh_bg_xy_.array(mfi);
    const amrex::Array4<amrex::Real> Lap2_xy = Laplac2_const_xy_.array(mfi);
    amrex::ParallelFor(loops::q_points(valid), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const amrex::Real grid_sp_q2 = (2.0 * dx2q(i, j, k) * dy2q(i, j, k)) /
                                     (dx2q(i, j, k) + dy2q(i, j, k));
      if (smag) Lap2_xy(i, j, k) = Smag_Lap_const * grid_sp_q2;
      Kh_bg_xy(i, j, k) = amrex::max(Kh, Kh_vel_scale * std::sqrt(grid_sp_q2));
    });

    if (bound) {
      const amrex::Array4<const amrex::Real> IdyCu = grid.IdyCu().const_array(mfi);
      const amrex::Array4<const amrex::Real> IdxCu = grid.IdxCu().const_array(mfi);
      const amrex::Array4<const amrex::Real> IdxCv = grid.IdxCv().const_array(mfi);
      const amrex::Array4<const amrex::Real> IdyCv = grid.IdyCv().const_array(mfi);
      const amrex::Array4<const amrex::Real> IareaCu = grid.IareaCu().const_array(mfi);
      const amrex::Array4<const amrex::Real> IareaCv = grid.IareaCv().const_array(mfi);
      const amrex::Array4<amrex::Real> KhMax_xx = Kh_Max_xx_.array(mfi);
      amrex::ParallelFor(amrex::grow(valid, amrex::IntVect(1, 1, 0)), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real denom = amrex::max(
            (dy2h(i, j, k) * DY_dxT(i, j, k) * (IdyCu(i + 1, j, k) + IdyCu(i, j, k)) *
             amrex::max(IdyCu(i + 1, j, k) * IareaCu(i + 1, j, k),
                        IdyCu(i, j, k) * IareaCu(i, j, k))),
            (dx2h(i, j, k) * DX_dyT(i, j, k) * (IdxCv(i, j + 1, k) + IdxCv(i, j, k)) *
             amrex::max(IdxCv(i, j + 1, k) * IareaCv(i, j + 1, k),
                        IdxCv(i, j, k) * IareaCv(i, j, k))));
        KhMax_xx(i, j, k) = (denom > 0.0) ? (bound_coef * 0.25 * Idt / denom) : 0.0;
      });
      const amrex::Array4<amrex::Real> KhMax_xy = Kh_Max_xy_.array(mfi);
      amrex::ParallelFor(loops::q_points(valid), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real denom = amrex::max(
            (dx2q(i, j, k) * DX_dyBu(i, j, k) * (IdxCu(i, j, k) + IdxCu(i, j - 1, k)) *
             amrex::max(IdxCu(i, j - 1, k) * IareaCu(i, j - 1, k),
                        IdxCu(i, j, k) * IareaCu(i, j, k))),
            (dy2q(i, j, k) * DY_dxBu(i, j, k) * (IdyCv(i, j, k) + IdyCv(i - 1, j, k)) *
             amrex::max(IdyCv(i - 1, j, k) * IareaCv(i - 1, j, k),
                        IdyCv(i, j, k) * IareaCv(i, j, k))));
        KhMax_xy(i, j, k) = (denom > 0.0) ? (bound_coef * 0.25 * Idt / denom) : 0.0;
      });
    }
  }

  Kh_bg_min_ = Kh_bg_min;
}

void HorizontalViscosity::calculate(amrex::MultiFab &diffu, amrex::MultiFab &diffv,
                                    const amrex::MultiFab &u, const amrex::MultiFab &v,
                                    const amrex::MultiFab &h, const Domain &domain,
                                    const Grid &grid, const VerticalGrid &vgrid) const {

  const int nk = vgrid.nk();
  const bool smag = Smagorinsky_Kh_;
  const bool bound = bound_Kh_;
  const amrex::Real Kh_bg_min = Kh_bg_min_;
  constexpr amrex::Real h_neglect3 = H_NEGLECT * H_NEGLECT * H_NEGLECT;

  // Per-layer work fields.
  amrex::MultiFab sh_xx = domain.make_field(Stagger::Cell, 1, 1);
  amrex::MultiFab sh_xy = domain.make_field(Stagger::Node, 1, 1);
  amrex::MultiFab h_u = domain.make_field(Stagger::XFace, 1, 1);
  amrex::MultiFab h_v = domain.make_field(Stagger::YFace, 1, 1);
  amrex::MultiFab str_xx = domain.make_field(Stagger::Cell, 1, 1);
  amrex::MultiFab str_xy = domain.make_field(Stagger::Node, 1, 1);
  sh_xx.setVal(0.0);
  sh_xy.setVal(0.0);
  h_u.setVal(0.0);
  h_v.setVal(0.0);
  str_xx.setVal(0.0);
  str_xy.setVal(0.0);

  for (int k = 0; k < nk; ++k) {
    for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
      const amrex::Box valid = mfi.validbox();

      const amrex::Array4<const amrex::Real> uu = u.const_array(mfi);
      const amrex::Array4<const amrex::Real> vv = v.const_array(mfi);
      const amrex::Array4<const amrex::Real> hh = h.const_array(mfi);
      const amrex::Array4<const amrex::Real> maskT = grid.mask2dT().const_array(mfi);
      const amrex::Array4<const amrex::Real> maskBu = grid.mask2dBu().const_array(mfi);
      const amrex::Array4<const amrex::Real> IdyCu = grid.IdyCu().const_array(mfi);
      const amrex::Array4<const amrex::Real> IdxCu = grid.IdxCu().const_array(mfi);
      const amrex::Array4<const amrex::Real> IdxCv = grid.IdxCv().const_array(mfi);
      const amrex::Array4<const amrex::Real> IdyCv = grid.IdyCv().const_array(mfi);
      const amrex::Array4<const amrex::Real> DY_dxT = DY_dxT_.const_array(mfi);
      const amrex::Array4<const amrex::Real> DX_dyT = DX_dyT_.const_array(mfi);
      const amrex::Array4<const amrex::Real> DY_dxBu = DY_dxBu_.const_array(mfi);
      const amrex::Array4<const amrex::Real> DX_dyBu = DX_dyBu_.const_array(mfi);

      // The horizontal tension at h points.
      const amrex::Array4<amrex::Real> shxx = sh_xx.array(mfi);
      amrex::ParallelFor(loops::flat(loops::h_points_grown(valid, 2)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        const amrex::Real dudx = DY_dxT(i, j, 0) * ((IdyCu(i + 1, j, 0) * uu(i + 1, j, k)) -
                                                    (IdyCu(i, j, 0) * uu(i, j, k)));
        const amrex::Real dvdy = DX_dyT(i, j, 0) * ((IdxCv(i, j + 1, 0) * vv(i, j + 1, k)) -
                                                    (IdxCv(i, j, 0) * vv(i, j, k)));
        shxx(i, j, 0) = dudx - dvdy;
      });

      // The shearing strain at q points, with the free-slip mask.
      const amrex::Array4<amrex::Real> shxy = sh_xy.array(mfi);
      amrex::ParallelFor(loops::flat(loops::q_points(valid, 2)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        const amrex::Real dvdx = DY_dxBu(i, j, 0) * ((vv(i, j, k) * IdyCv(i, j, 0)) -
                                                     (vv(i - 1, j, k) * IdyCv(i - 1, j, 0)));
        const amrex::Real dudy = DX_dyBu(i, j, 0) * ((uu(i, j, k) * IdxCu(i, j, 0)) -
                                                     (uu(i, j - 1, k) * IdxCu(i, j - 1, 0)));
        shxy(i, j, 0) = maskBu(i, j, 0) * (dvdx + dudy);
      });

      // The thicknesses at the velocity points, masked so that land values
      // never enter.
      const amrex::Array4<amrex::Real> hu = h_u.array(mfi);
      const amrex::Array4<amrex::Real> hv = h_v.array(mfi);
      amrex::ParallelFor(loops::flat(loops::u_points(valid, 2)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        hu(i, j, 0) = 0.5 * (maskT(i - 1, j, 0) * hh(i - 1, j, k) +
                             maskT(i, j, 0) * hh(i, j, k));
      });
      amrex::ParallelFor(loops::flat(loops::v_points(valid, 2)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        hv(i, j, 0) = 0.5 * (maskT(i, j - 1, 0) * hh(i, j - 1, k) +
                             maskT(i, j, 0) * hh(i, j, k));
      });

      // The viscosity and stress at h points.
      const amrex::Array4<amrex::Real> sxx = str_xx.array(mfi);
      const amrex::Array4<const amrex::Real> Kh_bg_xx = Kh_bg_xx_.const_array(mfi);
      const amrex::Array4<const amrex::Real> Lap2_xx = Laplac2_const_xx_.const_array(mfi);
      const amrex::Array4<const amrex::Real> KhMax_xx = Kh_Max_xx_.const_array(mfi);
      const amrex::Array4<const amrex::Real> red_xx = reduction_xx_.const_array(mfi);
      amrex::ParallelFor(loops::flat(loops::h_points_grown(valid, 1)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        amrex::Real Kh = Kh_bg_xx(i, j, 0);
        if (smag) {
          const amrex::Real sh_xx_sq = shxx(i, j, 0) * shxx(i, j, 0);
          const amrex::Real sh_xy_sq = 0.25 *
              (((shxy(i, j, 0) * shxy(i, j, 0)) + (shxy(i + 1, j + 1, 0) * shxy(i + 1, j + 1, 0))) +
               ((shxy(i, j + 1, 0) * shxy(i, j + 1, 0)) + (shxy(i + 1, j, 0) * shxy(i + 1, j, 0))));
          const amrex::Real Shear_mag = std::sqrt(sh_xx_sq + sh_xy_sq);
          Kh = amrex::max(Kh, Lap2_xx(i, j, 0) * Shear_mag);
        }
        Kh = amrex::max(Kh, Kh_bg_min);
        if (bound) {
          const amrex::Real h_min = amrex::min(amrex::min(hu(i + 1, j, 0), hu(i, j, 0)),
                                               amrex::min(hv(i, j + 1, 0), hv(i, j, 0)));
          const amrex::Real hrat_min =
              amrex::min(amrex::Real(1.0), h_min / (hh(i, j, k) + H_NEGLECT));
          Kh = amrex::min(Kh, hrat_min * KhMax_xx(i, j, 0));
        }
        sxx(i, j, 0) = (-Kh * shxx(i, j, 0)) * (hh(i, j, k) * red_xx(i, j, 0));
      });

      // The viscosity and stress at q points.
      const amrex::Array4<amrex::Real> sxy = str_xy.array(mfi);
      const amrex::Array4<const amrex::Real> Kh_bg_xy = Kh_bg_xy_.const_array(mfi);
      const amrex::Array4<const amrex::Real> Lap2_xy = Laplac2_const_xy_.const_array(mfi);
      const amrex::Array4<const amrex::Real> KhMax_xy = Kh_Max_xy_.const_array(mfi);
      const amrex::Array4<const amrex::Real> red_xy = reduction_xy_.const_array(mfi);
      amrex::ParallelFor(loops::flat(loops::q_points(valid)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        const amrex::Real h2uq = 4.0 * (hu(i, j - 1, 0) * hu(i, j, 0));
        const amrex::Real h2vq = 4.0 * (hv(i - 1, j, 0) * hv(i, j, 0));
        const amrex::Real hq = (2.0 * (h2uq * h2vq)) /
            (h_neglect3 + (h2uq + h2vq) * ((hu(i, j - 1, 0) + hu(i, j, 0)) +
                                           (hv(i - 1, j, 0) + hv(i, j, 0))));

        amrex::Real Kh = Kh_bg_xy(i, j, 0);
        if (smag) {
          const amrex::Real sh_xy_sq = shxy(i, j, 0) * shxy(i, j, 0);
          const amrex::Real sh_xx_sq = 0.25 *
              (((shxx(i - 1, j - 1, 0) * shxx(i - 1, j - 1, 0)) + (shxx(i, j, 0) * shxx(i, j, 0))) +
               ((shxx(i - 1, j, 0) * shxx(i - 1, j, 0)) + (shxx(i, j - 1, 0) * shxx(i, j - 1, 0))));
          const amrex::Real Shear_mag = std::sqrt(sh_xy_sq + sh_xx_sq);
          Kh = amrex::max(Kh, Lap2_xy(i, j, 0) * Shear_mag);
        }
        if (bound) {
          const amrex::Real h_min = amrex::min(amrex::min(hu(i, j - 1, 0), hu(i, j, 0)),
                                               amrex::min(hv(i - 1, j, 0), hv(i, j, 0)));
          const amrex::Real hrat_min = amrex::min(amrex::Real(1.0), h_min / (hq + H_NEGLECT));
          Kh = amrex::min(Kh, hrat_min * KhMax_xy(i, j, 0));
        }
        sxy(i, j, 0) = (-Kh * shxy(i, j, 0)) * (hq * maskBu(i, j, 0) * red_xy(i, j, 0));
      });

      // The divergence of the stress tensor.
      const amrex::Array4<amrex::Real> du = diffu.array(mfi);
      const amrex::Array4<const amrex::Real> dx2q = dx2q_.const_array(mfi);
      const amrex::Array4<const amrex::Real> dy2h = dy2h_.const_array(mfi);
      const amrex::Array4<const amrex::Real> IareaCu = grid.IareaCu().const_array(mfi);
      amrex::ParallelFor(loops::flat(loops::u_points(valid)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        du(i, j, k) = ((IdxCu(i, j, 0) * ((dx2q(i, j, 0) * sxy(i, j, 0)) -
                                          (dx2q(i, j + 1, 0) * sxy(i, j + 1, 0))) +
                        IdyCu(i, j, 0) * ((dy2h(i - 1, j, 0) * sxx(i - 1, j, 0)) -
                                          (dy2h(i, j, 0) * sxx(i, j, 0)))) *
                       IareaCu(i, j, 0)) / (hu(i, j, 0) + H_NEGLECT);
      });

      const amrex::Array4<amrex::Real> dv = diffv.array(mfi);
      const amrex::Array4<const amrex::Real> dy2q = dy2q_.const_array(mfi);
      const amrex::Array4<const amrex::Real> dx2h = dx2h_.const_array(mfi);
      const amrex::Array4<const amrex::Real> IareaCv = grid.IareaCv().const_array(mfi);
      amrex::ParallelFor(loops::flat(loops::v_points(valid)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        dv(i, j, k) = ((IdyCv(i, j, 0) * ((dy2q(i, j, 0) * sxy(i, j, 0)) -
                                          (dy2q(i + 1, j, 0) * sxy(i + 1, j, 0))) -
                        IdxCv(i, j, 0) * ((dx2h(i, j - 1, 0) * sxx(i, j - 1, 0)) -
                                          (dx2h(i, j, 0) * sxx(i, j, 0)))) *
                       IareaCv(i, j, 0)) / (hv(i, j, 0) + H_NEGLECT);
      });
    }
  }

  diffu.FillBoundary(domain.periodicity());
  diffv.FillBoundary(domain.periodicity());
}

} // namespace MOM
