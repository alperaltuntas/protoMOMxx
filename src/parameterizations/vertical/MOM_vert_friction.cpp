#include <cmath>

#include "MOM_vert_friction.h"

#include "MOM_fields.h"
#include "MOM_fp_contract.h"
#include "MOM_kernel_inline.h"
#include "MOM_logger.h"
#include "MOM_loop_boxes.h"

namespace MOM {

namespace {


// The coupling coefficient is capped in MOM6; with the current answer date
// the cap is effectively absent (I_amax = 0), so only a_cpl_max remains.
constexpr amrex::Real A_CPL_MAX = 1.0e37;

// 1 / (1 + 0.09 z^6), MOM6's near-boundary blending function.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real boundary_fn(const amrex::Real z) {
  return 1.0 / (1.0 + 0.09 * z * z * z * z * z * z);
}

void reject_if_set(RuntimeParams &params, const std::string &key, const bool default_value,
                   const std::string &desc) {
  bool value = default_value;
  params.get(key, value, {.default_value = default_value, .desc = desc});
  if (value != default_value) {
    logger::fatal("VertFriction: ", key, " is not implemented yet.");
  }
}

} // namespace

VertFriction::VertFriction(RuntimeParams &params, const Domain &domain,
                           const VerticalGrid &vgrid) {

  params.doc_module("MOM_vert_friction", "");

  params.get("DIRECT_STRESS", direct_stress_,
             {.default_value = false,
              .desc = "If true, the wind stress is distributed over the topmost HMIX_STRESS of "
                      "fluid (like in HYCOM), and KVML may be set to a very small value."});

  reject_if_set(params, "FIXED_DEPTH_LOTW_ML", false,
                "If true, use a Law-of-the-wall prescription for the mixed layer viscosity "
                "within a boundary layer of a fixed thickness.");
  reject_if_set(params, "LOTW_VISCOUS_ML_FLOOR", false,
                "If true, use a Law-of-the-wall prescription to set a lower bound on the "
                "viscous coupling between layers within the surface boundary layer.");

  bool harmonic_visc = false;
  params.get("HARMONIC_VISC", harmonic_visc,
             {.default_value = false,
              .desc = "If true, use the harmonic mean thicknesses for calculating the vertical "
                      "viscosity."});
  if (!harmonic_visc) {
    // defer: the arithmetic-mean thickness path and its z_clear bookkeeping.
    logger::fatal("VertFriction: HARMONIC_VISC = False is not implemented yet.");
  }

  params.get("HMIX_FIXED", Hmix_,
             {.desc = "The prescribed depth over which the near-surface viscosity and "
                      "diffusivity are elevated when the bulk mixed layer is not used.",
              .units = "m",
              .fail_if_missing = true});

  Hmix_stress_ = Hmix_;
  if (direct_stress_) {
    params.get("HMIX_STRESS", Hmix_stress_,
               {.default_value = Hmix_,
                .desc = "The depth over which the wind stress is applied if DIRECT_STRESS is "
                        "true.",
                .units = "m"});
  }

  reject_if_set(params, "USE_GL90_IN_SSW", false,
                "If true, use the GL90 parameterization for vertical momentum transport.");

  params.get("KV_ML_INVZ2", Kvml_invZ2_,
             {.default_value = 0.01,
              .desc = "An extra kinematic viscosity in a region near the top of the ocean, "
                      "typically of the mixed layer, that decays away with the square of the "
                      "distance from the surface.",
              .units = "m2 s-1"});

  params.get("CFL_TRUNCATE", CFL_trunc_,
             {.default_value = 0.5,
              .desc = "The value of the CFL number that will cause velocity components to be "
                      "truncated; instability can occur past 0.5.",
              .units = "nondim"});

  amrex::Real ramp_time = 0.0;
  params.get("CFL_TRUNCATE_RAMP_TIME", ramp_time,
             {.default_value = 0.0,
              .desc = "The time over which the CFL truncation value is ramped up at the "
                      "beginning of the run.",
              .units = "s"});
  if (ramp_time > 0.0) {
    // defer: the ramped truncation threshold.
    logger::fatal("VertFriction: CFL_TRUNCATE_RAMP_TIME > 0 is not implemented yet.");
  }

  params.get("VEL_UNDERFLOW", vel_underflow_,
             {.default_value = 0.0,
              .desc = "A negligibly small velocity magnitude below which velocity components "
                      "are set to 0.  A reasonable value might be 1e-30 m/s, which is 1e-47 "
                      "of the speed of light.",
              .units = "m s-1"});

  a_u_ = make_interface_field(domain, vgrid, Stagger::XFace);
  a_v_ = make_interface_field(domain, vgrid, Stagger::YFace);
  h_u_ = make_layer_field(domain, vgrid, Stagger::XFace);
  h_v_ = make_layer_field(domain, vgrid, Stagger::YFace);
  a_u_.setVal(0.0);
  a_v_.setVal(0.0);
  h_u_.setVal(0.0);
  h_v_.setVal(0.0);
}

void VertFriction::coefficients(const amrex::MultiFab &u, const amrex::MultiFab &v,
                                const amrex::MultiFab &h, const VertVisc &visc,
                                const amrex::Real Kv_interior, const Domain &domain,
                                const Grid &grid, const VerticalGrid &vgrid) {

  const int nk = vgrid.nk();
  const amrex::Real Kv = Kv_interior;
  // MOM6's GV%H_subroundoff; see VerticalGrid::H_subroundoff.
  const amrex::Real H_NEGLECT = vgrid.H_subroundoff();
  const amrex::Real Kvml_invZ2 = Kvml_invZ2_;
  const amrex::Real I_Hmix = 1.0 / (Hmix_ + H_NEGLECT);

  for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
    const amrex::Box valid = mfi.validbox();
    const amrex::Array4<const amrex::Real> hh = h.const_array(mfi);
    const amrex::Array4<const amrex::Real> uu = u.const_array(mfi);
    const amrex::Array4<const amrex::Real> vv = v.const_array(mfi);
    const amrex::Array4<const amrex::Real> maskCu = grid.mask2dCu().const_array(mfi);
    const amrex::Array4<const amrex::Real> maskCv = grid.mask2dCv().const_array(mfi);

    // --- u points ---
    const amrex::Array4<const amrex::Real> kv_bbl_u = visc.Kv_bbl_u.const_array(mfi);
    const amrex::Array4<const amrex::Real> bbl_u = visc.bbl_thick_u.const_array(mfi);
    const amrex::Array4<amrex::Real> au = a_u_.array(mfi);
    const amrex::Array4<amrex::Real> hu = h_u_.array(mfi);
    amrex::ParallelFor(loops::flat(loops::u_points(valid)), [=] AMREX_GPU_DEVICE(int i, int j, int) MOM_KERNEL_INLINE {
      if (maskCu(i, j, 0) <= 0.0) return;
      const amrex::Real kv_bbl = kv_bbl_u(i, j, 0);
      const amrex::Real bbl_thick = bbl_u(i, j, 0) + H_NEGLECT;
      const amrex::Real I_Hbbl = 1.0 / bbl_thick;

      // z_i is the height above the bottom in units of the boundary layer
      // thickness, built upward, so it is stored in a_u's spare interface
      // slots while the coefficients are being formed.
      amrex::Real z_i_kp1 = 0.0;
      for (int k = nk - 1; k >= 0; --k) {
        const amrex::Real h_harm = 2.0 * hh(i - 1, j, k) * hh(i, j, k) /
                                   (hh(i - 1, j, k) + hh(i, j, k) + H_NEGLECT);
        const amrex::Real h_arith = 0.5 * (hh(i, j, k) + hh(i - 1, j, k));
        const amrex::Real h_delta = hh(i, j, k) - hh(i - 1, j, k);
        amrex::Real hvel = h_harm;
        if (uu(i, j, k) * h_delta < 0.0) {
          const amrex::Real botfn = boundary_fn(z_i_kp1);
          hvel = (1.0 - botfn) * h_harm + botfn * h_arith;
        }
        hu(i, j, k) = hvel + H_NEGLECT;
        // Stash the harmonic thickness and z_i in the interface array; the
        // second pass below turns them into coupling coefficients.
        au(i, j, k) = h_harm;
        z_i_kp1 = z_i_kp1 + h_harm * I_Hbbl;
      }

      // Rebuild z_i on the way down and form the coupling coefficients.
      amrex::Real z_i[64];
      z_i[nk] = 0.0;
      for (int k = nk - 1; k >= 0; --k) {
        z_i[k] = z_i[k + 1] + au(i, j, k) * I_Hbbl;
      }
      amrex::Real h_harm_arr[64];
      for (int k = 0; k < nk; ++k) h_harm_arr[k] = au(i, j, k);

      amrex::Real z_t = H_NEGLECT * I_Hmix;
      au(i, j, 0) = 0.0;
      for (int k = 1; k < nk; ++k) {
        z_t = z_t + h_harm_arr[k - 1] * I_Hmix;
        amrex::Real Kv_tot = Kv + Kvml_invZ2 / ((z_t * z_t) * (1.0 + 0.09 * z_t * z_t * z_t *
                                                               z_t * z_t * z_t));
        const amrex::Real botfn = boundary_fn(z_i[k]);
        Kv_tot = Kv_tot + (kv_bbl - Kv) * botfn;
        const amrex::Real dhc = 0.5 * (hu(i, j, k) - H_NEGLECT + hu(i, j, k - 1) - H_NEGLECT);
        const amrex::Real h_shear = (dhc > bbl_thick)
                                        ? (((1.0 - botfn) * dhc + botfn * bbl_thick) + H_NEGLECT)
                                        : (dhc + H_NEGLECT);
        au(i, j, k) = amrex::min(A_CPL_MAX, Kv_tot / h_shear);
      }
      const amrex::Real dhc_bot = (hu(i, j, nk - 1) - H_NEGLECT) * 0.5;
      au(i, j, nk) = amrex::min(A_CPL_MAX,
                                kv_bbl / (amrex::min(dhc_bot, bbl_thick) + H_NEGLECT));
    });

    // --- v points ---
    const amrex::Array4<const amrex::Real> kv_bbl_v = visc.Kv_bbl_v.const_array(mfi);
    const amrex::Array4<const amrex::Real> bbl_v = visc.bbl_thick_v.const_array(mfi);
    const amrex::Array4<amrex::Real> av = a_v_.array(mfi);
    const amrex::Array4<amrex::Real> hv = h_v_.array(mfi);
    amrex::ParallelFor(loops::flat(loops::v_points(valid)), [=] AMREX_GPU_DEVICE(int i, int j, int) MOM_KERNEL_INLINE {
      if (maskCv(i, j, 0) <= 0.0) return;
      const amrex::Real kv_bbl = kv_bbl_v(i, j, 0);
      const amrex::Real bbl_thick = bbl_v(i, j, 0) + H_NEGLECT;
      const amrex::Real I_Hbbl = 1.0 / bbl_thick;

      amrex::Real z_i_kp1 = 0.0;
      for (int k = nk - 1; k >= 0; --k) {
        const amrex::Real h_harm = 2.0 * hh(i, j - 1, k) * hh(i, j, k) /
                                   (hh(i, j - 1, k) + hh(i, j, k) + H_NEGLECT);
        const amrex::Real h_arith = 0.5 * (hh(i, j, k) + hh(i, j - 1, k));
        const amrex::Real h_delta = hh(i, j, k) - hh(i, j - 1, k);
        amrex::Real hvel = h_harm;
        if (vv(i, j, k) * h_delta < 0.0) {
          const amrex::Real botfn = boundary_fn(z_i_kp1);
          hvel = (1.0 - botfn) * h_harm + botfn * h_arith;
        }
        hv(i, j, k) = hvel + H_NEGLECT;
        av(i, j, k) = h_harm;
        z_i_kp1 = z_i_kp1 + h_harm * I_Hbbl;
      }

      amrex::Real z_i[64];
      z_i[nk] = 0.0;
      for (int k = nk - 1; k >= 0; --k) {
        z_i[k] = z_i[k + 1] + av(i, j, k) * I_Hbbl;
      }
      amrex::Real h_harm_arr[64];
      for (int k = 0; k < nk; ++k) h_harm_arr[k] = av(i, j, k);

      amrex::Real z_t = H_NEGLECT * I_Hmix;
      av(i, j, 0) = 0.0;
      for (int k = 1; k < nk; ++k) {
        z_t = z_t + h_harm_arr[k - 1] * I_Hmix;
        amrex::Real Kv_tot = Kv + Kvml_invZ2 / ((z_t * z_t) * (1.0 + 0.09 * z_t * z_t * z_t *
                                                               z_t * z_t * z_t));
        const amrex::Real botfn = boundary_fn(z_i[k]);
        Kv_tot = Kv_tot + (kv_bbl - Kv) * botfn;
        const amrex::Real dhc = 0.5 * (hv(i, j, k) - H_NEGLECT + hv(i, j, k - 1) - H_NEGLECT);
        const amrex::Real h_shear = (dhc > bbl_thick)
                                        ? (((1.0 - botfn) * dhc + botfn * bbl_thick) + H_NEGLECT)
                                        : (dhc + H_NEGLECT);
        av(i, j, k) = amrex::min(A_CPL_MAX, Kv_tot / h_shear);
      }
      const amrex::Real dhc_bot = (hv(i, j, nk - 1) - H_NEGLECT) * 0.5;
      av(i, j, nk) = amrex::min(A_CPL_MAX,
                                kv_bbl / (amrex::min(dhc_bot, bbl_thick) + H_NEGLECT));
    });
  }

  (void)domain;
}

void VertFriction::apply(amrex::MultiFab &u, amrex::MultiFab &v, const amrex::MultiFab &h,
                         const MechForcing &forces, const amrex::Real dt,
                         const Domain &domain, const Grid &grid,
                         const VerticalGrid &vgrid) {

  const int nk = vgrid.nk();
  // H_to_RZ is Rho0 in Boussinesq mode without unit scaling.
  const amrex::Real dt_Rho0 = dt / vgrid.Rho0();
  const amrex::Real H_NEGLECT = vgrid.H_subroundoff();
  const amrex::Real Hmix = Hmix_stress_;
  const amrex::Real I_Hmix = 1.0 / Hmix;
  const bool direct_stress = direct_stress_;

  for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
    const amrex::Box valid = mfi.validbox();
    const amrex::Array4<const amrex::Real> hh = h.const_array(mfi);
    const amrex::Array4<const amrex::Real> maskCu = grid.mask2dCu().const_array(mfi);
    const amrex::Array4<const amrex::Real> maskCv = grid.mask2dCv().const_array(mfi);
    const amrex::Array4<const amrex::Real> taux = forces.taux().const_array(mfi);
    const amrex::Array4<const amrex::Real> tauy = forces.tauy().const_array(mfi);

    const amrex::Array4<amrex::Real> uu = u.array(mfi);
    const amrex::Array4<const amrex::Real> au = a_u_.const_array(mfi);
    const amrex::Array4<const amrex::Real> hu = h_u_.const_array(mfi);
    amrex::ParallelFor(loops::flat(loops::u_points(valid)), [=] AMREX_GPU_DEVICE(int i, int j, int) MOM_KERNEL_INLINE {
      if (maskCu(i, j, 0) <= 0.0) return;

      // The wind stress, either spread over the top HMIX_STRESS or applied at
      // the surface as a boundary condition on the tridiagonal solve.
      amrex::Real surface_stress = 0.0;
      if (direct_stress) {
        amrex::Real zDS = 0.0;
        const amrex::Real stress = dt_Rho0 * taux(i, j, 0);
        for (int k = 0; k < nk; ++k) {
          const amrex::Real h_a = 0.5 * (hh(i - 1, j, k) + hh(i, j, k)) + H_NEGLECT;
          amrex::Real hfr = 1.0;
          if ((zDS + h_a) > Hmix) hfr = (Hmix - zDS) / h_a;
          uu(i, j, k) = uu(i, j, k) + fp_rounded(I_Hmix * hfr * stress);
          zDS += h_a;
          if (zDS >= Hmix) break;
        }
      } else {
        surface_stress = dt_Rho0 * taux(i, j, 0);
      }

      // The tridiagonal solve, in MOM6's forward-elimination form.
      amrex::Real c1[64];
      amrex::Real b_denom_1 = hu(i, j, 0) + dt * au(i, j, 0);
      amrex::Real b1 = 1.0 / (b_denom_1 + dt * au(i, j, 1));
      amrex::Real d1 = b_denom_1 * b1;
      uu(i, j, 0) = b1 * (hu(i, j, 0) * uu(i, j, 0) + surface_stress);
      for (int k = 1; k < nk; ++k) {
        c1[k] = dt * au(i, j, k) * b1;
        b_denom_1 = hu(i, j, k) + dt * (au(i, j, k) * d1);
        b1 = 1.0 / (b_denom_1 + dt * au(i, j, k + 1));
        d1 = b_denom_1 * b1;
        uu(i, j, k) = (hu(i, j, k) * uu(i, j, k) + dt * au(i, j, k) * uu(i, j, k - 1)) * b1;
      }
      for (int k = nk - 2; k >= 0; --k) {
        uu(i, j, k) = uu(i, j, k) + c1[k + 1] * uu(i, j, k + 1);
      }
    });

    const amrex::Array4<amrex::Real> vv = v.array(mfi);
    const amrex::Array4<const amrex::Real> av = a_v_.const_array(mfi);
    const amrex::Array4<const amrex::Real> hv = h_v_.const_array(mfi);
    amrex::ParallelFor(loops::flat(loops::v_points(valid)), [=] AMREX_GPU_DEVICE(int i, int j, int) MOM_KERNEL_INLINE {
      if (maskCv(i, j, 0) <= 0.0) return;

      amrex::Real surface_stress = 0.0;
      if (direct_stress) {
        amrex::Real zDS = 0.0;
        const amrex::Real stress = dt_Rho0 * tauy(i, j, 0);
        for (int k = 0; k < nk; ++k) {
          const amrex::Real h_a = 0.5 * (hh(i, j - 1, k) + hh(i, j, k)) + H_NEGLECT;
          amrex::Real hfr = 1.0;
          if ((zDS + h_a) > Hmix) hfr = (Hmix - zDS) / h_a;
          vv(i, j, k) = vv(i, j, k) + fp_rounded(I_Hmix * hfr * stress);
          zDS += h_a;
          if (zDS >= Hmix) break;
        }
      } else {
        surface_stress = dt_Rho0 * tauy(i, j, 0);
      }

      amrex::Real c1[64];
      amrex::Real b_denom_1 = hv(i, j, 0) + dt * av(i, j, 0);
      amrex::Real b1 = 1.0 / (b_denom_1 + dt * av(i, j, 1));
      amrex::Real d1 = b_denom_1 * b1;
      vv(i, j, 0) = b1 * (hv(i, j, 0) * vv(i, j, 0) + surface_stress);
      for (int k = 1; k < nk; ++k) {
        c1[k] = dt * av(i, j, k) * b1;
        b_denom_1 = hv(i, j, k) + dt * (av(i, j, k) * d1);
        b1 = 1.0 / (b_denom_1 + dt * av(i, j, k + 1));
        d1 = b_denom_1 * b1;
        vv(i, j, k) = (hv(i, j, k) * vv(i, j, k) + dt * av(i, j, k) * vv(i, j, k - 1)) * b1;
      }
      for (int k = nk - 2; k >= 0; --k) {
        vv(i, j, k) = vv(i, j, k) + c1[k + 1] * vv(i, j, k + 1);
      }
    });
  }

  limit_velocity(u, v, dt, grid, vgrid);

  u.FillBoundary(domain.periodicity());
  v.FillBoundary(domain.periodicity());
}

void VertFriction::limit_velocity(amrex::MultiFab &u, amrex::MultiFab &v, const amrex::Real dt,
                                  const Grid &grid, const VerticalGrid &vgrid) {

  const int nk = vgrid.nk();
  const amrex::Real CFL_trunc = CFL_trunc_;
  const amrex::Real vel_underflow = vel_underflow_;
  // Truncations in layers thinner than this are not counted, because they
  // carry no momentum worth reporting.
  const amrex::Real H_report = 3.0 * vgrid.angstrom();

  for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
    const amrex::Box valid = mfi.validbox();
    const amrex::Array4<amrex::Real> uu = u.array(mfi);
    const amrex::Array4<const amrex::Real> hu = h_u_.const_array(mfi);
    const amrex::Array4<const amrex::Real> dy_Cu = grid.dy_Cu().const_array(mfi);
    const amrex::Array4<const amrex::Real> areaT = grid.areaT().const_array(mfi);
    const amrex::Array4<const amrex::Real> IareaT = grid.IareaT().const_array(mfi);
    const amrex::Box bx = loops::u_points(valid);
    for (int k = 0; k < nk; ++k) {
      for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
          if (std::abs(uu(i, j, k)) < vel_underflow) {
            uu(i, j, k) = 0.0;
          } else if ((uu(i, j, k) * (dt * dy_Cu(i, j, 0))) * IareaT(i, j, 0) < -CFL_trunc) {
            uu(i, j, k) = (-0.9 * CFL_trunc) * (areaT(i, j, 0) / (dt * dy_Cu(i, j, 0)));
            if (hu(i, j, k) > H_report) ++ntrunc_;
          } else if ((uu(i, j, k) * (dt * dy_Cu(i, j, 0))) * IareaT(i - 1, j, 0) > CFL_trunc) {
            uu(i, j, k) = (0.9 * CFL_trunc) * (areaT(i - 1, j, 0) / (dt * dy_Cu(i, j, 0)));
            if (hu(i, j, k) > H_report) ++ntrunc_;
          }
        }
      }
    }
  }

  for (amrex::MFIter mfi(v); mfi.isValid(); ++mfi) {
    const amrex::Box valid = mfi.validbox();
    const amrex::Array4<amrex::Real> vv = v.array(mfi);
    const amrex::Array4<const amrex::Real> hv = h_v_.const_array(mfi);
    const amrex::Array4<const amrex::Real> dx_Cv = grid.dx_Cv().const_array(mfi);
    const amrex::Array4<const amrex::Real> areaT = grid.areaT().const_array(mfi);
    const amrex::Array4<const amrex::Real> IareaT = grid.IareaT().const_array(mfi);
    const amrex::Box bx = loops::v_points(valid);
    for (int k = 0; k < nk; ++k) {
      for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
          if (std::abs(vv(i, j, k)) < vel_underflow) {
            vv(i, j, k) = 0.0;
          } else if ((vv(i, j, k) * (dt * dx_Cv(i, j, 0))) * IareaT(i, j, 0) < -CFL_trunc) {
            vv(i, j, k) = (-0.9 * CFL_trunc) * (areaT(i, j, 0) / (dt * dx_Cv(i, j, 0)));
            if (hv(i, j, k) > H_report) ++ntrunc_;
          } else if ((vv(i, j, k) * (dt * dx_Cv(i, j, 0))) * IareaT(i, j - 1, 0) > CFL_trunc) {
            vv(i, j, k) = (0.9 * CFL_trunc) * (areaT(i, j - 1, 0) / (dt * dx_Cv(i, j, 0)));
            if (hv(i, j, k) > H_report) ++ntrunc_;
          }
        }
      }
    }
  }
}

} // namespace MOM
