#include <cmath>

#include "MOM_set_viscosity.h"

#include "MOM_logger.h"
#include "MOM_loop_boxes.h"

namespace MOM {

namespace {

constexpr amrex::Real H_NEGLECT = 1.0e-30;

// Reject a parameter that selects an unimplemented branch.
void reject_if_set(RuntimeParams &params, const std::string &key, const bool default_value,
                   const std::string &desc) {
  bool value = default_value;
  params.get(key, value, {.default_value = default_value, .desc = desc});
  if (value != default_value) {
    logger::fatal("SetViscosity: ", key, " is not implemented yet.");
  }
}

} // namespace

VertVisc::VertVisc(const Domain &domain)
  : Kv_bbl_u(domain.make_field(Stagger::XFace, 1, 1)),
    Kv_bbl_v(domain.make_field(Stagger::YFace, 1, 1)),
    bbl_thick_u(domain.make_field(Stagger::XFace, 1, 1)),
    bbl_thick_v(domain.make_field(Stagger::YFace, 1, 1)) {
  Kv_bbl_u.setVal(0.0);
  Kv_bbl_v.setVal(0.0);
  bbl_thick_u.setVal(0.0);
  bbl_thick_v.setVal(0.0);
}

SetViscosity::SetViscosity(RuntimeParams &params) {

  params.doc_module("MOM_set_visc", "");

  bool bottomdraglaw = true;
  params.get("BOTTOMDRAGLAW", bottomdraglaw,
             {.default_value = true,
              .desc = "If true, the bottom stress is calculated with a drag law of the form "
                      "c_drag*|u|*u. The velocity magnitude may be an assumed value or it may "
                      "be based on the actual velocity in the bottommost HBBL, depending on "
                      "LINEAR_DRAG."});
  if (!bottomdraglaw) {
    // defer: the no-bottom-drag-law path.
    logger::fatal("SetViscosity: BOTTOMDRAGLAW = False is not implemented yet.");
  }

  reject_if_set(params, "DRAG_AS_BODY_FORCE", false,
                "If true, the bottom stress is imposed as an explicit body force applied over a "
                "fixed distance from the bottom, rather than as an implied stress.");
  reject_if_set(params, "CHANNEL_DRAG", false,
                "If true, the bottom drag is exerted directly on each layer proportional to the "
                "fraction of the bottom it overlies.");

  bool linear_drag = false;
  params.get("LINEAR_DRAG", linear_drag,
             {.default_value = false,
              .desc = "If LINEAR_DRAG and BOTTOMDRAGLAW are defined the drag law is "
                      "cdrag*DRAG_BG_VEL*u."});
  if (!linear_drag) {
    // defer: the velocity-dependent drag, which needs the near-bottom
    //        velocity average over the boundary layer.
    logger::fatal("SetViscosity: LINEAR_DRAG = False is not implemented yet.");
  }

  reject_if_set(params, "DYNAMIC_VISCOUS_ML", false,
                "If true, use a bulk Richardson number criterion to determine the mixed layer "
                "thickness for viscosity.");

  params.get("HBBL", Hbbl_,
             {.desc = "The thickness of a bottom boundary layer with a viscosity increased by "
                      "KV_EXTRA_BBL if BOTTOMDRAGLAW is not defined, or the thickness over "
                      "which near-bottom velocities are averaged for the drag law if "
                      "BOTTOMDRAGLAW is defined but LINEAR_DRAG is not.",
              .units = "m",
              .fail_if_missing = true});

  params.get("CDRAG", cdrag_,
             {.default_value = 0.003,
              .desc = "CDRAG is the drag coefficient relating the magnitude of the velocity "
                      "field to the bottom stress.",
              .units = "nondim"});

  params.get("DRAG_BG_VEL", drag_bg_vel_,
             {.default_value = 0.0,
              .desc = "DRAG_BG_VEL is either the assumed bottom velocity (with LINEAR_DRAG) or "
                      "an unresolved velocity that is combined with the resolved velocity to "
                      "estimate the velocity magnitude.",
              .units = "m s-1"});

  params.get("BBL_THICK_MIN", BBL_thick_min_,
             {.default_value = 0.0,
              .desc = "The minimum bottom boundary layer thickness that can be used with "
                      "BOTTOMDRAGLAW. This might be Kv/(cdrag*drag_bg_vel) to give Kv as the "
                      "minimum near-bottom viscosity.",
              .units = "m"});

  params.get("KV", Kv_,
             {.desc = "The background kinematic viscosity in the interior. The molecular value, "
                      "~1e-6 m2 s-1, may be used.",
              .units = "m2 s-1",
              .fail_if_missing = true});

  params.get("KV_BBL_MIN", Kv_BBL_min_,
             {.default_value = 1.0e-4,
              .desc = "The minimum viscosities in the bottom boundary layer.",
              .units = "m2 s-1"});

  reject_if_set(params, "CORRECT_BBL_BOUNDS", false,
                "If true, uses the correct bounds on the BBL thickness and viscosity so that "
                "the bottom layer feels the intended drag.");
}

void SetViscosity::set_viscous_BBL(VertVisc &visc, const amrex::MultiFab &u,
                                   const amrex::MultiFab &v, const amrex::MultiFab &h,
                                   const Domain &domain, const Grid &grid,
                                   const VerticalGrid &vgrid) const {

  const int nk = vgrid.nk();
  const amrex::Real cdrag_sqrt = std::sqrt(cdrag_);
  // With LINEAR_DRAG the friction velocity does not depend on the flow.
  const amrex::Real ustar = cdrag_sqrt * drag_bg_vel_;
  const amrex::Real Rho0x400_G = 400.0 * (vgrid.Rho0() / vgrid.g_Earth());
  const amrex::Real ustarsq = Rho0x400_G * ustar * ustar;
  const amrex::Real BBL_thick_min = BBL_thick_min_;
  const amrex::Real Kv_BBL_min = Kv_BBL_min_;
  // cdrag * u2_bg is positive here, which selects MOM6's second bbl_thick form.
  const bool drag_positive = (cdrag_ * drag_bg_vel_ * drag_bg_vel_ > 0.0);
  if (!drag_positive) {
    logger::fatal("SetViscosity: DRAG_BG_VEL must be positive with LINEAR_DRAG.");
  }

  amrex::Gpu::DeviceVector<amrex::Real> Rlay_d(vgrid.Rlay().size());
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, vgrid.Rlay().begin(), vgrid.Rlay().end(),
                   Rlay_d.begin());
  const amrex::Real *Rlay = Rlay_d.data();

  for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
    const amrex::Box valid = mfi.validbox();
    const amrex::Array4<const amrex::Real> hh = h.const_array(mfi);
    const amrex::Array4<const amrex::Real> uu = u.const_array(mfi);
    const amrex::Array4<const amrex::Real> vv = v.const_array(mfi);
    const amrex::Array4<const amrex::Real> maskCu = grid.mask2dCu().const_array(mfi);
    const amrex::Array4<const amrex::Real> maskCv = grid.mask2dCv().const_array(mfi);
    const amrex::Array4<const amrex::Real> f = grid.CoriolisBu().const_array(mfi);

    const amrex::Array4<amrex::Real> kv_u = visc.Kv_bbl_u.array(mfi);
    const amrex::Array4<amrex::Real> th_u = visc.bbl_thick_u.array(mfi);
    amrex::ParallelFor(loops::flat(loops::u_points(valid)), [=] AMREX_GPU_DEVICE(int i, int j, int) {
      if (maskCu(i, j, 0) <= 0.0) return;
      // Accumulate the thickness of the layers the bottom stress can stir,
      // stopping where the density difference across them balances the stress.
      amrex::Real Rhtot = 0.0;
      amrex::Real htot = 0.0;
      for (int k = nk - 1; k >= 1; --k) {
        amrex::Real h_at_vel;
        if (uu(i, j, k) * (hh(i, j, k) - hh(i - 1, j, k)) >= 0.0) {
          h_at_vel = 2.0 * hh(i - 1, j, k) * hh(i, j, k) /
                     (hh(i - 1, j, k) + hh(i, j, k) + H_NEGLECT);
        } else {
          h_at_vel = 0.5 * (hh(i - 1, j, k) + hh(i, j, k));
        }
        const amrex::Real oldfn = Rhtot - Rlay[k] * htot;
        if (oldfn >= ustarsq) continue;
        const amrex::Real Dfn = (Rlay[k] - Rlay[k - 1]) * (h_at_vel + htot);
        amrex::Real Dh;
        if ((oldfn + Dfn) <= ustarsq) {
          Dh = h_at_vel;
        } else {
          Dh = h_at_vel * std::sqrt((ustarsq - oldfn) / Dfn);
        }
        htot += Dh;
        Rhtot += Rlay[k] * Dh;
      }
      if (Rhtot - Rlay[0] * htot < ustarsq) {
        amrex::Real h_at_vel;
        if (uu(i, j, 0) * (hh(i, j, 0) - hh(i - 1, j, 0)) >= 0.0) {
          h_at_vel = 2.0 * hh(i - 1, j, 0) * hh(i, j, 0) /
                     (hh(i - 1, j, 0) + hh(i, j, 0) + H_NEGLECT);
        } else {
          h_at_vel = 0.5 * (hh(i - 1, j, 0) + hh(i, j, 0));
        }
        htot += h_at_vel;
      }

      // The Ekman-like thickness, bounded below.
      const amrex::Real C2f = f(i, j, 0) + f(i, j + 1, 0);
      amrex::Real bbl_thick =
          htot / (0.5 + std::sqrt(0.25 + htot * htot * C2f * C2f / (ustar * ustar)));
      if (bbl_thick < BBL_thick_min) bbl_thick = BBL_thick_min;

      const amrex::Real kv_bbl =
          amrex::max(Kv_BBL_min, (cdrag_sqrt * ustar) * bbl_thick);
      th_u(i, j, 0) = bbl_thick;
      kv_u(i, j, 0) = kv_bbl;
    });

    const amrex::Array4<amrex::Real> kv_v = visc.Kv_bbl_v.array(mfi);
    const amrex::Array4<amrex::Real> th_v = visc.bbl_thick_v.array(mfi);
    amrex::ParallelFor(loops::flat(loops::v_points(valid)), [=] AMREX_GPU_DEVICE(int i, int j, int) {
      if (maskCv(i, j, 0) <= 0.0) return;
      amrex::Real Rhtot = 0.0;
      amrex::Real htot = 0.0;
      for (int k = nk - 1; k >= 1; --k) {
        amrex::Real h_at_vel;
        if (vv(i, j, k) * (hh(i, j, k) - hh(i, j - 1, k)) >= 0.0) {
          h_at_vel = 2.0 * hh(i, j - 1, k) * hh(i, j, k) /
                     (hh(i, j - 1, k) + hh(i, j, k) + H_NEGLECT);
        } else {
          h_at_vel = 0.5 * (hh(i, j - 1, k) + hh(i, j, k));
        }
        const amrex::Real oldfn = Rhtot - Rlay[k] * htot;
        if (oldfn >= ustarsq) continue;
        const amrex::Real Dfn = (Rlay[k] - Rlay[k - 1]) * (h_at_vel + htot);
        amrex::Real Dh;
        if ((oldfn + Dfn) <= ustarsq) {
          Dh = h_at_vel;
        } else {
          Dh = h_at_vel * std::sqrt((ustarsq - oldfn) / Dfn);
        }
        htot += Dh;
        Rhtot += Rlay[k] * Dh;
      }
      if (Rhtot - Rlay[0] * htot < ustarsq) {
        amrex::Real h_at_vel;
        if (vv(i, j, 0) * (hh(i, j, 0) - hh(i, j - 1, 0)) >= 0.0) {
          h_at_vel = 2.0 * hh(i, j - 1, 0) * hh(i, j, 0) /
                     (hh(i, j - 1, 0) + hh(i, j, 0) + H_NEGLECT);
        } else {
          h_at_vel = 0.5 * (hh(i, j - 1, 0) + hh(i, j, 0));
        }
        htot += h_at_vel;
      }

      const amrex::Real C2f = f(i, j, 0) + f(i + 1, j, 0);
      amrex::Real bbl_thick =
          htot / (0.5 + std::sqrt(0.25 + htot * htot * C2f * C2f / (ustar * ustar)));
      if (bbl_thick < BBL_thick_min) bbl_thick = BBL_thick_min;

      const amrex::Real kv_bbl =
          amrex::max(Kv_BBL_min, (cdrag_sqrt * ustar) * bbl_thick);
      th_v(i, j, 0) = bbl_thick;
      kv_v(i, j, 0) = kv_bbl;
    });
  }

  visc.Kv_bbl_u.FillBoundary(domain.periodicity());
  visc.Kv_bbl_v.FillBoundary(domain.periodicity());
  visc.bbl_thick_u.FillBoundary(domain.periodicity());
  visc.bbl_thick_v.FillBoundary(domain.periodicity());
}

} // namespace MOM
