#include <string>

#include "MOM_CoriolisAdv.h"

#include "MOM_logger.h"
#include "MOM_loop_boxes.h"

namespace MOM {

namespace {

// A volume so small that it is expected to be lost in roundoff
// [H L2 ~> m3]: GV%H_subroundoff * (1e-4 m)^2 with H_subroundoff = 1e-30 at
// the default ANGSTROM.
constexpr amrex::Real VOL_NEGLECT = 1.0e-30 * (1.0e-4 * 1.0e-4);

} // namespace

CoriolisAdv::CoriolisAdv(RuntimeParams &params) {

  params.doc_module("MOM_CoriolisAdv", "");

  params.get("NOSLIP", no_slip_,
             {.default_value = false,
              .desc = "If true, no slip boundary conditions are used; otherwise free slip "
                      "boundary conditions are assumed. The implementation of the free slip "
                      "BCs on a C-grid is much cleaner than the no slip BCs. The use of the "
                      "free slip BCs is strongly encouraged, and no slip BCs are not used with "
                      "the biharmonic viscosity."});
  if (no_slip_) {
    // defer: the no-slip vorticity boundary condition.
    logger::fatal("CoriolisAdv: NOSLIP is not implemented yet.");
  }

  bool en_dis = false;
  params.get("CORIOLIS_EN_DIS", en_dis,
             {.default_value = false,
              .desc = "If true, two estimates of the thickness fluxes are used to estimate the "
                      "Coriolis term, and the one that dissipates energy relative to the other "
                      "one is used."});
  if (en_dis) {
    // defer: the energy-dissipating variant of the Sadourny energy scheme.
    logger::fatal("CoriolisAdv: CORIOLIS_EN_DIS is not implemented yet.");
  }

  std::string scheme = "SADOURNY75_ENERGY";
  params.get("CORIOLIS_SCHEME", scheme,
             {.default_value = std::string("SADOURNY75_ENERGY"),
              .desc = "CORIOLIS_SCHEME selects the discretization for the Coriolis terms. "
                      "Valid values are:\n"
                      "\t SADOURNY75_ENERGY - Sadourny, 1975; energy conserving\n"
                      "\t ARAKAWA_HSU90     - Arakawa & Hsu, 1990\n"
                      "\t SADOURNY75_ENSTRO - Sadourny, 1975; enstrophy conserving\n"
                      "\t ARAKAWA_LAMB81    - Arakawa & Lamb, 1981; both energy and enstrophy\n"
                      "\t ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with Arakawa & Hsu "
                      "and Sadourny energy"});
  if (scheme != "SADOURNY75_ENERGY") {
    // defer: the remaining vorticity schemes.
    logger::fatal("CoriolisAdv: CORIOLIS_SCHEME \"", scheme, "\" is not implemented yet.");
  }

  params.get("BOUND_CORIOLIS", bound_Coriolis_,
             {.default_value = false,
              .desc = "If true, the Coriolis terms at u-points are bounded by the four "
                      "estimates of (f+rv)v from the four neighboring v-points, and similarly "
                      "at v-points.  This option is always effectively false with "
                      "CORIOLIS_EN_DIS defined and CORIOLIS_SCHEME set to SADOURNY75_ENERGY."});

  std::string ke_scheme = "KE_ARAKAWA";
  params.get("KE_SCHEME", ke_scheme,
             {.default_value = std::string("KE_ARAKAWA"),
              .desc = "KE_SCHEME selects the discretization for acceleration due to the "
                      "kinetic energy gradient. Valid values are:\n"
                      "\t KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV, KE_UP3"});
  if (ke_scheme != "KE_ARAKAWA") {
    // defer: the Gudonov and third-order-upwind kinetic energy discretizations.
    logger::fatal("CoriolisAdv: KE_SCHEME \"", ke_scheme, "\" is not implemented yet.");
  }

  std::string pv_scheme = "PV_ADV_CENTERED";
  params.get("PV_ADV_SCHEME", pv_scheme,
             {.default_value = std::string("PV_ADV_CENTERED"),
              .desc = "PV_ADV_SCHEME selects the discretization for PV advection. Valid "
                      "values are:\n"
                      "\t PV_ADV_CENTERED - centered (aka Sadourny, 75)\n"
                      "\t PV_ADV_UPWIND1  - upwind, first order"});
  if (pv_scheme != "PV_ADV_CENTERED") {
    // defer: the first-order upwind PV advection.
    logger::fatal("CoriolisAdv: PV_ADV_SCHEME \"", pv_scheme, "\" is not implemented yet.");
  }
}

void CoriolisAdv::calculate(amrex::MultiFab &CAu, amrex::MultiFab &CAv,
                            const amrex::MultiFab &u, const amrex::MultiFab &v,
                            const amrex::MultiFab &h, const amrex::MultiFab &uh,
                            const amrex::MultiFab &vh, const Domain &domain,
                            const Grid &grid, const VerticalGrid &vgrid) const {

  const int nk = vgrid.nk();
  const bool bound_Coriolis = bound_Coriolis_;

  // The ocean area of an h cell [L2 ~> m2] and the sum of the four ocean areas
  // around a q point. Both are time-invariant, but MOM6 recomputes them each
  // call, so they stay here rather than moving onto the grid.
  amrex::MultiFab Area_h = domain.make_field(Stagger::Cell, 1, 1);
  amrex::MultiFab Area_q = domain.make_field(Stagger::Node, 1, 1);
  Area_h.setVal(0.0);
  Area_q.setVal(0.0);

  // Per-layer work fields at q points.
  amrex::MultiFab abs_vort = domain.make_field(Stagger::Node, 1, 1);
  amrex::MultiFab q = domain.make_field(Stagger::Node, 1, 1);
  amrex::MultiFab hArea_u = domain.make_field(Stagger::XFace, 1, 1);
  amrex::MultiFab hArea_v = domain.make_field(Stagger::YFace, 1, 1);
  amrex::MultiFab KE = domain.make_field(Stagger::Cell, 1, 1);
  abs_vort.setVal(0.0);
  q.setVal(0.0);
  hArea_u.setVal(0.0);
  hArea_v.setVal(0.0);
  KE.setVal(0.0);

  for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
    const amrex::Box valid = mfi.validbox();
    const amrex::Array4<amrex::Real> area_h = Area_h.array(mfi);
    const amrex::Array4<const amrex::Real> maskT = grid.mask2dT().const_array(mfi);
    const amrex::Array4<const amrex::Real> areaT = grid.areaT().const_array(mfi);
    amrex::ParallelFor(loops::flat(loops::h_points_grown(valid, 3)),
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      area_h(i, j, k) = maskT(i, j, k) * areaT(i, j, k);
    });

    const amrex::Array4<amrex::Real> area_q = Area_q.array(mfi);
    amrex::ParallelFor(loops::flat(loops::q_points(valid, 2)),
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      // MOM6's Area_q(I,J) sums the four h cells around the corner; AMReX
      // node (i,j) is MOM6's (I,J) = (i-1,j-1), so the four cells are
      // (i-1,j-1), (i,j), (i,j-1) and (i-1,j).
      area_q(i, j, k) = (area_h(i - 1, j - 1, k) + area_h(i, j, k)) +
                        (area_h(i, j - 1, k) + area_h(i - 1, j, k));
    });
  }

  for (int k = 0; k < nk; ++k) {
    for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
      const amrex::Box valid = mfi.validbox();

      const amrex::Array4<const amrex::Real> uu = u.const_array(mfi);
      const amrex::Array4<const amrex::Real> vv = v.const_array(mfi);
      const amrex::Array4<const amrex::Real> hh = h.const_array(mfi);
      const amrex::Array4<const amrex::Real> uuh = uh.const_array(mfi);
      const amrex::Array4<const amrex::Real> vvh = vh.const_array(mfi);
      const amrex::Array4<const amrex::Real> area_h = Area_h.const_array(mfi);
      const amrex::Array4<const amrex::Real> area_q = Area_q.const_array(mfi);
      const amrex::Array4<const amrex::Real> dyCv = grid.dyCv().const_array(mfi);
      const amrex::Array4<const amrex::Real> dxCu = grid.dxCu().const_array(mfi);
      const amrex::Array4<const amrex::Real> maskBu = grid.mask2dBu().const_array(mfi);
      const amrex::Array4<const amrex::Real> IareaBu = grid.IareaBu().const_array(mfi);
      const amrex::Array4<const amrex::Real> f = grid.CoriolisBu().const_array(mfi);

      const amrex::Array4<amrex::Real> hA_u = hArea_u.array(mfi);
      const amrex::Array4<amrex::Real> hA_v = hArea_v.array(mfi);
      const amrex::Array4<amrex::Real> av = abs_vort.array(mfi);
      const amrex::Array4<amrex::Real> qq = q.array(mfi);

      // The area-weighted thicknesses at the velocity points. MOM6's
      // hArea_u(I,j) is at AMReX u index I+1, between cells I and I+1.
      amrex::ParallelFor(loops::flat(loops::u_points(valid, 2)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        hA_u(i, j, 0) = 0.5 * ((area_h(i - 1, j, 0) * hh(i - 1, j, k)) +
                               (area_h(i, j, 0) * hh(i, j, k)));
      });
      amrex::ParallelFor(loops::flat(loops::v_points(valid, 2)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        hA_v(i, j, 0) = 0.5 * ((area_h(i, j - 1, 0) * hh(i, j - 1, k)) +
                               (area_h(i, j, 0) * hh(i, j, k)));
      });

      // The circulation around a q point, the relative and absolute
      // vorticities, and the potential vorticity q = (f + rv) / h.
      amrex::ParallelFor(loops::flat(loops::q_points(valid, 1)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        const amrex::Real dvdx = (vv(i, j, k) * dyCv(i, j, 0)) -
                                 (vv(i - 1, j, k) * dyCv(i - 1, j, 0));
        const amrex::Real dudy = (uu(i, j, k) * dxCu(i, j, 0)) -
                                 (uu(i, j - 1, k) * dxCu(i, j - 1, 0));
        const amrex::Real rel_vort = maskBu(i, j, 0) * (dvdx - dudy) * IareaBu(i, j, 0);
        av(i, j, 0) = f(i, j, 0) + rel_vort;

        const amrex::Real hArea_q = (hA_u(i, j - 1, 0) + hA_u(i, j, 0)) +
                                    (hA_v(i - 1, j, 0) + hA_v(i, j, 0));
        const amrex::Real Ih_q = area_q(i, j, 0) / (hArea_q + VOL_NEGLECT);
        qq(i, j, 0) = av(i, j, 0) * Ih_q;
      });

      // The Arakawa & Lamb kinetic energy, which carries the metric terms
      // that make the kinetic energy budget close.
      const amrex::Array4<amrex::Real> ke = KE.array(mfi);
      const amrex::Array4<const amrex::Real> areaCu = grid.areaCu().const_array(mfi);
      const amrex::Array4<const amrex::Real> areaCv = grid.areaCv().const_array(mfi);
      const amrex::Array4<const amrex::Real> IareaT = grid.IareaT().const_array(mfi);
      amrex::ParallelFor(loops::flat(loops::h_points_grown(valid, 1)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        ke(i, j, 0) = (((areaCu(i + 1, j, 0) * (uu(i + 1, j, k) * uu(i + 1, j, k))) +
                        (areaCu(i, j, 0) * (uu(i, j, k) * uu(i, j, k)))) +
                       ((areaCv(i, j + 1, 0) * (vv(i, j + 1, k) * vv(i, j + 1, k))) +
                        (areaCv(i, j, 0) * (vv(i, j, k) * vv(i, j, k))))) *
                      0.25 * IareaT(i, j, 0);
      });

      // CAu = q vh - d(KE)/dx, Sadourny 1975 energy-conserving form.
      const amrex::Array4<amrex::Real> cau = CAu.array(mfi);
      const amrex::Array4<const amrex::Real> IdxCu = grid.IdxCu().const_array(mfi);
      const amrex::Array4<const amrex::Real> IdxCu_OBC =
          grid.IdxCu_OBCmask().const_array(mfi);
      amrex::ParallelFor(loops::flat(loops::u_points(valid)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        amrex::Real ca = 0.25 *
            ((qq(i, j + 1, 0) * (vvh(i, j + 1, k) + vvh(i - 1, j + 1, k))) +
             (qq(i, j, 0) * (vvh(i - 1, j, k) + vvh(i, j, k)))) * IdxCu(i, j, 0);

        if (bound_Coriolis) {
          const amrex::Real fv1 = av(i, j + 1, 0) * vv(i, j + 1, k);
          const amrex::Real fv2 = av(i, j + 1, 0) * vv(i - 1, j + 1, k);
          const amrex::Real fv3 = av(i, j, 0) * vv(i, j, k);
          const amrex::Real fv4 = av(i, j, 0) * vv(i - 1, j, k);
          const amrex::Real max_fv = amrex::max(amrex::max(fv1, fv2), amrex::max(fv3, fv4));
          const amrex::Real min_fv = amrex::min(amrex::min(fv1, fv2), amrex::min(fv3, fv4));
          ca = amrex::min(ca, max_fv);
          ca = amrex::max(ca, min_fv);
        }

        const amrex::Real KEx = (ke(i, j, 0) - ke(i - 1, j, 0)) * IdxCu_OBC(i, j, 0);
        cau(i, j, k) = ca - KEx;
      });

      // CAv = -q uh - d(KE)/dy.
      const amrex::Array4<amrex::Real> cav = CAv.array(mfi);
      const amrex::Array4<const amrex::Real> IdyCv = grid.IdyCv().const_array(mfi);
      const amrex::Array4<const amrex::Real> IdyCv_OBC =
          grid.IdyCv_OBCmask().const_array(mfi);
      amrex::ParallelFor(loops::flat(loops::v_points(valid)),
                         [=] AMREX_GPU_DEVICE(int i, int j, int) {
        amrex::Real ca = -0.25 *
            ((qq(i, j, 0) * (uuh(i, j - 1, k) + uuh(i, j, k))) +
             (qq(i + 1, j, 0) * (uuh(i + 1, j - 1, k) + uuh(i + 1, j, k)))) * IdyCv(i, j, 0);

        if (bound_Coriolis) {
          const amrex::Real fu1 = -av(i + 1, j, 0) * uu(i + 1, j, k);
          const amrex::Real fu2 = -av(i + 1, j, 0) * uu(i + 1, j - 1, k);
          const amrex::Real fu3 = -av(i, j, 0) * uu(i, j, k);
          const amrex::Real fu4 = -av(i, j, 0) * uu(i, j - 1, k);
          const amrex::Real max_fu = amrex::max(amrex::max(fu1, fu2), amrex::max(fu3, fu4));
          const amrex::Real min_fu = amrex::min(amrex::min(fu1, fu2), amrex::min(fu3, fu4));
          ca = amrex::min(ca, max_fu);
          ca = amrex::max(ca, min_fu);
        }

        const amrex::Real KEy = (ke(i, j, 0) - ke(i, j - 1, 0)) * IdyCv_OBC(i, j, 0);
        cav(i, j, k) = ca - KEy;
      });
    }
  }

  CAu.FillBoundary(domain.periodicity());
  CAv.FillBoundary(domain.periodicity());
}

} // namespace MOM
