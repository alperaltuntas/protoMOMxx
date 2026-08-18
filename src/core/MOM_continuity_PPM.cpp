#include <cmath>

#include "MOM_continuity_PPM.h"

#include "MOM_logger.h"
#include "MOM_loop_boxes.h"

namespace MOM {

namespace {

// The piecewise parabolic edge values of a cell, limited to stay positive.
// The analogue of MOM6's PPM_limit_pos, applied in place.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ppm_limit_pos(amrex::Real &h_L, amrex::Real &h_R, const amrex::Real h_in,
                   const amrex::Real h_min) {
  const amrex::Real curv = 3.0 * ((h_L + h_R) - 2.0 * h_in);
  if (curv > 0.0) {  // Only minima are limited.
    const amrex::Real dh = h_R - h_L;
    if (std::abs(dh) < curv) {  // The parabola's minimum is inside the cell.
      if (h_in <= h_min) {
        h_L = h_in;
        h_R = h_in;
      } else {
        // This file is compiled without floating-point contraction, because
        // gfortran leaves MOM6's continuity solver unfused everywhere except
        // these two expressions, which are written out as fused
        // multiply-adds.
        const amrex::Real bound = std::fma(curv, curv, 3.0 * (dh * dh));
        if (12.0 * curv * (h_in - h_min) < bound) {
          const amrex::Real scale = 12.0 * curv * (h_in - h_min) / bound;
          h_L = std::fma(scale, h_L - h_in, h_in);
          h_R = std::fma(scale, h_R - h_in, h_in);
        }
      }
    }
  }
}

// The monotonized second-order slope of MOM6's PPM_reconstruction_*, Eq. B2 of
// Lin 1994. Returns zero where any of the three cells is land.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real ppm_slope(const amrex::Real h_m1, const amrex::Real h_0,
                      const amrex::Real h_p1, const amrex::Real mask_m1,
                      const amrex::Real mask_0, const amrex::Real mask_p1) {
  if ((mask_m1 * mask_0 * mask_p1) == 0.0) {
    return 0.0;
  }
  amrex::Real slp = 0.5 * (h_p1 - h_m1);
  const amrex::Real dMx = amrex::max(amrex::max(h_p1, h_m1), h_0) - h_0;
  const amrex::Real dMn = h_0 - amrex::min(amrex::min(h_p1, h_m1), h_0);
  const amrex::Real sgn = (slp >= 0.0) ? 1.0 : -1.0;
  return sgn * amrex::min(std::abs(slp), 2.0 * amrex::min(dMx, dMn));
}

// The transport through one face and the marginal face thickness. The analogue
// of MOM6's flux_elem with vol_CFL false and no porous barriers.
// h/h_L/h_R are the cell and its edge values on the low side of the face,
// and h_p1/h_L_p1/h_R_p1 those on the high side.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real flux_elem(const amrex::Real u, const amrex::Real h, const amrex::Real h_p1,
                      const amrex::Real h_L, const amrex::Real h_L_p1,
                      const amrex::Real h_R, const amrex::Real h_R_p1,
                      const amrex::Real face_length, const amrex::Real Idx,
                      const amrex::Real Idx_p1, const amrex::Real dt) {
  if (u > 0.0) {
    const amrex::Real CFL = u * dt * Idx;
    const amrex::Real curv_3 = (h_L + h_R) - 2.0 * h;
    const amrex::Real dh = h_L - h_R;
    return face_length * u * (h_R + CFL * (0.5 * dh + curv_3 * (CFL - 1.5)));
  }
  if (u < 0.0) {
    const amrex::Real CFL = -u * dt * Idx_p1;
    const amrex::Real curv_3 = (h_L_p1 + h_R_p1) - 2.0 * h_p1;
    const amrex::Real dh = h_R_p1 - h_L_p1;
    return face_length * u * (h_L_p1 + CFL * (0.5 * dh + curv_3 * (CFL - 1.5)));
  }
  return 0.0;
}

} // namespace

ContinuityPPM::ContinuityPPM(RuntimeParams &params) {

  params.doc_module("MOM_continuity_PPM", "");

  bool monotonic = false;
  params.get("MONOTONIC_CONTINUITY", monotonic,
             {.default_value = false,
              .desc = "If true, CONTINUITY_PPM uses the Colella and Woodward monotonic limiter. "
                      "The default (false) is to use a simple positive definite limiter."});
  if (monotonic) {
    // defer: the Colella and Woodward monotonic limiter (PPM_limit_CW84).
    logger::fatal("ContinuityPPM: MONOTONIC_CONTINUITY is not implemented yet.");
  }

  bool simple_2nd = false;
  params.get("SIMPLE_2ND_PPM_CONTINUITY", simple_2nd,
             {.default_value = false,
              .desc = "If true, CONTINUITY_PPM uses a simple 2nd order (arithmetic mean) "
                      "interpolation of the edge values. This may give better PV conservation "
                      "properties. While it formally reduces the accuracy of the continuity "
                      "solver itself in the strongly advective limit, it does not reduce the "
                      "overall order of accuracy of the dynamic core."});
  if (simple_2nd) {
    // defer: the arithmetic-mean edge reconstruction.
    logger::fatal("ContinuityPPM: SIMPLE_2ND_PPM_CONTINUITY is not implemented yet.");
  }

  bool upwind_1st = false;
  params.get("UPWIND_1ST_CONTINUITY", upwind_1st,
             {.default_value = false,
              .desc = "If true, CONTINUITY_PPM becomes a 1st-order upwind continuity solver. "
                      "This scheme is highly diffusive but may be useful for debugging or "
                      "in single-column mode where its minimal stencil is useful."});
  if (upwind_1st) {
    // defer: the first-order upwind reconstruction.
    logger::fatal("ContinuityPPM: UPWIND_1ST_CONTINUITY is not implemented yet.");
  }

  bool vol_CFL = false;
  params.get("CONT_PPM_VOLUME_BASED_CFL", vol_CFL,
             {.default_value = false,
              .desc = "If true, use the ratio of the open face lengths to the tracer cell "
                      "areas when estimating CFL numbers in the continuity solver."});
  if (vol_CFL) {
    // defer: the volume-based CFL estimate.
    logger::fatal("ContinuityPPM: CONT_PPM_VOLUME_BASED_CFL is not implemented yet.");
  }

  // The remaining CONT_PPM_* parameters only affect the barotropic-correction
  // iteration, which the unsplit scheme never enters; they are read by the
  // split scheme when it lands.
}

void ContinuityPPM::solve(amrex::MultiFab &h, amrex::MultiFab &uh, amrex::MultiFab &vh,
                          const amrex::MultiFab &u, const amrex::MultiFab &v,
                          const amrex::MultiFab &hin, const amrex::Real dt,
                          const Domain &domain, const Grid &grid,
                          const VerticalGrid &vgrid) const {

  const int nk = vgrid.nk();
  const int sten = stencil();
  const amrex::Real h_min = vgrid.angstrom();
  // MOM6 limits the PPM edge values against twice the Angstrom.
  const amrex::Real h_min_ppm = 2.0 * vgrid.angstrom();

  // The reconstruction edge values and slopes, on h's decomposition.
  amrex::MultiFab h_L(h.boxArray(), h.DistributionMap(), 1, h.nGrowVect());
  amrex::MultiFab h_R(h.boxArray(), h.DistributionMap(), 1, h.nGrowVect());
  amrex::MultiFab slp(h.boxArray(), h.DistributionMap(), 1, h.nGrowVect());
  h_L.setVal(0.0);
  h_R.setVal(0.0);
  slp.setVal(0.0);

  // --- Zonal pass, over the range the meridional pass will need. ---
  for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
    const amrex::Box valid = mfi.validbox();
    // MOM6's LB with the j stencil added: i over the computational domain,
    // j widened by the scheme's stencil.
    const amrex::Box lb = amrex::grow(valid, amrex::IntVect(0, sten, 0));

    const amrex::Array4<const amrex::Real> hi = hin.const_array(mfi);
    const amrex::Array4<const amrex::Real> maskT = grid.mask2dT().const_array(mfi);
    const amrex::Array4<amrex::Real> hL = h_L.array(mfi);
    const amrex::Array4<amrex::Real> hR = h_R.array(mfi);
    const amrex::Array4<amrex::Real> sl = slp.array(mfi);

    // PPM_reconstruction_x: slopes over isl-1..iel+1, edges over isl..iel,
    // where isl = ish-1 and iel = ieh+1.
    amrex::ParallelFor(amrex::grow(lb, amrex::IntVect(2, 0, 0)),
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      sl(i, j, k) = ppm_slope(hi(i - 1, j, k), hi(i, j, k), hi(i + 1, j, k),
                              maskT(i - 1, j, 0), maskT(i, j, 0), maskT(i + 1, j, 0));
    });
    amrex::ParallelFor(amrex::grow(lb, amrex::IntVect(1, 0, 0)),
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      // A land neighbour contributes the local thickness instead of its own.
      const amrex::Real h_im1 = maskT(i - 1, j, 0) * hi(i - 1, j, k) +
                                (1.0 - maskT(i - 1, j, 0)) * hi(i, j, k);
      const amrex::Real h_ip1 = maskT(i + 1, j, 0) * hi(i + 1, j, k) +
                                (1.0 - maskT(i + 1, j, 0)) * hi(i, j, k);
      constexpr amrex::Real oneSixth = 1.0 / 6.0;
      amrex::Real L = 0.5 * (h_im1 + hi(i, j, k)) + oneSixth * (sl(i - 1, j, k) - sl(i, j, k));
      amrex::Real R = 0.5 * (h_ip1 + hi(i, j, k)) + oneSixth * (sl(i, j, k) - sl(i + 1, j, k));
      ppm_limit_pos(L, R, hi(i, j, k), h_min_ppm);
      hL(i, j, k) = L;
      hR(i, j, k) = R;
    });

    // zonal_mass_flux: MOM6's I = ish-1..ieh, which is the u-face box of lb.
    const amrex::Array4<const amrex::Real> uu = u.const_array(mfi);
    const amrex::Array4<amrex::Real> uuh = uh.array(mfi);
    const amrex::Array4<const amrex::Real> dy_Cu = grid.dy_Cu().const_array(mfi);
    const amrex::Array4<const amrex::Real> IdxT = grid.IdxT().const_array(mfi);
    amrex::Box u_bx = loops::u_points(lb);
    u_bx.setSmall(2, 0);
    u_bx.setBig(2, nk - 1);
    amrex::ParallelFor(u_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      uuh(i, j, k) = flux_elem(uu(i, j, k), hi(i - 1, j, k), hi(i, j, k),
                               hL(i - 1, j, k), hL(i, j, k), hR(i - 1, j, k), hR(i, j, k),
                               dy_Cu(i, j, 0), IdxT(i - 1, j, 0), IdxT(i, j, 0), dt);
    });

    // continuity_zonal_convergence. MOM6 floors the updated thickness at one
    // Angstrom in both passes, not at zero.
    const amrex::Array4<amrex::Real> hh = h.array(mfi);
    const amrex::Array4<const amrex::Real> IareaT = grid.IareaT().const_array(mfi);
    amrex::Box cell_bx = lb;
    cell_bx.setSmall(2, 0);
    cell_bx.setBig(2, nk - 1);
    amrex::ParallelFor(cell_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      hh(i, j, k) = amrex::max(std::fma(-(dt * IareaT(i, j, 0)),
                                       uuh(i + 1, j, k) - uuh(i, j, k), hi(i, j, k)),
                               h_min);
    });
  }

  // --- Meridional pass, on the zonally updated thicknesses. ---
  for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
    const amrex::Box valid = mfi.validbox();
    const amrex::Box lb = valid;

    const amrex::Array4<const amrex::Real> hh_in = h.const_array(mfi);
    const amrex::Array4<const amrex::Real> maskT = grid.mask2dT().const_array(mfi);
    const amrex::Array4<amrex::Real> hL = h_L.array(mfi);
    const amrex::Array4<amrex::Real> hR = h_R.array(mfi);
    const amrex::Array4<amrex::Real> sl = slp.array(mfi);

    amrex::ParallelFor(amrex::grow(lb, amrex::IntVect(0, 2, 0)),
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      sl(i, j, k) = ppm_slope(hh_in(i, j - 1, k), hh_in(i, j, k), hh_in(i, j + 1, k),
                              maskT(i, j - 1, 0), maskT(i, j, 0), maskT(i, j + 1, 0));
    });
    amrex::ParallelFor(amrex::grow(lb, amrex::IntVect(0, 1, 0)),
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const amrex::Real h_jm1 = maskT(i, j - 1, 0) * hh_in(i, j - 1, k) +
                                (1.0 - maskT(i, j - 1, 0)) * hh_in(i, j, k);
      const amrex::Real h_jp1 = maskT(i, j + 1, 0) * hh_in(i, j + 1, k) +
                                (1.0 - maskT(i, j + 1, 0)) * hh_in(i, j, k);
      constexpr amrex::Real oneSixth = 1.0 / 6.0;
      amrex::Real S = 0.5 * (h_jm1 + hh_in(i, j, k)) + oneSixth * (sl(i, j - 1, k) - sl(i, j, k));
      amrex::Real N = 0.5 * (h_jp1 + hh_in(i, j, k)) + oneSixth * (sl(i, j, k) - sl(i, j + 1, k));
      ppm_limit_pos(S, N, hh_in(i, j, k), h_min_ppm);
      hL(i, j, k) = S;
      hR(i, j, k) = N;
    });

    const amrex::Array4<const amrex::Real> vv = v.const_array(mfi);
    const amrex::Array4<amrex::Real> vvh = vh.array(mfi);
    const amrex::Array4<const amrex::Real> dx_Cv = grid.dx_Cv().const_array(mfi);
    const amrex::Array4<const amrex::Real> IdyT = grid.IdyT().const_array(mfi);
    amrex::Box v_bx = loops::v_points(lb);
    v_bx.setSmall(2, 0);
    v_bx.setBig(2, nk - 1);
    amrex::ParallelFor(v_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      vvh(i, j, k) = flux_elem(vv(i, j, k), hh_in(i, j - 1, k), hh_in(i, j, k),
                               hL(i, j - 1, k), hL(i, j, k), hR(i, j - 1, k), hR(i, j, k),
                               dx_Cv(i, j, 0), IdyT(i, j - 1, 0), IdyT(i, j, 0), dt);
    });

    // continuity_merdional_convergence, on the zonally updated thicknesses.
    const amrex::Array4<amrex::Real> hh = h.array(mfi);
    const amrex::Array4<const amrex::Real> IareaT = grid.IareaT().const_array(mfi);
    amrex::Box cell_bx = lb;
    cell_bx.setSmall(2, 0);
    cell_bx.setBig(2, nk - 1);
    amrex::ParallelFor(cell_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      hh(i, j, k) = amrex::max(std::fma(-(dt * IareaT(i, j, 0)),
                                       vvh(i, j + 1, k) - vvh(i, j, k), hh(i, j, k)),
                               h_min);
    });
  }

  (void)domain;
}

} // namespace MOM
