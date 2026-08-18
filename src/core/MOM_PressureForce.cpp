#include <string>
#include <vector>

#include <AMReX_Gpu.H>

#include "MOM_PressureForce.h"

#include "MOM_logger.h"
#include "MOM_loop_boxes.h"

namespace MOM {

PressureForce::PressureForce(RuntimeParams &params) {

  params.doc_module("MOM_PressureForce", "");

  std::string form = "FV";
  params.get("PRESSUREFORCE", form,
             {.default_value = std::string("FV"),
              .desc = "Pressure force scheme.  Options are:\n"
                      "\t FV - finite volume, with analytic integration\n"
                      "\t Montgomery - Montgomery potential form"});

  if (form != "Montgomery") {
    // defer: the finite-volume pressure force (MOM_PressureForce_FV), which
    //        is MOM6's default and needs the density integrals.
    logger::fatal("PressureForce: PRESSUREFORCE \"", form,
                  "\" is not implemented yet; only \"Montgomery\" is.");
  }
}

void PressureForce::calculate(amrex::MultiFab &PFu, amrex::MultiFab &PFv,
                              const amrex::MultiFab &h, const Domain &domain,
                              const Grid &grid, const VerticalGrid &vgrid) const {

  const int nk = vgrid.nk();

  // g_prime lives on the host as a std::vector; copy it to a device-visible
  // buffer once per call so the kernel can index it.
  amrex::Gpu::DeviceVector<amrex::Real> g_prime_d(vgrid.g_prime().size());
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, vgrid.g_prime().begin(),
                   vgrid.g_prime().end(), g_prime_d.begin());
  const amrex::Real *g_prime = g_prime_d.data();

  // The Montgomery potential, M = p/rho + g z [L2 T-2 ~> m2 s-2], on the same
  // decomposition as h, and the interface heights e it is built from
  // [Z ~> m]. MOM6 keeps both as local arrays of the same shapes.
  amrex::MultiFab M(h.boxArray(), h.DistributionMap(), 1, h.nGrowVect());
  M.setVal(0.0);
  amrex::MultiFab e(domain.box_array(nk + 1), h.DistributionMap(), 1, h.nGrowVect());
  e.setVal(0.0);

  for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
    const amrex::Box valid = mfi.validbox();
    // MOM6 builds e and M over Isq..Ieq+1, Jsq..Jeq+1: one halo ring, which is
    // what the gradients below need.
    const amrex::Box cells = loops::h_points_grown(valid, 1);
    const amrex::Box columns(amrex::IntVect(cells.smallEnd(0), cells.smallEnd(1), 0),
                             amrex::IntVect(cells.bigEnd(0), cells.bigEnd(1), 0));

    const amrex::Array4<const amrex::Real> hh = h.const_array(mfi);
    const amrex::Array4<const amrex::Real> D = grid.bathyT().const_array(mfi);
    const amrex::Array4<amrex::Real> MM = M.array(mfi);
    const amrex::Array4<amrex::Real> ee = e.array(mfi);

    // Both recursions run in one thread per column: the interface heights
    // build upward from the bottom, and the Montgomery potential downward
    // from the surface. H_to_Z is 1 in Boussinesq mode without unit scaling.
    amrex::ParallelFor(columns, [=] AMREX_GPU_DEVICE(int i, int j, int) {
      ee(i, j, nk) = -D(i, j, 0);
      for (int k = nk - 1; k >= 0; --k) {
        ee(i, j, k) = ee(i, j, k + 1) + hh(i, j, k);
      }
      MM(i, j, 0) = g_prime[0] * ee(i, j, 0);
      for (int k = 1; k < nk; ++k) {
        MM(i, j, k) = MM(i, j, k - 1) + g_prime[k] * ee(i, j, k);
      }
    });
  }

  for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
    const amrex::Box valid = mfi.validbox();
    const amrex::Array4<const amrex::Real> MM = M.const_array(mfi);

    const amrex::Array4<amrex::Real> pfu = PFu.array(mfi);
    const amrex::Array4<const amrex::Real> IdxCu = grid.IdxCu().const_array(mfi);
    amrex::Box u_bx = loops::u_points(valid);
    u_bx.setSmall(2, 0);
    u_bx.setBig(2, nk - 1);
    amrex::ParallelFor(u_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      pfu(i, j, k) = -(MM(i, j, k) - MM(i - 1, j, k)) * IdxCu(i, j, 0);
    });

    const amrex::Array4<amrex::Real> pfv = PFv.array(mfi);
    const amrex::Array4<const amrex::Real> IdyCv = grid.IdyCv().const_array(mfi);
    amrex::Box v_bx = loops::v_points(valid);
    v_bx.setSmall(2, 0);
    v_bx.setBig(2, nk - 1);
    amrex::ParallelFor(v_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      pfv(i, j, k) = -(MM(i, j, k) - MM(i, j - 1, k)) * IdyCv(i, j, 0);
    });
  }

  PFu.FillBoundary(domain.periodicity());
  PFv.FillBoundary(domain.periodicity());
}

} // namespace MOM
