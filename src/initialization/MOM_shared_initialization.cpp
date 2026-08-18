#include <cmath>

#include "MOM_kernel_inline.h"
#include "MOM_logger.h"
#include "MOM_shared_initialization.h"

namespace MOM {

amrex::MultiFab planetary_rotation(const Domain &domain, const GridSpec &spec,
                                   const amrex::MultiFab &geoLatBu) {

  const int n_levels = 1;
  const int ncomp = 1;
  amrex::MultiFab CoriolisBu = domain.make_field(Stagger::Node, n_levels, ncomp);

  const amrex::Real PI = 4.0 * std::atan(1.0);
  const amrex::Real omega = spec.omega;

  for (amrex::MFIter mfi(CoriolisBu); mfi.isValid(); ++mfi) {
    const amrex::Box bx = mfi.growntilebox();
    const amrex::Array4<amrex::Real> f = CoriolisBu.array(mfi);
    const amrex::Array4<const amrex::Real> lat = geoLatBu.const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) MOM_KERNEL_INLINE {
      f(i, j, k) = (2.0 * omega) * std::sin((PI * lat(i, j, k)) / 180.0);
    });
  }

  return CoriolisBu;
}

amrex::MultiFab named_topography(const Domain &domain, const std::string &config,
                                 const GridSpec &grid_spec, const TopoSpec &topo_spec,
                                 const amrex::MultiFab &geoLonT,
                                 const amrex::MultiFab &geoLatT) {

  const int n_levels = 1;
  const int ncomp = 1;
  amrex::MultiFab D = domain.make_field(Stagger::Cell, n_levels, ncomp);
  // MOM6 allocates bathyT with source=0.0 and the named shapes only write the
  // computational domain, so the halos start at zero here too.
  D.setVal(0.0);

  const amrex::Real PI = 4.0 * std::atan(1.0);
  const amrex::Real max_depth = topo_spec.max_depth;
  const amrex::Real min_depth = topo_spec.min_depth;
  const amrex::Real Dedge = topo_spec.edge_depth;
  const amrex::Real expdecay = topo_spec.expdecay;
  const amrex::Real len_lat = grid_spec.len_lat;
  const amrex::Real len_lon = grid_spec.len_lon;
  const amrex::Real south_lat = grid_spec.south_lat;
  const amrex::Real west_lon = grid_spec.west_lon;
  const amrex::Real rad_earth = grid_spec.rad_earth;

  if (config == "flat") {
    for (amrex::MFIter mfi(D); mfi.isValid(); ++mfi) {
      const amrex::Box bx = mfi.tilebox();
      const amrex::Array4<amrex::Real> d = D.array(mfi);
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) MOM_KERNEL_INLINE {
        d(i, j, k) = max_depth;
      });
    }
  } else if (config == "spoon") {
    // A bowl-like basin with a vertical wall at the southern face. The
    // expression is transcribed from MOM6's initialize_topography_named,
    // including the squared denominator written as a product.
    const amrex::Real decay = 1.0 - std::exp(-0.5 * len_lat * rad_earth * PI / (180.0 * expdecay));
    const amrex::Real D0 = (max_depth - Dedge) / (decay * decay);
    for (amrex::MFIter mfi(D); mfi.isValid(); ++mfi) {
      const amrex::Box bx = mfi.tilebox();
      const amrex::Array4<amrex::Real> d = D.array(mfi);
      const amrex::Array4<const amrex::Real> lon = geoLonT.const_array(mfi);
      const amrex::Array4<const amrex::Real> lat = geoLatT.const_array(mfi);
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) MOM_KERNEL_INLINE {
        d(i, j, k) = Dedge + D0 *
            (std::sin(PI * (lon(i, j, k) - west_lon) / len_lon) *
             (1.0 - std::exp((lat(i, j, k) - (south_lat + len_lat)) * rad_earth * PI /
                             (180.0 * expdecay))));
      });
    }
  } else {
    // defer: the remaining named topographies (bowl, halfpipe) and the
    //        file-based and testcase-specific configurations.
    logger::fatal("named_topography: TOPO_CONFIG \"", config, "\" is not implemented yet.");
  }

  // limit_topography's default (MASKING_DEPTH unset) path, run over the grown
  // boxes as MOM6 runs it over the data domain.
  for (amrex::MFIter mfi(D); mfi.isValid(); ++mfi) {
    const amrex::Box bx = mfi.growntilebox();
    const amrex::Array4<amrex::Real> d = D.array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) MOM_KERNEL_INLINE {
      d(i, j, k) = amrex::min(amrex::max(d(i, j, k), 0.5 * min_depth), max_depth);
    });
  }

  // MOM6's pass_var(G%bathyT) before the masks are set. Halo points outside
  // the global domain are not touched and keep the clamped value above.
  D.FillBoundary(domain.periodicity());

  return D;
}

} // namespace MOM
