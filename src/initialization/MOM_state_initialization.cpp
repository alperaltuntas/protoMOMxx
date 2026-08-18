#include <string>

#include "MOM_logger.h"
#include "MOM_state_initialization.h"

namespace MOM {

namespace {

// Fill h with layers evenly distributed between the surface and MAXIMUM_DEPTH,
// squeezed into the local water column where it is shallower. The analogue of
// MOM6's initialize_thickness_uniform.
//
// MOM6 builds the column downward from the resting interface heights and the
// total depth, giving each layer at least an Angstrom. The recursion is
// sequential in k, so this launches one thread per column and loops over k
// inside it rather than one thread per cell.
void initialize_thickness_uniform(amrex::MultiFab &h, const Domain &domain,
                                  const StateSpec &spec,
                                  const amrex::MultiFab &bathyT) {

  const int nk = spec.nk;
  const amrex::Real max_depth = spec.max_depth;
  const amrex::Real angstrom = spec.angstrom;

  if (!(max_depth > 0.0)) {
    logger::fatal("initialize_thickness_uniform: MAXIMUM_DEPTH is not set.");
  }

  for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
    // MOM6 fills the computational domain only and leaves the halos to the
    // exchange below.
    const amrex::Box bx = mfi.tilebox();
    const amrex::Array4<amrex::Real> hh = h.array(mfi);
    const amrex::Array4<const amrex::Real> D = bathyT.const_array(mfi);

    const amrex::Box columns(amrex::IntVect(bx.smallEnd(0), bx.smallEnd(1), 0),
                             amrex::IntVect(bx.bigEnd(0), bx.bigEnd(1), 0));
    amrex::ParallelFor(columns, [=] AMREX_GPU_DEVICE(int i, int j, int) {
      // depth_tot = bathyT + Z_ref, with Z_ref at its default of zero.
      amrex::Real eta_below = -D(i, j, 0);
      for (int k = nk - 1; k >= 0; --k) {
        // The resting height of interface k, positive upward.
        const amrex::Real e0 = -max_depth * (amrex::Real(k) / amrex::Real(nk));
        amrex::Real eta = e0;
        if (eta < (eta_below + angstrom)) {
          eta = eta_below + angstrom;
          hh(i, j, k) = angstrom;
        } else {
          hh(i, j, k) = eta - eta_below;
        }
        eta_below = eta;
      }
    });
  }

  h.FillBoundary(domain.periodicity());
}

} // namespace

StateFields initialize_state(const Domain &domain, const StateSpec &spec,
                             const amrex::MultiFab &bathyT, RuntimeParams &params) {

  params.doc_module("MOM_state_initialization", "");

  StateFields fields;
  fields.h = domain.make_field(Stagger::Cell, spec.nk, 1);
  fields.u = domain.make_field(Stagger::XFace, spec.nk, 1);
  fields.v = domain.make_field(Stagger::YFace, spec.nk, 1);

  // MOM6 allocates the prognostic arrays with source=0.0; the halos outside
  // the global domain, which no initialization reaches, stay zero here too.
  fields.h.setVal(0.0);
  fields.u.setVal(0.0);
  fields.v.setVal(0.0);

  std::string h_config = "uniform";
  params.get("THICKNESS_CONFIG", h_config,
             {.default_value = std::string("uniform"),
              .desc = "A string that determines how the initial layer thicknesses are "
                      "specified for a new run:\n"
                      "\t file - read interface heights from the file specified\n"
                      "\t\t by (THICKNESS_FILE).\n"
                      "\t thickness_file - read thicknesses from the file specified\n"
                      "\t\t by (THICKNESS_FILE).\n"
                      "\t uniform - uniform thickness layers evenly distributed\n"
                      "\t\t between the surface and MAXIMUM_DEPTH.\n"
                      "\t list - read a list of positive interface depths.\n"
                      "\t param - use thicknesses from parameter THICKNESS_INIT_VALUES.\n"
                      "\t USER - call a user modified routine."});

  if (h_config == "uniform") {
    initialize_thickness_uniform(fields.h, domain, spec, bathyT);
  } else {
    // defer: the remaining THICKNESS_CONFIG options.
    logger::fatal("initialize_state: THICKNESS_CONFIG \"", h_config,
                  "\" is not implemented yet.");
  }

  std::string u_config = "zero";
  params.get("VELOCITY_CONFIG", u_config,
             {.default_value = std::string("zero"),
              .desc = "A string that determines how the initial velocities are specified "
                      "for a new run:\n"
                      "\t file - read velocities from the file specified\n"
                      "\t\t by (VELOCITY_FILE).\n"
                      "\t zero - the fluid is initially at rest.\n"
                      "\t uniform - the flow is uniform (determined by\n"
                      "\t\t parameters INITIAL_U_CONST and INITIAL_V_CONST).\n"
                      "\t USER - call a user modified routine."});

  if (u_config != "zero") {
    // defer: the remaining VELOCITY_CONFIG options. "zero" needs no work
    //        beyond the zeroing above, as in MOM6's initialize_velocity_zero.
    logger::fatal("initialize_state: VELOCITY_CONFIG \"", u_config,
                  "\" is not implemented yet.");
  }

  // defer: the surrounding MOM_initialize_state content -- restart reading,
  //        temperature and salinity (ENABLE_THERMODYNAMICS is false here),
  //        sponges, OBCs, ODA increments, and the initial ALE remap.

  return fields;
}

} // namespace MOM
