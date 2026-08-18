#include <string>

#include "MOM_logger.h"
#include "MOM_vertical_grid.h"

namespace MOM {

namespace {

// Set the layer target densities (Rlay) and the interface reduced gravities
// (g_prime) for COORD_CONFIG = "gprime". The analogue of MOM6's
// set_coord_from_gprime (MOM_coord_initialization.F90), Boussinesq branch.
void set_coord_from_gprime(RuntimeParams &params,
                           const amrex::Real g_Earth, const amrex::Real Rho0,
                           const bool Boussinesq,
                           std::vector<amrex::Real> &g_prime,
                           std::vector<amrex::Real> &Rlay) {

  const int nk = static_cast<int>(Rlay.size());

  // Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  amrex::Real g_fs = g_Earth;
  params.get("GFS", g_fs,
             {.default_value = g_Earth,
              .desc = "The reduced gravity at the free surface.",
              .units = "m s-2"});

  // Reduced gravity across the internal interfaces [L2 Z-1 T-2 ~> m s-2].
  amrex::Real g_int = 0.0;
  params.get("GINT", g_int,
             {.desc = "The reduced gravity across internal interfaces.",
              .units = "m s-2",
              .fail_if_missing = true});

  // Target density of the surface layer [R ~> kg m-3].
  amrex::Real Rlay_ref = Rho0;
  params.get("LIGHTEST_DENSITY", Rlay_ref,
             {.default_value = Rho0,
              .desc = "The reference potential density used for layer 1.",
              .units = "kg m-3"});

  g_prime[0] = g_fs;
  for (int k = 1; k < nk; ++k) {
    g_prime[k] = g_int;
  }

  Rlay[0] = Rlay_ref;
  if (Boussinesq) {
    for (int k = 1; k < nk; ++k) {
      Rlay[k] = Rlay[k - 1] + g_prime[k] * (Rho0 / g_Earth);
    }
  } else {
    // defer: the non-Boussinesq branch of this recursion,
    //          Rlay[k] = Rlay[k-1] * (g_Earth + 0.5*g_prime[k])
    //                              / (g_Earth - 0.5*g_prime[k]),
    //        taken only when BOUSSINESQ and SEMI_BOUSSINESQ are both false.
    logger::fatal("set_coord_from_gprime: the non-Boussinesq mode is not implemented.");
  }
}

} // namespace

VerticalGrid::VerticalGrid(RuntimeParams &params) {

  // The analogue of MOM6's verticalGridInit (MOM_verticalGrid.F90).
  params.doc_module("MOM_verticalGrid", "Parameters providing information about the vertical grid.");

  params.get("G_EARTH", g_Earth_,
             {.default_value = 9.80,
              .desc = "The gravitational acceleration of the Earth.",
              .units = "m s-2"});
  if (!(g_Earth_ > 0.0)) {
    logger::fatal("VerticalGrid: G_EARTH must be positive.");
  }

  params.get("RHO_0", Rho0_,
             {.default_value = 1035.0,
              .desc = "The mean ocean density used with BOUSSINESQ true to calculate accelerations "
                      "and the mass for conservation properties, or with BOUSSINESQ false to convert "
                      "some parameters from vertical units of m to kg m-2.",
              .units = "kg m-3"});
  if (!(Rho0_ > 0.0)) {
    logger::fatal("VerticalGrid: RHO_0 must be positive.");
  }

  params.get("BOUSSINESQ", Boussinesq_,
             {.default_value = true,
              .desc = "If true, make the Boussinesq approximation."});
  if (!Boussinesq_) {
    logger::fatal("VerticalGrid: the non-Boussinesq mode is not implemented.");
  }

  // defer: the remaining verticalGridInit content. SEMI_BOUSSINESQ and
  //        RHO_KV_CONVERT (both only meaningful when BOUSSINESQ is false,
  //        which aborts above), ANGSTROM and the subroundoff thicknesses,
  //        and H_RESCALE_POWER with the thickness-unit conversion factor
  //        family (H_to_m, Z_to_H, ...), all of which will be implemented
  //        with the thickness/units layer, and the mixed-layer layer counts
  //        (nkml, nk_rho_varies), which will be implemented with the bulk
  //        mixed layer.

  params.get("NK", nk_,
             {.desc = "The number of model layers.",
              .units = "nondim",
              .fail_if_missing = true});
  if (nk_ < 1) {
    logger::fatal("VerticalGrid: NK must be positive.");
  }

  // todo: read unlogged here; the logged, authoritative MAXIMUM_DEPTH read
  //       shoud happen in MOM_grid_init when it's introduced.
  params.get("MAXIMUM_DEPTH", max_depth_,
             {.desc = "The maximum depth of the ocean.",
              .units = "m",
              .fail_if_missing = true,
              .do_not_log = true});
  if (!(max_depth_ > 0.0)) {
    logger::fatal("VerticalGrid: MAXIMUM_DEPTH must be positive.");
  }

  // The analogue of MOM6's MOM_initialize_coord (MOM_coord_initialization.F90):
  // set up the layer densities, Rlay, and reduced gravities, g_prime.
  params.doc_module("MOM_coord_initialization", "");

  std::string coord_config = "none";
  params.get("COORD_CONFIG", coord_config,
             {.default_value = std::string("none"),
              .desc = "This specifies how layers are to be defined:\n"
                      "\t ALE or none - used to avoid defining layers in ALE mode\n"
                      "\t file - read coordinate information from the file\n"
                      "\t\t specified by (COORD_FILE).\n"
                      "\t BFB - Custom coords for buoyancy-forced basin case\n"
                      "\t\t based on SST_S, T_BOT and DRHO_DT.\n"
                      "\t linear - linear based on interfaces not layers\n"
                      "\t layer_ref - linear based on layer densities\n"
                      "\t ts_ref - use reference temperature and salinity\n"
                      "\t ts_range - use range of temperature and salinity\n"
                      "\t\t (T_REF and S_REF) to determine surface density\n"
                      "\t\t and GINT calculate internal densities.\n"
                      "\t gprime - use reference density (RHO_0) for surface\n"
                      "\t\t density and GINT calculate internal densities.\n"
                      "\t ts_profile - use temperature and salinity profiles\n"
                      "\t\t (read from COORD_FILE) to set layer densities.\n"
                      "\t USER - call a user modified routine."});

  g_prime_.assign(nk_ + 1, 0.0);
  Rlay_.assign(nk_, 0.0);

  if (coord_config == "gprime") {
    set_coord_from_gprime(params, g_Earth_, Rho0_, Boussinesq_, g_prime_, Rlay_);
  } else {
    // defer: the remaining COORD_CONFIG options (none/ALE, file, BFB, linear,
    //        layer_ref, ts_ref, ts_range, ts_profile, USER).
    logger::fatal("VerticalGrid: Unsupported COORD_CONFIG \"", coord_config,
                  "\". Only \"gprime\" is currently implemented.");
  }

  // There are nk+1 g_prime values because it is an interface field, but the
  // value at the bottom should not matter. This is here (as in MOM6) just to
  // avoid having an uninitialized value in some output.
  g_prime_[nk_] = 10.0 * g_Earth_;

  // defer: the DEBUG checksum output of Rlay and g_prime (MOM6's chksum
  //        calls in MOM_initialize_coord); the checksum oracle arrives with
  //        TIM::checksum at the State PR.
  // defer: setVerticalGridAxes (the sLayer/sInterface diagnostic coordinate
  //        axes and their names/units) until diagnostics are ready.
  // defer: the unit-string helpers (get_thickness_units, get_flux_units,
  //        get_tr_flux_units) until diagnostics and I/O registration.
}

} // namespace MOM
