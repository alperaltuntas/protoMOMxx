#include <cmath>

#include "MOM_logger.h"
#include "MOM_surface_forcing.h"

namespace MOM {

SurfaceForcing::SurfaceForcing(const Grid &grid, RuntimeParams &params)
  : grid_(grid) {

  params.doc_module("MOM_surface_forcing", "");

  params.get("VARIABLE_WINDS", variable_winds_,
             {.default_value = true,
              .desc = "If true, the winds vary in time after the initialization."});

  wind_config_ = "zero";
  params.get("WIND_CONFIG", wind_config_,
             {.default_value = std::string("zero"),
              .desc = "The character string that indicates how wind forcing is specified.  "
                      "Valid options include (file), (data_override), (2gyre), (1gyre), "
                      "(gyres), (zero), (const), (Neverworld), (scurves), (ideal_hurr), "
                      "(SCM_CVmix_tests) and (USER)."});

  if (wind_config_ == "2gyre" || wind_config_ == "1gyre") {
    params.get("TAUX_MAGNITUDE", taux_mag_,
               {.default_value = 0.1,
                .desc = "The magnitude of the wind stress.",
                .units = "Pa"});
  } else if (wind_config_ != "zero") {
    // defer: the remaining wind configurations (file, data_override, gyres,
    //        const, Neverworld, scurves, ideal_hurr, SCM_CVmix_tests, USER).
    logger::fatal("SurfaceForcing: WIND_CONFIG \"", wind_config_,
                  "\" is not implemented yet.");
  }

  if (variable_winds_ && wind_config_ != "zero") {
    // defer: time-varying winds. The analytic gyre profiles below are steady,
    //        so a run that asks for varying winds would silently get steady
    //        ones.
    logger::fatal("SurfaceForcing: VARIABLE_WINDS is not implemented yet; "
                  "set VARIABLE_WINDS = False.");
  }

  // defer: the buoyancy forcing (BUOY_CONFIG), the gustiness parameters, and
  //        the ustar derivation -- the driving testcase runs adiabatically
  //        with no buoyancy forcing.
}

void SurfaceForcing::set_forcing(MechForcing &forces, const amrex::Real time) const {

  (void)time;  // The implemented profiles are steady.

  if (wind_config_ == "zero") {
    forces.taux().setVal(0.0);
    forces.tauy().setVal(0.0);
    return;
  }

  const amrex::Real PI = 4.0 * std::atan(1.0);
  const amrex::Real taux_mag = taux_mag_;
  const amrex::Real south_lat = grid_.south_lat();
  const amrex::Real len_lat = grid_.len_lat();
  const bool two_gyre = (wind_config_ == "2gyre");

  for (amrex::MFIter mfi(forces.taux()); mfi.isValid(); ++mfi) {
    const amrex::Box bx = mfi.tilebox();
    const amrex::Array4<amrex::Real> taux = forces.taux().array(mfi);
    const amrex::Array4<const amrex::Real> lat = grid_.geoLatCu().const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      if (two_gyre) {
        taux(i, j, k) = taux_mag *
            (1.0 - std::cos(2.0 * PI * (lat(i, j, k) - south_lat) / len_lat));
      } else {
        taux(i, j, k) = taux_mag *
            std::cos(PI * (lat(i, j, k) - south_lat) / len_lat);
      }
    });
  }

  forces.tauy().setVal(0.0);
}

} // namespace MOM
