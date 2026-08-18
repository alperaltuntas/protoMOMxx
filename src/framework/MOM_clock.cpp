#include <cmath>

#include "MOM_clock.h"
#include "MOM_logger.h"

namespace MOM {

namespace {

// Convert a duration in seconds to an integer microsecond count, rounding to
// the nearest microsecond so that a value that is exact in seconds stays exact.
std::int64_t to_microseconds(const amrex::Real seconds) {
  return static_cast<std::int64_t>(std::llround(seconds * 1.0e6));
}

} // namespace

Clock::Clock(RuntimeParams &params) {

  params.doc_module("MOM_main (MOM_driver)", "");

  // DT is documented where MOM6 documents it, in the MOM module read by the
  // Model constructor; the driver re-reads it here without logging, as MOM6's
  // solo driver does.
  amrex::Real dt = 0.0;
  params.get("DT", dt,
             {.units = "s", .fail_if_missing = true, .do_not_log = true});
  if (!(dt > 0.0)) {
    logger::fatal("Clock: DT must be positive.");
  }

  dt_forcing_ = dt;
  params.get("DT_FORCING", dt_forcing_,
             {.default_value = dt,
              .desc = "The time step for changing forcing, coupling with other components, or "
                      "potentially writing certain diagnostics. The default value is given by DT.",
              .units = "s"});
  if (!(dt_forcing_ > 0.0)) {
    logger::fatal("Clock: DT_FORCING must be positive.");
  }

  // MOM6 shortens the dynamics step so that an integer number of them fills a
  // forcing interval exactly; the 0.001 slack absorbs the rounding of a ratio
  // that is meant to be integral.
  steps_per_forcing_ = 1;
  if (dt_forcing_ > dt) {
    steps_per_forcing_ = static_cast<int>(std::ceil(dt_forcing_ / dt - 0.001));
  }
  dt_ = dt_forcing_ / static_cast<amrex::Real>(steps_per_forcing_);

  // MOM6 documents TIMEUNIT under MOM_sum_output, which does not exist here
  // yet; it is documented under MOM_main until that module arrives.
  amrex::Real time_unit = 86400.0;
  params.get("TIMEUNIT", time_unit,
             {.default_value = 86400.0,
              .desc = "The time unit for DAYMAX, ENERGYSAVEDAYS, and RESTINT.",
              .units = "s"});
  if (!(time_unit > 0.0)) {
    logger::fatal("Clock: TIMEUNIT must be positive.");
  }

  amrex::Real daymax = 0.0;
  params.get("DAYMAX", daymax,
             {.desc = "The final time of the whole simulation, in units of TIMEUNIT seconds.  "
                      "This also sets the potential end time of the present run segment if the "
                      "end time is not set via ocean_solo_nml in input.nml.",
              .units = "days",
              .fail_if_missing = true});
  if (!(daymax > 0.0)) {
    logger::fatal("Clock: DAYMAX must be positive.");
  }

  dt_forcing_us_ = to_microseconds(dt_forcing_);
  end_us_ = to_microseconds(daymax * time_unit);
  now_us_ = 0;

  // defer: a start time other than zero (restarts), the calendar and date
  //        stamping, and the separate thermodynamic timestep DT_THERM --
  //        ADIABATIC is true in the driving testcase, so nothing consumes it.
}

void Clock::advance() {
  now_us_ += dt_forcing_us_;
  ++forcing_step_;
}

} // namespace MOM
