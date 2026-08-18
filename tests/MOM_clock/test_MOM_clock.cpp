// Unit tests for Clock (src/framework/MOM_clock.cpp): the timestep hierarchy
// and the run length read from the runtime parameters. Clock touches no AMReX
// state, so no infra bring-up is needed.

#include <filesystem>
#include <gtest/gtest.h>

#include "MOM_clock.h"
#include "MOM_logger.h"

namespace {

MOM::RuntimeParams params_from(const std::string &rel_path) {
  return MOM::RuntimeParams((std::filesystem::path(__FILE__).parent_path() / rel_path).string());
}

} // namespace

// The driving testcase: DT = 1200, DT_FORCING = 2400, DAYMAX = 10 days, so the
// run is 360 forcing intervals of two 1200 s dynamics steps each -- the 720
// steps legacy MOM6 reports for this configuration.
TEST(Clock, DoubleGyreTimestepHierarchy) {
  auto params = params_from("../double_gyre/MOM_input");
  MOM::Clock clock(params);

  EXPECT_DOUBLE_EQ(clock.dt_forcing(), 2400.0);
  EXPECT_DOUBLE_EQ(clock.dt(), 1200.0);
  EXPECT_EQ(clock.steps_per_forcing(), 2);
  EXPECT_DOUBLE_EQ(clock.end_time(), 10.0 * 86400.0);
  EXPECT_DOUBLE_EQ(clock.time(), 0.0);
  EXPECT_FALSE(clock.done());

  int intervals = 0;
  while (!clock.done()) {
    clock.advance();
    ++intervals;
  }
  EXPECT_EQ(intervals, 360);
  EXPECT_EQ(clock.forcing_step(), 360);
  // Integer microseconds, so the end time is hit exactly, not overshot.
  EXPECT_DOUBLE_EQ(clock.time(), clock.end_time());
}

// A forcing interval that is not an integer multiple of DT shortens the
// dynamics step rather than overshooting, as MOM6 does.
TEST(Clock, DynamicsStepDividesTheForcingInterval) {
  auto params = params_from("MOM_param_files/MOM_input_uneven");
  MOM::Clock clock(params);

  // DT = 1000, DT_FORCING = 2400: three steps of 800 s.
  EXPECT_EQ(clock.steps_per_forcing(), 3);
  EXPECT_DOUBLE_EQ(clock.dt(), 800.0);
}

TEST(Clock, NonPositiveIntervalIsFatal) {
  auto params = params_from("MOM_param_files/MOM_input_bad_dt");
  EXPECT_THROW({ MOM::Clock clock(params); }, MOM::logger::FatalError);
}
