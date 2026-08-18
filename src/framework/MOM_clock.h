#pragma once
/// @file MOM_clock.h
/// @brief The run's time axis: the current time, the end time, and the
///        timestep hierarchy. The analogue of the time_type variables and the
///        run-control parameters that MOM6's solo driver holds
///        (MOM_driver.F90), rather than of a MOM6 module.

#include <cstdint>

#include <AMReX_REAL.H>

#include "MOM_file_parser.h"

namespace MOM {

/// @class Clock
/// @brief The run's time axis: how long the run is, how long a forcing
/// interval is, and how many dynamics steps fit in one.
///
/// Time is kept as an integer count of microseconds since the start of the
/// run, so that advancing by a timestep neither drifts nor rounds. MOM6 does
/// the same job with FMS's time_type, a (days, ticks) integer pair. When TIM's
/// C++ time types (TIM::TimeStamp, TIM::Calendar) land, they replace this
/// representation; the accessors below are what the rest of the model uses,
/// and they do not change.
class Clock {
public:
  /// @brief Read the run-control parameters (DT, DT_FORCING, DAYMAX,
  /// TIMEUNIT) and set the time axis. The analogue of the run-control
  /// parameter reads in MOM6's solo driver.
  /// @param params Runtime parameters.
  /// @throws logger::FatalError on a non-positive or inconsistent interval.
  explicit Clock(RuntimeParams &params);

  /// @brief The baroclinic dynamics timestep actually used [T ~> s]. It is
  /// dt_forcing divided by an integer, which may be shorter than DT.
  /// @return The dynamics timestep.
  amrex::Real dt() const { return dt_; }

  /// @brief The forcing (coupling) interval [T ~> s].
  /// @return The forcing interval.
  amrex::Real dt_forcing() const { return dt_forcing_; }

  /// @brief The number of dynamics steps in one forcing interval.
  /// @return The dynamics step count per forcing interval.
  int steps_per_forcing() const { return steps_per_forcing_; }

  /// @brief The time elapsed since the start of the run [T ~> s].
  /// @return The current time.
  amrex::Real time() const { return static_cast<amrex::Real>(now_us_) * 1.0e-6; }

  /// @brief The time at which the run ends [T ~> s].
  /// @return The end time.
  amrex::Real end_time() const { return static_cast<amrex::Real>(end_us_) * 1.0e-6; }

  /// @brief Whether the run has reached its end time.
  /// @return True once the current time is at or past the end time.
  bool done() const { return now_us_ >= end_us_; }

  /// @brief Advance the clock by one forcing interval.
  void advance();

  /// @brief The number of forcing intervals completed so far.
  /// @return The forcing-interval count.
  int forcing_step() const { return forcing_step_; }

private:
  amrex::Real dt_ = 0.0;              ///< The dynamics timestep [T ~> s].
  amrex::Real dt_forcing_ = 0.0;      ///< The forcing interval [T ~> s].
  int steps_per_forcing_ = 1;         ///< Dynamics steps per forcing interval.
  int forcing_step_ = 0;              ///< Forcing intervals completed.
  std::int64_t now_us_ = 0;           ///< The current time [microseconds].
  std::int64_t end_us_ = 0;           ///< The end time [microseconds].
  std::int64_t dt_forcing_us_ = 0;    ///< The forcing interval [microseconds].
};

} // namespace MOM
