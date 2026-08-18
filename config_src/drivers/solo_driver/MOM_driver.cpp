/**
 * @file MOM_driver.cpp
 * @brief Driver program for protoMOMxx
 */

#include <cstdlib>
#include <exception>
#include <iostream>

#include "MOM.h"
#include <AMReX_Utility.H>

#include "MOM_clock.h"
#include "MOM_debug_dump.h"
#include "MOM_directories.h"
#include "MOM_infra.h"
#include "MOM_logger.h"
#include "MOM_surface_forcing.h"

/// @brief Main entry point for the protoMOMxx driver program.
/// @param argc Number of arguments including binary name.
/// @param argv Parameter command line as array pointers.
/// @return Exit code (0 for success, non-zero for failure)
int main(int argc, char* argv[]) {
  try {

    MOM::logger::info("Hello C++ world. This is protoMOMxx!");

    // Initialize the infrastructure layer via MOM::Infra, a wrapper
    // around TIM::Runtime. This initializes MPI and AMReX, which are
    // finalized automatically when this scope exits.
    const MOM::Infra infra(argc, argv);

    // todo: ensemble manager (from the infrastructure layer.)
    const int ensemble_num = -1;

    // Read input.nml to determine input/output directories and param file names
    const MOM::Directories directories(ensemble_num);

    // RuntimeParams reads the parameter files specified in input.nml.
    MOM::RuntimeParams params(directories.parameter_filenames(), "MOM_parameters_doc");

    // defer: set_calendar_type(); the run starts at time zero and carries no
    //        date, so nothing needs a calendar yet.
    // defer: time_interp_external_init()

    // Initialize the core MOM object (the analogue of MOM6's initialize_MOM).
    MOM::Model model(params);

    // The run's time axis. Constructed after the Model, as in MOM6's solo
    // driver, which reads DT_FORCING and DAYMAX after initialize_MOM.
    MOM::Clock clock(params);

    // defer: extract_surface_state() -- nothing consumes the surface state
    //        while the forcing is analytic and adiabatic.
    // The analytic surface forcing of the solo driver, and the fields it
    // writes into each forcing interval.
    const MOM::SurfaceForcing surface_forcing(model.grid(), params);
    MOM::MechForcing forces(model.domain());

    // defer: MOM_wave_interface_init(), data_override_init(), ice shelf hooks

    // Parity scaffolding: with DEBUG set, dump the fixed fields so they can be
    // diffed against legacy MOM6's ocean_geometry.nc. Retires with the I/O
    // layer.
    if (model.config().debug) {
      MOM::debug::set_reporting(true);
      MOM::debug::report_field(model.grid().areaT(), "areaT");
      MOM::debug::report_field(model.grid().dxT(), "dxT");
      MOM::debug::report_field(model.grid().CoriolisBu(), "CoriolisBu");
      MOM::debug::report_field(model.grid().bathyT(), "bathyT");
      MOM::debug::report_field(model.state().h(), "h init");
      MOM::debug::dump_field(model.grid().geoLatT(), "geoLatT");
      MOM::debug::dump_field(model.grid().geoLonT(), "geoLonT");
      MOM::debug::dump_field(model.grid().dxT(), "dxT");
      MOM::debug::dump_field(model.grid().dyT(), "dyT");
      MOM::debug::dump_field(model.grid().areaT(), "areaT");
      MOM::debug::dump_field(model.grid().dxCu(), "dxCu");
      MOM::debug::dump_field(model.grid().dyCu(), "dyCu");
      MOM::debug::dump_field(model.grid().dxCv(), "dxCv");
      MOM::debug::dump_field(model.grid().dyCv(), "dyCv");
      MOM::debug::dump_field(model.grid().dxBu(), "dxBu");
      MOM::debug::dump_field(model.grid().dyBu(), "dyBu");
      MOM::debug::dump_field(model.grid().CoriolisBu(), "CoriolisBu");
      MOM::debug::dump_field(model.grid().bathyT(), "bathyT");
      MOM::debug::dump_field(model.grid().mask2dT(), "mask2dT");
      MOM::debug::dump_field(model.state().h(), "h_init");
    }

    MOM::logger::note("Starting the time loop: ", clock.end_time() / 86400.0,
                      " days, ", clock.steps_per_forcing(),
                      " dynamics steps per forcing interval.");

    const double loop_start = amrex::second();
    int dyn_steps = 0;

    while (!clock.done()) {
      surface_forcing.set_forcing(forces, clock.time());
      model.step(forces, clock.dt_forcing(), clock.steps_per_forcing());
      dyn_steps += clock.steps_per_forcing();
      clock.advance();
    }

    const double loop_seconds = amrex::second() - loop_start;
    MOM::logger::note("Main loop: ", loop_seconds, " s for ", dyn_steps,
                      " dynamics steps (", loop_seconds / dyn_steps, " s/step).");

    // todo: finish_MOM_initialization() on the first iteration
    // defer: mech_forcing_diags(), forcing_diagnostics()
    // defer: save_MOM_restart(), write_ocean_solo_res(), diag_mediator_end()
    // todo: MOM_end()

    return 0;

  } catch (const MOM::logger::FatalError&) {
    // Already logged by logger::fatal.
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::cerr << "protoMOMxx terminated with an unhandled exception: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}
