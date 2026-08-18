#include "MOM.h"
#include "MOM_domains.h"
#include "MOM_fixed_initialization.h"
#include "MOM_state_initialization.h"
#include "MOM_logger.h"

namespace MOM {

Model::Model(RuntimeParams &params)
  : config_(read_config_switches(params)),
    domain_(make_domain(params)),
    grid_(initialize_fixed(domain_, params)),
    vgrid_(params),
    state_(initialize_state(domain_,
                            {.nk = vgrid_.nk(),
                             .max_depth = grid_.max_depth(),
                             .angstrom = vgrid_.angstrom()},
                            grid_.bathyT(), params)) {

  // Initialization phases, in the order of MOM6's initialize_MOM:
  initialize_dynamics(params);

  logger::note("MOM core initialization complete.");
}

Model::Config Model::read_config_switches(RuntimeParams &params) {

  logger::info("Initializing the MOM core...");

  params.doc_module("MOM", "Main MOM ocean model module"); // set current param module for documentation purposes

  int verbosity = 2;
  params.get("VERBOSITY", verbosity,
             {.default_value = 2,
              .desc = "Integer controlling level of messaging\n"
                      "\t0 = Only FATAL messages\n"
                      "\t2 = Only FATAL, WARNING, NOTE [default]\n"
                      "\t9 = All",
              .units = "",
              .fail_if_missing = false});

  logger::set_verbosity(verbosity);
  logger::info("Log verbosity: ", logger::get_verbosity());

  Config config;

  params.get("SPLIT", config.split, {.default_value = true, .desc = "Use the split time stepping if true."});

  params.get("SPLIT_RK4", config.split_rk4,
             {.default_value = false,
              .desc = "If true, use a version of the split explicit time stepping scheme that "
                      "exchanges velocities with step_MOM that have the average barotropic phase over "
                      "a baroclinic timestep rather than the instantaneous barotropic phase.",
              .do_not_log = !config.split});

  if (!config.split) {
    params.get("USE_RK2", config.use_RK2,
               {.default_value = false,
                .desc = "If true, use RK2 instead of RK3 in the unsplit time stepping."});
  }

  params.get("FPMIX", config.fpmix,
             {.default_value = false,
              .desc = "If true, use the FPMIX algorithm for tracer advection.",
              .do_not_log = true});

  if (config.fpmix && !config.split) {
    logger::fatal("FPMIX is only implemented for the split time stepping.");
  }

  params.get("DEBUG", config.debug,
             {.default_value = false,
              .desc = "If true, write out verbose debugging data.",
              .units = "nondim",
              .debugging_param = true});

  return config;
}

void Model::initialize_dynamics(RuntimeParams &params) {

  logger::note("initialize_dynamics: (stub)");

  // todo: dynamics subsystem construction, dispatching on config_.split /
  //       config_.use_RK2 as in MOM6's four-way initialize_dyn_* branch.
  // defer: restart registration (register_restarts_dyn_*), diagnostics.
  (void)params;
}

} // namespace MOM
