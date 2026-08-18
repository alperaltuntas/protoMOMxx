#include "MOM.h"
#include "MOM_domains.h"
#include "MOM_fixed_initialization.h"
#include "MOM_logger.h"

namespace MOM {

Model::Model(RuntimeParams &params)
  : config_(read_config_switches(params)),
    domain_(make_domain(params)),
    grid_(initialize_fixed(domain_, params)),
    vgrid_(params) {

  // Initialization phases, in the order of MOM6's initialize_MOM:
  initialize_state(params);
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

void Model::initialize_state(RuntimeParams &params) {

  logger::note("initialize_state: (stub)");

  // todo: State container (u, v, h) with proper staggering, filled from
  //       THICKNESS_CONFIG / VELOCITY_CONFIG. Analogue of MOM6's
  //       MOM_initialize_state.
  (void)params;

  // tmp: the original psi (stream function) demo, kept operational until the
  // real state initialization. Until then, this demo exercises the AMReX
  // machinery on the decomposition the Domain owns.
  const int n_components = 1;
  amrex::MultiFab psi(domain_.box_array(vgrid_.nk()), domain_.distribution_mapping(),
                      n_components, domain_.nghost());

  fill_psi_demo(psi);

  psi.FillBoundary(domain_.periodicity());
}

void Model::initialize_dynamics(RuntimeParams &params) {

  logger::note("initialize_dynamics: (stub)");

  // todo: dynamics subsystem construction, dispatching on config_.split /
  //       config_.use_RK2 as in MOM6's four-way initialize_dyn_* branch.
  // defer: restart registration (register_restarts_dyn_*), diagnostics.
  (void)params;
}

void Model::fill_psi_demo(amrex::MultiFab &psi) const
{
    // tmp: nominal cell sizes and physical extents. The demo deliberately
    // keeps its uniform Cartesian sizes (rather than adopting the Grid
    // metrics) until it retires with the real State initialization.
    const amrex::Real dx = 100000;
    const amrex::Real dy = 100000;

    const amrex::Real x_min = 0.0;
    const amrex::Real x_max = domain_.ni_global() * dx;
    const amrex::Real y_min = 0.0;
    const amrex::Real y_max = domain_.nj_global() * dy;

    //////////////////////////////////////////////////////////////////////////
    // Initialization of stream function (psi)
    //////////////////////////////////////////////////////////////////////////

    // coefficient for initialization psi
    const amrex::Real a = 1000000;
    const double pi = 4. * std::atan(1.);

    for (amrex::MFIter mfi(psi); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const amrex::Array4<amrex::Real>& phi_array = psi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const amrex::Real x_cell_center = (i+0.5) * dx;
            const amrex::Real y_cell_center = (j+0.5) * dy;

            const amrex::Real x_transformed = LinearMapCoordinates(x_cell_center, x_min, x_max, 0.0, 2*pi);
            const amrex::Real y_transformed = LinearMapCoordinates(y_cell_center, y_min, y_max, 0.0, 2*pi);

            phi_array(i,j,k) = a*std::sin(x_transformed)*std::sin(y_transformed);
        });
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real LinearMapCoordinates(const amrex::Real x,
                                 const amrex::Real x_min, const amrex::Real x_max,
                                 const amrex::Real xi_min, const amrex::Real xi_max)
{
    return x_min + ((xi_max-xi_min)/(x_max-x_min))*x;
}

} // namespace MOM
