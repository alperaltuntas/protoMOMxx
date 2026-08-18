#include "MOM_dynamics_unsplit_RK2.h"

#include "MOM_debug_dump.h"
#include "MOM_fields.h"
#include "MOM_logger.h"
#include "MOM_loop_boxes.h"

namespace MOM {

namespace {

// h_av = (h + hp) / 2, over the range the Coriolis stencil reads.
void average_thickness(amrex::MultiFab &h_av, const amrex::MultiFab &h,
                       const amrex::MultiFab &hp, const int nk, const int halo) {
  for (amrex::MFIter mfi(h_av); mfi.isValid(); ++mfi) {
    amrex::Box bx = loops::h_points_grown(mfi.validbox(), halo);
    bx.setSmall(2, 0);
    bx.setBig(2, nk - 1);
    const amrex::Array4<amrex::Real> a = h_av.array(mfi);
    const amrex::Array4<const amrex::Real> hh = h.const_array(mfi);
    const amrex::Array4<const amrex::Real> hhp = hp.const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      a(i, j, k) = (hh(i, j, k) + hhp(i, j, k)) * 0.5;
    });
  }
}

// vel_out = mask * (vel_in + weight * dt * ((PF + CA) + diff)), the velocity
// update MOM6 writes out twice in the corrector.
void accelerate(amrex::MultiFab &out, const amrex::MultiFab &in,
                const amrex::MultiFab &PF, const amrex::MultiFab &CA,
                const amrex::MultiFab &diff, const amrex::MultiFab &mask,
                const amrex::Real dt_wt, const int nk, const bool x_face) {
  for (amrex::MFIter mfi(out); mfi.isValid(); ++mfi) {
    amrex::Box bx = x_face ? loops::u_points(mfi.validbox())
                           : loops::v_points(mfi.validbox());
    bx.setSmall(2, 0);
    bx.setBig(2, nk - 1);
    const amrex::Array4<amrex::Real> o = out.array(mfi);
    const amrex::Array4<const amrex::Real> in_a = in.const_array(mfi);
    const amrex::Array4<const amrex::Real> pf = PF.const_array(mfi);
    const amrex::Array4<const amrex::Real> ca = CA.const_array(mfi);
    const amrex::Array4<const amrex::Real> df = diff.const_array(mfi);
    const amrex::Array4<const amrex::Real> m = mask.const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      o(i, j, k) = m(i, j, 0) * (in_a(i, j, k) + dt_wt *
                                 ((pf(i, j, k) + ca(i, j, k)) + df(i, j, k)));
    });
  }
}

} // namespace

DynamicsUnsplitRK2::DynamicsUnsplitRK2(RuntimeParams &params, const amrex::Real dt,
                                       const Domain &domain, const Grid &grid,
                                       const VerticalGrid &vgrid)
  : pressure_force_(params),
    coriolis_adv_(params),
    continuity_(params),
    hor_visc_(params, dt, domain, grid),
    set_visc_(params),
    vert_friction_(params, domain, vgrid),
    visc_(domain),
    CAu_(make_layer_field(domain, vgrid, Stagger::XFace)),
    CAv_(make_layer_field(domain, vgrid, Stagger::YFace)),
    PFu_(make_layer_field(domain, vgrid, Stagger::XFace)),
    PFv_(make_layer_field(domain, vgrid, Stagger::YFace)),
    diffu_(make_layer_field(domain, vgrid, Stagger::XFace)),
    diffv_(make_layer_field(domain, vgrid, Stagger::YFace)),
    up_(make_layer_field(domain, vgrid, Stagger::XFace)),
    vp_(make_layer_field(domain, vgrid, Stagger::YFace)),
    hp_(make_layer_field(domain, vgrid, Stagger::Cell)),
    h_av_(make_layer_field(domain, vgrid, Stagger::Cell)),
    uh_(make_layer_field(domain, vgrid, Stagger::XFace)),
    vh_(make_layer_field(domain, vgrid, Stagger::YFace)) {

  params.doc_module("MOM_dynamics_unsplit_RK2", "");

  params.get("BE", BE_,
             {.default_value = 0.6,
              .desc = "If SPLIT is true, BE determines the relative weighting of a 2nd-order "
                      "Runga-Kutta baroclinic time stepping scheme (0.5) and a backward Euler "
                      "scheme (1) that is used for the Coriolis and inertial terms.  BE may be "
                      "from 0.5 to 1, but instability may occur near 0.5.  BE is also "
                      "applicable if SPLIT is false and USE_RK2 is true.",
              .units = "nondim"});

  params.get("BEGW", begw_,
             {.default_value = 0.0,
              .desc = "If SPLIT is true, BEGW is a number from 0 to 1 that controls the extent "
                      "to which the treatment of gravity waves is forward-backward (0) or "
                      "simulated backward Euler (1).  0 is almost always used. If SPLIT is "
                      "false and USE_RK2 is true, BEGW can be between 0 and 0.5 to damp "
                      "gravity waves.",
              .units = "nondim"});

  bool dt_visc_bug = false;
  params.get("UNSPLIT_DT_VISC_BUG", dt_visc_bug,
             {.default_value = false,
              .desc = "If false, use the correct timestep in the viscous terms applied in the "
                      "first predictor step with the unsplit time stepping scheme, and in the "
                      "calculation of the turbulent mixed layer properties for viscosity with "
                      "unsplit or unsplit_RK2.  If true, an older incorrect value is used."});
  if (dt_visc_bug) {
    // defer: reproducing the older incorrect viscous timestep.
    logger::fatal("DynamicsUnsplitRK2: UNSPLIT_DT_VISC_BUG is not implemented.");
  }

  CAu_.setVal(0.0); CAv_.setVal(0.0);
  PFu_.setVal(0.0); PFv_.setVal(0.0);
  diffu_.setVal(0.0); diffv_.setVal(0.0);
  up_.setVal(0.0); vp_.setVal(0.0);
  hp_.setVal(0.0); h_av_.setVal(0.0);
  uh_.setVal(0.0); vh_.setVal(0.0);
}

void DynamicsUnsplitRK2::step(State &state, const MechForcing &forces, const amrex::Real dt,
                              const Domain &domain, const Grid &grid,
                              const VerticalGrid &vgrid) {

  const int nk = vgrid.nk();
  const amrex::Real dt_pred = dt * BE_;
  // MOM6's CoriolisAdv_stencil, 2 for every scheme implemented here.
  const int cor_stencil = 2;

  amrex::MultiFab &u = state.u();
  amrex::MultiFab &v = state.v();
  amrex::MultiFab &h = state.h();

  debug::report_field(u, "Start Predictor u");
  debug::report_field(v, "Start Predictor v");
  debug::report_field(h, "Start Predictor h");

  // The bottom boundary layer properties, from the state at the start of the
  // step. MOM6 calls this from step_MOM, before the time stepping scheme; it
  // sits here because the viscosity fields are owned by this class.
  set_visc_.set_viscous_BBL(visc_, u, v, h, domain, grid, vgrid);

  // diffu, diffv: the horizontal viscosity at the start of the step.
  hor_visc_.calculate(diffu_, diffv_, u, v, h, domain, grid, vgrid);

  // A continuity step whose only purpose is the Coriolis terms: it gives the
  // thickness fluxes and the half-step thicknesses the potential vorticity
  // denominator needs. MOM6 notes this duplicates the last continuity call of
  // the previous step and could be optimized out.
  continuity_.solve(hp_, uh_, vh_, u, v, h, dt_pred, domain, grid, vgrid);
  hp_.FillBoundary(domain.periodicity());
  uh_.FillBoundary(domain.periodicity());
  vh_.FillBoundary(domain.periodicity());

  average_thickness(h_av_, h, hp_, nk, cor_stencil);

  coriolis_adv_.calculate(CAu_, CAv_, u, v, h_av_, uh_, vh_, domain, grid, vgrid);
  pressure_force_.calculate(PFu_, PFv_, h, domain, grid, vgrid);

  debug::report_field(CAu_, "Predictor 1 accel CAu");
  debug::report_field(CAv_, "Predictor 1 accel CAv");
  debug::report_field(PFu_, "Predictor 1 accel PFu");
  debug::report_field(PFv_, "Predictor 1 accel PFv");
  debug::report_field(diffu_, "Predictor 1 accel diffu");
  debug::report_field(diffv_, "Predictor 1 accel diffv");
  // The predictor velocities, up = u + BE*dt*(PF + CA + diff).
  accelerate(up_, u, PFu_, CAu_, diffu_, grid.mask2dCu(), dt_pred, nk, true);
  accelerate(vp_, v, PFv_, CAv_, diffv_, grid.mask2dCv(), dt_pred, nk, false);

  vert_friction_.coefficients(up_, vp_, h_av_, visc_, set_visc_.Kv(), domain, grid, vgrid);
  vert_friction_.apply(up_, vp_, h_av_, forces, dt_pred, domain, grid, vgrid);

  // The half-step continuity, giving the fluxes the corrector's Coriolis
  // terms use and the thicknesses centred on the half step.
  continuity_.solve(hp_, uh_, vh_, up_, vp_, h, dt, domain, grid, vgrid);
  hp_.FillBoundary(domain.periodicity());
  uh_.FillBoundary(domain.periodicity());
  vh_.FillBoundary(domain.periodicity());

  average_thickness(h_av_, h, hp_, nk, cor_stencil);

  debug::report_field(up_, "Predictor 1 u");
  debug::report_field(vp_, "Predictor 1 v");
  debug::report_field(h_av_, "Predictor 1 h");
  debug::report_field(uh_, "Predictor 1 uh");
  debug::report_field(vh_, "Predictor 1 vh");

  coriolis_adv_.calculate(CAu_, CAv_, up_, vp_, h_av_, uh_, vh_, domain, grid, vgrid);

  debug::report_field(CAu_, "Corrector accel CAu");
  debug::report_field(CAv_, "Corrector accel CAv");

  // The corrector. up is extrapolated by (1+BEGW) to damp gravity waves and
  // is what the final continuity step transports with; u takes the plain
  // full step.
  accelerate(up_, u, PFu_, CAu_, diffu_, grid.mask2dCu(), dt * (1.0 + begw_), nk, true);
  accelerate(vp_, v, PFv_, CAv_, diffv_, grid.mask2dCv(), dt * (1.0 + begw_), nk, false);
  accelerate(u, u, PFu_, CAu_, diffu_, grid.mask2dCu(), dt, nk, true);
  accelerate(v, v, PFv_, CAv_, diffv_, grid.mask2dCv(), dt, nk, false);

  vert_friction_.coefficients(up_, vp_, h_av_, visc_, set_visc_.Kv(), domain, grid, vgrid);
  vert_friction_.apply(up_, vp_, h_av_, forces, dt, domain, grid, vgrid);
  vert_friction_.coefficients(u, v, h_av_, visc_, set_visc_.Kv(), domain, grid, vgrid);
  vert_friction_.apply(u, v, h_av_, forces, dt, domain, grid, vgrid);

  // The final continuity step advances the thicknesses with the extrapolated
  // velocities.
  continuity_.solve(h, uh_, vh_, up_, vp_, h, dt, domain, grid, vgrid);
  h.FillBoundary(domain.periodicity());
  uh_.FillBoundary(domain.periodicity());
  vh_.FillBoundary(domain.periodicity());

  debug::report_field(u, "Corrector u");
  debug::report_field(v, "Corrector v");
  debug::report_field(h, "Corrector h");
  debug::report_field(uh_, "Corrector uh");
  debug::report_field(vh_, "Corrector vh");
  debug::report_field(up_, "Corrector up");
}

} // namespace MOM
