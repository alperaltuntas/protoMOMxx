#include "MOM_sum_output.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include "MOM_fields.h"
#include "MOM_fp_contract.h"
#include "MOM_logger.h"
#include "MOM_loop_boxes.h"

namespace MOM {

namespace {

/// A directory name with exactly one trailing slash. MOM6's slasher.
std::string slasher(const std::string &dir) {
  if (dir.empty()) {
    return "./";
  }
  return (dir.back() == '/') ? dir : dir + "/";
}

/// The index list that sorts depths into ascending order. This is MOM6's
/// indexed heapsort, transcribed rather than replaced by std::sort, because
/// the order of equal depths decides the order the areas are accumulated in
/// and so shows up in the last digits of the volumes.
std::vector<int> sort_index_ascending(const std::vector<double> &d, const int n) {
  std::vector<int> indx(static_cast<std::size_t>(n) + 2);
  for (int j = 1; j <= n + 1; ++j) {
    indx[static_cast<std::size_t>(j)] = j;
  }
  if (n < 2) {
    return indx;
  }

  int k = n / 2 + 1;
  int ir = n;
  while (true) {
    int indxt = 0;
    double dnow = 0.0;
    if (k > 1) {
      --k;
      indxt = indx[static_cast<std::size_t>(k)];
      dnow = d[static_cast<std::size_t>(indxt)];
    } else {
      indxt = indx[static_cast<std::size_t>(ir)];
      dnow = d[static_cast<std::size_t>(indxt)];
      indx[static_cast<std::size_t>(ir)] = indx[1];
      --ir;
      if (ir == 1) {
        indx[1] = indxt;
        break;
      }
    }
    int i = k;
    int j = k * 2;
    while (j <= ir) {
      if (j < ir && d[static_cast<std::size_t>(indx[static_cast<std::size_t>(j)])] <
                        d[static_cast<std::size_t>(indx[static_cast<std::size_t>(j) + 1])]) {
        ++j;
      }
      if (dnow < d[static_cast<std::size_t>(indx[static_cast<std::size_t>(j)])]) {
        indx[static_cast<std::size_t>(i)] = indx[static_cast<std::size_t>(j)];
        i = j;
        j = j + i;
      } else {
        j = ir + 1;
      }
    }
    indx[static_cast<std::size_t>(i)] = indxt;
  }
  return indx;
}

} // namespace

SumOutput::SumOutput(RuntimeParams &params, const Domain &domain, const Grid &grid,
                     const VerticalGrid &vgrid, const std::string &directory) {

  params.doc_module("MOM_sum_output", "");

  params.get("CALCULATE_APE", do_APE_calc_,
             {.default_value = true,
              .desc = "If true, calculate the available potential energy of the interfaces. "
                      "Setting this to false reduces the memory footprint of high-PE-count "
                      "models dramatically."});

  bool write_stocks = true;
  params.get("WRITE_STOCKS", write_stocks,
             {.default_value = true,
              .desc = "If true, write the integrated tracer amounts to stdout when the energy "
                      "files are written."});

  params.get("DT", dt_,
             {.desc = "The (baroclinic) dynamics time step.",
              .units = "s",
              .fail_if_missing = true,
              .do_not_log = true});

  params.get("MAXTRUNC", maxtrunc_,
             {.default_value = 0,
              .desc = "The run will be stopped, and the day set to a very large value if the "
                      "velocity is truncated more than MAXTRUNC times between energy saves. "
                      "Set MAXTRUNC to 0 to stop if there is any truncation of velocities.",
              .units = "truncations save_interval-1"});

  params.get("MAX_ENERGY", max_energy_,
             {.default_value = 0.0,
              .desc = "The maximum permitted average energy per unit mass; the model will be "
                      "stopped if there is more energy than this.  If zero or negative, this "
                      "is set to 10*MAXVEL^2.",
              .units = "m2 s-2"});
  if (max_energy_ <= 0.0) {
    amrex::Real maxvel = 3.0e8;
    params.get("MAXVEL", maxvel,
               {.default_value = 3.0e8,
                .desc = "The maximum velocity allowed before the velocity components are "
                        "truncated.",
                .units = "m s-1"});
    max_energy_ = 10.0 * maxvel * maxvel;
  }

  std::string energyfile = "ocean.stats";
  params.get("ENERGYFILE", energyfile,
             {.default_value = std::string("ocean.stats"),
              .desc = "The file to use to write the energies and globally summed diagnostics."});
  energyfile_ = slasher(directory) + energyfile;

  params.get("TIMEUNIT", timeunit_,
             {.default_value = 86400.0,
              .desc = "The time unit in seconds a number of input fields",
              .units = "s"});
  if (timeunit_ < 0.0) {
    timeunit_ = 86400.0;
  }

  if (do_APE_calc_) {
    // defer: READ_DEPTH_LIST. The list is cheap to build and there is no I/O
    //        layer to read it back through; the values are the same either way.
    amrex::Real d_list_min_inc = 1.0e-10;
    params.get("DEPTH_LIST_MIN_INC", d_list_min_inc,
               {.default_value = 1.0e-10,
                .desc = "The minimum increment between the depths of the entries in the "
                        "depth-list file.",
                .units = "m"});
    lH_.assign(static_cast<std::size_t>(vgrid.nk()) + 1, 0);
    create_depth_list(domain, grid, d_list_min_inc);
    for (int k = 1; k <= vgrid.nk(); ++k) {
      lH_[static_cast<std::size_t>(k)] = listsize_ - 1;
    }
  } else {
    listsize_ = 1;
  }

  amrex::Real energysavedays = 1.0;
  params.get("ENERGYSAVEDAYS", energysavedays,
             {.default_value = 1.0,
              .desc = "The interval in units of TIMEUNIT between saves of the energies of the "
                      "run and other globally summed diagnostics.",
              .units = "days"});
  energysavedays_ = energysavedays * timeunit_;

  amrex::Real energysavedays_geometric = 0.0;
  params.get("ENERGYSAVEDAYS_GEOMETRIC", energysavedays_geometric,
             {.default_value = 0.0,
              .desc = "The starting interval in units of TIMEUNIT for the first call to save "
                      "the energies of the run and other globally summed diagnostics. The "
                      "interval increases by a factor of 2. after each call to write_energy.",
              .units = "days"});
  if (energysavedays_geometric > 0.0 && energysavedays_geometric < energysavedays) {
    // defer: the geometric save progression.
    logger::fatal("SumOutput: ENERGYSAVEDAYS_GEOMETRIC is not implemented yet.");
  }

  if (!vgrid.Boussinesq()) {
    // defer: the non-Boussinesq branches of the mass, volume and APE sums.
    logger::fatal("SumOutput: BOUSSINESQ = False is not implemented yet.");
  }
}

void SumOutput::create_depth_list(const Domain &domain, const Grid &grid,
                                  const amrex::Real min_depth_inc) {
  const int ni = domain.ni_global();
  const int nj = domain.nj_global();
  const int mls = ni * nj;

  // One entry per global cell, gathered the way MOM6 gathers it: each rank
  // writes only its own cells and the sum picks up everyone else's zeros.
  std::vector<double> Dlist(static_cast<std::size_t>(mls) + 2, 0.0);
  std::vector<double> Arealist(static_cast<std::size_t>(mls) + 2, 0.0);

  for (amrex::MFIter mfi(grid.bathyT()); mfi.isValid(); ++mfi) {
    const amrex::Box bx = loops::flat(mfi.validbox());
    const amrex::Array4<const amrex::Real> D = grid.bathyT().const_array(mfi);
    const amrex::Array4<const amrex::Real> mask = grid.mask2dT().const_array(mfi);
    const amrex::Array4<const amrex::Real> area = grid.areaT().const_array(mfi);
    for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
      for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
        // MOM6's list_pos, one-based.
        const std::size_t pos = static_cast<std::size_t>(j) * ni + i + 1;
        Dlist[pos] = D(i, j, 0);
        Arealist[pos] = mask(i, j, 0) * area(i, j, 0);
      }
    }
  }
  if (amrex::ParallelContext::NProcsSub() > 1) {
    amrex::ParallelAllReduce::Sum(Dlist.data(), mls + 1,
                                  amrex::ParallelContext::CommunicatorSub());
    amrex::ParallelAllReduce::Sum(Arealist.data(), mls + 1,
                                  amrex::ParallelContext::CommunicatorSub());
  }

  const std::vector<int> indx = sort_index_ascending(Dlist, mls);
  const auto at = [&indx](const int k) { return static_cast<std::size_t>(indx[static_cast<std::size_t>(k)]); };

  // Count the entries that survive culling: identical and very close depths
  // collapse onto one, and the shallowest and deepest are always kept.
  double d_list_prev = Dlist[at(mls)];
  int list_size = 2;
  for (int k = mls - 1; k >= 1; --k) {
    if (Dlist[at(k)] < d_list_prev - min_depth_inc) {
      ++list_size;
      d_list_prev = Dlist[at(k)];
    }
  }

  listsize_ = list_size + 1;
  depth_.assign(static_cast<std::size_t>(listsize_) + 1, 0.0);
  area_.assign(static_cast<std::size_t>(listsize_) + 1, 0.0);
  vol_below_.assign(static_cast<std::size_t>(listsize_) + 1, 0.0);

  double vol = 0.0;
  double area = 0.0;
  double dprev = Dlist[at(mls)];
  d_list_prev = dprev;

  int kl = 0;
  for (int k = mls; k >= 1; --k) {
    const std::size_t i = at(k);
    vol = vol + area * (dprev - Dlist[i]);
    area = area + Arealist[i];

    bool add_to_list = false;
    if (kl == 0 || k == 1) {
      add_to_list = true;
    } else if (Dlist[at(k - 1)] < d_list_prev - min_depth_inc) {
      add_to_list = true;
      d_list_prev = Dlist[at(k - 1)];
    }

    if (add_to_list) {
      ++kl;
      depth_[static_cast<std::size_t>(kl)] = Dlist[i];
      area_[static_cast<std::size_t>(kl)] = area;
      vol_below_[static_cast<std::size_t>(kl)] = vol;
    }
    dprev = Dlist[i];
  }

  // Pad out any entries the count reserved but the loop did not fill, then
  // cap the list with an entry no volume can reach. Both are MOM6's.
  while (kl + 1 < listsize_) {
    ++kl;
    vol_below_[static_cast<std::size_t>(kl)] = vol_below_[static_cast<std::size_t>(kl) - 1] * 1.000001;
    area_[static_cast<std::size_t>(kl)] = area_[static_cast<std::size_t>(kl) - 1];
    depth_[static_cast<std::size_t>(kl)] = depth_[static_cast<std::size_t>(kl) - 1];
  }
  const std::size_t last = static_cast<std::size_t>(listsize_);
  vol_below_[last] = vol_below_[last - 1] * 1000.0;
  area_[last] = area_[last - 1];
  depth_[last] = depth_[last - 1];
}

void SumOutput::write_energy(const State &state, const Domain &domain, const Grid &grid,
                             const VerticalGrid &vgrid, const amrex::Real time,
                             const int n_steps, const amrex::Real dt_forcing) {

  // The save interval, decided the way MOM6 decides it. The driver calls this
  // once per forcing interval rather than once per dynamics step, so a save
  // time that falls inside a forcing interval is written at its end.
  if (previous_calls_ == 0) {
    next_write_time_ = energysavedays_ * std::floor(1.0 + time / energysavedays_);
  } else if (time + 0.5 * dt_forcing < next_write_time_) {
    return;
  } else {
    next_write_time_ += energysavedays_;
  }

  const int nk = vgrid.nk();
  // Boussinesq with thicknesses in metres, so a thickness times this is a
  // mass per unit area. MOM6's GV%H_to_RZ.
  const amrex::Real H_to_RZ = vgrid.Rho0();
  const amrex::Real Rho0 = vgrid.Rho0();

  amrex::MultiFab work = make_layer_field(domain, vgrid, Stagger::Cell);
  work.setVal(0.0);

  // Total mass, and the mass in each layer.
  for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
    const amrex::Box bx = mfi.validbox();
    const amrex::Array4<amrex::Real> t = work.array(mfi);
    const amrex::Array4<const amrex::Real> h = state.h().const_array(mfi);
    const amrex::Array4<const amrex::Real> mask = grid.mask2dT().const_array(mfi);
    const amrex::Array4<const amrex::Real> area = grid.areaT().const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const amrex::Real areaTm = mask(i, j, 0) * area(i, j, 0);
      t(i, j, k) = h(i, j, k) * (H_to_RZ * areaTm);
    });
  }
  std::vector<double> mass_lay;
  EFP mass_EFP;
  const amrex::Real mass_tot = reproducing_sum(work, &mass_lay, &mass_EFP);

  // The potential energy, measured against the depth each layer's volume
  // would rest at.
  std::vector<amrex::Real> Z_0APE(static_cast<std::size_t>(nk) + 2, 0.0);
  amrex::Real PE_tot = 0.0;
  if (do_APE_calc_) {
    std::vector<amrex::Real> vol_lay(static_cast<std::size_t>(nk) + 1, 0.0);
    for (int k = 1; k <= nk; ++k) {
      vol_lay[static_cast<std::size_t>(k)] =
          (1.0 / Rho0) * mass_lay[static_cast<std::size_t>(k) - 1];
    }

    int lbelow = 1;
    double volbelow = 0.0;
    for (int k = nk; k >= 1; --k) {
      volbelow += vol_lay[static_cast<std::size_t>(k)];
      int li = lH_[static_cast<std::size_t>(k)];
      if (!(volbelow >= vol_below_[static_cast<std::size_t>(li)] &&
            volbelow < vol_below_[static_cast<std::size_t>(li) + 1])) {
        int labove = listsize_;
        li = (labove + lbelow) / 2;
        while (li > lbelow) {
          if (volbelow < vol_below_[static_cast<std::size_t>(li)]) {
            labove = li;
          } else {
            lbelow = li;
          }
          li = (labove + lbelow) / 2;
        }
        lH_[static_cast<std::size_t>(k)] = li;
      }
      lbelow = li;
      Z_0APE[static_cast<std::size_t>(k)] =
          depth_[static_cast<std::size_t>(li)] -
          (volbelow - vol_below_[static_cast<std::size_t>(li)]) / area_[static_cast<std::size_t>(li)];
    }
    Z_0APE[static_cast<std::size_t>(nk) + 1] = depth_[2];

    // The interface potential energies live on interfaces, but the layer
    // field is one level short of that; the deepest interface contributes
    // nothing, so the layer levels hold interfaces 1..nk.
    work.setVal(0.0);
    const std::vector<amrex::Real> &gprime = vgrid.g_prime();
    amrex::Gpu::DeviceVector<amrex::Real> d_Z0(Z_0APE.size());
    amrex::Gpu::DeviceVector<amrex::Real> d_gp(gprime.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, Z_0APE.begin(), Z_0APE.end(), d_Z0.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, gprime.begin(), gprime.end(), d_gp.begin());
    const amrex::Real *pZ0 = d_Z0.data();
    const amrex::Real *pgp = d_gp.data();

    for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
      const amrex::Box bx = loops::flat(mfi.validbox());
      const amrex::Array4<amrex::Real> t = work.array(mfi);
      const amrex::Array4<const amrex::Real> h = state.h().const_array(mfi);
      const amrex::Array4<const amrex::Real> mask = grid.mask2dT().const_array(mfi);
      const amrex::Array4<const amrex::Real> area = grid.areaT().const_array(mfi);
      const amrex::Array4<const amrex::Real> bathy = grid.bathyT().const_array(mfi);
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int) {
        const amrex::Real areaTm = mask(i, j, 0) * area(i, j, 0);
        amrex::Real hbelow = 0.0;
        for (int k = nk - 1; k >= 0; --k) {
          hbelow = hbelow + h(i, j, k);
          const amrex::Real hint = pZ0[k + 1] + (hbelow - bathy(i, j, 0));
          amrex::Real hbot = pZ0[k + 1] - bathy(i, j, 0);
          hbot = (hbot + std::abs(hbot)) * 0.5;
          t(i, j, k) = (0.5 * areaTm) * (Rho0 * pgp[k]) * (hint * hint - hbot * hbot);
        }
      });
    }
    PE_tot = reproducing_sum(work);
  }

  // The kinetic energy.
  work.setVal(0.0);
  for (amrex::MFIter mfi(work); mfi.isValid(); ++mfi) {
    const amrex::Box bx = mfi.validbox();
    const amrex::Array4<amrex::Real> t = work.array(mfi);
    const amrex::Array4<const amrex::Real> h = state.h().const_array(mfi);
    const amrex::Array4<const amrex::Real> u = state.u().const_array(mfi);
    const amrex::Array4<const amrex::Real> v = state.v().const_array(mfi);
    const amrex::Array4<const amrex::Real> mask = grid.mask2dT().const_array(mfi);
    const amrex::Array4<const amrex::Real> area = grid.areaT().const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const amrex::Real areaTm = mask(i, j, 0) * area(i, j, 0);
      // gfortran does not contract MOM6's u(I-1)**2 + u(I)**2, so the squares
      // are held back from the sum here too. The terms go into a fixed-point
      // sum, which turns a last-bit difference into a different reported total.
      t(i, j, k) = (0.25 * H_to_RZ * (areaTm * h(i, j, k))) *
                   ((fp_rounded(u(i, j, k) * u(i, j, k)) +
                     fp_rounded(u(i + 1, j, k) * u(i + 1, j, k))) +
                    (fp_rounded(v(i, j, k) * v(i, j, k)) +
                     fp_rounded(v(i, j + 1, k) * v(i, j + 1, k))));
    });
  }
  const amrex::Real KE_tot = reproducing_sum(work);

  // The maximum CFL numbers, transport-based and finite-difference.
  amrex::Real max_CFL_trans = 0.0;
  amrex::Real max_CFL_lin = 0.0;
  for (amrex::MFIter mfi(state.u()); mfi.isValid(); ++mfi) {
    const amrex::Box bxu = loops::u_points(mfi.validbox());
    const amrex::Array4<const amrex::Real> u = state.u().const_array(mfi);
    const amrex::Array4<const amrex::Real> IareaT = grid.IareaT().const_array(mfi);
    const amrex::Array4<const amrex::Real> dy_Cu = grid.dy_Cu().const_array(mfi);
    const amrex::Array4<const amrex::Real> IdxCu = grid.IdxCu().const_array(mfi);
    for (int k = 0; k < nk; ++k) {
      for (int j = bxu.smallEnd(1); j <= bxu.bigEnd(1); ++j) {
        for (int i = bxu.smallEnd(0); i <= bxu.bigEnd(0); ++i) {
          amrex::Real CFL_Iarea = IareaT(i - 1, j, 0);
          if (u(i, j, k) < 0.0) {
            CFL_Iarea = IareaT(i, j, 0);
          }
          const amrex::Real ut = std::abs(u(i, j, k) * dt_);
          max_CFL_trans = std::max(max_CFL_trans, ut * (dy_Cu(i, j, 0) * CFL_Iarea));
          max_CFL_lin = std::max(max_CFL_lin, ut * IdxCu(i, j, 0));
        }
      }
    }
  }
  for (amrex::MFIter mfi(state.v()); mfi.isValid(); ++mfi) {
    const amrex::Box bxv = loops::v_points(mfi.validbox());
    const amrex::Array4<const amrex::Real> v = state.v().const_array(mfi);
    const amrex::Array4<const amrex::Real> IareaT = grid.IareaT().const_array(mfi);
    const amrex::Array4<const amrex::Real> dx_Cv = grid.dx_Cv().const_array(mfi);
    const amrex::Array4<const amrex::Real> IdyCv = grid.IdyCv().const_array(mfi);
    for (int k = 0; k < nk; ++k) {
      for (int j = bxv.smallEnd(1); j <= bxv.bigEnd(1); ++j) {
        for (int i = bxv.smallEnd(0); i <= bxv.bigEnd(0); ++i) {
          amrex::Real CFL_Iarea = IareaT(i, j - 1, 0);
          if (v(i, j, k) < 0.0) {
            CFL_Iarea = IareaT(i, j, 0);
          }
          const amrex::Real vt = std::abs(v(i, j, k) * dt_);
          max_CFL_trans = std::max(max_CFL_trans, vt * (dx_Cv(i, j, 0) * CFL_Iarea));
          max_CFL_lin = std::max(max_CFL_lin, vt * IdyCv(i, j, 0));
        }
      }
    }
  }
  amrex::ParallelDescriptor::ReduceRealMax(max_CFL_trans);
  amrex::ParallelDescriptor::ReduceRealMax(max_CFL_lin);
  amrex::ParallelDescriptor::ReduceIntSum(ntrunc_);

  if (previous_calls_ == 0) {
    mass_prev_ = mass_EFP;
  }
  const amrex::Real mass_chg = (mass_EFP - mass_prev_).to_real();
  // No fresh water enters, so the anomalous change is the whole change.
  const amrex::Real mass_anom = mass_chg;

  const amrex::Real toten = KE_tot + PE_tot;
  const amrex::Real En_mass = toten / mass_tot;
  const amrex::Real reday = time / timeunit_;

  if (amrex::ParallelDescriptor::IOProcessor()) {
    char day_str[32];
    char n_str[32];
    std::snprintf(day_str, sizeof(day_str), "%12.3f", reday);
    std::snprintf(n_str, sizeof(n_str), "%6d", n_steps);

    std::FILE *f = std::fopen(energyfile_.c_str(), previous_calls_ == 0 ? "w" : "a");
    if (f == nullptr) {
      logger::fatal("SumOutput: cannot open ", energyfile_, " for writing.");
    }
    if (previous_calls_ == 0) {
      std::fprintf(f,
                   "  Step,       Day,  Truncs,      Energy/Mass,      Maximum CFL,"
                   "  Mean sea level,   Total Mass,    Frac Mass Err\n");
      std::fprintf(f,
                   "            [days]                 [m2 s-2]           [Nondim]"
                   "        [m]             [kg]           [Nondim]\n");
    }
    std::fprintf(f, "%s,%s,%6d, En %22.16E, CFL %8.5f, SL %11.4E, Mass %11.5E, Me %9.2E\n",
                 n_str, day_str, ntrunc_, En_mass, max_CFL_trans, -Z_0APE[1], mass_tot,
                 mass_anom / mass_tot);
    std::fclose(f);

    logger::note("MOM Day", day_str, " ", n_str, ": En ", En_mass, ", MaxCFL ", max_CFL_trans,
                 ", Mass ", mass_tot);
  }

  if (!std::isfinite(En_mass)) {
    logger::fatal("SumOutput: NaNs in the total model energy.");
  }
  if (En_mass > max_energy_) {
    logger::fatal("SumOutput: energy per unit mass of ", En_mass, " exceeds ", max_energy_, ".");
  }
  if (ntrunc_ > maxtrunc_) {
    logger::fatal("SumOutput: the velocity has been truncated ", ntrunc_,
                  " times, more than MAXTRUNC.");
  }

  ntrunc_ = 0;
  ++previous_calls_;
  mass_prev_ = mass_EFP;
}

} // namespace MOM
