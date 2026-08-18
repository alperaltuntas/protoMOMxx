# Throwaway prototype: the double-gyre backbone and dynamical core

Branch `throwaway_dg`, branched off `horGrid`: the remaining PRs of
`docs/PR_pipeline.md` (5-9), then a working subset of the dynamical core, then
the verification scaffolding the comparison needed, then `ocean.stats` and the
work that closing the bit-for-bit gap took.

The point of the exercise was not the code. It was to find out where the
design bends under the weight of the dynamics, how far bit-for-bit parity
with legacy MOM6 actually goes, and what the remaining work costs.

## 1. What was built

| Commit | Content |
|---|---|
| `fe68f66` | PR 6: VerticalGrid + gprime coordinate (cherry-picked from `verticalGrid`) |
| `443fb58` | PR 5: topography (`flat`, `spoon`) + land/sea masks |
| `b207cee` | PR 7: State (u/v/h), `THICKNESS_CONFIG=uniform`, `VELOCITY_CONFIG=zero` |
| `2089868` | PR 8: Clock + the driver time loop |
| `f92e6b6` | PR 9: MechForcing, the solo driver's 2-gyre winds, `Model::step` dispatch |
| `fc8dd23` | ASCII field dump for parity checking |
| `ba12efc` | The derived grid metrics deferred at PR 4 |
| `44cb557` | The dynamical core: unsplit RK2 + continuity/Coriolis/pressure/viscosity |
| `ae258ed` | `NBOXES`, so layout independence is testable without MPI |
| (new) | `MOM_coms`: the order-invariant global sum (Hallberg & Adcroft 2014) |
| (new) | `MOM_sum_output`: `ocean.stats`, byte-compatible with MOM6's |
| (new) | The bit-for-bit fixes of section 2 |

New code, this branch: about 3,700 lines. The dynamics is 1,903 of them,
against 18,555 lines in the seven MOM6 modules it corresponds to. Most of
that ratio is options protoMOMxx rejects rather than implements. `ocean.stats`
and the reproducing sum are a further 900.

The configuration that runs end to end is the double-gyre testcase with
`SPLIT = False`, `USE_RK2 = True`, `ANALYTIC_FV_PGF = False`,
`BIHARMONIC = False`, `DT = DT_FORCING = 300`. The last three deviations are
forced: the biharmonic branch is not ported; the unsplit scheme resolves the
external gravity wave explicitly, so MOM6 itself is unstable here at the
testcase's 1200 s step; and MOM6 recalculates the bottom boundary layer once
per forcing cycle rather than once per dynamics step, so its answer depends on
`DT_FORCING` (see section 2).

### `ocean.stats`

The ASCII energy file is written by `src/diagnostics/MOM_sum_output`, the
analogue of MOM6's `write_energy`, and it is byte-compatible with MOM6's:
the same columns, the same Fortran edit descriptors, the same save cadence.
Three pieces had to come with it, and all three are reusable:

- **`MOM_coms`**, the extended-fixed-point sum of Hallberg & Adcroft (2014).
  A real number is held as six 64-bit integers, so adding terms is integer
  addition, which is exact and associative. This is `DESIGN.md` section 5's
  reduction service, and it is what makes a global sum independent of the
  decomposition. About 190 lines, and unit-tested for exactness, order
  invariance and layout independence.
- **The hypsometry**, MOM6's depth list: the globally sorted bottom depths
  with the open area and open volume below each, used to find the depth each
  layer's volume would rest at and hence the available potential energy. It
  needs MOM6's indexed heapsort transcribed rather than replaced, because the
  order of equal depths decides the order the areas are accumulated in and so
  shows up in the last digits.
- **The velocity truncation** (MOM6's `vertvisc_limit_vel`), which caps
  velocities at `CFL_TRUNCATE` and counts what it capped. It never fires in
  this testcase, but without it the two models would part company as soon as a
  run went unstable, and the `Truncs` column would be a lie.

The first line of both models' `ocean.stats` -- energy, CFL, mean sea level and
total mass at step 0 -- is byte-identical, which exercises the whole of the
above against MOM6 in one comparison.

## 2. Bit-for-bit parity with legacy MOM6

Reference: `MOM6_using_TIM` built with gcc 14.3 on derecho, same compiler as
protoMOMxx, same 44x40x2 configuration, one rank.

**The fixed initialization is bit-identical.** Every field in MOM6's
`ocean_geometry.nc` -- `geolat`/`geolon` at h/u/v/q, `dxT`, `dyT`, `dxCu`,
`dyCu`, `dxCv`, `dyCv`, `dxBu`, `dyBu`, `Ah`, `D`, `wet`, `f` -- matches
exactly at every point. MOM6's own bit-count checksums agree too: `areaT`
53108, `dxT` 53284, initial `h` 82985.

**The prognostic state is bit-identical for 108 dynamics steps** -- nine
simulated hours. Every value of `u`, `v` and `h`, at every point, matches
MOM6's per-step diagnostic output exactly through step 108; the first
difference is one thickness value, one ulp, at step 109.

**Ten simulated days agree to roundoff.** `ocean.stats`, written by both
models, is byte-identical at step 0, and thereafter the energy per unit mass
agrees to about 1e-13 relative while the CFL, mean sea level and total mass
columns stay byte-identical:

| dynamics steps | max rel. diff, u | max rel. diff, h | rel. diff in reported energy |
|---:|---:|---:|---:|
| 1 - 108 | 0 | 0 | 0 at step 0 |
| 200 | 7.5e-12 | 8.0e-16 | -- |
| 288 (1 day) | 7.7e-12 | 1.0e-15 | 2.4e-13 |
| 2880 (10 days) | -- | -- | 4.9e-14 |

Getting there took several fixes, and the interesting thing about them is that
only two are numerical errors.

### The reference was running a different scheme

protoMOMxx read a parameter called `PRESSUREFORCE`, a string with values
`"FV"` and `"Montgomery"`. MOM6 has no such parameter: it selects the pressure
gradient with `ANALYTIC_FV_PGF`, a Boolean that defaults to true. So
`#override PRESSUREFORCE = "Montgomery"` in the MOM6 side of the comparison
was read by nothing, MOM6 ran the finite-volume form, and protoMOMxx ran the
Montgomery form. The two agree analytically for a layered adiabatic Boussinesq
ocean and differ in their discretization, which is why the symptom looked like
a small kernel error -- about 1e-13 in the pressure force at the first step --
rather than the wrong equation.

That cost most of the debugging time in this exercise, and the lesson is
narrow and practical: **an invented parameter name is worse than an
unimplemented option, because MOM6 does not complain about an override it
never reads.** The parameter is now `ANALYTIC_FV_PGF`, with the same name,
type, default and description as MOM6's, and protoMOMxx aborts when it is
true. Every parameter protoMOMxx reads should be checked against MOM6's, and
that check is mechanical: the names appear side by side in the two
`MOM_parameter_doc.all` files.

### Two real errors

`set_viscous_BBL` was being called with `h_av`, the half-step thickness, and
from inside the time stepping scheme. MOM6 calls it from `step_MOM` before the
scheme runs, with the thickness at the start of the step. The two agree at the
first step, when `h_av` is `h`, and diverge afterwards.

The continuity solver floored the updated thickness at zero. MOM6 floors it at
one Angstrom, in both the zonal and the meridional pass. It does not bind in
this configuration, so it cost nothing here and would have cost a great deal
somewhere else.

A third, smaller one is worth its own line because of what it says about
transcription. `GV%H_subroundoff` -- the thickness MOM6 adds to a denominator
so that a vanishing layer does not divide by zero -- had been transcribed as
the literal `1.0e-30`. MOM6 computes it as `1e-20 * max(Angstrom, 1e-17)`,
which at the default Angstrom is **not** the same double: it is one ulp below
`1e-30`. The value appears in the horizontal viscosity, the vertical friction,
the bottom boundary layer and the potential vorticity denominator, and in the
last of those it is what a vanishing layer's thickness is compared against. A
constant that "obviously" equals a decimal literal is not a constant, it is an
expression.

### Floating-point contraction has to be matched expression by expression

`DESIGN.md` section 5 says "matched FMA contraction", and the earlier draft of
this document read that as a build flag. It is not. Both compilers fuse
`a*b + c` into an FMA at `-O`, but contraction only ever happens inside a
single expression, and the two models write the same formula with different
expression boundaries:

- Where MOM6 writes an intermediate into an array before using it, gfortran
  cannot fuse across the store. `MOM_CoriolisAdv`, `MOM_continuity_PPM` and
  `MOM_hor_visc` do this for nearly everything they compute -- `dvdx`, `dudy`,
  `rel_vort`, `KEx`, the strain rates, the PPM edge values -- so gfortran
  contracts almost nothing in them, while the same quantities written as C++
  locals are fused into whatever comes next. Those three files are compiled
  with `-ffp-contract=off`.
- Inside each of them there are a few single statements gfortran *does*
  contract, and those are written out as `std::fma`: the outer sum of the
  stress divergence in `diffu`/`diffv`, and two expressions in the continuity
  solver's positive-definite limiter. Three call sites in about four thousand
  lines.
- Where MOM6 reads and writes the same array element in one statement --
  `u(I,j,k) = u(I,j,k) + I_Hmix*hfr*stress` in `vertvisc` -- gfortran also
  declines to fuse, and g++ does. `MOM_fp_contract.h` adds an empty-asm barrier
  for that one.
- Everywhere else contraction must stay **on**: a build-wide
  `-ffp-contract=off` moves the fixed initialization *away* from MOM6.

The cost is about 4% of the main loop, which is what disabling contraction in
the three kernels buys back in rounding fidelity.

This is the least portable thing in the branch, and it should be recorded as a
cost of the parity requirement rather than hidden. It also gives a design
rule with teeth: **transcribe MOM6's statement boundaries, not just its
formulas.** An intermediate that MOM6 stores in an array is part of the
numerical contract.

Finding which expressions those are is not guesswork. Given a stage whose
inputs are known bit-for-bit -- and after the first step, every input to every
kernel is -- the fusion pattern can be *solved for*: reimplement the kernel in
Python with exact (rational) fused multiply-add, enumerate the placements, and
score each against MOM6's own output. That is how the horizontal viscosity was
closed (one FMA, in the outer sum of the stress divergence, and nothing else in
the routine) and how the continuity solver was closed (two FMAs, both inside
the positive-definite limiter, and nothing else). Each search took seconds and
returned an exact match on every point.

What survives is smaller than the tool can see. The Coriolis term at v points
differs from MOM6 by two or three population counts at the second step, at
points in the vanishing bottom layer over the shelf, and its minima, maxima and
every other stage of the step agree exactly. It is not localizable further,
because MOM6 posts no diagnostic for `CAv` and none for the intermediate
transports the predictor stage uses, so there is no pointwise reference to
solve against -- only a checksum, which is too weak to discriminate. The
prognostic state absorbs it: `u`, `v` and `h` stay bit-identical for another
hundred steps.

### MOM6's answer depends on DT_FORCING

With `ADIABATIC = True`, MOM6 recalculates the bottom boundary layer only when
`bbl_time_int > 0`, which for the solo driver means once per forcing cycle,
not once per dynamics step. Running the same configuration with
`DT_FORCING = 600` and `DT_FORCING = 300` changes MOM6's own answer by 2e-7
over a day. `Model::step` has no notion of a cycle and recalculates every
step, which is MOM6's behaviour when the two timesteps are equal, so the
testcase now sets them equal.

This is a genuine gap, not a bookkeeping detail: MOM6's `step_MOM` carries
`start_cycle`/`end_cycle`/`cycle_length` through the whole call and several
subsystems change behaviour on them. protoMOMxx's `step(forces, dt_forcing,
n_steps)` has thrown that away, and the cycle will have to come back before
the thermodynamic timestep does.

### Three specification details, unchanged from the earlier draft

- **Halo fill values are part of the numerical contract.** MOM6 allocates `h`
  at `GV%Angstrom_H`, not zero, and the thickness configurations write only
  the computational domain, so that value is what the halo outside the global
  domain keeps and what the boundary pressure gradients are computed from.
- **MOM6 reports a velocity field over its non-symmetric index range**, which
  drops the low boundary face. Reproducing that range is what makes the
  stage-by-stage checksum comparison usable.
- **The Adcroft reciprocal** (`1/x`, or 0 where `x` is 0) is a convention, not
  an optimization, and it appears in every derived metric.

**Conclusion on parity:** it is achievable and it was achieved -- exactly for
the first 108 dynamics steps, and to roundoff for ten simulated days. What it
costs is a parity harness used continuously, a discipline about parameter
names, and per-expression attention to contraction. None of that is
discoverable by reading the Fortran, and all of it is mechanical once the
harness exists.

## 3. Layout independence

`NBOXES` splits a rank's share of the domain into several AMReX boxes, which
exercises the halo paths without MPI. After one day of the unsplit RK2
dynamics, `NBOXES` = 1, 2, 4, 8 and 16 give **identical checksums and
identical minima and maxima** for u, v and h.

That is the strongest single result of the exercise. All of the index-
convention translation -- MOM6's `I = i-1` face labelling, the loop bounds,
the halo widths each kernel reads -- came out right on the first try, and
AMReX's decomposition and halo exchange needed no special handling.

The reduction service of `DESIGN.md` §5 is now built (`MOM_coms`), and with it
`ocean.stats` is **byte-identical across `NBOXES` = 1, 2, 4, 8 and 16** for the
whole ten-day run, global means included. The earlier draft recorded the
missing reproducing sum as the one place layout independence did not hold; it
now holds everywhere that is testable without MPI.

## 4. Performance

44x40x2, 2,880 dynamics steps, one rank, gcc 14.3, MOM6 at its own release
flags and protoMOMxx at `-O3`:

| | main loop | per step |
|---|---:|---:|
| MOM6 | 4.038 s | 1.402 ms |
| protoMOMxx | 2.362 s | 0.820 ms |

Both runs write `ocean.stats` and no other output. protoMOMxx is **1.7x
faster**, after paying about 4% for the disabled contraction of section 2.
Read that carefully:

- MOM6's diagnostics cost only 2.7% here (measured by running with an empty
  `diag_table`), so their absence is not the explanation.
- Neither model runs the barotropic solver, tracers or thermodynamics in this
  configuration, so the comparison is over the same physics.
- MOM6's kernels carry runtime branches for the many options protoMOMxx
  rejects outright, and MOM6 does halo exchanges that a single-box serial run
  does not.
- 44x40x2 fits in cache, so this measures instruction count, not memory
  behaviour. Do not extrapolate to a real grid.

The honest reading: C++ and AMReX are not a performance liability at this
scale, and the per-call `MultiFab` allocation in the kernels (see §5) is not
yet visible in the profile.

## 5. Design assessment

### What held up

**Constructor injection and full initialization (§14).** Adding
`VerticalGrid`, `State`, `MechForcing` and `Dynamics` to `Model` was
mechanical each time, and the initializer list reads as the phase order of
`initialize_MOM`. No two-phase initialization appeared anywhere, including in
the dynamics, which in MOM6 is the worst offender.

**The GridFields to Grid handoff.** The grid grew four times in this branch
(topography, masks, derived metrics, OBC-masked reciprocals) and each growth
cost a struct field and an accessor. The constructor's "every field is
created" check caught two omissions the moment a test ran.

**Deferred branches abort (§9).** This is the highest-value convention in the
document, by a wide margin. Running the stock double-gyre configuration
walked the prototype from fatal to fatal: `SPLIT`, then `ANALYTIC_FV_PGF`, then
`BIHARMONIC`, then `HARMONIC_VISC`, then `LINEAR_DRAG`. Every one of those was
a parameter that would otherwise have silently taken a default the code does
not implement. The set of live aborts is, in effect, a machine-checked list of
the remaining work.

**`TIM::Stagger` plus the field factory.** Every field creation reads
`make_field(Stagger::XFace, nk, 1)`, and AMReX's `IndexType` congruence turns
most stagger confusions into runtime failures instead of silent aliasing.

### What bent or broke

**The per-directory library graph did not survive the dynamics, and has been
fixed.** The one-way `framework <- initialization <- core` graph held through
PR 5. PR 7 bent it: `src/initialization` could not see `Grid` or
`VerticalGrid`, so `initialize_state` took a `StateSpec` of plain values
copied out of them. The dynamics broke it outright: `src/parameterizations`
uses `MOM_grid`, and `src/core`'s stepper calls the parameterizations, so the
two shared one CMake library. MOM6 has the same cycles and hides them behind a
single link unit; protoMOMxx had inherited the directory layout without
inheriting that.

The fix, and the one place this document's recommendation has already been
acted on, is a `src/types` layer below both: `Grid`, `VerticalGrid`, `State`,
`MechForcing` and their construction precursors, with the MOM6 file names
kept. The graph is now

    infra <- framework <- types <- { initialization, parameterizations,
                                     diagnostics } <- core

with no cycles and nothing sharing a library to hide one. Three things fell
out of it immediately: `StateSpec` is gone, and `initialize_state` takes the
`Grid` and the `VerticalGrid` it always wanted; `src/parameterizations` is its
own library; and `src/diagnostics` depends on the types rather than on the
whole core. `MOM_loop_boxes.h` moved to `src/framework` on the way, since the
parameterizations need it and it is a helper over `amrex::Box` with no model
dependency.

The cost is the one deviation from MOM6's directory layout in the branch:
MOM6 keeps `MOM_grid.F90`, `MOM_verticalGrid.F90` and `MOM_forcing_type.F90`
in its `src/core`. The file names are unchanged, so a MOM6 developer still
recognizes them, and `src/types/README.md` states the rule for what belongs
there. That is the right trade: the alternative is a directory layout that
matches MOM6's and a link graph that cannot be expressed.

**Two-dimensional work fields under a three-dimensional MFIter.** The only
crash of the exercise, and it is instructive. `mfi.validbox()` taken from an
`MFIter` over a 3-D `MultiFab` carries `k` in `[0, nk-1]`. Using it to drive a
`ParallelFor` that writes a 2-D work field writes past the end of the FAB.
AMReX does not catch that in a Release build: it corrupted the heap and
surfaced at exit inside AMReX's own `FabArrayBase` registry, with a backtrace
pointing nowhere near the bug. This is `DESIGN.md` §4's "C++ UB silently
corrupts results" threat, realized on the first serious kernel.

Three mitigations, in increasing order of value:
1. `loops::flat()` and `loops::layers()` (added here) make the vertical extent
   of a loop explicit at the call site.
2. A debug build with AMReX assertions in CI. There is none today.
3. A field type that carries its own vertical residency, so the mismatch is a
   compile error. This is the layer-versus-interface half of the typed-field
   decision parked in `PR_pipeline.md`, and the prototype is a concrete
   argument for doing it.

**Named iteration boxes turned out to be a design element, not a helper.**
`MOM_loop_boxes.h` was written under duress and became the single place where
MOM6's `I = i-1` face labelling is translated. Every index bug in the exercise
was in that translation. It should be promoted, and extended so a box carries
its stagger and vertical extent.

**Work fields are allocated per call.** `CoriolisAdv` creates seven
`MultiFab`s per call, `hor_visc` six, `continuity` three. MOM6 uses stack
automatics. It costs nothing measurable at this size but it is the obvious
first optimisation and it argues for the term modules owning their scratch,
the way the stepper already owns the accelerations.

**Fixed-size stack arrays in the column solves.** `MOM_vert_friction` uses
`double z_i[64]` and `c1[64]` for the tridiagonal solve. That is an unguarded
limit on `nk`. MOM6 uses automatic arrays sized `SZK_(GV)`. This needs either
a real check or per-column scratch.

**The argument-train claim needs qualifying.** `DESIGN.md` §3 says the
`(G, GV, US, CS, ...)` trains are what real classes replace. Half of that
happened: every term module owns its own state, so `CS` is gone, and `US` was
never born. But the kernels still take `(domain, grid, vgrid)` because they
genuinely need all three. The spec-struct pattern that worked for `Domain` and
`Grid` construction does not scale to the dynamics, and pretending otherwise
would just rename the train. The honest statement is that the control
structure is replaced, not the grid arguments -- and the control structure was
the larger problem.

**Two parameter-system gaps.** `RuntimeParams` has no log-only entry point, so
`MAXIMUM_DEPTH` is read twice to reproduce MOM6's documentation order. And the
override semantics differ from MOM6's: `#override` requires the key to already
exist, and a plain reassignment across files is an error, so a fixture that
selects a non-default path has to know which of the two forms to use. Both are
small and both cost time repeatedly.

### The AMReX boundary reassessment (§7)

The threshold is passed: this branch adds about 1,900 lines of non-trivial
AMReX code. The reassessment answer, on this evidence:

- **Do not abstract `MultiFab`.** It carried the decomposition, the halo
  exchange and the layout independence with no friction at all.
- **Do abstract the index space.** Every bug was in `Box` arithmetic, not in
  field access. The kernel signature `(Array4s..., Box)` that §7 wants is
  already what the code looks like inside the `MFIter` loops; what is missing
  is that the `Box` is untyped -- it does not know its stagger or its vertical
  extent, and both mismatches are silent.
- `amrex::Real` as the only scalar type has not pushed back yet, because the
  transcription is unit-for-unit with MOM6. The units layer will be the first
  real test.

So: a thin layer over the *iteration box* and the *field residency*, and
nothing over the field itself. That is small, it is debuggable, and it keeps
the exit path.

## 6. Refined plan for the full double-gyre

To run the stock configuration (`SPLIT = True`, `ANALYTIC_FV_PGF = True`,
biharmonic viscosity, `DT = 1200`), in dependency order:

1. ~~**Types layer split**~~ -- done: `src/types` holds `Grid`,
   `VerticalGrid`, `State` and `MechForcing` below `src/initialization`,
   `src/parameterizations` and `src/diagnostics`. It was as cheap as predicted
   and it is now a precondition rather than a task.
2. **Parity harness** -- `report_field` plus a script that diffs a protoMOMxx
   log against a MOM6 `DEBUG` log stage by stage, in CI. About 200 lines, and
   the diff script is another 60. This made the exercise tractable and belongs
   before the dynamics, not after. `PR_pipeline.md` currently schedules the
   checksum oracle at PR 7; it should be its own PR, it should reproduce MOM6's
   ranges exactly, and it should also diff the two `MOM_parameter_doc.all`
   files, because the largest error found here was a parameter name that only
   one of the two models had.
3. **Typed field residency** -- layer versus interface versus 2-D, and stagger,
   in the type system. The one crash of this exercise is a compile error under
   it.
4. **Units layer (§6)** -- unchanged as the named prerequisite for kernels.
   Everything below re-types, so it must precede them.
5. **Continuity's barotropic machinery** -- `uhbt`/`vhbt` targets, the Newton
   iteration in `zonal_flux_adjust`, `BT_cont`. About 1,200 lines of MOM6 that
   this prototype skipped, and it is testable against MOM6 in isolation.
6. **`MOM_barotropic` + `MOM_dynamics_split_RK2`** -- 8,977 lines of MOM6, the
   single largest item, and larger than everything in this branch by three to
   five times. It should be its own multi-PR track, not one PR.
7. **`MOM_PressureForce_FV`** -- needs `MOM_density_integrals`; a reduced form
   is possible for the layered no-EOS case.
8. **Biharmonic horizontal viscosity** -- about 400 lines of `MOM_hor_visc`.
9. **Diagnostics** -- the NetCDF side, through TIM. `ocean.stats` and the
   reproducing sum it needs are now built (`MOM_sum_output`, `MOM_coms`) and
   should move to the mainline early: it is the cheapest whole-model
   comparison there is, it proved layout independence of the global means, and
   it is what surfaced the `DT_FORCING` behaviour above.
10. **Restart** -- needed for the restart-exactness invariant of §5.

Not needed for double-gyre and correctly deferred: non-linear bottom drag,
open boundaries, ALE, thermodynamics, tracers.

Items 1-3 are new relative to `PR_pipeline.md` and all three are cheap. Item 6
is the schedule.

## 7. Lessons, in the "build one to throw away" sense

**The prototype found the load-bearing walls the plan did not predict.** The
PR pipeline expected the risk to be AMReX fit. The actual friction was the
library dependency graph and the index-convention translation, neither of
which appears in `DESIGN.md`. AMReX itself was the easy part. The graph one
has since been fixed on the mainline branch and merged back here, which is the
throwaway working as intended: the finding outlived the code that produced
it.

**A throwaway is only trustworthy if it refuses to guess.** The
abort-on-deferred-branch policy is what made this safe. A prototype that runs
the stock input file and produces plausible numbers while silently defaulting
five unimplemented options is worse than no prototype, because it launders
missing work into apparent progress.

**Parity against the reference is a design tool, not just a verification
tool.** The bit-level comparison surfaced specification details -- the halo
fill value, per-expression contraction, MOM6's reporting ranges, the Adcroft
reciprocal, the forcing-cycle dependence of the bottom boundary layer -- that
no amount of reading the Fortran produced. Each one would have been discovered
eventually, at much higher cost, inside a science bug.

**The reference has to be pinned as hard as the port.** The single largest
error in this exercise was on the MOM6 side: an override MOM6 silently ignored,
so the two models ran different pressure gradients for the whole comparison. A
harness that only compares outputs cannot see that. What would have caught it
in a minute is a diff of the two `MOM_parameter_doc.all` files.

**"Matched contraction" is a per-expression property, and it is solvable.** It
was written down as a build flag and it is not one. What has to match is where
each model rounds, and that is set by where MOM6 stores an intermediate in an
array. A transcription that preserves MOM6's formulas but collapses its
statements into one expression will not reproduce it, and no build flag fixes
that. But it does not have to be guessed either: with the inputs to a kernel
known bit-for-bit, reimplementing it in Python with exact fused multiply-add
and enumerating the placements finds the pattern that reproduces MOM6 on every
point, in seconds. That technique should be part of the parity harness, not a
one-off.

**The cost is now measurable.** Backbone plus a working subset dynamical core:
about 3,700 lines. That is small enough that discarding it is
genuinely affordable, and small enough that the estimate for the full
double-gyre -- dominated by the barotropic solver at three to five times this
size -- is credible rather than a guess.

**What to keep from the throwaway:** the `src/types` layering, which is
already back on `horGrid`; `MOM_loop_boxes.h`; `report_field` and its
MOM6-compatible checksum; `MOM_coms` and `MOM_sum_output` (both faithful to
MOM6 and both cheap to test); the `NBOXES` layout test; the contraction
findings; and the inventory of live aborts. **What to throw away:** the kernels
themselves. They should be re-derived on top of the units layer and the typed
fields, not retrofitted -- retrofitting numerics is exactly what §6 warns
against.

## 8. What was not tested

- **MPI.** No host list is available on the login node, so the `mpi`-labelled
  tests could not run. Layout independence was tested with multiple boxes on
  one rank instead, which exercises the same halo code but not the exchange.
- **GPU.** Not attempted.
- **Restart exactness and rotational symmetry**, two of the four §5
  invariants; neither has an implementation to test.
- **The bit-for-bit result on another machine or compiler.** The contraction
  matching of §2 is specific to this pair of builds. A different `-march`, a
  different GCC, or a non-GNU compiler will move the fusion decisions again.
  Whether MOM6 parity should be claimed per build pair, or whether both models
  should be built with contraction off for the comparison, is a decision the
  project has not made and should.
