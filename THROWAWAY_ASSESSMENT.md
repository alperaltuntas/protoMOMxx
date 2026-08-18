# Throwaway prototype: the double-gyre backbone and dynamical core

Branch `throwaway_dg`, branched off `horGrid`. Nine commits: the remaining
PRs of `docs/PR_pipeline.md` (5-9), then a working subset of the dynamical
core, then the verification scaffolding the comparison needed.

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

New code, this branch: about 2,700 lines. The dynamics is 1,903 of them,
against 18,555 lines in the seven MOM6 modules it corresponds to. Most of
that ratio is options protoMOMxx rejects rather than implements.

The configuration that runs end to end is the double-gyre testcase with
`SPLIT = False`, `USE_RK2 = True`, `PRESSUREFORCE = "Montgomery"`,
`BIHARMONIC = False`, `DT = 300`. The last two deviations are forced:
the biharmonic branch is not ported, and the unsplit scheme resolves the
external gravity wave explicitly, so MOM6 itself is unstable here at the
testcase's 1200 s step.

## 2. Bit-for-bit parity with legacy MOM6

Reference: `MOM6_using_TIM` built with gcc 14.3 on derecho, same compiler as
protoMOMxx, same 44x40x2 configuration, one rank.

**The fixed initialization is bit-identical.** Every field in MOM6's
`ocean_geometry.nc` -- `geolat`/`geolon` at h/u/v/q, `dxT`, `dyT`, `dxCu`,
`dyCu`, `dxCv`, `dyCv`, `dxBu`, `dyBu`, `Ah`, `D`, `wet`, `f` -- matches
exactly at every point (1,760 to 1,845 points each, zero mismatches). MOM6's
own bit-count checksums agree too: `areaT` 53108, `dxT` 53284.

**The initial state is bit-identical**: `h` checksum 82985.

**The dynamics agrees to roundoff for a few steps, then diverges
exponentially**, which is what a chaotic system does with a roundoff-level
seed:

| dynamics steps | max(u) rel. diff | min(v) rel. diff | max(h) rel. diff |
|---:|---:|---:|---:|
| 2 | 1.2e-13 | 1.5e-16 | 3.0e-14 |
| 10 | 7.2e-15 | 0 | 1.0e-13 |
| 50 | 2.9e-10 | 2.9e-10 | 1.1e-12 |
| 288 (1 day) | 1.8e-07 | 2.9e-07 | 4.3e-13 |

The first stage that differs is the pressure force, by about 1e-13 relative.
The difference is a few ulp at interior points and larger at the domain-edge
faces, where the Montgomery potential differences involve heavy cancellation.
It is not yet localized further.

Three specification details had to be matched before even the initialization
agreed, and none of them is visible in the equations:

- **Floating-point contraction must *match* the Fortran side, not be turned
  off.** Adding `-ffp-contract=off` to protoMOMxx alone moved `bathyT` away
  from MOM6, because the reference gfortran build contracts. `DESIGN.md` §5
  already says "matched FMA contraction"; the prototype confirms the word
  that matters is *matched*. The build option is now `PROTOMOM_NO_FP_CONTRACT`,
  default off.
- **Halo fill values are part of the numerical contract.** MOM6 allocates `h`
  at `GV%Angstrom_H`, not zero. Since the thickness configurations write only
  the computational domain, that value is what the halo outside the global
  domain keeps, and the boundary pressure gradients are computed from it.
  Filling with zero instead changed boundary values by 4e-10 relative.
- **MOM6 reports a velocity field over its non-symmetric index range**, which
  drops the low boundary face. Comparing statistics without reproducing that
  range produces differences that look like bugs and are not.

**Conclusion on parity:** it is achievable, and it already holds for
everything up to the first kernel. Closing the last 1e-13 in the kernels is
ordinary work, but it needs the stage-by-stage comparison as a standing tool.
`report_field` (MOM6's statistics ranges and bit-count checksum, reproduced
exactly) is that tool and should move into the mainline before the dynamics
PRs, not after.

## 3. Layout independence

`NBOXES` splits a rank's share of the domain into several AMReX boxes, which
exercises the halo paths without MPI. After one day of the unsplit RK2
dynamics, `NBOXES` = 1, 2, 4, 8 and 16 give **identical checksums and
identical minima and maxima** for u, v and h.

That is the strongest single result of the exercise. All of the index-
convention translation -- MOM6's `I = i-1` face labelling, the loop bounds,
the halo widths each kernel reads -- came out right on the first try, and
AMReX's decomposition and halo exchange needed no special handling.

The domain *means* differ in their last digits across box counts, because the
global sum is not reproducing. That is the reduction service of `DESIGN.md`
§5, and it is not built.

## 4. Performance

44x40x2, 2,880 dynamics steps, one rank, gcc 14.3 `-O2`:

| | main loop | per step |
|---|---:|---:|
| MOM6 | 3.955 s | 1.373 ms |
| protoMOMxx | 2.197 s | 0.763 ms |

protoMOMxx is **1.8x faster**. Read that carefully:

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
walked the prototype from fatal to fatal: `SPLIT`, then `PRESSUREFORCE`, then
`BIHARMONIC`, then `HARMONIC_VISC`, then `LINEAR_DRAG`. Every one of those was
a parameter that would otherwise have silently taken a default the code does
not implement. The set of live aborts is, in effect, a machine-checked list of
the remaining work.

**`TIM::Stagger` plus the field factory.** Every field creation reads
`make_field(Stagger::XFace, nk, 1)`, and AMReX's `IndexType` congruence turns
most stagger confusions into runtime failures instead of silent aliasing.

### What bent or broke

**The per-directory library graph does not survive the dynamics.** The
one-way `framework <- initialization <- core` graph held through PR 5. PR 7
bent it: `src/initialization` cannot see `Grid` or `VerticalGrid`, so
`initialize_state` takes a `StateSpec` of plain values instead. The dynamics
broke it outright: `src/parameterizations` uses `MOM_grid`, and `src/core`'s
stepper calls the parameterizations, so the two now share one CMake library.
MOM6 has the same cycles and hides them behind a single link unit; protoMOMxx
inherited the directory layout without inheriting that.

The fix is to split the *types* -- `Grid`, `VerticalGrid`, `State`,
`MechForcing` -- into a layer below both `src/initialization` and the
dynamics, keeping their MOM6 file names. It is cheap now and expensive later.

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

To run the stock configuration (`SPLIT = True`, `PRESSUREFORCE = "FV"`,
biharmonic viscosity, `DT = 1200`), in dependency order:

1. **Types layer split** -- `Grid`, `VerticalGrid`, `State`, `MechForcing`
   below `src/initialization` and `src/parameterizations`. Small; do it first,
   before more code is written against the current graph.
2. **Parity harness** -- `report_field` plus a script that diffs a protoMOMxx
   log against a MOM6 `DEBUG` log stage by stage, in CI. About 200 lines. This
   made the exercise tractable and belongs before the dynamics, not after.
   `PR_pipeline.md` currently schedules the checksum oracle at PR 7; it should
   be its own PR and it should reproduce MOM6's ranges exactly.
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
9. **Diagnostics and `ocean.stats`** -- needed for energy-level comparison,
   and needs the reproducing-sum service, which is also what makes the domain
   means layout-independent.
10. **Restart** -- needed for the restart-exactness invariant of §5.

Not needed for double-gyre and correctly deferred: non-linear bottom drag,
open boundaries, ALE, thermodynamics, tracers.

Items 1-3 are new relative to `PR_pipeline.md` and all three are cheap. Item 6
is the schedule.

## 7. Lessons, in the "build one to throw away" sense

**The prototype found the load-bearing walls the plan did not predict.** The
PR pipeline expected the risk to be AMReX fit. The actual friction was the
library dependency graph and the index-convention translation, neither of
which appears in `DESIGN.md`. AMReX itself was the easy part.

**A throwaway is only trustworthy if it refuses to guess.** The
abort-on-deferred-branch policy is what made this safe. A prototype that runs
the stock input file and produces plausible numbers while silently defaulting
five unimplemented options is worse than no prototype, because it launders
missing work into apparent progress.

**Parity against the reference is a design tool, not just a verification
tool.** The bit-level comparison surfaced four specification details -- the
halo fill value, contraction matching, MOM6's reporting ranges, the Adcroft
reciprocal convention -- that no amount of reading the Fortran produced. Each
one would have been discovered eventually, at much higher cost, inside a
science bug.

**The cost is now measurable.** Backbone plus a working subset dynamical core:
about 2,700 lines over nine commits. That is small enough that discarding it
is genuinely affordable, and small enough that the estimate for the full
double-gyre -- dominated by the barotropic solver at three to five times this
size -- is credible rather than a guess.

**What to keep from the throwaway:** `MOM_loop_boxes.h`, `report_field` and
its MOM6-compatible checksum, the `NBOXES` layout test, the contraction
finding, and the inventory of live aborts. **What to throw away:** the kernels
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
