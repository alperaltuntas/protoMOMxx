# Loop order: why k innermost or outermost decides the speed

A three-dimensional kernel has to pick which index runs in the innermost loop.
The choice is invisible in the physics and it changes the run time by a factor
of three to eleven. It is also the one place where a CPU and a GPU want
opposite things, which is why it is worth writing down before more kernels are
ported.

This note explains the mechanism, shows the measurement that isolates it, and
records where protoMOMxx and MOM6 each stand today.

## Background: what the memory looks like

Both models store a field as a single contiguous block with `i` fastest. In
AMReX an `Array4` indexes it as

```
    a(i, j, k)  ->  base[ (i - lo.x)
                        + (j - lo.y) * jstride
                        + (k - lo.z) * kstride ]

    jstride = nx                 (the grown x extent: NIGLOBAL + 2*NIHALO)
    kstride = nx * ny
```

Fortran's `a(i,j,k)` with `i` declared first is the same layout. So:

- stepping `i` by one moves **8 bytes** -- the next element in the same cache line;
- stepping `j` by one moves **`nx` * 8 bytes** -- a different line, usually a
  different page once the grid is large;
- stepping `k` by one moves **`nx*ny` * 8 bytes** -- a whole horizontal plane.

For the double_gyre grid at `NIGLOBAL=44, NJGLOBAL=40` with four halo points,
`nx*ny = 52*48 = 2496` doubles, so one step in `k` jumps **20 KB**. At
`132x120` it jumps **143 KB**.

Two hardware facts follow from that layout, and they are the whole story:

1. **A cache line is 64 bytes -- eight doubles.** Touching one element brings
   in eight. If the loop then moves along `i`, the next seven iterations are
   free. If it moves along `k`, the other seven are wasted and the line is
   likely evicted before anything comes back for them.
2. **A vector instruction works on four doubles at once** (256-bit AVX2; eight
   with AVX-512). The compiler can only use it when consecutive iterations
   touch consecutive addresses -- which means the vectorised loop must be the
   one running along `i`.

## The two shapes

Take the simplest kernel with a vertical dependence: the interface heights,
built upward from the bottom.

    e(i,j,nk) = -D(i,j)
    e(i,j,k)  = e(i,j,k+1) + h(i,j,k)     for k = nk-1 ... 0

Each column is independent; within a column the recursion is strictly
sequential. There are two ways to write it.

**Shape A -- column-parallel, k innermost.** One task per column, the
recursion inside it.

```cpp
    ParallelFor(columns_2d, [=] (int i, int j, int) {
      ee(i, j, nk) = -D(i, j, 0);
      for (int k = nk - 1; k >= 0; --k)
        ee(i, j, k) = ee(i, j, k + 1) + hh(i, j, k);   // stride nx*ny
    });
```

**Shape B -- plane sweeps, k outermost.** One pass per interface, each pass a
full horizontal plane.

```cpp
    ParallelFor(cells_2d, [=] (int i, int j, int) { ee(i, j, nk) = -D(i, j, 0); });
    for (int k = nk - 1; k >= 0; --k)
      ParallelFor(cells_2d, [=] (int i, int j, int) {
        ee(i, j, k) = ee(i, j, k + 1) + hh(i, j, k);   // stride 8 bytes in i
      });
```

The arithmetic is identical and so is its order, so the two produce
bit-identical answers. Only the traversal differs.

## Why the CPU prefers k outermost

In shape A the innermost loop steps by 20 KB per iteration. Consecutive
iterations are in different cache lines and usually different pages, nothing
prefetches usefully, and -- decisively -- the loop **cannot be vectorised**,
because vectorising means doing four consecutive iterations at once and those
four addresses are 80 KB apart.

In shape B the innermost loop steps by 8 bytes. Eight iterations share one
cache line, the hardware prefetcher recognises the stream, and four iterations
fit in one AVX2 register. The `k` dependence has been lifted out to a loop that
runs `nk` times and does no work of its own.

This is not a theoretical argument. Disassembling `PressureForce::calculate`
from two builds that differ only in this rewrite, both icpx 2025.2 at `-O2`:

```
  shape A (today)                  shape B (plane sweeps)

  vmovsd (%r9,%r14,1),%xmm0        vmovupd (%r11,%r15,8),%ymm1
  vaddsd (%r15,%r12,1),%xmm0,%xmm0 vaddpd  (%rbx,%r15,8),%ymm1,%ymm1
  vmovsd %xmm0,(%r12,%r14,1)       vmovupd %ymm1,(%r14,%r15,8)
```

`sd` is scalar double, one element per instruction; `pd` on `ymm` is packed
double, four. Shape A is entirely scalar. Shape B is 4-wide with a scalar
remainder loop, and the Montgomery potential's `g_prime[k] * e` picks up a
`vbroadcastsd` + `vmulpd` pair.

## Why the GPU prefers k innermost

A GPU wants two things a CPU does not care about.

**Enough independent work.** A kernel launch has to fill tens of thousands of
threads. Shape B launches one kernel per interface, each with `nx*ny` threads
and a barrier between them -- for `44x40` that is 1760 threads per launch and
`nk` launches, which is mostly launch latency. Shape A launches once with
`nx*ny` threads, each doing `nk` work: one launch, enough parallelism, and the
sequential dependence lives inside a thread where it costs nothing.

**Coalescing.** Adjacent threads in a warp should touch adjacent addresses.
In shape A the thread index maps to `(i,j)` and `i` is contiguous, so the 32
threads of a warp read 32 consecutive doubles -- fully coalesced -- and then
all step to the next `k` plane together. The stride that hurts a CPU's
*innermost loop* costs a GPU nothing, because on the GPU that stride is walked
by the *whole warp in lockstep*, not by successive iterations of one scalar
loop.

So the rule is not "one order is better". It is:

> The CPU wants the contiguous index in the innermost **loop**.
> The GPU wants the contiguous index across the **threads**, and is happy for
> the loop to walk the strided index.

For a kernel with no vertical dependence both can be satisfied at once -- an
AMReX `ParallelFor` over a 3D box gives `k` outermost and `i` innermost on the
CPU, and maps `i` to the fast thread index on the GPU, from the same source.
The conflict appears exactly when there is a **k-recursion**, because then `k`
cannot be the parallel index and someone has to walk it serially.

## Where each model stands

Verified by reading the source, not inferred:

| routine | MOM6 (dev/ncar, what we benchmark) | protoMOMxx |
|---|---|---|
| pressure force (Montgomery) | `do j; do k=nz,1,-1; do i` -- shape B | shape A |
| vertical viscosity, tridiagonal | `do j; do I; ... do k=2,nz` -- shape A | shape A |
| continuity, Coriolis, hor. visc. | `do k; do j; do i` (no k-recursion) | `ParallelFor(3d box)`, same nest |

The important row is the second one: **MOM6's tridiagonal solve is already
column-parallel with `k` innermost**, the same shape protoMOMxx uses. MOM6 has
not chosen shape B everywhere; it has chosen it where the recursion is cheap to
lift out (the pressure force) and shape A where it is not (the tridiagonal).

## What it costs, measured

From the CPU size sweep (`turbo-prof`, intel 2025.2 at `-O2` on both sides,
`ocean.stats` byte-identical at every point), in nanoseconds per gridpoint per
timestep:

| pressure force | 44x40x2 | 440x400x2 | 44x40x64 | 132x120x40 |
|---|---:|---:|---:|---:|
| MOM6 (shape B) | 2.0 | 2.7 | 2.4 | 2.6 |
| protoMOMxx (shape A) | 10.5 | 13.9 | 19.0 | 30.1 |

The gap widens with both grid size and depth, because both make the strided
walk worse: a larger grid makes each `k` step a bigger jump, and a deeper
column makes more of them. At the largest point the routine is 9% of
protoMOMxx's main loop and 0.7% of MOM6's.

Rewriting just this kernel as shape B, with no other change:

| | shape A | shape B | MOM6 |
|---|---:|---:|---:|
| 44x40x2 | 10.22 | **3.23** | 2.0 |
| 44x40x64 | 19.08 | **8.29** | 2.4 |

(A separate 480-step pair of runs, which is why the shape-A column differs by a
few percent from the sweep table above.) 2.3x to 3.2x, and `ocean.stats` is
byte-for-byte unchanged, as the identical arithmetic requires. The residue over MOM6 is not loop order: it is two
`MultiFab` allocations plus `setVal(0.0)` and a `Gpu::DeviceVector` copy of
`g_prime`, all performed **per timestep**, which want hoisting into the object.

## What loop order does *not* explain

protoMOMxx's advantage over MOM6 shrinks from 1.81x at `NK=2` to 1.18x at
`NK=64`, and the routine driving that is the vertical viscosity. It is
tempting to blame shape A. That is wrong, and the table above says why: **both
models use shape A there.** Per gridpoint per step,

| vertical viscosity | NK=2 | NK=64 | 132x120x40 |
|---|---:|---:|---:|
| MOM6 | 265.7 | 183.6 | 205.3 |
| protoMOMxx | 141.5 | 167.8 | 206.3 |

MOM6 gets *cheaper* per cell as the column deepens because it carries a large
fixed per-column setup cost that amortises; protoMOMxx starts with much less of
that and has less to gain. They converge to within a few percent. So
protoMOMxx's large advantage at `NK=2` is an advantage in per-column overhead,
not in the layer-by-layer work, and the honest expectation at a production
column depth is a small advantage rather than a large one.

## What MOM6 is doing about this

`marshallward/MOM6` `dev/gpu` (not in `mom-ocean/main`, which has the
preparatory refactor but no offload) resolves the conflict by making the loop
shape a **runtime parameter** instead of a source-level decision. The
continuity control structure gained `niblock`, `njblock` and `nkblock`, where a
block size of zero means "the whole computational domain" -- the GPU setting --
and a small value gives CPU-style blocking. The same treatment was applied to
the horizontal viscosity and to `CorAdCalc`.

Its vertical viscosity, meanwhile, is unchanged in shape: `!$omp target teams
loop collapse(2)` wrapped around the `do j; do I` nest that was already there,
with a private `c1(k)`. That is protoMOMxx's kernel written in Fortran.

Two incidental confirmations from that branch: it needed `!NVF$ INLINE` (later
Intel's `forceinline`) on `flux_elem` and `ratio_max`, and `bind(parallel,
teams)` to work around an nvhpc memory error -- the same inlining wall
`MOM_kernel_inline.h` exists to get around, reached from the Fortran side.

## What this implies here

1. **The pressure force should be shape B.** It is a measured 2.3-3.2x on the
   CPU for a rewrite that cannot change the answer. It costs GPU-friendliness,
   which is the next point.
2. **Loop shape should be a policy, not a literal.** MOM6 can pick its
   traversal at run time; protoMOMxx cannot, because whether the `k` loop sits
   inside or outside the lambda is written into each kernel. AMReX supplies the
   equivalent levers -- `MFIter` tiling, `max_grid_size`, and the rank of the
   box handed to `ParallelFor` -- but a kernel has to be *written* so they
   apply. Kernels with a vertical recursion should be authored so both
   traversals can be generated from one body.
3. **Only k-recursive kernels need the treatment.** Continuity, Coriolis and
   the horizontal viscosity have no vertical dependence; a `ParallelFor` over a
   3D box already gives the CPU `i`-innermost and the GPU a coalesced thread
   map. Those three are at parity or better and need nothing here.
4. **A CPU-only comparison against `dev/gpu` will read differently.** Where
   MOM6 adopts shape A for offload, its CPU cost should move toward
   protoMOMxx's. The 11x pressure-force gap is measured against a loop that
   MOM6 is in the process of giving up.

## Reproducing the measurement

The sweep and its report generator live in `turbo-prof`:

```bash
qsub /path/to/turbo-prof/scripts/jobs/job-protomomxx-sweep.sh
python3 /path/to/turbo-prof/scripts/generators/gen_protomomxx_sweep_report.py \
    --run-dir <run_dir> --stack protoMOMxx=/path/to/protoMOMxx
```

See `turbo-prof/docs/PROTOMOMXX_SWEEP.md` for the two axes and what is held
fixed. The shape-A/shape-B comparison is the sweep run twice against two
builds that differ only in `src/core/MOM_PressureForce.cpp`; the check that it
is the same computation is that `ocean.stats` does not move.
