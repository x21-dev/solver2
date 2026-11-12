# solver2

The next version of solver (`../solver`), the deterministic black-box
optimisation core extracted from pleaser: CMA-ES stages and a
quality-diversity archive over a problem interface. solver2 adds a problem
as a class that need not be known at compile time, and a C interface and
Python module over it. It started as solver's code at commit `ae68dbb`,
with the package and modules renamed to `solver2`.

**Read `docs/work-plan.md` before starting**: it is the plan and the
hand-over (design, steps, where the work stands, open questions for Dev).
Update it in the same change when a step is done, a decision is made or
the plan changes. See README.md for a summary and docs/manual.md for the
developer manual, which is still solver's until the plan's step 7.

## Working here

- **Never change `../solver`, `../pleaser` or `../satisfaction` from
  here.** The projects that use solver stay on it, untouched; that is why
  solver2 is a separate project. To build or run solver for a comparison,
  export it elsewhere (`git -C ../solver archive <commit> | tar -x -C
  <dir>`) rather than building in its checkout.
- LDC only, developed with LDC 1.42.0. In cloud sessions the SessionStart
  hook installs it. solver2 is a dub source library with no dependencies.
- Before every commit: `dub test --compiler=ldc2` in solver2, and in
  `docs/example` both `dub test --compiler=ldc2` and
  `dub run --compiler=ldc2 --build=release` (the example checks a resumed
  run against the uninterrupted one bit for bit). The work plan adds the
  binding's and the Python tests as those steps land.
- The compile-time route must keep solver's bits: the worked example
  prints the same output as solver's at `ae68dbb` on one platform (but
  for a change the plan lists under "Changes a migrating project will
  meet"). Check it when a change touches the searches.
- Commit only when Dev asks. The repository has no remote and no CI yet
  (`.github/workflows/ci.yml` is solver's, waiting for one).
- `docs/example/source/checkpoint.d` is printed whole in the manual (To a
  file) as the checkpoint template. Change both together: a
  `static assert` in `docs/example/source/app.d` fails the build when
  they differ.
- The manual's Gaps and caveats list what the library doesn't do. Keep
  it, and its date, true when a change closes or opens one.

## Determinism

- Results are bit-identical on any thread count and CPU of one platform.
  Seeded noise only through `solver2.rng.SplitMix64`, never `std.random`.
- Runs differ between platforms, and that is accepted: LDC lowers `exp`,
  `log`, `pow`, `sin`, `cos` to the C maths library, whose last bits differ
  between glibc and the MSVC runtime. Keep golden runs and reference
  outputs per platform. A test that needs something to happen during a run
  (an explorer aimed, a restart) must force it, so it happens on Linux and
  Windows alike.
- `solver2.detmath` gives the same bits everywhere; its unittests hold it
  to reference bit patterns. Don't add flags that change floating-point
  results (`-fp-contract=fast`, `@fastmath`).

## D style (target LDC)

- Write float selects as `cond ? a : b`, not multiply-by-bool.
- Prefer a functional style wherever it costs neither performance nor
  readability: pure functions (marked `pure`), `const`/`immutable` by
  default, a struct changed only once everything it needs has been
  computed and checked. Not at the cost of vectorisation or allocations in
  hot loops: there keep indexed loops and array operations (`zip`/`lockstep`
  usually stop LDC vectorising).
- Aim for loops LDC can vectorise, and check the asm when it matters. A
  hot loop that writes through one slice and reads a struct's other
  slices: copy them into locals first, or LLVM reloads their pointers every
  iteration and leaves the loop scalar.
- Branchless binary searches need `-x86-cmov-converter=false`.
- No `@fastmath`. Vectorise reductions by splitting into lanes under strict
  IEEE: `enum lanes = 4;` by default, `static foreach` over independent
  accumulators and a fixed-order combine at the end. Floating-point sums
  keep their fixed order: never `std.algorithm.sum` (pairwise; changes the
  bits).
- Biquads: double-precision coefficients and state.
- Keep large values off the stack. With `N > 0` an `Archive` holds its
  cells by value (pleaser's 75-cell one is about 40 KB, and pleaser links
  Windows with an 8 MB stack); Windows' default main-thread stack is 1 MB,
  and an exit with no output and code -1073741571 (0xC00000FD) is a stack
  overflow.
