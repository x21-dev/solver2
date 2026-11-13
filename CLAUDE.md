# solver2

The next version of solver (`../solver`), the deterministic black-box
optimisation core extracted from pleaser: CMA-ES stages and a
quality-diversity archive over a problem interface. solver2 adds a problem
as a class that need not be known at compile time, and a C interface and
Python module over it. It started as solver's code at commit `ae68dbb`,
with the package and modules renamed to `solver2`.

See README.md for a summary and docs/manual.md for the developer manual,
which documents all three routes (a compile-time problem, the class,
Python and C), the migration from solver, and, under Gaps and caveats,
what the library doesn't do. The work that built solver2 is finished and
released (v0.1.0, 2026-11-13); this file and the manual are the
hand-over.

## Decided by Dev

- The core gets no serialiser: a D project copies the manual's
  checkpoint template, and the C interface has its own copy
  (`bindings/c/source/solver2c/checkpoint.d`).
- The Python module is a folder inside a release archive, not a wheel,
  until someone needs pip.
- Python doesn't own the loop (no ask and tell); not planned.
- Bounds checks are on in the shared library only; the core takes its
  user's flags.
- The Python example is held to `docs/example-class`, not to
  `docs/example`: a problem behind the C interface has the class's
  learning spec, and the interface isn't shaped around one example.
- macOS and arm64 are not planned.

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
  run against the uninterrupted one bit for bit). Then the C interface:
  in `bindings/c`, `dub test --compiler=ldc2` and
  `dub build --compiler=ldc2 --build=release`, and in `bindings/c/test`
  `dub run --compiler=ldc2` (the driver loads the library just built and
  must end with `all passed`). Then Python, in `bindings/python`:
  `python -m unittest discover -s tests` (it loads the same library, and
  builds and runs `docs/example-class` to compare the Python example
  with it).
- `docs/example-class/source/app.d` and
  `bindings/python/examples/beam.py` are one program twice, in D and in
  Python, and must print the same lines. Change them together.
- The C interface's header (`bindings/c/include/solver2.h`), its D
  mirror (`source/solver2c/c.d`, `api.d`) and its Python one
  (`bindings/python/solver2/_c.py`) change together: the binding's tests
  read the header with ImportC and fail to build when the D side
  differs, and the Python tests parse it and fail when `_c.py` does. A
  change of a struct or a function there is a change of `SOLVER2_ABI`,
  in all three.
- In the binding, a thread unknown to the D runtime is attached once, on
  its first call, and detached by a hook when it ends. Don't attach and
  detach per call: with this druntime that corrupts memory (the comment
  above `enter` in `bindings/c/source/solver2c/api.d` says how).
- The shared library links the D runtime statically on Windows and uses
  LDC's shared runtime on Linux (`bindings/c/dub.json`). Never link it
  statically on Linux: the collector then scans the wrong thread-local
  storage and crashes about once in 50 runs of the driver (the manual's
  "The C interface, for other languages" says why). CI checks it with
  `ldd`.
- A crash that shows on Linux only can be hunted on CI: a temporary
  branch with a workflow of its own (triggered by a push to that
  branch) that builds with `DFLAGS=-g`, turns core dumps on
  (`ulimit -c unlimited`, `sudo sysctl -w kernel.core_pattern=<file>`),
  runs the failing program a few hundred times in several jobs, and on
  the first failure prints what it printed and
  `gdb -batch -ex "thread apply all bt" <program> <core>`. Ask before
  pushing it, and delete it afterwards.
- The compile-time route must keep solver's bits. The worked example's
  stage prints what solver's at `ae68dbb` prints on one platform; its
  archive's lines differ only because its `minScore` is −50 (with
  `MIN_SCORE` set to 0 in a copy it prints all of solver's lines). Check
  both when a change touches the searches: export solver
  (`git -C ../solver archive ae68dbb | tar -x -C <dir>`) and a copy of
  solver2 (`git archive HEAD | tar -x -C <dir>`, with `MIN_SCORE` 0),
  run `dub run --compiler=ldc2 --build=release` in both `docs/example`
  folders, and compare what the programs print (dub's own lines name
  the package and differ). A deliberate change of bits goes into the
  manual's "Migrating from solver", under "What changes".
- solver is frozen by default: a bug in code both hold is fixed here and
  listed there, and Dev decides for each whether solver gets the fix, in
  a session of solver's own.
- The manual documents the Python module and the C interface
  ("Interfacing with Python"): a change of `bindings/` that a user can
  see changes that section too, its C program included (it was compiled
  and run by hand, with `ldc2 demo.c solver2.lib`).
- Commit only when Dev asks, and ask before pushing. The repository is
  github.com/x21-dev/solver2 (public): a change is a branch and a pull
  request, and merging is Dev's. CI (`.github/workflows/ci.yml`) runs
  every check above on Linux and Windows for each push to `main` and
  each pull request, and once more on Linux on the oldest Python tested.
- Releases: `.github/workflows/release.yml` builds an archive per
  platform with `release/package.py` (the library, what it needs beside
  it, the Python module, the header, the example, `release/README.md`
  and `LICENSE.md`), tests each where there is no LDC with
  `release/check.py` and the driver, and for a tag `v<version>` makes a
  **draft** release. It also runs, without the release, on every pull
  request. Tagging is Dev's word, and publishing the draft is Dev's.
- The Linux archive is built in an `ubuntu:22.04` container, not on a
  runner's own system: Ubuntu 22.04 stays supported whatever
  `ubuntu-latest` becomes. Don't move the build to a newer system
  without Dev's word, and keep the test on a newer one beside it.
- solver2 is proprietary, with licences granted individually
  (`LICENSE.md`, and `"license": "proprietary"` in the dub files). Don't
  add an open-source licence, and keep `release/README.md` true of what
  the archive holds and needs.
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
