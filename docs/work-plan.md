# solver2: work plan

This is the plan for solver2 and the hand-over for the sessions that carry
it out. Read it all before starting, and keep it true: when a step is
done, a decision is made or the plan changes, update it in the same
change.

Written 2026-11-12.

## Where the work stands

| Step | State |
|---|---|
| 0. Set up | done 2026-11-12 |
| 1. Spike: Python calls D | done on Windows 2026-11-12: it works; Linux owed |
| 2. The class, completed | next |
| 3. Entry points and run-time `minScore` | |
| 4. The C interface | |
| 5. The Python module | |
| 6. CI | waits for Dev to create the GitHub repository; Linux runs of steps 1 to 5 are owed to it |
| 7. The manual | |

solver2 is solver's code with the package and modules renamed, and
nothing else changed: its tests pass and the worked example prints what
solver's does (step 0 has the details).

## What solver2 is

solver2 is the next version of [solver](https://github.com/x21-dev/solver)
(`../solver`), the deterministic black-box optimiser extracted from
pleaser: CMA-ES stages and a quality-diversity archive (CMA-MAE) over a
problem interface, with a learned surrogate model.

It adds two things solver doesn't have:

1. **A problem as a class that need not be known at compile time.** In
   solver a problem is a compile-time concept: the searches are templates
   over the problem's type. solver2 makes the class route complete, so a
   problem chosen, built or loaded at run time gets everything a
   compile-time one does.
2. **Python.** A C interface over that class, built as a shared library
   (a DLL on Windows, a `.so` on Linux), and a Python module that loads
   it, so a problem written in Python runs through the D searches.

## Why a separate project

Decided by Dev on 2026-11-12:

- Several projects use solver as it is today: `../pleaser`,
  `../pleaser-bus-chain` and `../pleaser-console-tilt` name it in their
  `dub.json`. They must not need rewriting or rebuilding for this work.
- So solver stays as it is, and the work happens here, in `../solver2`,
  starting from a copy of solver's code.
- No project has to move. One that wants to follows the manual's
  migration section (step 7).

**Never change `../solver`, `../pleaser` or `../satisfaction` from here.**

## Starting point

solver `main` at commit `ae68dbba0e9f04e78bb778fdacd8266f3c983474` (the
merge of solver#17), developed with LDC 1.42.0. The machine this was
written on has LDC 1.42.0, dub 1.41.0 and Python 3.15.0, on Windows 11.

Before starting, read in `../solver`: `CLAUDE.md`, `README.md`,
`docs/manual.md` and `source/solver/problem.d`. The line numbers below are
solver's at that commit.

## What stops the class route in solver

solver already has `DynamicProblem` (`source/solver/problem.d:461`), an
abstract class that satisfies the problem concept, with `Outcome` as its
evaluation type. The searches instantiate once over it. It falls short in
five ways:

1. **Optional members are looked up on the static type.** `screenPenalty`
   (`stage.d:224`, `archive.d:687`), `noFeasibleSeed` (`archive.d:775`),
   `adoptCauses` (`archive.d:905`) and the cause counts (`hasCauses`,
   `stage.d:68`) are found with `__traits(hasMember, P, ...)`. With
   `P = DynamicProblem`, a subclass that defines one is silently ignored.
2. **No learned model.** `DynamicProblem` has none of the surrogate
   members (`featureCount`, `features`, `targetCount`, `targets`,
   `predicted`), so `LearnedModel!DynamicProblem` fails `isLearning`
   (`model.d:44`). `Outcome.targets` is declared and nothing reads it.
3. **`minScore` is a template argument** of `Archive`, `ArchiveRun` and
   `Archiving` (`archive.d:127`, `413`, `699`). A problem met at run time
   can't say where the thresholds start.
4. **The caller has to know the problem.** The grid's axes, the ranking's
   offsets and the start genes are written by the caller
   (`docs/example/source/app.d:140`, `203`). A caller handed an unknown
   class has nowhere to get them.
5. **Nothing checks at run time what `isProblem` checks at compile time**,
   and the class has never been through a search: its only test saves and
   loads an evaluation (`problem.d:669`). `descriptors` copies with
   `d[] = e.descriptors[]`, whatever the lengths.

## Design

### The searches stay templates; the class is one problem among them

The compile-time route is kept, unchanged in behaviour: a struct with its
own `Eval` type still runs with direct, inlined calls. The class is a `P`
like any other, and one instantiation of each search serves every
subclass.

Why not a class-only interface:

- a project's typed evaluation (pleaser's `Evaluation` in its elites, its
  `Best` and every callback) would have to be erased, which is a rewrite
  for anyone migrating;
- the checkpoint template walks the evaluation's fields at compile time;
- an archive held by value (`N > 0`) needs the type.

The cost of the class route is a few virtual calls and three small arrays
per candidate: nothing next to a real evaluation, noticeable on a
microsecond-scale problem.

### The class, completed

- **Defaults for the optional members**, so the static look-ups always
  find them and a default behaves as "absent" does today:
  `screenPenalty` returns 0, `noFeasibleSeed` returns null, `adoptCauses`
  does nothing.
- **Cause counts**: a `Causes` type with one count per `ConstraintSpec`,
  filled from `Outcome.violations`.
- **Surrogate members** with defaults: the features are the genes; the
  targets are the objective (NaN when not scored), the hard plus terminal
  breach, and the descriptors; `predicted` reads that layout back with a
  headroom. A subclass may override them. `features` and `predicted` are
  called from several threads at once.
- **Self-description**, so a caller needs to know nothing else: the grid's
  axes, the ranking, `minScore`, and start genes from `GeneSpec.neutral`.
- **Run-time checks** in the final `evaluate`: array lengths, breaches
  not negative and not NaN. And a `checkProblem(p)` to run on a problem
  when it is loaded: rows evaluated alone and in a batch, and twice, give
  the same bits; genes outside [0, 1] are accepted.

### Entry points without template arguments

A module (working name `dynamic`) with the stage and the archive over the
class and `BinnedGrid`, taking the grid, ranking and `minScore` from the
object. These are plain functions, so they can be exported from a
compiled library.

`minScore` becomes a run-time field of `Archive`, initialised from the
template argument, so compile-time users see no change. `el.value -
minScore` gives the same bits read from a field as from a constant (no
`@fastmath`); the equivalence test below must show it.

### The C interface

A dub package of its own in this repository (working name `bindings/c`,
`targetType: dynamicLibrary`) that depends on the core by path. The core
stays a source library.

- **The problem as data.** The caller passes genes, constraints, axes,
  ranking and `minScore` once; D builds a `ForeignProblem` (a subclass)
  from them.
- **One callback per batch.** D calls the caller's `evaluate` with all
  the genes as one flat array; the caller fills flat arrays for scored,
  objective, violations and descriptors. D sums each breach with
  `breachOf`, so the order and the bits are fixed.
- **Opaque handles and copies.** The caller holds a handle to a problem,
  a run or a model, and D keeps the object alive. Results (elites,
  thresholds, progress) are copied into buffers the caller owns. The
  caller never holds memory of D's garbage collector.
- **No exception crosses the boundary.** Every exported function catches
  everything and returns an error code, with the message from a
  `last error` call. A callback reports failure by its return value, and
  the run stops.
- **Saving**: a run and its model exported and imported as bytes.

What D needs inside another process:

- an exported `init` that starts the D runtime before any other call,
  and a `term`;
- every thread that calls in known to the D runtime: the one that called
  `init` is, and any other must attach first (`thread_attachThis`). The
  spike (step 1) showed this working, and that the exported functions
  can do it themselves;
- the model's worker threads never call back: for a foreign problem,
  features and predictions use the D-side defaults;
- the archive on the heap (`N = 0`), as the class route has it anyway. A
  stack overflow on Windows shows as an exit with no output and code
  -1073741571.

### The Python module

Pure Python over the shared library: `ctypes` and numpy, no compiled
extension, so it isn't tied to a Python version. A problem is a Python
class with its description and an `evaluate(genes)` that returns arrays.
The wrapper catches a Python exception in the callback, stops the run and
raises it again from the call that started the run (`ctypes` on its own
prints it and carries on).

Not used: `pyd` and `autowrap`, which build Python extension modules from
D. They tie the build to Python versions and add a dependency. Their
current state wasn't checked.

### Determinism

As in solver: bit-identical on any thread count and CPU of one platform,
different between platforms (the C maths library). The D side gives the
same bits for the same outcomes. Whether a Python `evaluate` is
deterministic is its author's job: numpy's threaded maths can change the
last bits with the thread count. The manual must say so.

## Steps

One branch and one pull request per step once the repository has a
remote; until then, one commit per step. Each step lists what makes it
done.

### 0. Set up (done 2026-11-12)

- solver's tree at the commit above, copied with `git archive`: every
  tracked file but `docs/next-session.md` (solver's hand-over; this plan
  is solver2's).
- The package is `solver2` and the modules `solver2.*`, in
  `source/solver2/` (Dev's choice, 2026-11-12). The worked example is
  `solver2-example`, and its checkpoint's magic string and temporary file
  are named so too, so that it never reads or overwrites solver's.
- `CLAUDE.md` and `README.md` are solver2's. The manual is solver's with
  the module names changed and a note at its top saying so; step 7
  rewrites it.
- `.github/workflows/ci.yml` and the SessionStart hook are solver's; the
  workflow runs nothing until there is a remote.
- Checked on Windows 11 with LDC 1.42.0: `dub test --compiler=ldc2`
  passes (9 modules), and in `docs/example` `dub test --compiler=ldc2`
  (10 modules) and `dub run --compiler=ldc2 --build=release`. The
  example's 12 lines of output are byte-identical to those of solver's
  example at `ae68dbb`, built from an export of that commit outside
  `../solver`. Not run on Linux.

To repeat the comparison after a change to the searches: export solver
(`git -C ../solver archive ae68dbb | tar -x -C <dir>`), run
`dub run --compiler=ldc2 --build=release` in both `docs/example` folders,
and compare what the programs print (dub's own lines name the package and
differ).

### 1. Spike: Python calls D (done on Windows 2026-11-12)

The risky parts, tested before anything is designed around them. The
verdict: **it works**, on Windows 11 with LDC 1.42.0 and Python 3.15.0
(64-bit). Its Linux run waits for a Linux machine (step 6, "Linux before
the remote").

The code is under `spike/`, throwaway, kept as the record:

- `source/spike.d`: a shared library that runs `cmaStage` over a problem
  whose `evaluate` is a C callback, plain or pre-screened by a
  `LearnedModel` on its own `TaskPool`.
- `source/native.d`: 50 plain stages and 5 pre-screened ones with the
  callback in D, printing each stage's bits (a hash of every generation's
  numbers and genes, the best and its genes).
- `spike.py`: the same 55 stages with the callback in Python through
  `ctypes`, then the checks below.

To run it, in `spike/`:

```
dub build --compiler=ldc2 --build=release --config=native
dub build --compiler=ldc2 --build=release --config=library
./spike-native > native.out
python spike.py > python.out      # then compare the two files
python spike.py --checks          # also --stress, --probe, --numpy, --unattached
```

#### What it showed

- **The same bits.** The 55 stages through the Python callback print
  what the native ones print, byte for byte, the pre-screened ones with
  three worker threads included. D's collector ran 80 times on the way.
- **The runtime must be linked statically.** dub builds a
  `dynamicLibrary` with `-fvisibility=public -dllimport=all`, which makes
  the DLL need LDC's `druntime-ldc-shared.dll` and
  `phobos2-ldc-shared.dll`. Python doesn't search `PATH` for a DLL's
  dependencies, so loading fails ("Could not find module ... or one of
  its dependencies"). With `"dflags": ["-fvisibility=hidden",
  "-dllimport=none", "-link-defaultlib-shared=false"]` the DLL holds the
  runtime (919 KB), exports only the functions marked `export`, and
  loads. On Linux LDC also links the shared runtime into a shared library
  by default, so the same flag is expected to be needed there (untested).
- **The runtime doesn't start by itself.** A call before `init`
  (`rt_init`) is an access violation. No `DllMain` is needed.
- **Errors stay on their side.** A D exception is caught at the exported
  function and returned as a code with its message. A Python exception in
  the callback is kept by the wrapper, the callback returns non-zero, D
  throws and unwinds, and the wrapper raises the Python exception again;
  the next stage runs as before.
- **Other threads.** A stage run from a second Python thread after
  `thread_attachThis` gives the same bits. The D runtime doesn't know a
  thread that hasn't attached (`Thread.getThis()` is null). Unattached,
  51 stages still gave the same bits, but a collection started there
  reclaimed nothing (192 MB in use where an attached thread had 64 MB).
  So: attach, as druntime asks; unattached is not shown to be safe.
- **A busy Python thread slows the callbacks.** `ctypes` releases the
  interpreter lock during the D call and takes it again for each
  callback. Beside a Python thread that never sleeps, each callback
  waited about 5 ms for it (2600 callbacks in 14 s); the bits were the
  same.
- **The callback's cost**: about 17 µs for a batch of 8 candidates of 6
  genes, the Python arithmetic included (55 stages, 25,500 callbacks:
  2.68 s against 2.25 s native).
- **numpy views work** (`np.ctypeslib.as_array` on D's arrays, nothing
  copied), but `ndarray.sum` doesn't add in D's order: 49 of the 55
  stages then differ from the native ones. A numpy problem is
  deterministic by itself; it only can't be compared bit for bit with a
  D reference unless it sums in the same order.
- **A D bug in a `release` build is not caught.** An index one past the
  end wrote out of bounds and returned success; a failed `assert` did
  nothing (asserts are compiled out); a wild index and a null class
  reference raised a Windows exception, which `ctypes` turns into an
  `OSError`, with D's state left as it was. With `-boundscheck=on` both
  bad indexes come back as errors ("index [3] is out of bounds for array
  of length 3"), the bits are the same, and the 55 native stages take
  2.84 s instead of 2.25 s (26% more, on a problem whose evaluation
  costs nothing).

#### What it changes in the design

- **Build the binding with the three flags above**, on both platforms.
- **Bounds checks on in the binding, and only there** (`-boundscheck=on`;
  Dev's decision, 2026-11-12): a real problem's evaluation dwarfs the
  26%, and a bug then reaches Python as an error instead of corrupting
  the process. The core adds no flag: a D project builds it with its own.
  The spike's library is built so (`spike/dub.json`), with the same
  bits. Asserts stay off (`-enable-asserts` would keep them; not asked
  for).
- **How the flag can be switched** (checked with dub 1.41.0): a
  package's own `dflags` come last on LDC's command line and win, so
  neither the `DFLAGS` environment variable nor `--build=release-nobounds`
  turns off a `-boundscheck=on` written in the binding's `dub.json`. The
  binding therefore gets a second configuration without the flag
  (working name `unchecked`), and the manual says how to build either. A
  project without the flag in its `dub.json` is switched by `DFLAGS` or
  the build type.
- **Every exported function attaches the calling thread when the runtime
  doesn't know it, and detaches it before returning**, so a caller on any
  thread is safe and no thread ends while still attached. The rule "one
  thread only" is then not needed; "one call at a time per object" still
  is.
- **The callback's status and the wrapper's kept exception**, as the
  spike has them, are the error path for step 4 and step 5.
- **The manual** (step 7) gets: the 5 ms wait beside a busy Python
  thread, numpy's summation order, and the build flags.

Not tried in the spike: the archive (the same callback path as a stage),
saving and loading, Linux.

### 2. The class, completed

Everything under "The class, completed", with tests that run a stage, an
archive, a pre-screened stage and a resumed run over a subclass used only
through the base class. Done when those pass and a subclass's
`screenPenalty` and cause counts are shown to take effect.

### 3. Entry points and run-time `minScore`

The `dynamic` module and the `Archive` field. Done when the class runs
with no template argument at the call site, and the equivalence test
passes: the worked example's compile-time problem gives the same bytes as
in solver.

### 4. The C interface

The `bindings/c` package, with a C header written by hand and a test in D
that drives the library through the exported functions only. Built with
the runtime linked statically and bounds checks on, and with a second
configuration without the checks (step 1, "What it changes in the
design"). Done when a
stage, an archive, a learned model, and a save, load and resume all run
through it, and a resumed run equals the uninterrupted one bit for bit.

### 5. The Python module

The module, its tests, and a Python version of the worked example
(`Beam`), whose result equals the D example's bit for bit on one
platform. Done when the tests pass against the library built in step 4.

### 6. CI

Dev creates the GitHub repository once solver2 is mature enough
(2026-11-12); this step waits for it. Until then every check runs on
Dev's Windows machine only (see "Linux before the remote" below).

`.github/workflows/ci.yml` is solver's and already runs, on
`ubuntu-latest` and `windows-latest` with LDC 1.42.0, the core's tests,
the example's tests and the example. To add:

- **The shared library**: build `bindings/c` and run its D test, which
  drives it through the exported functions only.
- **Python**: `actions/setup-python`, then numpy and the test runner with
  pip, then the Python tests against the library just built.
- **The same bits through Python**: the Python worked example's result
  compared with the D example's, within one job (so on one platform).
- **Python versions**: the module is `ctypes` only, so the version should
  not matter. One job per platform on the newest Python, and one more on
  Linux on the oldest supported, shows it without multiplying the matrix.
- **The equivalence with solver** (the example's output against solver's
  at `ae68dbb`) needs solver's code in the job. Keep it a local check
  unless Dev wants the workflow to fetch solver.

Keep the matrix small: on a private repository Windows minutes cost
more than Linux's and macOS minutes far more (check GitHub's current
rates when the repository exists). macOS is not planned;
runs differ between platforms anyway, so adding it would mean a third
set of reference outputs.

#### Linux before the remote

The machine has no WSL and no Docker (checked 2026-11-12), so nothing
here has run on Linux, and nothing will until the repository exists (CI,
or a cloud session, whose SessionStart hook installs LDC) or Dev offers
another Linux machine. Steps 1 to 5 are therefore done on Windows first,
and each one's Linux run is owed to step 6. Write Linux-specific parts
(the `.so` name, how the library is found, `-fPIC`) as the documentation
says, and mark them untested in this plan until they have run.

### 7. The manual

`docs/manual.md`, brought to solver2's state throughout, with its date and
Gaps and caveats true. Two new sections:

**Interfacing with Python**

- building the shared library, and where the Python module finds it;
- writing a problem in Python: genes, limits, descriptors, `evaluate`;
- running a stage and an archive; reading the results; progress and
  stopping;
- the learned model from Python;
- saving and resuming;
- what to know: threads, exceptions, determinism of a Python
  `evaluate` (numpy's summation order), cost per batch, the wait beside a
  busy Python thread;
- building the library with or without bounds checks: completes the
  manual's "Bounds checks" (under Performance), which already covers a D
  project and was written on 2026-11-12;
- the C interface itself, for other languages.

**Migrating from solver**

- who needs to: nobody; solver stays as it is;
- what is the same: the problem concept and its members, the options, the
  callbacks;
- what changes: the package and module names, `minScore`, and whatever
  else the steps above change, kept as a list while working (below);
- the same bits: a migrated project gives byte-identical runs on one
  platform, and how to prove it with a golden run (the same run before
  and after, files compared);
- saved states: whether a checkpoint written with solver reads in
  solver2;
- moving a compile-time problem to the class, for a project that wants
  the run-time route.

Also the `README.md`.

## Changes a migrating project will meet

Keep this list as the work goes, for the manual's migration section.

- **The names** (step 0): the dub package is `solver2` and every module
  `solver2.*`. A project changes its dependency (`"solver2": { "path":
  "../solver2" }`) and its imports (`import solver.problem;` becomes
  `import solver2.problem;`, and so for `stage`, `archive`, `model`,
  `ensemble`, `cmaes`, `bytes`, `rng`, `detmath`). Nothing else differs
  yet, and the bits are the same.

## Open questions for Dev

1. **Names.** Decided 2026-11-12: package `solver2`, modules `solver2.*`.
   A build can never pick up the wrong library, and migration is a search
   and replace of the imports.
2. **Checkpoints.** solver has no checkpoint API by decision: the save
   format is the manual's template, copied by each project. Python can't
   walk D structs, so the binding has to export a run as bytes.
   Recommended: the binding copies the template, as any project does, and
   the core still gets no serialiser.
3. **Python packaging.** A plain folder whose module loads the library by
   path (recommended for now), or a wheel.
4. **solver from here on.** Whether a bug found in code both hold
   (`cmaes.d`, `ensemble.d`) is fixed in both, or solver is frozen.
5. **A remote and CI.** Decided 2026-11-12: Dev creates the GitHub
   repository once solver2 is mature enough, and step 6 then extends
   solver's workflow. Don't create a repository or push from a session.
6. **Python owning the loop** (ask for candidates, evaluate, tell the
   results, as pycma and pyribs do). Possible later without rewriting the
   searches, by running them in a D fiber that pauses at each `evaluate`.
   Not planned.
7. **Bounds checks in the binding.** Decided 2026-11-12: on in the
   Python binding only, and the manual says how to turn them on or off in
   any case, for the binding and for a D project. The core, as a source
   library, takes whatever flags its user builds with.

## How to work

- LDC only, 1.42.0. Checks before every commit: `dub test --compiler=ldc2`
  in the core, and in `docs/example` both `dub test --compiler=ldc2` and
  `dub run --compiler=ldc2 --build=release`; from step 4 the binding's
  tests, from step 5 the Python tests.
- Determinism and D style: `CLAUDE.md`, as in solver. In
  short: seeded noise only through `SplitMix64`; no `@fastmath` or
  `-fp-contract=fast`; fixed-order sums; goldens per platform; a test
  that needs an event during a run forces it; large values off the stack.
- `docs/example/source/checkpoint.d` is printed whole in the manual; a
  `static assert` fails the build when they differ. Change both together.
- Commit only when Dev asks.
