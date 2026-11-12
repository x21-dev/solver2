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
| 1. Spike: Python calls D | next |
| 2. The class, completed | |
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
- calls from one thread, the one that called `init`. D's collector only
  knows threads its runtime has seen; another thread must attach first
  (`thread_attachThis`);
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

### 1. Spike: Python calls D

The risky parts, tested before anything is designed around them. Throwaway
code under `spike/`.

- A shared library with `init` and one function that runs `cmaStage` over
  a problem whose `evaluate` is a C callback.
- A Python script that loads it with `ctypes`, supplies the callback and
  reads the result.
- The same problem run natively in D, and the two results compared bit
  for bit.
- Done when it passes on Windows, including a run long enough for the
  collector to run many times, and the findings are written into this
  plan (build flags, how the runtime starts, anything that surprised).
  Its Linux run waits for a Linux machine (step 6, "Linux before the
  remote"): the spike's verdict covers Windows only until then.
- If the D runtime can't be made to behave inside a Python process, stop
  and bring it to Dev. The fallback is a D program driven over a pipe.

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
that drives the library through the exported functions only. Done when a
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
- what to know: one thread, exceptions, determinism of a Python
  `evaluate`, cost per batch;
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
