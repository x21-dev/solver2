# solver2: work plan

This is the plan for solver2 and the hand-over for the sessions that carry
it out. Read it all before starting, and keep it true: when a step is
done, a decision is made or the plan changes, update it in the same
change.

Written 2026-11-12. All its steps were done by 2026-11-13: from here it
is the record of how solver2 was built, and the place for what comes
next ("What is left open", below).

## Where the work stands

| Step | State |
|---|---|
| 0. Set up | done 2026-11-12 |
| 1. Spike: Python calls D | done on Windows 2026-11-12: it works (never run on Linux; the driver of step 4 is, in CI) |
| 2. The class, completed | done 2026-11-12 |
| 3. Entry points and run-time `minScore` | done 2026-11-12 |
| 4. The C interface | done 2026-11-12; on Linux since step 6 |
| 5. The Python module | done 2026-11-12; on Linux since step 6 |
| 6. CI | done 2026-11-12: Linux and Windows, every push and pull request |
| 7. The manual | done 2026-11-13 |
| 8. Releases | done 2026-11-13: the workflow drafts a release for a tag; none made yet |

solver2 is solver's code with the package and modules renamed, plus the
completed `DynamicProblem` (step 2), the archive's `minScore` as a value
known at run time, the searches over the class as plain functions in
`solver2.dynamic` (step 3), their C interface as a shared library in
`bindings/c` (step 4), and the Python module over it in
`bindings/python` (step 5), and a manual that documents all of it
(step 7). The searches otherwise are unchanged: the worked example's
stage prints what solver's does, and so did its archive, to the last bit,
until step 7 gave the example a `minScore` of −50 (with 0 it still
does).

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
  axes, the ranking, `minScore` (with step 3), and start genes from
  `GeneSpec.neutral`.
- **Run-time checks** in the final `evaluate`: array lengths, breaches
  not negative and not NaN. And a `checkProblem(p)` to run on a problem
  when it is loaded: rows evaluated alone and in a batch, and twice, give
  the same bits; genes outside [0, 1] are accepted.

### Entry points without template arguments

The module `solver2.dynamic` has the stage and the archive over the class
and `BinnedGrid` (`runStage`, `runArchive`), taking the grid, ranking,
start genes and `minScore` from the object. These are plain functions, so
they can be exported from a compiled library.

`minScore` is a run-time field of `Archive`, initialised from the
template argument, so compile-time users see no change. `el.value -
minScore` gives the same bits read from a field as from a constant (no
`@fastmath`); step 3's tests and its equivalence check show it.

How the run-time value reaches the archive (decided in step 3): a
problem's own `minScore` member is one more optional member the archive
looks for on the static type, as it does `screenPenalty` or
`noFeasibleSeed`. `DynamicProblem` declares it, default 0. The other way,
a parameter of `runArchive`, would have left a subclass's `minScore`
silently unread by anyone who calls the template directly, which is the
first of the five faults above.

### The C interface

A dub package of its own in this repository (`bindings/c`,
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
  exported functions do it themselves: once per thread, on its first
  call, with a hook that detaches the thread when it ends. Not per call:
  step 4 found that unsafe (see there);
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

One branch and one pull request per step now that the repository has a
remote (steps 0 to 5 were one commit each on `main`, before it had).
Each step lists what makes it done.

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

Since step 7 the example's archive starts its thresholds at −50, so only
the stage's five lines equal solver's. To compare the archive too, do it
in a copy of solver2 (`git archive HEAD | tar -x -C <dir>`) with
`MIN_SCORE` set to 0 in `docs/example/source/app.d`: it then prints all
twelve of solver's lines.

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
  (Wrong for Linux: there the static runtime is a fault, and the library
  uses the shared one. See step 7.)
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
  (On Windows only, since step 7.)
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
- **(Superseded in step 4: attaching per call is unsafe; a thread is
  attached once and detached when it ends.)** Every exported function
  attaches the calling thread when the runtime doesn't know it, and
  detaches it before returning, so a caller on any
  thread is safe and no thread ends while still attached. The rule "one
  thread only" is then not needed; "one call at a time per object" still
  is.
- **The callback's status and the wrapper's kept exception**, as the
  spike has them, are the error path for step 4 and step 5.
- **The manual** (step 7) gets: the 5 ms wait beside a busy Python
  thread, numpy's summation order, and the build flags.

Not tried in the spike: the archive (the same callback path as a stage),
saving and loading, Linux.

### 2. The class, completed (done 2026-11-12)

`DynamicProblem` in `source/solver2/problem.d` now has everything under
"The class, completed" but `minScore`, which comes with the code that
reads it (step 3). `stage.d`, `archive.d` and `model.d` are untouched.

- **Defaults**: `screenPenalty` 0, `noFeasibleSeed` null, `adoptCauses`
  sizing the counts.
- **Cause counts**: `CauseCounts`, one count per limit, from
  `Outcome.violations`.
- **The learning spec**: `LearnedModel!DynamicProblem` compiles and runs.
  The default targets are the value (not the objective), the hard and
  the terminal breach, and the descriptors; `predictionHeadroom` is the
  one knob.
- **Self-description**: `axes` and `grid`, `ranking`, `startGenes`;
  `descriptorCount` defaults to `axes.length`.
- **Run-time checks**: `evaluate` checks every batch; `checkProblem`
  checks a problem once.
- The required members are now `abstract`, so a subclass that forgets one
  fails where it is constructed, not at link time.

The tests are in the new `source/solver2/dynamic.d` (which the entry
points of step 3 join): a subclass used only through a `DynamicProblem`
variable, against the same problem as a struct with none of the optional
members.

- A stage and an archive give the struct's bits, with cause counts on
  top.
- A pre-screened stage gives the struct's bits with the default
  `screenPenalty`; a subclass's is called through the base class
  (28 generations × 32 candidates) and changes the picks.
- A subclass's `noFeasibleSeed` hint reaches `NoFeasibleSeed`.
- An archive stopped, saved as plain data and resumed equals the
  uninterrupted one, cause counts included.
- A pre-screened stage and model-guided emitters with
  `LearnedModel!DynamicProblem` run, the same twice.
- `checkProblem` passes a sound problem and names what is wrong with five
  unsound ones; `evaluate` refuses three kinds of malformed outcome.

Checked on Windows: `dub test` (10 modules), the example's tests and its
run, whose output is still solver's; the spike rebuilt against the
changed core gives the same lines.

Found on the way, for step 7's manual and its Gaps and caveats:

- **A prediction is ranked by its objective, ignoring a predicted soft
  breach** (`predictedFitness` in `archive.d`, `predictScore` in
  `stage.d`), while an evaluation is ranked by its value. A learning spec
  for a problem with soft limits should therefore learn the value. The
  class's defaults do; the worked example's `Beam` learns the objective.
- `checkProblem` can only catch a problem that isn't deterministic on the
  rows it tries.

### 3. Entry points and run-time `minScore` (done 2026-11-12)

The class runs with no template argument at the call site, and the
worked example's compile-time problem gives the same bytes as in solver.
`stage.d`, `model.d` and the searches' loops are untouched; `archive.d`
changed only where an archive is made and where a resumed run is checked.

**`minScore` at run time**

- `Archive` has a field `minScore`, last in the struct, which starts as
  the template's third argument (renamed `defaultMinScore`; the name
  isn't visible to a caller). `qdScore` reads the field. With `N` 0 the
  constructor takes it: `Archive!E(cells, minScore)`.
- `hasMinScore!P`: the problem has a member `minScore`. `runArchive` then
  starts the thresholds there (with `N > 0` too, by filling them) and
  refuses one that isn't finite. A resumed run on such a problem must
  have an archive whose `minScore` equals the problem's.
- `Archiving!(P, ..., minScore)` with a non-zero argument on a problem
  that has the member fails to compile: one source only.
- `DynamicProblem.minScore()`, default 0. `checkProblem` wants it finite
  and, for a problem with descriptors, no checked row within the limits
  with a value below it.

**The entry points** (`source/solver2/dynamic.d`)

- `runStage(p, best, x0, o, model, generationDone, run, afterGeneration)`
  and `runArchive(p, seeds, o, progress, resume, afterIteration,
  evaluated, model)`: `cmaStage` and `Archiving.runArchive` with the same
  parameters in the same order, less the grid.
- From the object: the start genes when `x0` or `seeds` is empty, the
  `BinnedGrid` over its axes, the ranking and `minScore`. **The options'
  `ranking` is replaced by the problem's**, so a caller (the C interface)
  never has to carry one; a run records it, and resumes through the same
  entry point.
- They check their arguments and throw with the class's name: no
  problem, a start or seed of another length or with a gene not finite,
  no axes or not one per descriptor, a resumed archive of another cell
  count.
- Aliases for the types, so nothing is instantiated by the caller:
  `DynamicBest`, `DynamicStageGeneration`, `DynamicArchiving`,
  `DynamicRun`, `DynamicArchiveGeneration`, `DynamicBatch`,
  `DynamicModel`.
- Saving is the templates': `StageRun`, `runState`, `restoreRun`, the
  archive and the causes. Nothing was added for it (open question 2).

**The tests**

- `archive.d`: `minScore` as the template's argument, as a constructor's
  argument and by value give the same thresholds and QD score bits; a
  struct problem with the member runs as the same problem with the
  template's argument, with `N` 0 and `N` 16; a resumed run with another
  `minScore` is refused, and so is one not finite.
- `dynamic.d`: `runStage` and `runArchive` give the bits of the templates
  called with the object's start, grid and ranking, whatever ranking the
  options name; an archive from a subclass's `minScore` 5 has the bits of
  the compile-time problem with 5.0 as the template's argument, QD score
  included; a stage and an archive stopped and resumed through the entry
  points equal the uninterrupted ones; the learned model through them;
  and each thing they refuse. `checkProblem` names two more unsound
  problems (seven now).

**The equivalence check**, on Windows with LDC 1.42.0: the example's 12
lines are byte-identical to those of solver's example at `ae68dbb`. And
stronger than the printed lines, which round: copies of both examples
(outside both checkouts) changed to print the final archive as
hexadecimal floats gave identical output: all 16 thresholds, every
elite's value, genes and `found`, and the QD score of the resumed and of
the uninterrupted run.

Checked on Windows: `dub test` (10 modules), the example's tests
(11 modules) and its run. Not rerun: the spike (it uses `cmaStage` only,
and `stage.d` didn't change). Not run on Linux.

Found on the way:

- **The worked example's values go below its `minScore`.** `Beam`'s
  values run from about −48 (width 2, rate 100) to 10, and its archive's
  thresholds start at 0: the final archive holds elites at −3.06 and
  −0.40, in cells whose thresholds never moved. The manual says to keep
  `minScore` at or below the lowest value. It is solver's example and
  its bytes are the equivalence check, so it is left as it is; see open
  question 8.
- **dub leaves the example's `app.d` out of its test build** ("Excluding
  main source file"), so a `unittest` there never runs. The example's
  tests are `checkpoint.d`'s and the core's.
- `StageEnd.rise` is NaN unless the stage stalled, so two `StageEnd`
  never compare equal with `==`: compare their fields.

For step 4:

- The binding's `ForeignProblem` overrides `axes`, `ranking` and
  `minScore` from the caller's data, and the exported functions call
  `runStage` and `runArchive`. Their checks already turn a wrong length
  into an exception with a message, before any index is used.
- An archive saved as bytes must include `Archive.minScore`; the
  checkpoint template does it by walking the fields.

### 4. The C interface (done on Windows 2026-11-12)

`bindings/c`: a stage, an archive, a learned model, and a save, load and
resume all run through the exported functions, and a resumed run equals
the uninterrupted one to the last byte of its saved state. Its Linux run
waits for a Linux machine (step 6).

**What is there** (`bindings/c/README.md` has the commands and the file
list):

- `include/solver2.h`, by hand: 37 functions (38 since step 5), 14
  structs with fields, 4 handle types. Its opening comment is the
  contract.
- The library, `solver2.dll` / `libsolver2.so`, in two configurations:
  `library` (the D runtime linked in, bounds checks on) and `unchecked`
  (the same without the checks). Both write the same file;
  `solver2_checked()` tells them apart.
- `ForeignProblem` (`foreign.d`): the problem as data, one callback per
  batch with flat arrays, each breach summed by `breachOf` on the D side.
  Its description is checked when it is made (`checkDescription`, the
  half of `checkProblem` that evaluates nothing, split out in the core
  for this), and its fingerprint is the caller's pairs plus the
  description itself, so a run saved on other genes, limits, axes,
  ranking or `min_score` is refused with what differs.
- A stage, an archive and a model between calls (`runs.d`), over
  `runStage` and `runArchive`. One report per generation or iteration,
  whose return value goes on, stops (to resume) or reports a failure.
- Saving (`checkpoint.d`): the binding copies the manual's checkpoint
  template, as open question 2 recommends, and the core gets no
  serialiser. A saved run holds its model's state, and how the model was
  made, so that it is restored into one made the same way.
- Handles are a slot and a count, never an address: a handle freed, of
  another kind, or never one is an error.

**Decisions made here**

- **A failed evaluation leaves a run whole.** When the caller's
  evaluation fails (any Python exception, a Ctrl-C), the run and its
  model are put back as after their last whole generation or iteration,
  and the same call made again ends as the run that never failed, bit
  for bit. A stage already keeps its state after every generation; an
  archive is copied after every iteration (its emitters deep, its cells
  shallow), and a model saved as bytes, because the searches train the
  model and sample from the emitters before they evaluate. The cost is
  below the training's own, and nothing when a run has no model but for
  the archive's copy. The alternative, a run that is broken after a
  failure, would have lost a notebook user's run to one Ctrl-C.
- **The options' structs have no ranking**: it is the problem's, as
  step 3 made it.
- **The default functions** (`solver2_*_default`) fill the structs from
  the D defaults, so a wrapper never repeats a number.
- **`checkpoint.d` is a copy, not the example's file used in place**: the
  bytes are the binding's format from here on, with its own magic
  strings and version, and a change of the manual's template doesn't
  change it unseen.

**The thread rule changed: found by the driver.** The plan had every
exported function attach an unknown thread and detach it before
returning. That is not safe with druntime 2.112 (LDC 1.42), in two ways:

- `thread_detachThis` takes the thread off the runtime's list but leaves
  it as the runtime's "this thread". The thread's next call finds itself
  known, doesn't attach, and runs where the collector doesn't see its
  stack. With a model's threads allocating, the collector then frees
  memory the thread is using: of 8 runs of the driver, its second
  thread got other bits in 5, and in one more the process died without a
  word.
- With that cleared (`thread_setThis(null)`), the second attach makes a
  new `Thread` object, while the thread's block cache
  (`core.internal.gc.blkcache`) stays recorded on the first, which is
  gone. The collector no longer clears the cache of blocks it frees, and
  the heap was corrupted in every run (25 of 25, exit code 0xC0000374).

So a thread is attached **once, on its first call, and stays attached**,
which is what the spike had shown working, and is detached when it ends
by a hook the operating system calls on the thread (fiber-local storage
on Windows, a pthread key elsewhere). Checked on Windows: the hook runs
(a one-off build printed the runtime's thread count going from 2 to 1 as
the thread ended), and the driver passed 31 of 31 runs with two foreign
threads at once and twelve one after another. What it means for a
caller:

- a thread that has called the library is stopped for a moment, like
  any D thread, whenever the D collector runs, for as long as it lives;
- `solver2_init`'s first call is to come from a thread that outlives the
  library's use: it becomes the D runtime's main thread;
- a process ended from a foreign thread while another holds the
  runtime's thread lock could hang in the hook. Not seen; Python ends
  from its main thread.

The spike's attach test passed because it attached each thread once, by
hand. Its conclusion that the exported functions "can do it themselves"
was right only for that.

**The tests**

- `dub test` in `bindings/c` (13 modules): the interface in one process
  with the searches. A stage, an archive and both with a model give the
  bits of `runStage` and `runArchive` on the same problem as a D class:
  every report, the best, every cell and threshold, the counts, the
  diverse picks. A run stopped by its report, saved, freed, loaded and
  resumed has the saved bytes of the uninterrupted one, with a model
  too, and on any number of the model's threads. The evaluation made to
  fail at five different batches, the seeds' and the first iteration's
  among them, and the report made to fail: the call says so, the run is
  where the header says, and gone on with it has the uninterrupted run's
  bytes. Each thing a load refuses. With the archive's roll-back
  disabled, or the model's in a stage or in an archive, these tests fail
  (checked by hand, one at a time).
- **The header is held to the library** (`conformance.d`): `solver2.h`
  is read by ImportC and compared with the D declarations, each struct's
  fields, offsets and size, each function's parameters, each constant,
  in both directions. A wrong field, constant or parameter planted on
  either side failed the build with its name.
- **The driver** (`bindings/c/test`): a program built with `-betterC`,
  so with no D runtime, that loads the library at run time and knows
  only the header. It is the Python situation without Python, and it is
  what found the thread fault. It ran against both configurations.

Checked on Windows with LDC 1.42.0: the core's tests (10 modules), the
example's (11) and its run, whose output is still solver's; the
binding's tests (13); the driver against `library` and `unchecked`.

And from Python itself, with a throwaway `ctypes` script (not kept;
step 5 writes the module): a stage with its evaluation in Python gave
the same bits on the main thread, after a Python exception in batch 9
(the call returned `SOLVER2_CALLBACK` at generation 8 and the stage went
on), and from 22 Python threads, four at a time and one after another.
Six runs of it, Python 3.15.0.

**Not done, or not exposed**

- Linux: nothing here has run there. Written from the documentation and
  untested: the pthread half of the exit hook, `dlopen` in the driver,
  `libs-posix` in the driver's `dub.json`, and whether `-fvisibility=hidden`
  with the runtime linked in gives a `.so` that loads beside nothing
  else.
- A bug's `Error` caught at the boundary (an index out of range with the
  checks on) has no test: the spike showed it, and nothing in the
  interface can be made to raise one on purpose.
- Not in the interface: the archive's `evaluated` callback and its
  per-emitter reports, focused runs, a grid other than `BinnedGrid`, a
  subclass's own learning spec, `Best` carried from one stage to the
  next (a stage has its own; start the next from the last one's best
  genes).
- `solver2_term` and then `solver2_init` again in one process: druntime
  doesn't promise it, and the header says not to.

**For step 5**

- The wrapper's job is small: build the spec from a Python object, keep
  the callback object alive as long as the problem, keep the exception
  an `evaluate` raised and raise it again when the call returns
  `SOLVER2_CALLBACK` (the run can then be continued).
- `ctypes` structs must match the header. A test that compares their
  sizes and offsets with numbers from the library would need those
  exported; simpler, check `solver2_abi()` and keep the Python structs
  next to a note of the header's version.
- A Python `evaluate` gets flat arrays: numpy views with
  `np.ctypeslib.as_array`, nothing copied (the spike).

### 5. The Python module (done on Windows 2026-11-12)

`bindings/python`: the module, its tests, and the worked example in
Python, which prints what its D counterpart prints, bit patterns
included. Its Linux run waits for a Linux machine (step 6).

**What is there** (`bindings/python/README.md` has the setting up and
what a user needs to know):

- `solver2/_c.py`: the header as `ctypes` declarations, by hand, and the
  loading of the library: the file named by `SOLVER2_LIBRARY`, else
  beside the module, else `bindings/c` of the checkout (open question 3:
  a plain folder, as recommended).
- `solver2/__init__.py`: `Problem` (a class with its description as
  class attributes and an `evaluate(genes, feasible_only)` that returns
  the objectives, the violations and the descriptors as arrays), `Model`,
  `Stage` and `Archive` (`run`, `result`, `best`, `elites`, `thresholds`,
  `causes`, `diverse`, `save`, `load`), and `Error`. Options are
  keywords over the library's defaults (`lambda_` for `lambda`).
- An exception raised in `evaluate` or in a report is kept, the library
  is told the callback failed, and `run` raises the exception itself, a
  `KeyboardInterrupt` included. The run can then be gone on with (step
  4's roll-back).
- `examples/beam.py`, and `docs/example-class`, its D counterpart.
- The tests, with the standard library's `unittest` (pytest isn't
  installed here, and the module then needs nothing but numpy): 30 of
  them, in about 30 s.

**The plan's goal had to change: the Python example equals a D example
on the class route, not `docs/example`.** `docs/example`'s `Beam` is a
compile-time problem with a learning spec of its own (it learns the
objective, and one breach for the hard and the terminal tier). A problem
behind the C interface has the class's default learning spec (the value,
a breach per tier), by design: the model's threads never call back. The
two models differ, so no Python `Beam` can give `docs/example`'s bits
once a model is in the run, and every part of the example has one. So:

- `docs/example-class` is the worked example as a `DynamicProblem`,
  through `runStage` and `runArchive`: the same problem, flow and
  options, the class's learning spec, and `minScore` −50, below the
  problem's lowest value, so that `checkProblem` passes (open question 8
  doesn't arise for it). It is a second worked example in its own
  right, for the manual's class route (step 7).
- `examples/beam.py` is that program in Python. The two print the same
  16 lines; two of them are bit patterns: the stage's best value and a
  hash of its genes, and hashes of every threshold, every elite's value
  and genes, and the QD score. `tests/test_example.py` builds the D one
  with dub and compares.
- `docs/example` is untouched, and still prints what solver's does.

**Added to the C interface: `solver2_problem_physical`** (`SOLVER2_ABI`
is 2; 38 functions). It gives genes as physical values by the library's
own `GeneSpec.toPhysical`. The Python example needs it for its bits: the
rate's log scale goes through `exp` and `log`, the library's are the C
library LDC links, and numpy has its own, which round otherwise. It is
also what a Python problem wants anyway, instead of writing the clamping
and the scales again. It may be called from inside the evaluation.

**Also found**

- dub doesn't see that the driver and the binding's tests depend on the
  header, which they reach through a `#include`: after the header
  changed, `dub run` ran the old driver, whose check of the interface
  version then failed (as it should). Both now list `include` as a
  string import path, which dub tracks.
- The worked example takes 43 s, the original too: nearly all of it is
  its model training on one thread (5 networks of 64 units, trained
  before each of about 340 generations and iterations, on up to 4000
  rows). The two new
  examples give their models 4 threads, which changes no bit (checked:
  the same output as on none): 8 s in D, 13 s in Python.
- What Python costs: about 26 µs a batch when `evaluate` computes
  nothing, and 55 µs with a numpy evaluation of 8 candidates. In the
  example, 12.5 s against the D program's 7.6 s: 1.7 s of the
  difference is the bounds checks (10.8 s with the `unchecked` library)
  and the rest mostly the model saved after every iteration for the
  roll-back. If that matters, a model could be rolled back by forgetting
  the rows it learned since, which needs a way into `LearnedModel`.
- `ctypes` wants a null function pointer, not `None`, for an absent
  callback (`ctypes.cast(None, kind)`).

**The tests** (`python -m unittest discover -s tests`, in
`bindings/python`)

- `test_conformance.py` reads `solver2.h` with a small parser and holds
  `_c.py` to it: every struct's fields and types in order, every
  function's and callback's result and parameters, every constant, in
  both directions. A wrong field and a wrong parameter planted in
  `_c.py` each failed it.
- `test_solver2.py`: a problem's description and what is refused of it;
  `physical`; `check`; a stage and an archive, with and without a model,
  the same on any number of the model's threads; stopped, saved, closed,
  loaded and resumed to the uninterrupted run's bytes; an exception in
  `evaluate` at several batches and in a report, raised as itself, and
  the run gone on with to the same bytes; what `evaluate` returns
  checked; a load's refusals; seven Python threads, four at once.
- `test_example.py`, above.

Checked on Windows with LDC 1.42.0 and Python 3.15.0 (numpy 2.5.4): the
core's tests (10 modules), the example's (11) and its run, whose output
is still solver's; the binding's tests (13) and the driver; the Python
tests (30).

**Not done**

- Linux, as for step 4.
- A wheel, or any packaging beyond the folder (open question 3).
- Older Python versions: only 3.15 was run. The module uses nothing
  newer than `namedtuple` defaults (3.7).
- The interface's omissions are the module's too (step 4, "Not done, or
  not exposed").

### 6. CI (done 2026-11-12)

Dev published the repository on 2026-11-12:
[x21-dev/solver2](https://github.com/x21-dev/solver2), public. From here
on a step is a branch and a pull request, and `main` is what has passed
CI.

`.github/workflows/ci.yml` runs on every push to `main` and every pull
request, three jobs: `ubuntu-latest` and `windows-latest` on the newest
Python (3.15.2 when it was written), and `ubuntu-latest` once more on
Python 3.10, the oldest tested. Each job, with LDC 1.42.0:

1. the core's tests, the example's tests and the example (solver's
   workflow, as it was);
2. the C interface's tests (`bindings/c`), which read the header with
   ImportC;
3. the shared library built without bounds checks and the driver
   against it, then the same with the checks, which is the build left
   in place;
4. numpy from pip, then the Python tests (`bindings/python`), a
   `RuntimeWarning` failing them. `tests/test_example.py` builds
   `docs/example-class` and compares the Python example with it: the
   same bits through Python, within the job, so on one platform.

A job takes two to three minutes. The repository is public, so the
minutes aren't billed.

**Linux, at last.** Nothing of steps 2 to 5 had run on Linux. All of it
passed there on its first run (pull request 1): the core's tests, both
examples, the binding's tests with ImportC's `cc`, `libsolver2.so` with
the D runtime linked in and hidden visibility, loaded by `dlopen` from
the driver with no D runtime of its own, the threads' exit hook on a
pthread key (two foreign threads at once, twelve one after another),
and the Python module on 3.15.2 and 3.10.22, its example printing what
the D class-route example prints. What was written from the
documentation and marked untested in steps 4 and 5 needed no change.
(Not so: the runtime linked statically was a fault on Linux that these
runs happened to miss. Step 7 found it and has the fix.)

One thing the Linux run on Python 3.10 showed, in the module: reading an
archive's elites went through numpy's view of a `ctypes` array of
structs, and Python before 3.11 describes a padded struct wrongly to
numpy, which then guesses the layout and warns. The guess was right, but
the module now gives numpy the layout from the `ctypes` declaration
itself, and CI fails on a `RuntimeWarning`.

Not in CI:

- **The equivalence with solver** (the example's output against solver's
  at `ae68dbb`). It is still a local check (see step 0). solver's
  repository is public too, so the workflow could fetch it and compare;
  say if you want that.
- **The spike** (step 1) was never run on Linux and won't be: the driver
  is what it led to, and that runs there.
- **macOS**, as planned.

The old plan for this step, kept for what it says of the matrix:

Keep the matrix small: on a private repository Windows minutes cost
more than Linux's and macOS minutes far more. macOS is not planned;
runs differ between platforms anyway, so adding it would mean a third
set of reference outputs.

### 7. The manual (done 2026-11-13)

`docs/manual.md` is solver2's: every section brought to solver2's state,
two new sections, the date and Gaps and caveats true on 2026-11-13.

**What changed in the manual**

- **The head**: what solver2 is, the three routes a problem can take (a
  compile-time struct, the class, Python or C) with a worked example
  each, and where solver stands. The note that it was still solver's
  manual is gone, and "solver" now means the older project wherever it
  appears.
- **Interfacing with Python** (new, section 8): setting up, a problem in
  Python, a stage, the archive, a learned model, stopping and saving,
  errors and exceptions, threads, determinism with numpy, what Python
  costs, the library's bounds checks, and the C interface for other
  languages. `bindings/python/README.md`'s "What to know" moved here;
  that README is now the quick start, the setting up, the tests and a
  pointer.
- **The C interface's small program** in that section (a stage over a
  one-gene problem) was compiled and run by hand against the built
  library (`ldc2 demo.c solver2.lib`, ImportC; it printed `best
  10.000000 after 400 evaluations`). It is not a tracked file and no
  test builds it.
- **Migrating from solver** (new, section 12): who needs to, what stays,
  what changes (from the list below), how to prove the same bits with a
  golden run, saved states, and moving a compile-time problem to the
  class, with a table of what becomes what and what the move costs.
- **Gaps and caveats**: 18 items, in three groups (the searches, the
  class route, Python and the C interface). New: a prediction ranked
  without its soft breach (found in step 2), `StageEnd.rise` (step 3),
  `checkProblem`'s limits, what the entry points fix, the library owning
  the loop (question 6), the folder that isn't a wheel (question 3),
  what the interface leaves out (step 4), the library not started
  twice, the Pythons tested, the saved bytes' format, and the model
  saved for the roll-back (step 5).
- Smaller: the class route's cost under Performance, the model's pool in
  the library, the module reference's table of what lies outside the
  core, a pointer from Determinism to numpy's rules.

**`docs/example` has a `minScore` of −50** (question 8), as the template's
argument: `Archiving!(Beam, BinnedGrid, ArchiveOptions, 0, MIN_SCORE)`.
Its stage's five lines are unchanged, and are solver's. From the archive
on the output differs, as expected: the QD score counts from −50 (586.33
after 10 iterations, where it was 39.17), and the run takes another path
(4 explorer targets reached where there were 8; 14 of 16 cells and a
best of 9.9880 both ways). There were no reference outputs kept in the
repository to renew: the example checks itself (the resumed run against
the uninterrupted one), and the comparison with solver is by hand (step
0 says how, now).

**A checkpoint written with solver reads in solver2: checked.** The plan
had it from the template's own test only. On Windows: solver's example
at `ae68dbb`, exported outside its checkout and changed to keep its
checkpoint (written after 10 iterations, 598,631 bytes), and a copy of
solver2's example at `b673f92` (so `minScore` 0) given solver's magic
string and made to load that file. It loaded, resumed, and passed the
example's own check against its uninterrupted run, bit for bit, printing
solver's twelve lines. The manual's "Saved states" says so.

**Also changed**: `README.md` (what solver2 is, not what it is becoming),
`bindings/c/README.md` (a pointer to the manual), `CLAUDE.md` (the
manual is solver2's; the equivalence check as it now stands; solver
frozen by default; the manual's Python section changes with
`bindings/`).

Checked on Windows with LDC 1.42.0 and Python 3.15.0: the core's tests,
the example's tests and its run, the binding's tests, the library and
the driver, the Python tests; and on Linux by CI. The manual's anchors and file links were
checked by a throwaway script (every `](#...)` has its heading, no
heading twice).

**Found by CI on this step's pull request: the library crashed on Linux,
about once in 50 runs.** The driver segfaulted in the job "ubuntu-latest,
Python 3.x" (exit code -11, nothing printed), while the other Linux job
passed the same step. The step changes no code of the library, so the
fault was there since step 4, and five earlier CI runs had missed it.

- **Hunted on CI**, there being no Linux on this machine: a temporary
  branch, `debug-driver-linux`, with a workflow of its own (on a push to
  `debug-**`): the library built with `DFLAGS=-g`, core dumps on
  (`ulimit -c unlimited`, `kernel.core_pattern` set to a file), the
  driver run up to 250 times in each of six jobs, and on the first
  failure what it had printed and `gdb -batch -ex "thread apply all bt"`
  on the core. All six jobs crashed, at runs 21 to 154, all at the same
  place.
- **Where**: the main thread, in the archive with a model, in a
  collection: `Gcx.collectRoots` under `rt.sections_elf_shared.scanTLSRanges`
  under `thread_scanAll`. Not the foreign threads, which a first guess
  had blamed (druntime's `thread_attachThis` allocates before it
  registers the thread; a stress test of 7,200 attaches on Windows found
  nothing, and that guess stays unproven either way).
- **Why**: druntime's `getTLSRange`, with a static runtime on Linux x86
  (`Static_Linux_X86_Any` in LDC 1.42's `rt/sections_elf_shared.d`),
  takes a thread's TLS block to be the `size` bytes just below its
  thread control block. That holds when the D code is the executable.
  In a shared library the block is elsewhere, so every collection
  scanned memory that wasn't the library's, and the library's own TLS
  was never scanned. On the main thread the wrong range can run off the
  start of a mapping, depending on where the loader put things: hence
  one run in 50.
- **The fix**: on Linux the library uses LDC's shared runtime, which
  finds TLS through the loader (`bindings/c/dub.json`:
  `-link-defaultlib-shared=true` and a run path of `$ORIGIN` for posix;
  Windows as before). The same hunt on the fixed library: six jobs of
  250 runs, 1500 in all, no crash. `ldd` shows the two libraries, and
  the run path is `$ORIGIN` and then LDC's `lib`.
- **What it costs**: on Linux the library is no longer one file. It needs
  `libdruntime-ldc-shared.so.112` and `libphobos2-ldc-shared.so.112`,
  found in the LDC that built it, or copied beside it. The manual, both
  READMEs and its gap 18 say so. Step 1's "the runtime must be linked
  statically" was shown on Windows and holds there; its guess that Linux
  would need the same flag was wrong, and step 6's first green runs on
  Linux were luck.
- **So that it stays found**: CI checks on Linux that the library uses
  the shared runtime (`ldd`), and the driver's output is unbuffered, so
  a crash no longer takes the lines before it along.
- **Not checked**: whether, with the shared runtime, the TLS of a thread
  attached by `thread_attachThis` is scanned (druntime keeps the list of
  loaded libraries per thread, and fills it for threads it starts). The
  binding keeps nothing of the collector's in thread-local storage (the
  last error is a plain buffer), so nothing depends on it today; keep it
  so.

**Decided here, for Dev to overrule**

- **`docs/example`'s learning spec still learns the objective**, though
  the manual now says a problem with soft limits should learn the value
  (gap 3, and a bullet under "The learning spec"). Changing it would
  change the stage's lines, the last of the example's output that equals
  solver's. The manual says so where it shows the spec. To change it:
  `t[0] = e.scored ? e.objective - e.breach.soft : double.nan;` in
  `Beam.targets`, in `app.d` and in the manual's copy, and the
  comparison with solver then needs the old line in a copy.
- **`checkProblem` refusing a neutral value outside its range is not in
  the manual's migration list**: solver has no `checkProblem`, so no
  migrating project meets it. It stays on the list below as a change
  within solver2.

### 8. Releases (2026-11-13)

Dev's decision, after step 7 made the Linux library need LDC's runtime
beside it: distribute built archives as GitHub releases, the library
bundled with the runtime libraries it was built against. Ubuntu 22.04 is
to stay supported even when `ubuntu-latest` moves on. The licence stays
proprietary, with licences granted individually.

**What is there**

- **`release/package.py`** makes the archive for the platform it runs
  on, from the library built in `bindings/c`:
  `solver2-<version>-linux-x86_64.tar.gz` or
  `solver2-<version>-windows-x86_64.zip`. It unpacks to one folder:
  `solver2/` (the Python module with the library inside it; on Linux
  LDC's two runtime libraries under the names the library asks for, on
  Windows the import library `solver2.lib`), `include/solver2.h`,
  `examples/beam.py`, `README.md` (from `release/README.md`, with the
  version and the interface version filled in) and `LICENSE.md`. One
  copy of the library serves Python and C. It refuses a library built
  without bounds checks, and on Linux one that doesn't use the shared
  runtime.
- **`release/check.py`** checks an unpacked archive where it lies, with
  nothing of the repository: the module loads the archive's own
  library, the worked example runs to its end, and with `--expected` it
  prints exactly the lines `docs/example-class` printed where the
  archive was built.
- **`.github/workflows/release.yml`**: for a tag `v<version>`, for every
  pull request, and by hand.
  - *Linux archive*: built in an `ubuntu:22.04` container on
    `ubuntu-latest`, with LDC from its own release tarball. So the
    system built on doesn't depend on which runner images GitHub keeps.
    The job also runs the driver and the Python tests there, on 22.04's
    Python 3.10.
  - *Linux tests*: in `ubuntu:22.04` and `ubuntu:24.04` containers with
    no LDC and no checkout, with the system's own Python and numpy: the
    library finds its runtime in the archive's folder (`ldd`), the
    driver passes against it, and the example runs; on 22.04 to the
    same bits as where it was built. On 24.04 only to its end: another
    glibc may round otherwise.
  - *Windows archive* and its test on a runner without LDC, the example
    to the same bits.
  - *Draft release*, for a tag only, once all tests pass: the two
    archives and notes. A draft, so that publishing is Dev's.
- **`LICENSE.md`**: proprietary, licences granted individually; and what
  the releases include of LDC's runtime, under its own licences. The
  wording is the session's: Dev to read it.

**To make a release**: on `main`, with CI green, `git tag v0.1.0` and
`git push origin v0.1.0` (Dev's word for both); the workflow drafts the
release; Dev reads and publishes it. No release has been made yet.

**Checked on Windows** (LDC 1.42.0, Python 3.15.0): the archive made by
`package.py`, unpacked elsewhere; `check.py` on it with the lines of
`docs/example-class`, the same 16; the driver against the archive's
`solver2.dll`. The workflow itself can only be shown on CI.

**Not done**

- A wheel (question 3 stands): the archive is most of the way there.
- arm64, macOS.
- The module has no version of its own: a release's version is its tag,
  in the archive's name and its `README.md`.
- The archives aren't signed, and have no checksums beyond GitHub's.

### What is left open

Nothing more is planned. What the steps left open, should it be wanted:

- **The equivalence with solver in CI** (step 6): fetch solver at
  `ae68dbb`, run both examples, compare the stage's lines, and the
  archive's with `MIN_SCORE` 0.
- **The manual's C program as a tested file**, built and run by CI
  against the library, so that it can't rot.
- **A wheel** (question 3) and **Python owning the loop** (question 6):
  decided against for now. Step 8's archives are what a wheel would
  hold.
- **A cheaper roll-back for the model** (step 5): forget the rows learned
  since the last whole iteration instead of saving the model each time.
- **The interface's omissions** (step 4, "Not done, or not exposed"; the
  manual's gap 13).
- **macOS**: untested, and it would need a third set of goldens.

## Changes a migrating project will meet

Kept as the work went, for the manual's migration section.

- **The names** (step 0): the dub package is `solver2` and every module
  `solver2.*`. A project changes its dependency (`"solver2": { "path":
  "../solver2" }`) and its imports (`import solver.problem;` becomes
  `import solver2.problem;`, and so for `stage`, `archive`, `model`,
  `ensemble`, `cmaes`, `bytes`, `rng`, `detmath`). The bits are the same.
- **A subclass of `DynamicProblem`** (step 2; no project is known to have
  one): `genes`, `constraints`, `fingerprint` and `evaluateInto` are
  `abstract`. The class now declares `axes`, `ranking`, `screenPenalty`,
  `noFeasibleSeed`, `adoptCauses`, `featureCount`, `features`,
  `targetCount`, `targets`, `predicted`, `predictionHeadroom`, `grid`,
  `startGenes`, `newCauses` and `countCauses`: a subclass with a member of
  one of those names must `override` it (the last four are final). Its
  outcomes are now checked: one within the limits needs `scored`, an
  objective and its descriptors. Step 3 adds `minScore` to those names.
- **`Archive` has one more field, `minScore`** (step 3), after `alpha`.
  The bits of a run are the same. What a project may meet:
  - a checkpoint written by the manual's template gains the field. One
    written with solver lacks it and still reads: the field keeps the
    template argument's value, which is where that archive's thresholds
    started. (The template's own test covers a field added since; no
    solver file was read in solver2 to show it.)
  - a save format of the project's own (pleaser's `elites.state`) is
    untouched: an archive it builds has the template argument's value.
  - an `Archive` built as a struct literal with every field in order
    would need the new one. With `N` 0 it has a constructor, so only an
    `N > 0` archive could be written so; none is known.
- **`checkProblem` refuses a gene whose neutral value is outside its
  range** (step 4). It used to pass: the start genes are clamped, so the
  check on them could never fail. `checkDescription` is the half of
  `checkProblem` that evaluates nothing. (A change within solver2:
  solver has no `checkProblem`, so the manual's migration section leaves
  it out.)
- **A problem's member named `minScore` is now read by the archive**
  (step 3), as where its thresholds start, and `Archiving`'s `minScore`
  argument must then stay 0. None of pleaser, pleaser-bus-chain and
  pleaser-console-tilt has such a member (they pass `MIN_SCORE` as the
  argument, which goes on working).
- **The worked example's archive starts its thresholds at −50** (step 7,
  question 8), where solver's starts them at 0. Not a change a project
  meets in its own code; one that compares its build's example output
  with solver's does. The stage's lines are the same.

The manual's "Migrating from solver" has this list, in its own order.
Add to both.

## Open questions for Dev

None is open: all nine were decided, the last six on 2026-11-13. They
keep their numbers, which the steps above refer to.

1. **Names.** Decided 2026-11-12: package `solver2`, modules `solver2.*`.
   A build can never pick up the wrong library, and migration is a search
   and replace of the imports.
2. **Checkpoints.** Decided 2026-11-13: the core gets no serialiser.
   solver has no checkpoint API by decision: the save format is the
   manual's template, copied by each project. Python can't walk D
   structs, so the binding exports a run as bytes, and for that it
   copies the template as any project does
   (`bindings/c/source/solver2c/checkpoint.d`, step 4).
3. **Python packaging.** Decided 2026-11-13: the plain folder of step 5,
   whose module loads the library by path. A wheel (per-platform builds
   with LDC, versions, releases) waits until someone needs to install
   the module without a checkout.
4. **solver from here on.** Decided 2026-11-13: solver is frozen by
   default. A bug found in code both hold (`cmaes.d`, `ensemble.d`) is
   fixed in solver2 and listed under "Changes a migrating project will
   meet"; tell Dev, who decides for each bug whether solver gets the fix
   too, in a session of solver's own (never from here). A fix that
   changes bits would invalidate the golden runs of the projects on
   solver, which is why it is not automatic.
5. **A remote and CI.** Decided 2026-11-12: Dev creates the GitHub
   repository once solver2 is mature enough, and step 6 then extends
   solver's workflow. Done the same day: Dev published
   x21-dev/solver2, and step 6 went in as a branch and a pull request,
   pushed from the session with Dev's word. Ask before pushing; merging
   is Dev's.
6. **Python owning the loop** (ask for candidates, evaluate, tell the
   results, as pycma and pyribs do). Decided 2026-11-13: not planned.
   Possible later without rewriting the searches, by running them in a D
   fiber that pauses at each `evaluate`. Step 7 lists it in the manual's
   Gaps and caveats.
7. **Bounds checks in the binding.** Decided 2026-11-12: on in the
   Python binding only, and the manual says how to turn them on or off in
   any case, for the binding and for a D project. The core, as a source
   library, takes whatever flags its user builds with.
8. **The worked example's `minScore`** (found in step 3). `Beam`'s values
   go down to about −48 and its thresholds start at 0, against the
   manual's own rule, so `checkProblem` would refuse it. Decided
   2026-11-13: `docs/example` keeps its bytes until step 7 (they are the
   equivalence check with solver, one comparison of two outputs), and
   step 7 gives it a `minScore` of −50. Its output then differs from
   solver's from the archive on, the check narrows to the stage's lines,
   and the change joins "Changes a migrating project will meet".
9. **Which D example the Python one must equal** (step 5). Decided
   2026-11-13: `docs/example-class`, as step 5 built it. The plan had
   said `docs/example`; that can't be, because its learning spec is its
   own and a problem behind the C interface has the class's. The other
   way, letting a foreign problem choose its learning spec's layout
   (which targets, one breach or two) so that a Python `Beam` could be
   given `docs/example`'s, is not done: it would shape the interface
   around one example.

## How to work

- CI runs the checks below on Linux and Windows for every pull request;
  run them locally first all the same, they take two minutes.
- LDC only, 1.42.0. Checks before every commit: `dub test --compiler=ldc2`
  in the core, and in `docs/example` both `dub test --compiler=ldc2` and
  `dub run --compiler=ldc2 --build=release`; in `bindings/c`
  `dub test --compiler=ldc2`, `dub build --compiler=ldc2 --build=release`
  and then in `bindings/c/test` `dub run --compiler=ldc2`; and in
  `bindings/python` `python -m unittest discover -s tests`, which also
  builds and runs `docs/example-class`.
- Determinism and D style: `CLAUDE.md`, as in solver. In
  short: seeded noise only through `SplitMix64`; no `@fastmath` or
  `-fp-contract=fast`; fixed-order sums; goldens per platform; a test
  that needs an event during a run forces it; large values off the stack.
- `docs/example/source/checkpoint.d` is printed whole in the manual; a
  `static assert` fails the build when they differ. Change both together.
- Commit only when Dev asks.
