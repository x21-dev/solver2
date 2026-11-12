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
| 2. The class, completed | done 2026-11-12 |
| 3. Entry points and run-time `minScore` | done 2026-11-12 |
| 4. The C interface | done on Windows 2026-11-12; Linux owed |
| 5. The Python module | next |
| 6. CI | waits for Dev to create the GitHub repository; Linux runs of steps 1 to 5 are owed to it |
| 7. The manual | |

solver2 is solver's code with the package and modules renamed, plus the
completed `DynamicProblem` (step 2), the archive's `minScore` as a value
known at run time, the searches over the class as plain functions in
`solver2.dynamic` (step 3), and their C interface as a shared library in
`bindings/c` (step 4). The searches otherwise are unchanged: the worked
example still prints what solver's does, and its final archive has the
same bits.

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

- `include/solver2.h`, by hand: 37 functions, 14 structs with fields, 4
  handle types. Its opening comment is the contract.
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

- **The shared library**: in `bindings/c`, `dub test --compiler=ldc2`,
  then `dub build --compiler=ldc2 --build=release`, then in
  `bindings/c/test` `dub run --compiler=ldc2`, which drives the library
  through the exported functions only (it needs the library built
  first, and a C compiler for ImportC's preprocessor, which both runners
  have).
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
  `checkProblem` that evaluates nothing.
- **A problem's member named `minScore` is now read by the archive**
  (step 3), as where its thresholds start, and `Archiving`'s `minScore`
  argument must then stay 0. None of pleaser, pleaser-bus-chain and
  pleaser-console-tilt has such a member (they pass `MIN_SCORE` as the
  argument, which goes on working).

## Open questions for Dev

1. **Names.** Decided 2026-11-12: package `solver2`, modules `solver2.*`.
   A build can never pick up the wrong library, and migration is a search
   and replace of the imports.
2. **Checkpoints.** solver has no checkpoint API by decision: the save
   format is the manual's template, copied by each project. Python can't
   walk D structs, so the binding has to export a run as bytes.
   Recommended: the binding copies the template, as any project does, and
   the core still gets no serialiser. Step 4 did so
   (`bindings/c/source/solver2c/checkpoint.d`); say if the core should
   have one after all.
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
8. **The worked example's `minScore`** (found in step 3). `Beam`'s values
   go down to about −48 and its thresholds start at 0, against the
   manual's own rule. Step 5 wants a Python `Beam` whose result equals
   the D example's bit for bit: as a `DynamicProblem` with `minScore` 0
   it runs the same, but `checkProblem` refuses it (a candidate within
   the limits below `minScore`). Either:
   - keep the example's bytes (they are the equivalence check with
     solver) until step 7, and have the Python example say why it
     doesn't call `checkProblem`, or call it on another problem
     (recommended: the check with solver stays one comparison of two
     outputs for as long as the searches are being changed); or
   - give the example a `minScore` of −50 now. Its output then differs
     from solver's from the archive on, the check narrows to the stage's
     lines, and the change joins the list above.

## How to work

- LDC only, 1.42.0. Checks before every commit: `dub test --compiler=ldc2`
  in the core, and in `docs/example` both `dub test --compiler=ldc2` and
  `dub run --compiler=ldc2 --build=release`; in `bindings/c`
  `dub test --compiler=ldc2`, `dub build --compiler=ldc2 --build=release`
  and then in `bindings/c/test` `dub run --compiler=ldc2`; from step 5
  the Python tests.
- Determinism and D style: `CLAUDE.md`, as in solver. In
  short: seeded noise only through `SplitMix64`; no `@fastmath` or
  `-fp-contract=fast`; fixed-order sums; goldens per platform; a test
  that needs an event during a run forces it; large values off the stack.
- `docs/example/source/checkpoint.d` is printed whole in the manual; a
  `static assert` fails the build when they differ. Change both together.
- Commit only when Dev asks.
