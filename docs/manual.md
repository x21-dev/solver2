# solver2: developer manual

solver2 is a deterministic black-box optimiser written in D, built with
LDC. You describe a problem: genes in [0, 1], an objective to maximise,
and limits in three tiers. solver2 searches it with two tools:

- **A CMA-ES stage** (`solver2.stage`) finds the best candidate from a
  start point. It stops early when progress stalls, and a learned model
  can pre-screen its candidates.
- **An archive** (`solver2.archive`) runs CMA-MAE, aimed at reaching cells
  no candidate has found yet. It keeps the best candidate in every cell of
  a grid over your problem's descriptors, and its explorers head straight
  for empty cells. It has optional model-guided and pre-screened
  emitters. A run can stop between iterations and resume bit for bit.

Both tools are deterministic. The same problem, seed and options give the
same bits on any thread count of one platform, provided the problem keeps
the rules in [Determinism](#determinism). Across platforms runs differ.

A problem can be written in three ways, and all three run the same
searches:

| Route | The problem is | Use it when | Worked example |
|---|---|---|---|
| Compile time | A D struct with an evaluation type of your own; the searches are templates over it ([Writing a problem](#writing-a-problem)) | The project is in D and knows its problem when it is built. Calls inline, and evaluations keep your type. | [`docs/example`](example/source/app.d) |
| The class | A subclass of `DynamicProblem`, run through plain functions ([Dynamic problems](#dynamic-problems)) | The problem is chosen, built or loaded at run time | [`docs/example-class`](example-class/source/app.d) |
| Python, or C | A Python class, or data and one function, over solver2's shared library ([Interfacing with Python](#interfacing-with-python)) | The evaluation is written in another language | [`bindings/python/examples/beam.py`](../bindings/python/examples/beam.py) |

solver2 is the next version of
[solver](https://github.com/x21-dev/solver), the optimiser extracted from
pleaser, and adds the second and third routes to it. A project on solver
needs no change; [Migrating from solver](#migrating-from-solver) says
what moving one takes.

This manual describes solver2 `main` as of 2026-11-13. Everything here
was checked against that code. The three worked examples build and run
against it, and CI runs them, with every test, on Linux and Windows.
[Gaps and caveats](#gaps-and-caveats) lists what solver2 doesn't do yet
and the traps in the current API. Read it before you rely on a feature.

Contents:

1. [Adding solver2 to a project](#adding-solver2-to-a-project)
2. [The model: genes, evaluations, limits, ranking](#the-model)
3. [Writing a problem](#writing-a-problem)
4. [A CMA-ES stage](#a-cma-es-stage)
5. [The archive](#the-archive)
6. [A learned model](#a-learned-model)
7. [Stopping, saving and resuming](#stopping-saving-and-resuming)
8. [Interfacing with Python](#interfacing-with-python)
9. [Lower-level pieces](#lower-level-pieces)
10. [Determinism](#determinism)
11. [Performance](#performance)
12. [Migrating from solver](#migrating-from-solver)
13. [Gaps and caveats](#gaps-and-caveats)
14. [Module reference](#module-reference)
15. [References](#references)

---

## Adding solver2 to a project

This section is for a D project. From Python, see
[Setting up](#setting-up).

solver2 is a dub **source library** (`"targetType": "sourceLibrary"`). Your
project compiles solver2's modules together with its own, so calls into
solver2 can inline like calls within one package. Add it as a dependency:

```json
"dependencies": {
    "solver2": { "path": "../solver2" }
}
```

Use a path dependency (a checkout beside your project, at `../solver2`)
or a git dependency with a pinned commit. solver2 isn't on the dub
registry.

Requirements:

- **LDC only.** `solver2.ensemble` imports `ldc.simd`, so DMD and GDC can't
  build any module that imports it: `solver2.ensemble`, `solver2.model`, or
  anything that imports those. solver2 is developed with LDC 1.42.0.
- No other dependencies. solver2 uses Phobos and `ldc.simd` only.
- solver2 has no `package.d`. Import the modules you need by name:
  `import solver2.problem; import solver2.stage; import solver2.archive;`
  and so on.

Build flags are your project's. solver2 adds none. pleaser builds with
`--build=release-native`. If you change the target CPU, check that the
floating-point results stay identical across the machines you need
(see [Determinism](#determinism)). A release build leaves array bounds
checks out of most of solver2's code; see
[Bounds checks](#bounds-checks) to turn them on or off.

Run solver2's own tests with:

```
dub test --compiler=ldc2
```

Build and run the worked example with:

```
cd docs/example && dub run --compiler=ldc2 --build=release
```

`docs/example-class` is the same example on the class route, and runs the
same way.

---

## The model

### Gene space

Every search works on `dimension` genes. Each gene is meant to lie in
[0, 1]. Your problem maps genes to physical values (millimetres, Hz, a
count) itself. `GeneSpec` helps with that (see [Genes](#genes-genespec)).

CMA-ES samples a normal distribution, so its candidates can fall outside
[0, 1]. **`evaluate` gets genes exactly as CMA-ES sampled them, which can
be outside [0, 1].** Your problem must clamp them (`clamp01`, or
`GeneSpec.toPhysical`, which clamps). The search ranks a candidate with
genes out of range lower, by `Ranking.boundPenalty` times `boundExcess(x)`
(the sum of the squared overshoots). Everything after the evaluation gets
the genes clamped to [0, 1]: the best, the archive's elites, the model,
and the progress reports.

### Evaluations

You choose the evaluation type, `P.Eval`. solver2 never looks inside it.
It reads it only through your problem's members: `scored`, `objective`,
`breach` and `descriptors`. `P.Eval` must be a **value type**: a copy must
not alias anything your problem changes later. The archive copies an elite
with `copyOf`, which duplicates arrays one level deep. So a struct of
scalars and arrays of scalars is safe, but an array of arrays, a class
reference or a pointer is shared between copies.

### Limits: three tiers

`Breach` holds how far an evaluation is over its limits in each tier. Each
field is 0 within the limits and never negative.

| Tier | Meaning | How it ranks |
|---|---|---|
| `hard` | Structure and survival | Any hard breach ranks below every candidate without one |
| `terminal` | Conditions on the outcome | Ranks after hard, before the objective |
| `soft` | Limits that may give, at a cost | Subtracted from the objective: `value = objective − soft` |

A candidate is **feasible** when `hard == 0 && terminal == 0`. Only
feasible candidates enter the archive or become a stage's best. Soft limits
never make a candidate infeasible. They only lower its value.

Your problem computes `Breach` itself, adding each tier up in the order it
chooses, so it controls the bits. If it has no order of its own to keep,
`breachOf(specs, violations)` adds per-limit violations in spec order, with
soft ones weighted by `ConstraintSpec.weight`.

### Ranking

CMA-ES minimises a fitness. `Ranking.fitness(breach, value, excess)` gives:

```
hard > 0      →  hardOffset + hard
terminal > 0  →  terminalOffset + terminal
otherwise     →  −value
              +  boundPenalty × excess
```

The defaults are `hardOffset = 100`, `terminalOffset = 50` and
`boundPenalty = 10`. They are pleaser's settings, for values from 0 to 28.

The archive ranks an emitter's candidates the same way, but in place of
the value it uses the CMA-MAE improvement (the value minus its cell's
threshold) or, for an explorer heading for a target, minus the distance
to it.

#### Setting the offsets

**The offsets must rank every breach below every candidate within the
limits**: a fitness of `terminalOffset + terminal` must be larger than any
`−value` (CMA-ES minimises). Nothing checks this. With offsets too small
for your objective, the search runs and quietly prefers candidates over
the limits: in the archive, a candidate within the limits but 200 below its
cell's threshold gets fitness 200, and one over a hard limit by 1 gets
101 and ranks better.

Work them out from your objective, where *value* is what the ranking
reads: the objective less the soft tier's costs (see
[Limits: three tiers](#limits-three-tiers)).

| | A stage | The archive |
|---|---|---|
| The most negative quantity ranked | the lowest value, *v*<sub>min</sub> | an improvement, at least *v*<sub>min</sub> − *v*<sub>max</sub> (a cell's threshold can reach the highest value); an explorer's −distance, at least −*d*<sub>max</sub> |
| `terminalOffset` must exceed | −*v*<sub>min</sub> | *v*<sub>max</sub> − *v*<sub>min</sub>, and *d*<sub>max</sub> |
| `hardOffset` must exceed | `terminalOffset` + the largest terminal breach | the same |

- **Leave a margin**, about twice the bound, since the extremes are hard
  to know in advance. A larger offset costs nothing: it only has to be
  above the values, not close to them.
- **A problem without a terminal tier** still needs `terminalOffset`
  above the values. It is the line every breach is ranked from.
- **The breaches are in your own units**: `hard` and `terminal` are the
  sums of your violations (`breachOf`). A terminal breach that can grow
  without bound needs `hardOffset` far above `terminalOffset`, or a capped
  violation.
- **`boundPenalty` scales with the objective too.** It is the cost per
  squared gene unit outside [0, 1]: pleaser's 10 makes a gene 0.1 outside
  cost 0.1, about a useful gain on its scale. For an objective 1000 times
  smaller, the default 10 would outweigh every difference in value, and
  CMA-ES would chase the bounds instead. Set it so a gene 0.1 outside costs
  about one useful step of your objective: 100 times that step.
- **The archive's `minScore`** (a template argument of `Archiving`, 0 by
  default, or the problem's own `minScore` member) is where every cell's
  threshold starts: keep it at or below *v*<sub>min</sub>. See
  [Running](#running).

Example: an objective from 0 to 5000 in the archive, a terminal breach of
at most 100, explorer distances below 10, and a useful step of about 20:

```d
Ranking r = {terminalOffset: 10_000,     // 2 × (5000 − 0), above the distances too
             hardOffset: 20_000,         // well above 10_000 + 100
             boundPenalty: 2_000};       // 100 × 20
```

pleaser's defaults pass these tests: its values run from 0 to about 28, so
its stages need `terminalOffset` above 0 and its archive above 28 (and its
explorers' distances, a few cells), and its terminal breaches stay small.

**Where to set it**: `StageOptions.ranking` for a stage, and
`ArchiveOptions.ranking` for the archive (in an options type of your own,
a `Ranking ranking` field; see [Options](#options)).

### Scoring

`evaluate` takes a `Scoring`. With `Scoring.feasibleOnly`, a candidate
over a hard limit may skip computing its objective: the ranking never reads
it. Report that through `scored(e) == false`. Both searches always call
`evaluate` with `feasibleOnly`. `Scoring.all` is there for your own calls,
such as a model's training data.

---

## Writing a problem

A problem is a compile-time concept: the searches are templates over
`P`, so every call into your problem is direct and can inline.
`isProblem!P` checks the required members. A problem that must be chosen
at run time is a class instead (see
[Dynamic problems](#dynamic-problems)), and one written in Python a
Python class (see [Interfacing with Python](#interfacing-with-python));
both run through the same templates.

### Required members

```d
struct MyProblem
{
    struct Eval { ... }            // or alias Eval = SomeType;

    size_t dimension() const;      // the number of genes
    Eval[] evaluate(const(double[])[] genes, Scoring scoring);
    bool   scored(const ref Eval e) const;      // objective computed?
    double objective(const ref Eval e) const;   // higher is better; read only when scored
    Breach breach(const ref Eval e) const;      // per tier, 0 within
}
static assert(isProblem!MyProblem);
```

`evaluate` gets a whole batch (a stage's population, or every emitter's
population in an archive iteration) and returns results **in input
order**. It may run the batch in parallel. It must be deterministic: the
same genes must give the same bits whatever the batch size, whatever the
other candidates in the batch, and on any thread count. It is the only
non-`const` member the searches call. The problem is passed by `ref`, so
it may count calls or cache.

### Optional members

The searches detect each one with `__traits(hasMember)` or a concept
template. This table lists only the members that the code reads.

| Member | Read by | Purpose |
|---|---|---|
| `size_t descriptorCount() const`, `void descriptors(const ref Eval, double[] d) const` | Archive (required there: `hasDescriptors!P`) | Where an evaluation lands. The grid maps these values to a cell. Read only for feasible candidates. |
| `alias Causes`, `Causes newCauses() const`, `void countCauses(ref Causes, const ref Eval) const` | Stage and archive (`hasCauses!P`, in `solver2.stage`) | Counts why candidates broke their limits, for progress reports. `Causes` needs `opOpAssign!"+"`. |
| `double screenPenalty(const(double)[] genes)` | Pre-screening in the stage and the archive | A known cost subtracted from a predicted objective (pleaser's processing cost, which needs no render). |
| `double minScore() const` | Archive (`hasMinScore!P`) | Where every cell's threshold starts, for a problem that knows it only when it runs. Without it, `Archiving`'s template argument. See [Running](#running). |
| `string noFeasibleSeed() const` | Archive | A hint appended to `NoFeasibleSeed`'s message. |
| `void adoptCauses(ref Causes)` | Archive, on resume | Lets the problem fix up resumed cause counts. |

`solver2.problem` also documents the surrogate members `featureCount`,
`features`, `targetCount`, `targets` and `predicted` (`hasSurrogate!P`):
with them a problem is its own learning spec for the learned model (see
[A learned model](#a-learned-model)). And it documents `fingerprint` and
`saveEval`/`loadEval` (`isSaveable!P`), a convention for a project's own
save format: no search reads them, and saving is up to you (see
[To a file](#to-a-file)).

### Genes: `GeneSpec`

```d
static immutable GeneSpec[2] spec = [GeneSpec("width", 1, 10, 5),
                                     GeneSpec("rate", 0.1, 100, 1, Scale.log)];
immutable double w = spec[0].toPhysical(g[0]);   // clamps, then maps
immutable double x0 = spec[0].toGene(spec[0].neutral);
```

`GeneSpec(name, lo, hi, neutral, scale)`:

- `Scale.linear` maps 0 to `lo` and 1 to `hi`.
- `Scale.log` needs `lo > 0`.
- `Scale.integer` rounds to the nearest whole value in [lo, hi].
- `fixed` is true when `lo == hi`. Leave a fixed gene out of the search.
  `toGene` asserts on one.

`toPhysical` clamps its gene to [0, 1], so a problem that maps every gene
through it is safe from out-of-range genes.

### Limits: `ConstraintSpec` and `breachOf`

```d
static immutable ConstraintSpec[3] limits = [ConstraintSpec("load", Tier.hard),
                                             ConstraintSpec("min_width", Tier.terminal),
                                             ConstraintSpec("rate_cap", Tier.soft, 0.05)];
immutable double[3] v = [w * r > 200 ? w * r - 200 : 0, w < 2 ? 2 - w : 0, r > 40 ? r - 40 : 0];
e.breach = breachOf(limits[], v[]);
```

`ConstraintSpec.group` is a free label for your own cause counts. solver2
doesn't read it.

### The worked example's problem

[`docs/example/source/app.d`](example/source/app.d) has a two-gene problem,
`Beam`. It has one limit in each tier, cause counts, descriptors, a
fingerprint and a `noFeasibleSeed` hint. Start from a copy of it. Its
archive starts its thresholds at −50, below the lowest value a beam
within its limits can have (−48.25, at width 2 and rate 100).

[`docs/example-class`](example-class/source/app.d) is the same problem as
a class, and
[`bindings/python/examples/beam.py`](../bindings/python/examples/beam.py)
the same in Python.

### Dynamic problems

`DynamicProblem` is an abstract class: a problem that need not be known
at compile time, such as one chosen or loaded at run time, a binding to
another language, or a quick test. A search instantiated over
`DynamicProblem` itself runs every subclass, and `solver2.dynamic` has
them instantiated, as plain functions: `runStage` and `runArchive` (see
[The entry points](#the-entry-points) below). Each call is virtual, and
there is one `evaluateInto` per batch, so the dispatch costs nothing next
to an evaluation. A subclass gets everything a compile-time problem gets,
and the same bits: solver2's tests hold a subclass, used only through
the base class, to the bits of the same problem written as a struct.

Its `Eval` is `Outcome`: `scored`, `objective`, `violations` (one per
limit, 0 within), `breach`, `descriptors` and `targets`.

**A subclass must give** four members:

```d
final class Field : DynamicProblem
{
    override const(GeneSpec)[] genes() const { return spec[]; }
    override const(ConstraintSpec)[] constraints() const { return limits[]; }
    override const(string[2])[] fingerprint() const { return [["problem", "field 1"]]; }

    override void evaluateInto(const(double[])[] genes, Scoring scoring, Outcome[] outcomes)
    {
        foreach (k, g; genes)                    // genes as sampled: toPhysical clamps
        {
            immutable double x = spec[0].toPhysical(g[0]), y = spec[1].toPhysical(g[1]);
            outcomes[k].violations = [x + y > 1.9 ? x + y - 1.9 : 0];
            outcomes[k].breach = breachOf(limits[], outcomes[k].violations);
            outcomes[k].scored = scoring == Scoring.all || outcomes[k].breach.hard == 0;
            if (outcomes[k].scored)
                outcomes[k].objective = 10 - (x - 1) * (x - 1) - (y - 0.8) * (y - 0.8);
            outcomes[k].descriptors = [x, y];
        }
    }
}
```

Each outcome needs its `breach`; when it is within the hard and terminal
limits, also `scored`, its `objective` and its `descriptors`; and its
`violations` if you want cause counts.

**A subclass may override** the rest. The searches look for a problem's
optional members on its static type, which here is the base class, so
each one is declared there with a default that behaves as its absence
does:

| Member | Default | Purpose |
|---|---|---|
| `const(Axis)[] axes() const` | none | The archive's axes, one per descriptor. `grid()` is the `BinnedGrid` over them. |
| `size_t descriptorCount() const` | `axes.length` | Override it for a grid of another kind. |
| `Ranking ranking() const` | `Ranking.init` | Offsets that fit the problem's values (see [Setting the offsets](#setting-the-offsets)). The entry points rank by it. |
| `double minScore() const` | 0 | The value's lower bound: where the archive's thresholds start, and the zero of its QD score. Finite, and at or below the lowest value of a candidate within the limits. |
| `double screenPenalty(const(double)[] genes)` | 0 | A known cost subtracted from a predicted objective. |
| `string noFeasibleSeed() const` | null | The hint in `NoFeasibleSeed`'s message. |
| `void adoptCauses(ref CauseCounts) const` | sizes them | A resumed archive's cause counts, when they weren't restored. |
| `featureCount`, `features`, `targetCount`, `targets`, `predicted` | see below | The learning spec for `LearnedModel!DynamicProblem`. |
| `double predictionHeadroom() const` | 0 | How far over a limit a prediction may be before the default `predicted` counts it as a breach. |

Remember `override`: a member with one of these names and another
signature doesn't replace the default.

**What a caller needs comes from the object**, so the caller can run a
problem it knows nothing else about, with no template argument at the
call site:

```d
import solver2.dynamic;

DynamicProblem p = load(...);                   // any subclass
checkProblem(p);                                // see below

DynamicBest best;
StageOptions so = {lambda: 8, sigma: 0.2, budget: 800, seed: 42};
runStage(p, best, null, so);                    // from p.startGenes, ranked by p.ranking

ArchiveOptions ao = {evals: 4000, seed: 9};
auto run = runArchive(p, [best.genes], ao);     // p's grid, ranking and minScore
```

`startGenes()` is each gene's `neutral` value, or the middle of its range
where it has none. [`docs/example-class`](example-class/source/app.d) is
the worked example on this route: a pre-screened stage, the archive with
a model, saved, loaded and resumed.

#### The entry points

`solver2.dynamic` has the two searches over the class as plain functions.
They are `cmaStage` and `Archiving.runArchive` with the same parameters in
the same order, less the grid, and they give the same bits as the
templates called with what the object says of itself:

```d
StageEnd runStage(DynamicProblem p, ref DynamicBest best, const(double)[] x0, StageOptions o,
                  StageModel!Outcome model = null,
                  scope void delegate(const ref DynamicStageGeneration) generationDone = null,
                  StageRun* run = null,
                  scope bool delegate(const ref StageRun) afterGeneration = null);

DynamicRun runArchive(DynamicProblem p, const(double[])[] seeds, ArchiveOptions o,
                      scope void delegate(const ref DynamicArchiveGeneration) progress = null,
                      DynamicRun* resume = null,
                      scope bool delegate(const ref DynamicRun) afterIteration = null,
                      scope void delegate(const ref DynamicBatch) evaluated = null,
                      ArchiveSurrogate!Outcome model = null);
```

What they take from the object:

| | `runStage` | `runArchive` |
|---|---|---|
| Where to start | `p.startGenes`, when `x0` is empty | `[p.startGenes]`, when `seeds` is empty and the run isn't resumed |
| The ranking | `p.ranking` | `p.ranking` |
| The grid | | `BinnedGrid` over `p.axes` |
| Where the thresholds start | | `p.minScore` |

- **`o.ranking` is replaced by `p.ranking`**, whatever it held: a problem
  met at run time knows its own values. To rank it otherwise, override
  its `ranking`, or call the templates. A run records the ranking it
  used in its options, so it resumes through the same entry point.
- **They check what they are given** and throw an `Exception` naming the
  class and what is wrong: a start or a seed with another number of genes
  than the problem, or a gene that isn't finite (outside [0, 1] is fine);
  for the archive, a problem without axes or with another number of them
  than descriptors, and a resumed archive with another number of cells
  than the grid.
- **The types** are aliases, so nothing is instantiated at the call site:
  `DynamicBest` (`Best!Outcome`), `DynamicStageGeneration`,
  `DynamicArchiving` (`Archiving!(DynamicProblem, BinnedGrid)`),
  `DynamicRun`, `DynamicArchiveGeneration`, `DynamicBatch`, and
  `DynamicModel` (`LearnedModel!DynamicProblem`):

  ```d
  auto model = new DynamicModel(p, 3, pool);    // on the problem's learning spec
  StageOptions so = {lambda: 8, sigma: 0.2, budget: 800, seed: 42, surrogatePool: 4};
  runStage(p, best, null, so, model);
  ```

- **Stopping, saving and resuming** are the templates': a `StageRun` for a
  stage; for an archive `runState`, `restoreRun`, the run's archive and
  its causes (see
  [Stopping, saving and resuming](#stopping-saving-and-resuming)). A
  resumed archive must still have its `minScore`, which is a field of the
  archive, and the problem must still give the same one.
- A problem whose cells aren't a `BinnedGrid` uses the template with its
  own grid: `Archiving!(DynamicProblem, MyGrid)`. The archive still reads
  its `minScore` from the object.

**Cause counts** come for free: `Causes` is `CauseCounts`, whose `over[i]`
counts the candidates over limit `i` of `constraints`, from each
outcome's `violations`.

**The learning spec's defaults** learn from the genes: the value (the
objective less the soft breach, which is what a prediction is ranked
by), the hard breach, the terminal breach, and the descriptors. Set
`predictionHeadroom` to about the model's error on a breach, in the
breaches' units; at 0 a prediction a hair over a limit counts as over it.
Override `targets` and `predicted` together to learn something else
(`Outcome.targets` is there to carry it). `features` and `predicted` are
called from several threads at once.

**Checks at run time.** The compiler checks a compile-time problem; a
class is checked when it runs:

- `evaluate` (final; it calls your `evaluateInto`) checks every batch:
  no breach negative or NaN, as many violations as limits, and for a
  candidate within the limits an objective and as many descriptors as
  the problem has. It throws an `Exception` naming the class, the
  candidate and what it lacks.
- `checkProblem(p, rows = 8, seed = 1)` checks a problem once, when you
  meet it. Its description (`checkDescription(p)` checks that half
  alone, with no evaluation): finite gene ranges, neutral values within
  them, axes that fit `BinnedGrid`, a finite ranking and `minScore`.
  And its evaluation, on the start genes, both corners and `rows` random
  rows (some outside [0, 1]): the same bits when evaluated twice, in a
  batch in another order, and alone, and `Scoring.feasibleOnly` changing
  no breach and no objective it computes. For a problem with
  descriptors, no row within the limits may have a value below
  `minScore`. It costs about 5 × (`rows` + 3) evaluations and throws an
  `Exception` saying what is wrong. It can't prove a problem
  deterministic, or its `minScore` low enough, only catch one that isn't
  on these rows.

---

## A CMA-ES stage

```d
import solver2.stage;

Beam p;
Best!(Beam.Eval) best;                          // carried across stages
double[] x0 = [...];                            // start genes
StageOptions so = {lambda: 8, sigma: 0.2, budget: 800, seed: 42};
const end = cmaStage(p, best, x0, so, null, (const ref StageGeneration!Beam g) {
    // after each generation
});
// best.eval, best.genes (clamped), best.feasible, best.evals
```

`cmaStage(ref P p, ref Best!(P.Eval) best, const(double)[] x0, const ref
StageOptions o, M screenModel, generationDone = null, StageRun* run = null,
afterGeneration = null)` runs CMA-ES from `x0` for `o.budget` evaluations,
rounded up to whole generations. `screenModel` is `null` or a
`StageModel!(P.Eval)`. `run` and `afterGeneration` stop and resume it (see
[Stopping and resuming a stage](#stopping-and-resuming-a-stage)).

**`StageOptions`:**

| Field | Default | Meaning |
|---|---|---|
| `lambda` | none | Population size. Must be > 0 (`cmaStage` throws otherwise). |
| `sigma` | none | Initial step size, in gene units |
| `budget` | none | Evaluations, rounded up to whole generations |
| `seed` | none | CMA-ES and pre-screening seed |
| `stallGenerations`, `stallTolerance` | 100, 0.01 | Stop when the best value rose by less than the tolerance over the last N generations (0 = never) |
| `surrogatePool` | 0 | With a model: sample this many times `lambda`, evaluate `lambda`. Must be ≥ 1 when you pass a model (`cmaStage` throws otherwise). |
| `surrogateShare` | 0.5 | Share of the evaluated candidates that the model picks (the rest are random) |
| `surrogateFrom` | 5 | Generation where pre-screening starts |
| `surrogateEpochs` | 10 | Epochs trained before each pre-screened generation |
| `ranking` | `Ranking.init` | See [Ranking](#ranking) |

**`Best!E`** is the best feasible candidate so far. Start a run with
`Best.init`, or set `genes` and `eval` to your own start. Pass the same
`Best` to later stages: each stage updates it only with a better feasible
candidate. Its fields:

- `feasible` stays false until some candidate is feasible.
- `evals` counts the evaluations of every stage.

To chain stages, pass `best.genes` as the next stage's `x0`, perhaps with a
smaller `sigma`.

**`StageGeneration!P`** is passed to the callback after each generation,
with these fields:

- `generation` counts from 1.
- `before` is the number of evaluations before this generation.
- `genes` and `evals` are this generation's candidates, genes clamped.
- `best` is the best value so far, −∞ until a candidate is feasible.
- `populationBest` is this generation's best value, NaN if none is
  feasible.
- `infeasible`, `sigma`.
- `causes`, when `hasCauses!P`.
- The pre-screening fields `picked`, `byModel`, `sampled` and
  `predictedFeasible`. `picked` is null when a generation isn't
  pre-screened.

**`StageEnd`** returns `generations`, `spent`, `stalled`, `rise` (the
best's rise over the window it stalled in) and `stopped` (by
`afterGeneration`).

**Pre-screening.** From generation `surrogateFrom` on, the stage trains the
model, samples `surrogatePool × lambda` candidates and ranks them on the
model's prediction (minus `screenPenalty`). It evaluates `lambda` of them:
`surrogateShare` picked by the model and the rest at random. CMA-ES is then
told the real fitness of the evaluated candidates and the predicted fitness
of the others. A `StageModel!E` has three members:

```d
interface StageModel(E)
{
    void learn(const(double[])[] genes, const(E)[] evals);
    void train(size_t epochs);
    void predictScore(const(double[])[] genes, double[] score, double[] violation);
}
```

`violation` is ranked as a hard breach. `LearnedModel` implements this
interface (see [A learned model](#a-learned-model)).

### Stopping and resuming a stage

Pass a `StageRun` and the stage keeps its state there after every
generation. `afterGeneration` is asked after each one whether to stop; the
stage then returns with `end.stopped` set. Pass the same `StageRun` back,
with the same `Best` and model, and the stage goes on from there:

```d
StageRun run;
auto end = cmaStage(p, best, x0, so, model, null, &run, (const ref StageRun r) => r.generation == 20);
// ... save run, best and model.save(...), load them ...
end = cmaStage(p, best, x0, so, model, null, &run);    // generations 21 on
```

- A stage stopped and resumed ends exactly as the stage that never
  stopped, bit for bit (on one platform; see [Determinism](#determinism)).
- **Adding evaluations to a finished stage**: keep its `StageRun` (it is
  filled at the end too) and resume it with a larger `budget`. It ends as
  the stage given that budget from the start. The budget is always the
  whole stage's.
- **A stage that stalled** stays stopped when resumed with the same stall
  settings: the test is applied to its saved history at once. To go on,
  resume with a longer `stallGenerations` (or a smaller `stallTolerance`):
  the stage then ends as one given those settings from the start, since a
  longer window or a smaller tolerance never stops a stage sooner.
- **Checks**: a resumed stage throws unless its options equal the saved
  ones (all but `budget`, `stallGenerations` and `stallTolerance`), `x0`
  has the saved number of genes (its values are ignored), and a model is
  given once its pre-screening has begun.
- **`StageRun`** is plain data: the options, the generations and
  evaluations done, the best value after each generation (`history`, for
  the stall test), CMA-ES's `CmaEs.State` (with the larger population once
  pre-screening has begun) and the pre-screening generator's state. The
  caller's `Best` and model aren't in it: save them alongside (see
  [To a file](#to-a-file)).

---

## The archive

The archive is CMA-MAE [1]: CMA-ES [4, 5] emitters fill a MAP-Elites
archive [6], the best candidate found in each cell of a grid over the
descriptors, with candidates ranked by how far they beat their cell's
threshold. CMA-MAE builds on CMA-ME [2], which it becomes at `alpha` 1.

**What solver2 adds: going after cells not discovered yet.** CMA-MAE
explores indirectly. A new cell scores high once a candidate lands in it
(its improvement counts from `minScore`), and as thresholds rise in
well-visited cells the emitters drift to less-visited ones: once the
thresholds have caught up, or on a flat objective, CMA-MAE is a descent
on the archive's density histogram [1, Theorem 5.4]. Both depend on
candidates landing near the edge of what has been found; nothing points
an emitter at a cell that no candidate has reached. solver2's explorers
do: each takes the empty cell nearest the occupied ones as its target and
ranks its candidates by their distance to it (`grid.distance`), until it
reaches the cell or gives it up. Where the improvers refine what the
archive holds, the explorers widen it. pleaser wants every character of a
mix it can reach, so reaching new cells is the point of its runs; a
problem that only wants the best value in each cell already found can
set `explorers` to 0.

The explorers, the constraint tiers and the infeasible restarts are
solver2's additions, from pleaser (see
[Where the rules come from](#where-the-rules-come-from)).

### Grids

A grid maps descriptor values to cells. `isGrid!G` requires these members:

| Member | Purpose |
|---|---|
| `size_t cells() const` | The number of cells |
| `size_t cellOf(const(double)[] d) const` | The cell of descriptors `d` |
| `double distance(const(double)[] d, size_t cell) const` | How far `d` is from `cell`: 0 inside, in bin widths. Explorers move down it. |
| `long binDistance2(size_t a, size_t b) const` | Squared distance between cells in bins. Used for explorer starts and `selectDiverse`. |
| `string code(size_t) const`, `size_t parse(string) const` | A cell's name in files and on a command line |

`BinnedGrid` is the built-in grid. It's a product of `Axis` values, with
the first axis varying fastest in the cell index:

```d
static immutable Axis[2] axes = [Axis("width", [3, 5, 7], "nmwW", 0, 2),
                                 Axis("rate", [1, 10, 40], "slfF", 0, 10)];
const grid = BinnedGrid(axes[]);      // 4 × 4 = 16 cells, codes like "wf"
```

`Axis(name, edges, letters, nanBin, width)` has these parts:

- `edges` are ascending. A value below `edges[0]` goes in bin 0, and a
  value from the last edge up goes in the last bin.
- `letters` has one character per bin (`edges.length + 1`).
- `nanBin` is the bin a NaN value goes in.
- `width` is the unit of `distance` on this axis. Set it to a typical bin
  width so that the axes weigh alike.

If your cells aren't a regular product of bins, bring your own `G`.
pleaser's character grid is an example.

### Running

```d
import solver2.archive;

alias A = Archiving!(Beam, BinnedGrid, ArchiveOptions, 0, -50);     // P, G, O, N, minScore
ArchiveOptions ao = {evals: 4000, emitters: 4, explorers: 2, lambda: 8, seed: 9};
A.Run run = A.runArchive(p, grid, seeds, ao);
```

`Archiving!(P, G, O, N, minScore)` is a template of functions; `O`
defaults to `ArchiveOptions`, and `N` and `minScore` to 0.
`P` must satisfy `hasDescriptors!P` and `G` must satisfy `isGrid!G`. It
has these members:

- `Run` is `ArchiveRun!(P, O, N, minScore)`.
- `Generation` is `ArchiveGeneration!(CausesOf!P)`.
- `Batch` is `ArchiveBatch!(P.Eval)`.
- `Model` is `ArchiveSurrogate!(P.Eval)`.
- `runArchive`:

```d
Run runArchive(ref P p, const ref G grid, const(double[])[] seeds, O opts,
               scope void delegate(const ref Generation) progress = null,
               Run* resume = null,
               scope bool delegate(const ref Run) afterIteration = null,
               scope void delegate(const ref Batch) evaluated = null,
               Model model = null);
```

A run goes like this:

1. The seeds are evaluated, and the feasible ones enter the archive. If
   none is feasible, the run throws `NoFeasibleSeed` (with your
   `noFeasibleSeed()` hint).
2. All emitters start from the best seed. The last `explorers` of them are
   explorers.
3. Each iteration, every emitter asks for `lambda` candidates. All
   `emitters × lambda` are evaluated in **one** `evaluate` call. Then the
   candidates are offered to the archive emitter by emitter, in population
   order.
4. Improvers rank candidates by CMA-MAE improvement: value minus the cell's
   threshold. An improver restarts from a random elite when it has been
   idle, when nearly all its candidates are infeasible, or when its step
   size collapses (see [Restart and explorer rules](#restart-and-explorer-rules)).
5. An explorer heads for the empty cell nearest the occupied ones, down
   `grid.distance`. It retargets when it reaches the cell or gets stuck,
   and a cell it gives up on is blocked for a while.
6. The run stops once `opts.evals` evaluations beyond the seeds are spent
   (in whole iterations), or when `afterIteration` returns true.

**Template parameters `N` and `minScore`.**

- `N = 0` (the default) holds the cells in arrays sized by `grid.cells`.
- `N > 0` fixes the cell count at compile time, and the archive holds
  `Elite!E[N]` by value. That makes `Run` large; read
  [Performance](#performance) before you use it.
- `minScore` is where every cell's threshold starts, and the zero of
  `qdScore`. Set it at or below the lowest feasible value: a candidate
  whose value is below it can't beat its cell's threshold. The worked
  example's −50 is below its lowest value, −48.25.

**`minScore` known only at run time.** The archive keeps it as a field,
`run.archive.minScore`, which starts as the template's argument. A
problem that learns its lower bound when it runs gives it as a member
instead, and leaves the template's argument at 0 (the build fails
otherwise):

```d
double minScore() const { return lowest; }      // in the problem; must be finite
```

`runArchive` then starts every threshold there, with `N > 0` too.
`DynamicProblem` has the member, 0 unless a subclass overrides it. The
two routes give the same bits: a threshold or a QD score computed from
the field is the one computed from the constant. The field is part of
the archive, so it is saved and restored with the thresholds; a resumed
run on a problem with the member checks that the two still agree. An
archive of your own takes it in its constructor:
`Archive!E(cells, minScore)`.

### Options

`ArchiveOptions` fields (defaults in brackets):

| Field | Meaning |
|---|---|
| `evals` [8000] | Budget beyond the seeds, rounded up to whole iterations. On resume it is still the whole run's budget. |
| `emitters` [5] | Emitters, explorers included |
| `lambda` [0] | Population per emitter. 0 uses CMA-ES's default for the dimension (4 + ⌊3 ln n⌋). |
| `sigma` [0.05] | Improvers' initial step size |
| `alpha` [0.1] | Archive learning rate. 0 optimises value alone; 1 is CMA-ME. See [The archive learning rate](#the-archive-learning-rate-alpha). |
| `seed` [1] | Run seed |
| `explorers` [2] | How many of the emitters are explorers (the last ones) |
| `explorerSigma` [0.1] | Explorers' step size |
| `surrogatePool` [0] | With a model, model-guided explorers sample this many times their population. 0 = off. **A model must be passed if and only if this is > 0.** |
| `surrogateShare` [0.5], `surrogateFrom` [5], `surrogateEpochs` [10] | As in the stage. `surrogateFrom` counts iterations. |
| `surrogateImprovers` [false] | Pre-screen the improvers too |
| `surrogateAim` [false] | Start a retargeted explorer where the model predicts nearest its target |
| `ranking` [`Ranking.init`] | How candidates rank. Change it when your values or distances reach its offsets: see [Setting the offsets](#setting-the-offsets). |
| `rules` [`ArchiveRules.init`] | When emitters restart and explorers give a target up: see [Restart and explorer rules](#restart-and-explorer-rules). |

You may pass your own options type `O` instead. It needs every field above
but `ranking` and `rules`, with the same names and types, and may have
more. Without `ranking` or `rules`, the archive uses their defaults. `O`
must support `==`, because a resumed run checks that its options equal the
saved ones (all but `evals`): a run keeps its ranking and rules when
resumed.

### The archive learning rate (`alpha`)

Each candidate that beats its cell's threshold moves the threshold
towards its value by `alpha` of the gap. After about 1 / `alpha` such
candidates a cell's threshold has nearly caught up with what the cell
holds, and the cell stops drawing improvers away from the rest of the
archive. At 0 the thresholds stay at `minScore` and the emitters optimise
value alone; at 1 a threshold is its cell's best value, CMA-ME [2].

- **The paper's values** [1]: 0.01 in its main experiments (a 100 × 100
  archive), 0.1 and 0.02 in its image-generation ones. It reports CMA-MAE
  robust to the exact value.
- **Another resolution** [1, appendix G]: for archives of *M*<sub>1</sub>
  and *M*<sub>2</sub> cells, α<sub>2</sub> = 1 − (1 − α<sub>1</sub>)<sup>*M*<sub>2</sub>/*M*<sub>1</sub></sup>
  keeps the annealing of each cell the same for the same number of
  candidates offered. A coarser archive needs a smaller `alpha`, since each
  of its cells gets more of them. The derivation assumes they spread
  evenly over the cells and a constant objective; the paper conjectures it
  holds for convex ones.
- **Your run's length** (solver2's rule of thumb, not the paper's): a cell
  can be offered at most about `evals` ÷ cells candidates in a run. Keep
  1 / `alpha` a fraction of that, or the thresholds hardly move and the
  run behaves as at 0. pleaser's 0.1 with 75 cells and 8000 evaluations:
  about 100 candidates per cell, thresholds caught up after about 10. The
  conversion above would give a far smaller `alpha` for 75 cells than the
  paper's 0.01 for 10,000, but the paper's runs are hundreds of times
  longer.

### Restart and explorer rules

An emitter runs CMA-ES from one start until a rule says it has stopped
paying off; it then starts again elsewhere (an improver from a random
elite, at `sigma`; an explorer at a new target, from the elite nearest
it). `ArchiveOptions.rules`, an `ArchiveRules`, holds those rules. Each
trades depth for breadth: a patient emitter refines where it is, an
impatient one covers more of the archive. The defaults are pleaser's: a
few hundred genes, λ at CMA-ES's default (17 to 22 there), 5 emitters,
8000 evaluations, 75 cells; a run of 80 to 95 iterations.

A *generation* is one emitter's λ candidates; an *iteration* is one
generation of every emitter. A run has about evals ÷ (emitters × λ)
iterations.

| Rule | Default | An emitter restarts or retargets when | Unit |
|---|---|---|---|
| `idleGenerations` | 5 | an improver had no candidate above its cell's threshold for this many generations in a row | generations |
| `infeasibleGenerations`, `infeasibleShare` | 5, 0.9 | any emitter's last `infeasibleGenerations` generations had at least this share of their candidates over a hard or terminal limit | generations, share of λ |
| `minSigma` | 1e-4 | any emitter's CMA-ES step size fell below this | gene units |
| `staleGenerations`, `minProgress` | 10, 0.01 | an explorer's closest candidate within the limits got no closer to its target by `minProgress` for this many generations (its first generation on a target counts as progress when it has such a candidate) | generations, grid distance |
| `blockIterations` | 30 | (after) no explorer aims at a cell given up on for this many iterations | iterations |
| `aimSamples` | 200 | (with `surrogateAim`) a retargeted explorer starts at the best, on the model, of each elite and this many points around it | points per elite |

The rules are checked when a run starts: windows at least 1, a share
above 0 and at most 1, `minSigma` and `minProgress` finite and not
negative. Under their former names, `RESTART_GENERATIONS`,
`RESTART_INFEASIBLE`, `RESTART_SIGMA`, `EXPLORE_STALE`, `EXPLORE_PROGRESS`,
`EXPLORE_BLOCK` and `AIM_SAMPLES` are the defaults.

#### Where the rules come from

- **The restart itself is CMA-MAE's**: a new CMA-ES from a random elite,
  its covariance reset [1, Algorithm 1], as in CMA-ME [2].
- **Its triggers are pleaser's.** CMA-MAE's main experiments restart an
  emitter only when CMA-ES converges [1]: pyribs' "basic" rule, which stops
  when σ·√(largest eigenvalue of C) falls below 1e-11, the condition number
  of C exceeds 1e14, or a generation's ranking is flat within 1e-12, after
  pycma [3]. CMA-ME's improvement emitters, and pyribs' default
  ("no_improvement"), restart as soon as one generation adds nothing to the
  archive [2, 3]. `idleGenerations` waits 5 such generations; `minSigma`
  stops at 1e-4, far above 1e-11, because pleaser's evaluations are
  expensive and changes below that are inaudible; the infeasible restart
  has no counterpart, since the papers have no constraint tiers.
- **The explorers are pleaser's own.** Neither paper has emitters that
  head for an empty cell, give it up or block it. CMA-MAE reaches new
  cells through its thresholds and density descent [1]; the nearest
  directed idea is CMA-ME's random direction emitters [2], which travel in
  a chosen direction of descriptor space and restart when they stop
  improving the archive. An explorer instead picks the empty cell to
  reach (see [The archive](#the-archive)).
- So the guidance below follows from how solver2's code works and from
  pleaser's runs, not from published results.

#### The papers' configuration

To run solver2 as close to CMA-MAE as it gets [1]:

```d
ArchiveOptions ao;
ao.alpha = 0.01;                      // at the paper's resolution; see the conversion above
ao.explorers = 0;                     // improvers only
ao.rules.idleGenerations = size_t.max;   // no idle restarts: CMA-ES convergence only
ao.rules.minSigma = 1e-11;            // pyribs' "basic" step-size limit
ao.rules.infeasibleGenerations = size_t.max;   // a window that never fills: no infeasible restarts
```

pyribs' default rule instead is `idleGenerations = 1`. solver2's
convergence test is the step size alone, without the condition-number and
flat-ranking tests, and a problem with limits ranks breaches as
[Ranking](#ranking) says, where the paper folds penalties into the
objective.

#### Set from your units: `minProgress` and `minSigma`

These two are measured in your problem's units, so the defaults are only
right for problems shaped like pleaser's. Set them first.

- **`minProgress`** is in `grid.distance` units: for a `BinnedGrid`, the
  axis's `width` (1 unless you set it). Make it about 1% of the distance
  across one bin. With `width` set to the bin size, distances are in bins
  and 0.01 is right. With `width` left at 1 on an axis whose bins are
  0.002 wide, 0.01 is five bins, and every explorer gives its target up
  after `staleGenerations` generations however well it closes in; with
  bins 1000 wide, any step counts as progress.
- **`minSigma`** is in gene units, [0, 1] across each gene's range. Set it
  to about a tenth of the smallest gene change that still matters. A gain
  over 30 dB that matters to 0.03 dB changes by 0.001 of its range, so
  1e-4 is right; a parameter that matters to 1% of its range can stop at
  1e-3 and save the generations spent below it.

#### Patience: `idleGenerations`, `infeasibleGenerations`, `staleGenerations`

These windows ask how many generations without progress mean CMA-ES has
stalled. That depends on how fast CMA-ES makes progress on your problem,
not on the budget:

| Your problem | Effect on CMA-ES | Change the windows |
|---|---|---|
| Many more genes than pleaser, or a small λ | fewer improvements per generation; slower to adapt its step size and shape | longer |
| Few genes, a large λ, a smooth objective | improves within a generation or two when it can | shorter |
| A noisy objective | improvements come and go by chance | `idleGenerations` longer, so a run of unlucky generations doesn't restart an emitter that improves on average |
| A narrow feasible region (most candidates over a limit near the good ones) | needs several generations to get back inside | `infeasibleShare` up to 1, or `infeasibleGenerations` longer |
| Large regions over the limits that lead nowhere | wastes whole windows in them | `infeasibleGenerations` shorter |
| Cells far apart in descriptor space, or slow to reach | an explorer closes in slowly | `staleGenerations` longer |

Keep every window well under the run's length in generations, so each
emitter restarts or retargets many times: pleaser's 5 and 10 are 5 to 12%
of its 80 to 95 iterations. A run of 20 iterations needs shorter windows,
or more evaluations.

#### Run length and cost: `blockIterations`, `aimSamples`

- **`blockIterations`** keeps explorers off cells they failed to reach.
  pleaser's 30 is about a third of its run: a cell given up on is tried
  again once the archive has grown towards it. Scale it with your
  run's iterations; shorter when the archive keeps spreading (new
  neighbours make old targets reachable), longer when the same cells fail
  again and again.
- **`aimSamples`** costs occupied cells × (`aimSamples` + 1) model
  predictions each time an explorer is aimed. More genes leave the points
  around an elite sparser: raise it with the dimension when the model is
  cheap, lower it when predictions are slow next to an evaluation.

#### Checking them

The result says whether a rule is off. Read `restartsIdle`,
`restartsInfeasible`, `restartsConverged`, `reached` and `abandoned`, and,
in the `progress` callback, each generation's `restartReason`, `improved`
and `distance`:

| You see | Change |
|---|---|
| Improvers restart idle while each restart still lifts the QD score | `idleGenerations` up |
| Most of the improvers' generations improve nothing (`improved` 0), in a short run | `idleGenerations` down |
| `restartsInfeasible` high, few elites near the limits | `infeasibleShare` up, or `infeasibleGenerations` up |
| Generations with nearly all candidates infeasible go on and on | `infeasibleGenerations` down |
| `restartsConverged` high, elites still improving | `minSigma` down |
| Long runs of generations at a tiny `sigma` with `improved` 0 | `minSigma` up |
| `abandoned` far above `reached`, and `distance` still falling when an explorer gave up | `staleGenerations` up, or `minProgress` down |
| Explorers return to the same cells and give them up again | `blockIterations` up |

Change one rule at a time, and compare runs over a few seeds rather than
one: a rule moves every later restart, so one run says little.

### Callbacks

- `progress(const ref Generation g)` is called after every emitter's
  generation. `ArchiveGeneration` has these fields:
  - `iteration`, `emitter`, `evals`
  - `cells` (coverage), `qdScore`, `bestScore`
  - `improved`, `newCells`, `sigma`
  - `restartReason` (`reasonText` gives its name)
  - `causes`, `explorer`, `target`, `distance`, `focus`
- `afterIteration(const ref Run r)` is called after each iteration except
  the last, with the run exactly as it would be saved. Return true to stop
  there. `run.stopped` is then set.
- `evaluated(const ref Batch b)` gets every batch as soon as it's
  evaluated: the seeds (`iteration == 0`, `emitters` empty), then each
  iteration. Candidate `i` comes from emitter `i / lambda`. Use it to log
  every evaluation.

### Results

`Run` (an `ArchiveRun`) holds the result and everything needed to continue:

- `archive`: an `Archive` with these members:
  - `elites[c]`: `occupied`, `genes` in [0, 1], `eval`, `found`
    (evaluation count when it entered) and `value`
  - `threshold[c]`, `alpha`
  - `coverage()`, `qdScore()`, `best()` (cell index, `size_t.max` when
    empty), `bestScore()`, `occupiedCells()`
- `evals`, `iterations`, `seedEvals`, `lambda`
- `causes`
- `restartsIdle`, `restartsInfeasible`, `restartsConverged`
- `reached` and `abandoned` explorer targets
- `seeds` (the seeds' evaluations)
- `emitters`, `blockedUntil`
- `guide` and `screen`: how the model's picks turned out
- `stopped`

`selectDiverse(archive, grid, within, count)` picks elites to show: the
best first, then repeatedly the elite farthest in bins from those already
picked, among elites at most `within` below the best.

### Focused runs

To spend the rest of a stopped run improving chosen cells:

1. Set `run.focus` to the cells and `run.focusFrom = run.iterations`.
2. Call `refocus(run, k, modelReady)` for each emitter `k`.
3. Resume the run.

Every emitter then improves, restarting from the focused elites in turn.
`focusScores` and `focusImprovers` are bookkeeping for the caller: the
archive saves them but doesn't read them.

---

## A learned model

`solver2.model.LearnedModel!L` trains a small neural-network ensemble on
every evaluation of a run, and predicts candidates before they are
evaluated. It is both a `StageModel!(L.Eval)` and an
`ArchiveSurrogate!(L.Eval)`, so you can pass it to `cmaStage` or to
`runArchive`.

### The learning spec

`L` describes what the network learns of your problem: `alias Eval`,
`dimension()`, and the surrogate members of the problem interface.
`isLearning!L` checks it. A problem with those members (`hasSurrogate!P`)
is a learning spec itself, as the worked example's `Beam` is; a separate
struct with them works too, when the features need state the problem
doesn't hold (pleaser's spec keeps its feature map apart from its
evaluator).

```d
struct Beam
{
    // ... the problem's members, then:
    size_t featureCount() const { return 2; }     // network inputs
    void features(const(double)[] genes, double[] f) const { f[] = genes[]; }
    size_t targetCount() const { return 4; }      // network outputs

    void targets(const ref Eval e, double[] t) const
    {
        t[0] = e.scored ? e.objective : double.nan;     // NaN: not learned from this row
        t[1] = e.breach.hard + e.breach.terminal;
        t[2] = e.width;
        t[3] = e.rate;
    }

    Predicted predicted(const(double)[] m) const
    {
        // 0.5 of headroom on a prediction.
        return Predicted(m[0], Breach(m[1] > 0.5 ? m[1] : 0), m[2 .. 4].dup);
    }
}
```

- `features` maps genes (clamped to [0, 1]) to the network's inputs. It's
  called **from several threads at once** when the model has a pool, so it
  must not write shared state.
- `targets` gives one value per output. NaN, or any value that isn't
  finite, counts as missing and adds nothing to the loss. An unscored
  objective is the usual case.
- `predicted` turns outputs into a `Predicted`: an objective, a `Breach`
  and the descriptors (as many as the problem's `descriptorCount`; the
  archive asks for them, a stage ignores them). It's also called from
  several threads at once. Add headroom here to make predictions err
  towards infeasible.
- **For a problem with soft limits, learn the value**, the objective less
  the soft breach: a prediction is ranked by its objective alone, with
  no soft tier to subtract (see [gap 3](#gaps-and-caveats)).
  `DynamicProblem`'s default spec learns the value. `Beam` above learns
  the objective, so its model doesn't see what its soft limit costs.

### Using it

```d
alias Model = LearnedModel!Beam;
auto model = new Model(Beam(), seed, pool);   // pool: a TaskPool, or null

// A pre-screened stage:
StageOptions so = {lambda: 8, sigma: 0.2, budget: 800, seed: 42, surrogatePool: 4};
cmaStage(p, best, x0, so, model);

// Model-guided explorers, pre-screened improvers, model-aimed explorers:
ArchiveOptions ao = {..., surrogatePool: 4, surrogateFrom: 3, surrogateEpochs: 5,
                     surrogateImprovers: true, surrogateAim: true};
A.runArchive(p, grid, seeds, ao, null, null, null, null, model);
```

The constructor standardises the inputs on the features of 4000 uniformly
random gene vectors drawn from `seed`, so `features` must accept any genes
in [0, 1]. The ensemble uses `EnsembleOptions` defaults (5 members, two
hidden layers of 64, Adam at 1e-3) with your seed, unless you pass your
own as a fourth argument (its `seed` is replaced by the model's):
`new Model(Beam(), seed, pool, EnsembleOptions(3, 32))`.

Other members:

- `rows()` is the number of evaluations learned.
- `predictEach(genes, take)` gives raw outputs per candidate, for your own
  studies.
- `trainSeconds`, `predictSeconds` and `trainings` are timings.
- `save(evals)` and `restore(bytes, evals)`; see the next section.

`columnStats(rows)` is the column mean and population standard deviation
helper the model uses.

---

## Stopping, saving and resuming

### Within one process

A stage stops and resumes through a `StageRun` (see
[Stopping and resuming a stage](#stopping-and-resuming-a-stage)). An
archive run: stop it with `afterIteration`, then pass it back as `resume`. Use the
same problem, grid, options (`evals` may be larger) and model:

```d
auto run = A.runArchive(p, grid, seeds, ao, null, null, (const ref A.Run r) => r.iterations == 10);
A.runArchive(p, grid, null, ao, null, &run);     // changes run in place, returns a copy
```

A run stopped and resumed ends exactly as the run that never stopped. A
finished run resumed with a larger `evals` ends as a run given that budget
from the start. A resumed run checks three things, and throws if any
fails:

- its options equal the saved ones (all but `evals`);
- it still has its emitters;
- a model is given when its emitters sample for one.

### To a file

solver2 gives you the state as plain data and leaves **the file format to
you**. It has no serialiser, by design (see [gap 1](#gaps-and-caveats)):
a file format is a project's own choice. This section is a template to
copy into your project and adapt instead, the code the worked example
runs. (From Python a run saves itself as bytes: see
[Stopping, saving and resuming from Python](#stopping-saving-and-resuming-from-python).)

#### What to save

| What | From | Back with |
|---|---|---|
| A fingerprint of what the run was computed from | your problem | check with `fingerprintMismatch` |
| The run's state between iterations | `runState(run)` → `RunState!O` | `restoreRun(run, state)` |
| The archive: elites, **thresholds, alpha and minScore** | `run.archive` | assign `run.archive` |
| The cause counts | `run.causes` | assign `run.causes` |
| The model, if any | `model.save(run.evals)` → `ubyte[]` | `model.restore(bytes, run.evals)` |

`RunState!O` holds every field of the run except the archive, the causes
and the seeds' evaluations. Each emitter's CMA-ES is a `CmaEs.State`, and
its RNG is a `ulong`. `restoreRun` doesn't touch `run.archive`, so restore
the archive yourself, thresholds, `alpha` and `minScore` included. A
resumed run doesn't reset `alpha` from the options. A file written
before `Archive` had its `minScore` field (by solver) lacks it; read with
the template below, the field keeps the template argument's value, which
is where that archive's thresholds started (see
[Saved states](#saved-states)). `LearnedModel.save` is versioned
and holds the training data and the ensemble. `restore` refuses bytes
from another run size or shape, and leaves the model as it was when it
does.

The **fingerprint** is name/value pairs describing what the run was
computed from, such as problem inputs, settings and versions. Check it on
load with `solver2.problem.fingerprintMismatch(saved, now)`. It returns
null when every saved pair is in `now` with the same value, in any order,
or the first that isn't as a message. Pairs only in `now` are new and
accepted, so a project can add a setting to its fingerprint without
refusing older files; when a new setting changes results, give it a
default value its older files can be read as, or check it yourself.

#### The format

A file is a header, then one value:

- **Header**: your project's magic string (so another file is refused
  outright) and a format version (`uint`).
- **A struct**: its field count, then each field as its name, its byte
  length and its bytes. The reader matches fields **by name**, so their
  order doesn't matter, and the lengths let it step over a field.
- **A static array** or **a slice**: its length, then its elements.
- **A scalar**: its bytes, through `solver2.bytes.Bytes` (native layout, so
  little-endian on little-endian machines, and a `size_t` field is as wide
  as the build's). Lengths are `ulong`.

Writing a struct straight through `Bytes.put` is wrong for any struct
holding arrays, because it writes their pointers. The template walks the
fields instead, and refuses at compile time a type holding a pointer, a
class, a delegate or an associative array.

#### The template

Copy this module into your project as it is; adapt `MAGIC`, `FORMAT` and
`Saved` (below) rather than the module. It is
[`docs/example/source/checkpoint.d`](example/source/checkpoint.d), and the
example checks at compile time that the two are the same text.

<!-- checkpoint.d -->
```d
/// A checkpoint file format to copy into a project and adapt (the developer
/// manual's "Stopping, saving and resuming → To a file"). solver2 gives a
/// run's state as plain data (RunState, the archive, the causes) and the
/// model's as bytes, and leaves the file format to the caller.
///
/// A file is a header (the project's magic string and its format version)
/// and one value. A struct is written as its named fields, each with its
/// byte length, so the reader matches fields by name: a field the file
/// lacks keeps the default the reader starts from (older files stay
/// readable when a field is added), and a field the reader doesn't know is
/// refused (a newer file isn't half read) unless it is listed as retired.
/// A static array is its elements, a slice its length then its elements, a
/// scalar its bytes; lengths are ulong. Everything goes through
/// solver2.bytes, so the bytes are little-endian on little-endian machines.
module checkpoint;

import std.exception : enforce;
import std.traits : isPointer, isAssociativeArray, isDelegate, isFunctionPointer;

import solver2.bytes : Bytes;

/// Appends the header: `magic`, then `version_`.
void writeHeader(ref Bytes o, string magic, uint version_)
{
    write(o, magic);
    o.put(version_);
}

/// Reads the header `writeHeader` wrote and returns its version. Throws
/// when the magic isn't `magic` (not this project's file) or the version
/// is above `current` (a newer program's file).
uint readHeader(ref Bytes r, string magic, uint current)
{
    string m;
    read(r, m);
    enforce(m == magic, "not a " ~ magic ~ " file");
    immutable uint v = r.get!uint;
    enforce(v <= current, "the file is of a newer format than this program reads");
    return v;
}

/// Appends `v` to `o`: a struct its fields by name, a static array element
/// by element, a slice its length then its elements, a scalar its bytes.
void write(T)(ref Bytes o, const ref T v)
{
    static assert(!isPointer!T && !isAssociativeArray!T && !isDelegate!T && !isFunctionPointer!T
                  && !is(T == class) && !is(T == interface) && !is(T == union),
                  T.stringof ~ " holds a reference; save what it refers to instead");
    static if (is(T == struct))
    {
        o.put(ulong(v.tupleof.length));
        static foreach (i; 0 .. v.tupleof.length)
        {{
            immutable string name = __traits(identifier, T.tupleof[i]);
            write(o, name);
            // The field's byte length, filled in once it is written.
            immutable size_t at = o.data.length;
            o.put(ulong(0));
            write(o, v.tupleof[i]);
            immutable ulong n = o.data.length - at - ulong.sizeof;
            o.data[at .. at + ulong.sizeof] = (cast(const(ubyte)*) &n)[0 .. ulong.sizeof];
        }}
    }
    else static if (is(T == E[n], E, size_t n))
    {
        o.put(ulong(n));
        foreach (ref e; v)
            write(o, e);
    }
    else static if (is(T == U[], U))
    {
        o.put(ulong(v.length));
        foreach (ref e; v)
            write(o, e);
    }
    else
        o.put(cast() v);
}

/// Reads into `v` what `write` wrote. A struct field the file lacks keeps
/// the value `v` has, so start from a default (`T v;`). A field in the file
/// that `T` has no longer is skipped when `retired` names it as
/// "Struct.field" (`Struct` the type's name without template arguments),
/// and throws otherwise.
void read(T)(ref Bytes r, ref T v, const(string)[] retired = null)
{
    static if (is(T == struct))
    {
        immutable ulong fields = r.get!ulong;
        foreach (_; 0 .. fields)
        {
            string name;
            read(r, name);
            immutable ulong n = r.get!ulong;
            enforce(n <= r.data.length - r.at, "checkpoint cut short");
            auto field = Bytes(r.data[r.at .. r.at + cast(size_t) n]);
            r.at += cast(size_t) n;
            bool known;
            static foreach (i; 0 .. v.tupleof.length)
                if (name == __traits(identifier, T.tupleof[i]))
                {
                    read(field, v.tupleof[i], retired);
                    enforce(field.done, "checkpoint field " ~ name ~ " of another type");
                    known = true;
                }
            if (!known)
            {
                import std.algorithm : canFind;
                enforce(retired.canFind(__traits(identifier, T) ~ "." ~ name),
                        "checkpoint field " ~ __traits(identifier, T) ~ "." ~ name ~ " unknown to this program");
            }
        }
    }
    else static if (is(T == E[n], E, size_t n))
    {
        enforce(r.get!ulong == n, "checkpoint array of another length");
        foreach (ref e; v)
            read(r, e, retired);
    }
    else static if (is(T == immutable(char)[]))
    {
        char[] s;
        read(r, s);
        v = s.idup;
    }
    else static if (is(T == U[], U))
    {
        immutable ulong n = r.get!ulong;
        enforce(n <= r.data.length - r.at, "checkpoint cut short");
        v = new U[cast(size_t) n];
        foreach (ref e; v)
            read(r, e, retired);
    }
    else
        v = r.get!T;
}

unittest
{
    import std.exception : assertThrown;

    // Values come back as they went in, whatever they nest.
    static struct S { double a; int[2] b; string s; double[][] m; bool f; }
    const S x = {1.5, [3, 4], "beam", [[1.0], [], [2.0, 3.0]], true};
    Bytes o;
    writeHeader(o, "test", 1);
    write(o, x);
    auto r = Bytes(o.data);
    assert(readHeader(r, "test", 1) == 1);
    S y;
    read(r, y);
    assert(y == x && r.done);

    // Another project's file, or a newer format's, is refused.
    auto other = Bytes(o.data);
    assertThrown(readHeader(other, "other", 1));
    Bytes newer;
    writeHeader(newer, "test", 2);
    auto rn = Bytes(newer.data);
    assertThrown(readHeader(rn, "test", 1));

    // A field added since: the older file still reads, the new field keeps
    // its default. Fields are matched by name, so their order may change.
    static struct Added { bool f; string s; double a; int[2] b; double[][] m; int count = 7; }
    Added z;
    auto ra = Bytes(o.data);
    readHeader(ra, "test", 1);
    read(ra, z);
    assert(z.a == 1.5 && z.b == [3, 4] && z.s == "beam" && z.m == x.m && z.f && z.count == 7 && ra.done);

    // A field dropped since: refused unless retired.
    static struct Dropped { double a; int[2] b; string s; bool f; }
    Dropped d;
    auto rd = Bytes(o.data);
    readHeader(rd, "test", 1);
    assertThrown(read(rd, d));
    rd = Bytes(o.data);
    readHeader(rd, "test", 1);
    read(rd, d, ["Dropped.m"]);
    assert(d.a == 1.5 && d.s == "beam" && d.f && rd.done);

    // A field whose type changed is refused.
    static struct Changed { double a; int[2] b; string s; double[][] m; double f; }
    Changed c;
    auto rc = Bytes(o.data);
    readHeader(rc, "test", 1);
    assertThrown(read(rc, c));
}
```

The worked example uses it like this: one struct, `Saved`, holds
everything a checkpoint holds, so every part of the file is a named field.

```d
/// The checkpoint file's magic string and format version. Bump the
/// version when a saved field changes meaning or type, and teach load to
/// read the older one; adding a field needs no bump (see checkpoint.d).
enum string MAGIC = "solver2-example checkpoint";
enum uint FORMAT = 1;

/// What a checkpoint holds: the problem's fingerprint, the run's state, its
/// archive and causes, and the model.
struct Saved
{
    string[2][] fingerprint;
    RunState!ArchiveOptions state;
    typeof(A.Run.archive) archive;
    A.Run.Causes causes;
    ubyte[] model;
}

/// Saves a run stopped between iterations.
void save(string path, const ref Beam p, const ref A.Run run, Model model)
{
    const Saved s = {fingerprint: p.fingerprint, state: runState(run), archive: run.archive, causes: run.causes,
                     model: model.save(run.evals)};
    Bytes o;
    checkpoint.writeHeader(o, MAGIC, FORMAT);
    checkpoint.write(o, s);
    write(path, o.data);
}

/// Loads what save wrote into a new run, ready to pass as runArchive's
/// `resume`, and restores the model.
A.Run load(string path, const ref Beam p, Model model)
{
    auto r = Bytes(cast(ubyte[]) read(path));
    checkpoint.readHeader(r, MAGIC, FORMAT);
    Saved s;
    checkpoint.read(r, s);
    if (!r.done)
        throw new Exception("the checkpoint has bytes left over");
    if (const why = fingerprintMismatch(s.fingerprint, p.fingerprint))
        throw new Exception("the checkpoint is of another problem: " ~ why);
    A.Run run;
    restoreRun(run, s.state);
    run.archive = s.archive;
    run.causes = s.causes;
    model.restore(s.model, run.evals);
    return run;
}
```

The example then checks that the resumed run's archive is the
uninterrupted run's, bit for bit.

#### Adapting it

- **Name the file**: give `MAGIC` your project's name, and put what a
  checkpoint holds in `Saved`, one field per thing to save. Anything your
  problem keeps between iterations (its own counters, a schedule) goes in
  `Saved` too, as plain data.
- **Adding a field** (to `Saved`, or when a solver2 upgrade adds one to
  `RunState`, `Archive` or `ArchiveOptions`) needs no new version. An
  older file lacks it, and the field keeps the default it has in the
  `Saved s;` the reader starts from. Give a field you add a default
  initialiser that means "as before" (pleaser does the same in its text
  formats, below). When no default can mean that, bump `FORMAT`.
- **Removing a field**: an older file still has it, and `read` refuses
  any field it doesn't know, so that a newer program's file is never half
  read. Name the removed field in `read`'s `retired` list, as
  `"Struct.field"` (the struct's name without template arguments):
  `checkpoint.read(r, s, ["Saved.schedule"]);`.
- **Changing a field's type or meaning, or renaming it**: bump `FORMAT`
  and convert in `load`. Keep the old layout as a struct of its own and
  read it when the header says so:

  ```d
  enum uint FORMAT = 2;            // for example, 2: the fingerprint became a struct

  struct SavedV1 { string[2][] fingerprint; RunState!ArchiveOptions state; /* ... */ ubyte[] model; }

  A.Run load(string path, const ref Beam p, Model model)
  {
      auto r = Bytes(cast(ubyte[]) read(path));
      Saved s;
      if (checkpoint.readHeader(r, MAGIC, FORMAT) == 1)
      {
          SavedV1 old;
          checkpoint.read(r, old);
          s = upgrade(old);        // your conversion
      }
      else
          checkpoint.read(r, s);
      // ... as before
  }
  ```

  A field whose type changed without a bump is refused when its bytes
  don't fit the new type (`"of another type"`). One of the same size
  (`int` to `float`) isn't caught: bump the version.
- **solver2 upgrades**: a field solver2 adds reads as its default, and
  nothing promises that default means "as before". A resumed run is only
  guaranteed to be the uninterrupted run, bit for bit, on the solver2 it
  was saved with. If that matters, put your solver2 revision in the
  fingerprint.
- **Platforms**: a checkpoint loads on any platform, but resumes bit for
  bit only on the one it was saved on; elsewhere it goes on as that
  platform's run would (runs differ between platforms). Record the
  platform in the file (a field of `Saved`, or the fingerprint's
  information) if you need to tell; pleaser's states name theirs.
- **Write through a temporary file and a rename**, so a stop in mid-write
  leaves the last checkpoint whole: `write(path ~ ".tmp", o.data);
  rename(path ~ ".tmp", path);`.

pleaser's `elites.state` (`elitesText` and `parseElites` in
`pleaser/source/pleaser/archive.d`) is a real versioned format, as text: a
header line, then one tagged line per thing (`check`, `alpha`,
`threshold`, `elite`, `run ...`), with numbers in hexadecimal floating
point so they read back bit for bit. Fields added since its first version
are written only when they differ from their default (`cost` on an elite
line, the model-guided explorers' `run surrogate` line), so older files
read as before and files of runs without them are unchanged. A text
format costs more code per field and is worth it when people read or edit
the files; the template above is the shorter route.

A stage saves the same way: a `StageRun` and the caller's `Best` (its
`eval`, `genes`, `feasible` and `evals`), plus the model's bytes when it
pre-screens. Both are plain data, so the template writes them as they are:
add them to `Saved` as fields.

---

## Interfacing with Python

A problem written in Python runs through the same searches. They run in
D, in solver2's shared library ([`bindings/c`](../bindings/c/README.md));
the Python module ([`bindings/python/solver2`](../bindings/python/solver2/__init__.py))
is `ctypes` and numpy over it, with no compiled part of its own, so it
isn't tied to a Python version. CI runs its tests on Linux and Windows
with the newest Python, and on Linux with Python 3.10, the oldest tested.

```python
import numpy as np
import solver2

class Beam(solver2.Problem):
    genes = [solver2.Gene("width", 1, 10, neutral=5),
             solver2.Gene("rate", 0.1, 100, neutral=1, scale="log")]
    limits = [solver2.Limit("load", "hard"),
              solver2.Limit("rate_cap", "soft", weight=0.05)]
    axes = [solver2.Axis("width", [3, 5, 7], "nmwW", width=2),
            solver2.Axis("rate", [1, 10, 40], "slfF", width=10)]
    min_score = -50

    def evaluate(self, genes, feasible_only):
        w, r = self.physical(genes).T
        violations = np.column_stack([np.where(w * r > 200, w * r - 200, 0.0),
                                      np.where(r > 40, r - 40, 0.0)])
        objective = 10 - (w - 7) * (w - 7) / 4 - (r - 30) * (r - 30) / 100
        return objective, violations, np.column_stack([w, r])

problem = Beam()
problem.check()

stage = solver2.Stage(problem, lambda_=8, sigma=0.2, budget=800, seed=42)
stage.run()
print(stage.best.value, stage.best.descriptors)

archive = solver2.Archive(problem, seeds=[stage.best.genes], evals=4000, seed=9)
archive.run(report=lambda it: print(it.iteration, it.coverage, it.qd_score))
elites = archive.elites()
```

[`bindings/python/examples/beam.py`](../bindings/python/examples/beam.py)
is the whole worked example: a pre-screened stage, the archive with a
model, saved, loaded and resumed. It is
[`docs/example-class`](example-class/source/app.d) in Python, and prints
the same lines, bit patterns included; a test holds it to that. The
module's docstrings are the reference; this section says how the pieces
fit and what to know. The rest of the manual applies as it stands: the
model of [genes, limits and ranking](#the-model), the options of
[a stage](#a-cma-es-stage) and of [the archive](#the-archive), and the
rules of [Determinism](#determinism).

### Setting up

1. Build the library, with LDC 1.42.0, in `bindings/c`:

   ```
   dub build --compiler=ldc2 --build=release
   ```

2. Put `bindings/python` on Python's path (`sys.path`, or `PYTHONPATH`)
   and `import solver2`. It needs numpy and a 64-bit Python.

There is no package to install with pip (see
[Gaps and caveats](#gaps-and-caveats)). The module looks for
`solver2.dll` (`libsolver2.so` on Linux) in this order:

1. the file named by the environment variable `SOLVER2_LIBRARY`;
2. beside the module, in `solver2/`;
3. in `bindings/c` of the repository the module is in, which is where
   dub builds it.

So in a checkout nothing needs copying; elsewhere, copy the library into
`solver2/`. The library has the D runtime linked in and needs no other
file beside it. The module and the library must be of the same interface
version, and the import fails with a message that says so when they
aren't. `solver2.library_path` is the file that was loaded,
`solver2.abi` the interface version, and `solver2.checked` whether the
library checks array bounds (see
[Bounds checks in the library](#bounds-checks-in-the-library)).

### A problem in Python

A problem is a subclass of `solver2.Problem`: its description as class
attributes, and an `evaluate`. The description is read once, when the
problem is first used.

| Attribute | Default | Meaning |
|---|---|---|
| `genes` | none | The searched parameters: `Gene(name, lo, hi, neutral=nan, scale="linear")`. Gene 0 is `lo` and gene 1 is `hi`; `neutral` is where a search starts (NaN: the middle); `scale` is `"linear"`, `"log"` (needs `lo > 0`) or `"integer"`. |
| `limits` | none | `Limit(name, tier="hard", weight=1.0, group=None)`, in the order of `evaluate`'s violations. `tier` is `"hard"`, `"terminal"` or `"soft"` (see [Limits: three tiers](#limits-three-tiers)); `weight` is a soft limit's cost per unit over. |
| `axes` | none | The archive's grid, one `Axis(name, edges, letters, nan_bin=0, width=1.0)` per descriptor (see [Grids](#grids)). A problem without axes runs stages only. |
| `ranking` | the library's (100, 50, 10) | `Ranking(hard_offset, terminal_offset, bound_penalty)`. Set it when your values reach the offsets: see [Setting the offsets](#setting-the-offsets). |
| `min_score` | 0 | The value's lower bound, where the archive's thresholds start: at or below the lowest value (the objective less the soft breach) a candidate within the limits can have. |
| `prediction_headroom` | 0 | How far over a limit a model's prediction may be before it counts as a breach. |
| `fingerprint` | `{}` | Names and values of your own (a version, a setting) that a saved run must match to be loaded on this problem. The description itself is checked as well. |

**`evaluate(self, genes, feasible_only)`** gets a whole batch (a stage's
population, or every emitter's in an archive iteration): an array with
one row of genes per candidate, its own copy. It returns three arrays,
one row per candidate:

- the objectives, higher is better. NaN means not scored, which a
  candidate over a hard limit may be when `feasible_only` is true (the
  searches always pass true);
- the violations, one column per limit: how far over, 0 within, never
  negative. `None` for a problem without limits;
- the descriptors, one column per axis. `None` without axes.

The library adds each tier's violations up itself (as `breachOf` does),
so the order of the sum and its bits are fixed.

- **Genes come as sampled** and may lie outside [0, 1].
  `self.physical(genes)` clamps them and maps them over their ranges,
  for one gene vector or a whole batch.
- **The same genes must give the same bits whatever the batch** (see
  [Determinism in Python](#determinism-in-python)).
- **`problem.check(rows=8, seed=1)`** is `checkProblem` (see
  [Dynamic problems](#dynamic-problems)): it evaluates the start genes,
  both corners and `rows` random rows, twice, in another order and
  alone, and raises a `solver2.Error` that says what is wrong, a
  candidate within the limits below `min_score` included. Run it once on
  a new problem.

Also on a problem: `dimension` (the number of genes), `cells` (the
grid's), `start` (the genes a search starts from when given none),
`cell_code(cell)` (a cell's name, one letter per axis) and `close()`.

Underneath, the problem is a [`DynamicProblem`](#dynamic-problems) on
the D side, so what that section says of the class route holds here. Its
learning spec is the class's default; a Python problem can't give its
own, because the model's threads never call back into Python.

### A stage from Python

```python
stage = solver2.Stage(problem, lambda_=8, sigma=0.2, budget=800, seed=42)
end = stage.run(report=lambda g: print(g.generation, g.best))
best = stage.best                         # None while no candidate is within the limits
print(end.generations, end.stalled, best.value, best.genes, stage.causes)
```

`Stage(problem, start=None, model=None, **options)`. The options are
[`StageOptions`](#a-cma-es-stage)' fields in snake case, as keywords
(`lambda_` for `lambda`, which Python reserves), less `ranking`, which
is the problem's. `lambda_`, `sigma` and `budget` have no defaults.

- **`run(budget=None, report=None)`** runs until the budget is spent,
  the stage stalls, or the report stops it, and returns a `StageResult`:
  `generations`, `evaluations`, `finished`, `stalled`, `stopped`,
  `feasible`, `best_value`. Called again, it goes on.
- **`report`** is called after each generation with a `StageGeneration`:
  `generation`, `evaluations`, `best`, `population_best`, `sigma`,
  `infeasible`, `sampled`, `predicted_feasible`. A true return stops the
  stage there.
- **`stage.best`** is a `Candidate`: `genes` (in [0, 1]), `value`,
  `objective`, `scored`, `hard`, `terminal`, `soft`, `violations` and
  `descriptors`.
- **`stage.result`** is where the stage stands, and **`stage.causes`**
  the candidates over each limit so far, by the limit's name.
- **To chain stages**, start the next from `stage.best.genes`. Each stage
  has a best of its own: the D route's `Best` carried across stages
  isn't in the interface.

### The archive from Python

```python
archive = solver2.Archive(problem, seeds=[best.genes], evals=4000, emitters=4, lambda_=8, seed=9)
end = archive.run(report=lambda it: print(it.iteration, it.coverage, it.qd_score))
e = archive.elites()
for cell in archive.diverse(3, 4):
    print(problem.cell_code(cell), e.value[cell], e.descriptors[cell])
```

`Archive(problem, seeds=None, model=None, **options)`. `seeds` is one
gene vector or several (the problem's start when None). The options are
[`ArchiveOptions`](#options)' fields in snake case, less `ranking`;
`rules` is a dict of [`ArchiveRules`](#restart-and-explorer-rules)'
fields in snake case (`rules={"idle_generations": 8}`).

- **`run(evals=None, report=None)`** runs until `evals` evaluations
  beyond the seeds are spent, or the report stops it, and returns an
  `ArchiveResult`: `cells`, `coverage`, `evaluations`,
  `seed_evaluations`, `iterations`, `lambda_`, `started`, `finished`,
  `stopped`, `qd_score`, `best_score`, `best_cell` (None while the
  archive is empty), `restarts_idle`, `restarts_infeasible`,
  `restarts_converged`, `reached`, `abandoned`. The first call evaluates
  the seeds, and raises a `solver2.Error` when none is within the
  limits.
- **`report`** is called after each iteration (every emitter's
  generation, not each emitter's as on the D route) with an
  `ArchiveIteration`: `iteration`, `evaluations`, `coverage`,
  `qd_score`, `best_score`, `improved`, `new_cells`, `restarts`. A true
  return stops the run there.
- **`elites()`** gives every cell at once, as an `Elites` of arrays with
  cell `c` at index `c`: `occupied`, `value`, `found`, `genes`,
  `objective`, `scored`, `hard`, `terminal`, `soft`, `violations`,
  `descriptors`. An empty cell's numbers are NaN.
- **`thresholds`** is each cell's threshold, **`causes`** the candidates
  over each limit, **`result`** where the run stands, and
  **`diverse(within, count)`** is [`selectDiverse`](#results).

### A learned model from Python

```python
model = solver2.Model(problem, 5, threads=4)
archive = solver2.Archive(problem, seeds=[best.genes], model=model, evals=4000, seed=9,
                          surrogate_pool=4, surrogate_from=3, surrogate_epochs=5, surrogate_aim=True)
```

`Model(problem, seed, threads=0, **options)` is a
[`LearnedModel`](#a-learned-model) on the class's default learning spec:
from the genes it learns the value, the hard and the terminal breach,
and the descriptors. Its options are `members`, `hidden`, `batch`,
`learning_rate` and `l2`.

- A run with a model needs `surrogate_pool` above 0, and an archive
  without one must leave it at 0.
- `threads` are worker threads of the model's own, for training and
  predicting. They never call `evaluate`, and their number doesn't
  change the bits.
- A model serves one stage or archive at a time, and keeps what it
  learned in one for the next. `model.rows` is the number of evaluations
  it has learned from.
- `model.close()` stops its threads; a model also works in a `with`
  statement, as a stage and an archive do.

### Stopping, saving and resuming from Python

A report that returns true stops a run after that generation or
iteration, and `run` called again goes on. A run stopped and gone on
with ends as the one that never stopped, bit for bit, and a larger
`budget` or `evals` extends a finished one.

`save()` gives a stage or an archive as bytes, its model's state
included, and `load` makes a run from them:

```python
first = solver2.Archive(problem, seeds=[best.genes], model=solver2.Model(problem, 5), **options)
first.run(report=lambda it: it.iteration == 10)
with open("run.state", "wb") as f:
    f.write(first.save())

with open("run.state", "rb") as f:
    run = solver2.Archive.load(problem, f.read(), model=solver2.Model(problem, 5))
run.run()                                  # ends as the run that never stopped
```

- **The model** given to `load` must be made as the saved run's was: the
  same seed and options (its `threads` may differ). `load` restores its
  state. A run saved without a model loads without one.
- **The problem** must have the saved run's fingerprint and description
  (genes, limits, axes, ranking, `min_score`); `load` raises a
  `solver2.Error` that says what differs.
- **The bytes** are the library's own format, versioned, and not a file
  the D checkpoint template of [To a file](#to-a-file) reads. What that
  section says of upgrades and platforms holds: a run resumes bit for
  bit on the platform and the library it was saved with.
- **Writing the file** is yours: write through a temporary file and a
  rename, so a stop in mid-write leaves the last one whole.

### Errors and exceptions

- **`solver2.Error`** (a `RuntimeError`) is the library refusing a call:
  its message says why. A wrong option's name is a `TypeError`, and
  arrays of the wrong shape, from you or from `evaluate`, a `ValueError`.
- **An exception raised in `evaluate` or in a report** stops the run and
  comes out of `run` as itself, a `KeyboardInterrupt` (Ctrl-C) included.
- **The run is then left whole**: it and its model are as they were
  after their last whole generation or iteration. Call `run` again and
  it ends as the run that never failed, bit for bit. A long run in a
  notebook survives an interrupt.

### Threads

- **Import the module from the main thread** the first time, or from one
  that lives as long as the module is used: the D runtime in the library
  takes that thread for its main one.
- **Any thread may then use it**, several at once. One call at a time
  per object: a run, its model and its problem are one object for this,
  so two runs on one problem must not run at the same time.
- **`evaluate` and the reports are called on the thread that called
  `run`.** `ctypes` releases Python's interpreter lock while a search
  runs and takes it again for each callback, so other Python threads run
  meanwhile.
- **A busy Python thread slows the callbacks.** Beside a Python thread
  that never sleeps, each callback waits for the lock: about 5 ms each
  (measured: 2600 callbacks in 14 s), which is Python's default switch
  interval (`sys.setswitchinterval`). The bits are the same.
- **A thread that has called the library** is paused for a moment, like
  any D thread, whenever the D garbage collector runs, for as long as it
  lives.

### Determinism in Python

The searches give the same bits for the same outcomes, on any number of
the model's threads. Keeping `evaluate` deterministic is yours:

- **Element-wise numpy arithmetic** gives the same bits for the same
  genes whatever the batch.
- **A sum over a row** (`ndarray.sum`) is the same from run to run, but
  isn't added in D's order: a problem that sums that way is
  deterministic by itself and can't be compared bit for bit with the
  same problem in D.
- **BLAS-backed operations** (`dot`, `@`) may change their last bits
  with their thread count. Avoid them, or pin their threads.
- **To match a D problem bit for bit**, read the genes through
  `physical`: it is the library's own mapping, and a log scale goes
  through the library's `exp` and `log`, where numpy's round otherwise.
- Runs differ between platforms, as on the D route.

`problem.check()` catches an evaluation that depends on its batch, on
the rows it tries.

### What Python costs

Measured on Windows with Python 3.15 and LDC 1.42.0 (2026-11-12):

- **A batch's call** costs about 25 µs when `evaluate` computes nothing,
  and about 55 µs with a numpy evaluation of 8 candidates.
- **A learned model** is trained before every generation or iteration,
  in the library: that, not Python, is where a model-guided run spends
  its time. Give the model `threads`.
- **The worked example** takes 12.5 s in Python against 7.6 s for the
  same program in D. 1.7 s of the difference is the library's bounds
  checks, and most of the rest is the model's state being saved after
  every iteration, so that a failed evaluation can leave the run whole.

### Bounds checks in the library

The library is built with every array index checked, so a bug on the D
side comes back as a `solver2.Error` instead of corrupting the Python
process. That is the `library` configuration, the default. To build
without the checks:

```
dub build --compiler=ldc2 --build=release --config=unchecked
```

Both configurations write the same file, so the last one built is the
one in place; `solver2.checked` (`solver2_checked()` in C) tells which
it is. The results are the same bits either way. The checks cost about a
quarter more time in solver2's own code (see
[Bounds checks](#bounds-checks)): the Python example takes 12.5 s with
them and 10.8 s without. Keep them unless a profile says they matter.

### The C interface, for other languages

The Python module is one caller of the shared library; any language that
calls C can be another.
[`bindings/c/include/solver2.h`](../bindings/c/include/solver2.h) is the
interface and its documentation: its opening comment has the rules, and
each function its own. In short:

- `solver2_init` first, from a thread that outlives the library's use,
  and `solver2_term` at the end. The library isn't to be started again
  in the same process.
- Every other function returns `SOLVER2_OK` or an error code, and
  `solver2_last_error` gives the message. Nothing is thrown across the
  boundary.
- A problem, a model, a stage and an archive are handles. A handle that
  was freed, or never was one, is an error, not a crash.
- The problem is data (`solver2_problem_spec`) and one function,
  `solver2_evaluate_fn`, called once per batch with flat arrays.
  Returning non-zero from it fails the call with `SOLVER2_CALLBACK`, and
  the run is left whole.
- Results are copied into buffers the caller owns.
- Any thread may call, several at once; one call at a time per object.
- `solver2_abi()` is the interface's version (`SOLVER2_ABI`, 2 today).
  Check it against the header a wrapper was written for.

A stage over a one-gene problem, in C:

```c
#include <stdio.h>
#include "solver2.h"

/* Maximise 10 - (x - 3)^2, for x in 0 .. 10. */
static int evaluate(void* context, size_t rows, const double* genes, int feasible_only,
                    unsigned char* scored, double* objective, double* violations, double* descriptors)
{
    solver2_problem* problem = *(solver2_problem**) context;
    for (size_t k = 0; k < rows; ++k)
    {
        double x;
        solver2_problem_physical(problem, 1, genes + k, &x);    /* clamps, then maps */
        scored[k] = 1;
        objective[k] = 10 - (x - 3) * (x - 3);
    }
    return 0;
}

int main(void)
{
    const solver2_gene gene = {"x", 0, 10, 5, SOLVER2_LINEAR};
    solver2_problem* problem = NULL;
    solver2_stage* stage = NULL;
    solver2_problem_spec spec;
    solver2_stage_options options;
    solver2_stage_result result;

    if (solver2_init() != SOLVER2_OK)
        return 1;
    solver2_problem_spec_default(&spec);
    spec.genes = &gene;
    spec.gene_count = 1;
    spec.evaluate = evaluate;
    spec.context = &problem;
    solver2_stage_options_default(&options);
    options.lambda = 8;
    options.sigma = 0.2;
    options.budget = 400;
    options.seed = 1;
    if (solver2_problem_new(&spec, &problem) != SOLVER2_OK
        || solver2_stage_new(problem, &options, NULL, NULL, &stage) != SOLVER2_OK
        || solver2_stage_run(stage, 0, NULL, NULL, &result) != SOLVER2_OK)
    {
        fprintf(stderr, "%s\n", solver2_last_error());
        return 1;
    }
    printf("best %.6f after %zu evaluations\n", result.best_value, result.evaluations);
    solver2_stage_free(stage);
    solver2_problem_free(problem);
    return solver2_term();
}
```

An archive goes the same way (`solver2_archive_new`,
`solver2_archive_run`, `solver2_archive_elites`), and so do a model and
the saved bytes.

**The header, the library and a wrapper must agree.** The header is
written by hand, and solver2's tests hold both sides to it: the
binding's own read it with ImportC and fail to build when the D
declarations differ, and the Python tests parse it and fail when the
module's `ctypes` declarations do. A wrapper in another language should
check `solver2_abi()` when it loads the library.

**How the library is built** (`bindings/c/dub.json`): the D runtime is
linked in statically and only the interface's functions are exported
(`-link-defaultlib-shared=false -fvisibility=hidden -dllimport=none`),
so the library loads with nothing beside it and no D installation.
Building it needs LDC and dub only. Its tests
([`bindings/c/README.md`](../bindings/c/README.md)) also need a C
compiler, for ImportC's preprocessor.

**What the interface leaves out** of the D route is in
[Gaps and caveats](#gaps-and-caveats).

---

## Lower-level pieces

- **`solver2.cmaes.CmaEs`**: deterministic CMA-ES for minimisation.
  - `CmaEs(x0, sigma0, lambda, seed)`, where `lambda = 0` gives
    4 + ⌊3 ln n⌋.
  - `ask()` returns `lambda` rows, valid until the next `ask`.
    `tell(fitness)` takes their fitness.
  - `save()` returns a `CmaEs.State`. `CmaEs.restore(state)` continues
    bit-identically when the state was taken between `tell` and `ask`.
  - Accessors: `stepSize`, `currentMean`, `populationSize`, `generation`,
    `bestFitness`, `bestPoint`.
  - `jacobiEigen` is its eigensolver.
- **`solver2.rng.SplitMix64`**: the only RNG solver2 uses, and the one your
  problem should use for any seeded noise (see
  [Determinism](#determinism)).
  - `next()`, `uniform01()` (53 bits), `bipolar()`.
  - Its state is the public `ulong state`.
  - It is the same generator as pleaser's `pleaser.rng`, bit for bit.
- **`solver2.detmath`**: `log2` and `exp2` that give the same bits on every
  platform, for a problem whose results must match across platforms.
  - `log2d(x)` for finite `x > 0` (no denormals or zero: floor the input
    first), `exp2d(y)` for `y` in [−1000, 1000], exactly 1 at `y = 0`.
    About 1e-15 relative error.
  - `log2v` and `exp2v` do the same on a 4-lane `__vector(double[4])`
    (`V`); every lane equals the scalar result bit for bit, so a loop can
    run vectors with a scalar tail. `select(mask, a, b)` is the lane-wise
    `mask ? a : b`.
  - Plain IEEE arithmetic, no `libm`. Its unittests hold it to reference
    bit patterns, and CI runs them on Linux and Windows.
  - It makes the maths it is used for platform-independent, not a whole
    run: solver2 itself still uses `libm` (see
    [Determinism](#determinism), item 5).
  - pleaser and its plugins keep their own identical copy,
    `pleaser.detmath`.
- **`solver2.ensemble.Ensemble`**: the network ensemble under
  `LearnedModel`, usable on its own.
  - `new Ensemble(inputs, outputs, EnsembleOptions)`
  - `scaleInputs(mean, sd)`
  - `train(x, y, epochs, pool)`, with NaN targets treated as missing
  - `predict(x, mean, spread)`
  - `save()` and `restore(bytes)`
  - `targetStats` and `Stats`
- **`solver2.bytes.Bytes`**: little-endian put and get of scalars and
  arrays, with an array's length first. `getArray!T(expect)` checks the
  length, and `done()` is true when everything has been read. It writes
  native layout with no endianness conversion, so the bytes are
  little-endian only on little-endian machines.
- **`solver2.stage.pickRendered`, `stalled`**: the pre-screening pick and
  the stall test, exposed for custom loops.
- **`solver2.archive`'s smaller helpers**:
  - `emitterSeed`, `explorerGain`, `predictedFitness`
  - `nearlyAllInfeasible`, `explorerRestart`, `improverRestart` (each
    takes an `ArchiveRules` last, the defaults when left out)
  - `pickTarget`, `refocus`, `copyOf`
  - `NoCauses`, `CausesOf!P`

---

## Determinism

solver2 keeps these rules itself:

- Every random draw comes from `SplitMix64`, seeded from the options.
  Emitter `k` draws from its own seed.
- Every population is evaluated in one batch, and results are consumed in
  input order.
- Sorts are stable, with ties broken by index.
- Sums run in a fixed order.
- `Ensemble` and `LearnedModel` give the same bits with any `TaskPool`
  size, or none.

Your problem must keep them too:

1. `evaluate` returns the same bits for the same genes whatever the batch,
   the order of the other candidates in it, and the number of threads. If
   you parallelise, compute each candidate independently and write it to
   its own slot.
2. Draw any seeded noise from `solver2.rng.SplitMix64`, seeded from the
   genes or your own fixed seed. Don't use `std.random`, whose algorithms
   aren't guaranteed to stay the same across Phobos versions.
3. Keep floating-point reductions in a fixed order. `std.algorithm.sum`
   sums pairwise and changes the bits. Don't use `@fastmath`. To vectorise
   a reduction under strict IEEE, split it into lanes with independent
   accumulators and combine them in a fixed order.
4. `features`, `targets` and `predicted` in a learning spec must be pure
   functions of their arguments.
5. Runs are the same on any thread count and CPU of one platform, but
   **not across platforms**. LDC turns `exp`, `log`, `pow`, `sin`, `cos`
   and the like into calls to the C maths library, and glibc and the MSVC
   runtime can round the last bit differently. solver2 itself calls it
   (CMA-ES's normal samples and step size, the explorers' noise,
   log-scaled genes), and so does almost any objective. A last-bit
   difference flips a ranking sooner or later, and from there the runs
   take different paths: pleaser's Linux and Windows runs part within a
   few generations. So:
   - keep golden runs and reference outputs **per platform**, and compare
     across platforms by results over seeds, not run for run;
   - a test that needs something to happen during a run (an explorer
     aimed, a restart) must force it, so it happens on every platform;
   - a saved run resumes bit for bit on the platform it was saved on. On
     another it resumes, but goes on as that platform's run would.

   `solver2.detmath`'s `log2d`/`exp2d` (and their 4-lane forms) give the
   same bits everywhere. Use them where a problem's own result must not
   depend on the platform, say a value it stores or compares against a
   fixed reference; they don't make a whole run platform-independent. With
   LDC 1.42.0, plain arithmetic gives the same bits with or without FMA
   instructions (`detmath`'s reference bits are the same built with
   `-mcpu=native` or `-mattr=+fma`); don't add flags such as
   `-fp-contract=fast` or `@fastmath` that change that.

A problem written in Python keeps the same rules with numpy: see
[Determinism in Python](#determinism-in-python).

---

## Performance

- **Inlining**: solver2 is a source library and the searches are templates
  over your problem, so `objective`, `breach` and `descriptors` inline.
  Keep them cheap and `const`. On the class route those calls are
  virtual and each outcome holds three small arrays: nothing next to a
  real evaluation, noticeable on a microsecond-scale problem.
- **Threads**: solver2 never starts threads of its own. Parallelism comes
  from your `evaluate` and from the `TaskPool` you give `LearnedModel`
  (the shared library makes that pool for a model given `threads`).
  Make `lambda × emitters` a multiple of your thread count, so each archive
  batch fills the pool.
- **Stack size**: with `N > 0`, `Archive` holds `Elite!E[N]` by value, and
  `ArchiveRun` holds the archive. Each `Run` local, or a `Run` returned by
  value, takes `N × (sizeof(E) + about 40)` bytes of its caller's frame.
  pleaser's 75-cell archive is about 40 KB, and pleaser needs an 8 MB stack
  on Windows. With `N = 0` (the default) the cells live on the heap, and a
  `Run` is a few hundred bytes. Use `N > 0` only when you have measured a
  gain, and then keep runs off the stack, for example `auto run = new
  A.Run;`.
- **Allocation**: each generation allocates the batch's results and the
  clamped genes. That costs nothing next to a real evaluation. For a
  microsecond-scale problem, allocation can dominate.

### Bounds checks

D can check every array index and slice against the array's length, and
throw a `RangeError` when one is out of range. Whether the check is
compiled in is a build flag. solver2 is a source library, so your
project's flag covers solver2's code as well as your own.

| Build | What dub passes to LDC | Bounds checks |
|---|---|---|
| `dub build` (debug, the default) | nothing | on everywhere |
| `--build=release` | `-release` | in `@safe` functions only |
| `--build=release-nobounds` | `-release -boundscheck=off` | off everywhere |

Most of solver2's loops are not `@safe`, so a release build runs them
unchecked: an index past the end reads or writes whatever is there, and
nothing reports it.

**To turn them on in a release build**, give LDC `-boundscheck=on`:

- for every build of your project, in its `dub.json`:

  ```json
  "dflags": ["-boundscheck=on"]
  ```

- for one build, in the `DFLAGS` environment variable, which dub adds to
  its own flags:

  ```
  DFLAGS=-boundscheck=on dub build --compiler=ldc2 --build=release
  ```

  In PowerShell: `$env:DFLAGS = '-boundscheck=on'`, then the `dub`
  command, then `Remove-Item Env:DFLAGS`.

**To turn them off**, use `--build=release-nobounds`, or give
`-boundscheck=off` in either of those two ways. `-boundscheck=safeonly`
is the release default.

A `-boundscheck` in a package's own `dflags` comes last on LDC's command
line and wins: neither `DFLAGS` nor the build type changes it. To switch
a project that sets it there, change its `dub.json`.

**What they cost.** On 55 CMA-ES stages, five of them pre-screened by a
learned model, over a problem whose evaluation costs nothing: 2.84 s
with the checks against 2.25 s without, 26% more (LDC 1.42.0, Windows,
2026-11-12). That is the cost in solver2's own code; a real evaluation
makes it a small part of a run. The results are the same bits with or
without them.

**Asserts are a separate switch.** A release build compiles `assert`
out, and `-enable-asserts` keeps it, passed the same way. solver2's
asserts check what its callers hand it (matching array lengths, a grid's
cell count).

**The C interface** (`bindings/c`, which the Python module loads) is
built with bounds checks on: there a bug in D would otherwise corrupt
the caller's process instead of coming back as an error. It sets the
flag in its own `dub.json`, so it has a second configuration without it:
see [Bounds checks in the library](#bounds-checks-in-the-library).

---

## Migrating from solver

solver2 is the next version of
[solver](https://github.com/x21-dev/solver). It started as solver's code
at commit `ae68dbb`, with the package and its modules renamed, and adds
the class route, the C interface and the Python module. The searches are
the same code.

### Who needs to

Nobody. solver stays as it is, and a project that uses it needs neither
rewriting nor rebuilding. Move a project when it wants what solver2
adds: a problem chosen at run time, or Python. solver is frozen by
default: a bug found in the code both hold is fixed in solver2, and goes
into solver only when its owner decides so for that bug, because a fix
that changes bits would invalidate the golden runs of the projects on
it.

### What stays the same

- The problem concept and its members, required and optional.
- `StageOptions`, `ArchiveOptions`, `ArchiveRules`, `Ranking` and their
  defaults.
- `cmaStage`, `Archiving` and `runArchive`, their parameters and their
  callbacks; `Best`, `StageRun`, `ArchiveRun`, `RunState`.
- `LearnedModel`, `Ensemble`, `CmaEs`, `Bytes`, `SplitMix64`, `detmath`.
- The checkpoint template of [To a file](#to-a-file).
- The bits: a migrated project's runs are the same, on one platform.

### What changes

1. **The names.** The dub package is `solver2` and every module
   `solver2.*`. Change the dependency:

   ```json
   "dependencies": {
       "solver2": { "path": "../solver2" }
   }
   ```

   and the imports: `import solver.problem;` becomes
   `import solver2.problem;`, and so for `stage`, `archive`, `model`,
   `ensemble`, `cmaes`, `bytes`, `rng` and `detmath`. That is a search
   and replace of `solver.` after `import`. Because the names differ, a
   build can never pick up the wrong library.
2. **A problem's member named `minScore` is now read by the archive**,
   as where its thresholds start (see [Running](#running)), and
   `Archiving`'s `minScore` argument must then stay 0: the build fails
   otherwise. A project that passes its lower bound as the template's
   argument, as pleaser does, goes on working unchanged.
3. **`Archive` has one more field, `minScore`**, after `alpha`. A run's
   bits are the same. What a project may meet:
   - A checkpoint written by the template gains the field. One written
     with solver lacks it and still reads; see
     [Saved states](#saved-states).
   - A save format of the project's own (pleaser's `elites.state`) is
     untouched: an archive it builds has the template argument's value.
   - An `Archive` built as a struct literal with every field in order
     would need the new one. With `N` 0 it has a constructor, so only an
     `N > 0` archive could be written so.
4. **A subclass of `DynamicProblem`** (solver had the class, without
   what makes it complete):
   - `genes`, `constraints`, `fingerprint` and `evaluateInto` are
     `abstract`, so a subclass that lacks one fails where it is
     constructed.
   - The class now declares `axes`, `ranking`, `minScore`,
     `screenPenalty`, `noFeasibleSeed`, `adoptCauses`, `featureCount`,
     `features`, `targetCount`, `targets`, `predicted` and
     `predictionHeadroom`, and, final, `grid`, `startGenes`, `newCauses`
     and `countCauses`. A subclass with a member of one of those names
     must `override` it, or rename it when it is final.
   - Its outcomes are checked: one within the limits needs `scored`, an
     objective and its descriptors.
5. **The worked example's archive starts its thresholds at −50**, below
   its lowest value, where solver's starts them at 0 against the
   manual's own rule. Its stage prints what solver's example prints; its
   archive's lines differ.

### The same bits

A migrated project gives byte-identical runs on one platform: the names
changed, the arithmetic didn't. To prove it for your project, make a
golden run before and after:

1. On solver, run with a fixed seed and write the result to a file with
   every bit: a checkpoint, or each elite's genes, value and threshold
   as hexadecimal floats (`%a`). Rounded decimals hide a last-bit
   difference.
2. Migrate, build with the same compiler and flags, run the same.
3. Compare the two files byte for byte.

Do it on each platform you ship: runs differ between platforms on either
library (see [Determinism](#determinism)), so a golden run is one
platform's.

What was checked for solver2 itself, on Windows with LDC 1.42.0: before
change 5, the worked example printed solver's lines, and copies of both
examples changed to print the final archive as hexadecimal floats gave
identical output (every threshold, every elite's value, genes and
`found`, and the QD score).

### Saved states

- **A checkpoint written by the template with solver reads in
  solver2.** The magic string is your project's own, so it doesn't
  change. The file lacks `Archive.minScore`, and the field keeps the
  template argument's value, which is where that archive's thresholds
  started. Checked on Windows (2026-11-13): a checkpoint written by
  solver's example at `ae68dbb` after 10 iterations was loaded by
  solver2's example (given solver's magic string) and resumed to the
  uninterrupted run's archive, bit for bit.
- **A `LearnedModel`'s bytes** (`save`, `restore`) have the same format.
- **A format of your own** reads as before: it builds its runs from
  plain data, and that data's meaning didn't change.
- **Don't add a `minScore` member to a problem with saved runs** unless
  it gives the value those runs started from: a resumed run checks that
  the archive's `minScore` equals the problem's.

### Moving a compile-time problem to the class

A compile-time problem needs no change to run in solver2. Move one to
[the class](#dynamic-problems) when it must be chosen, built or loaded
at run time, or handed to code that knows nothing about it.

| The compile-time problem | The class |
|---|---|
| `struct P` with its own `Eval` | `final class P : DynamicProblem`; the evaluation is `Outcome` |
| `dimension()` | `genes()`, a `GeneSpec` per gene: the dimension is their number |
| `evaluate`, returning `Eval[]` | `evaluateInto(genes, scoring, outcomes)`, filling `outcomes` |
| `scored`, `objective`, `breach` | fields of each `Outcome` |
| `descriptorCount`, `descriptors` | `Outcome.descriptors`, and `axes()` |
| `Causes`, `newCauses`, `countCauses` | `constraints()` and `Outcome.violations`: the counts are `CauseCounts`, per limit |
| `fingerprint` (a convention) | `fingerprint()` (required) |
| `featureCount`, `features`, `targetCount`, `targets`, `predicted` | defaults that learn the value, the breaches and the descriptors; override `targets` and `predicted` together for your own |
| The grid, the ranking and the start genes, written at the call site | `axes()`, `ranking()` and each gene's `neutral` |
| `Archiving`'s `minScore` argument | `minScore()` |
| `cmaStage`, `Archiving!(P, G).runArchive` | `runStage`, `runArchive` (`solver2.dynamic`) |
| `Best!(P.Eval)`, `A.Run`, `LearnedModel!P` | `DynamicBest`, `DynamicRun`, `DynamicModel` |

[`docs/example`](example/source/app.d) and
[`docs/example-class`](example-class/source/app.d) are one problem both
ways; compare them.

What the move costs:

- **Your own evaluation type.** `Outcome` has the fields the searches
  read and `targets`; anything else an evaluation carried (pleaser's
  per-candidate measurements) has no place in it.
- **An archive held by value** (`N > 0`): the entry points use `N` 0.
- **Inlined calls**: a few virtual calls and three small arrays per
  candidate, which only a microsecond-scale evaluation notices.
- **Saved runs**: a checkpoint holds the evaluations by type, so one
  written on the compile-time route doesn't load on the class route.

The bits carry over when the arithmetic and the learning spec do: a
subclass gives the bits of the same problem as a struct, and solver2's
tests hold it to that. With a learned model, the class's default spec
must be what the struct's was, or the model's picks differ: that is why
the two worked examples print different numbers.

---

## Gaps and caveats

These are true of solver2 `main` on 2026-11-13. Each one says what it
costs you and how to work around it today.

The searches:

1. **No checkpoint API, by design.** solver2 has `RunState`,
   `CmaEs.State`, `LearnedModel.save` and `Bytes`, and nothing in the
   core that writes a run to a file. The intended route is the template
   in [To a file](#to-a-file): copy `checkpoint.d` into your project and
   adapt it. It versions the file and stays readable when a field is
   added. `isSaveable`'s `saveEval`/`loadEval` are a convention for a
   project's own format; nothing in solver2 calls them. From Python and
   C a run is saved as bytes in the library's own format.
2. **The terminal tier on stage pre-screening.** `StageModel.predictScore`
   returns one violation, which the stage ranks as hard. `LearnedModel`
   reports `hard + terminal` there, so a predicted terminal breach ranks as
   a hard one in a pre-screened stage. The archive's model interface keeps
   the tiers.
3. **A prediction is ranked by its objective, without a predicted soft
   breach.** An evaluation is ranked by its value (the objective less
   the soft breach); a `Predicted` has no soft tier to subtract. So a
   learning spec for a problem with soft limits should learn the value
   as its objective. The class's default spec does. The worked example's
   `Beam` in `docs/example` learns the objective, as solver's example
   did, so its model doesn't see the cost of its soft limit.
4. **LDC only, and no dub registry package.** See
   [Adding solver2 to a project](#adding-solver2-to-a-project).
5. **Runs differ between platforms.** solver2's own CMA-ES and explorers
   use the C maths library, whose last bits differ between glibc and the
   MSVC runtime, so the same run takes different paths on Linux and
   Windows. This is accepted, not planned to change: keep goldens per
   platform (see [Determinism](#determinism)). `solver2.detmath` has
   platform-identical `log2` and `exp2` for a problem's own use. macOS
   isn't tested.
6. **Evaluation log, survey and model test tools stay in pleaser.**
   pleaser's `evals.bin` log, random survey and offline surrogate tests
   are mix-specific and stay there. Use the `evaluated` callback to build
   your own log.
7. **Nothing checks the ranking's offsets against your values**, on any
   route: see [Setting the offsets](#setting-the-offsets).
8. **`StageEnd.rise` is NaN unless the stage stalled**, so two `StageEnd`
   never compare equal with `==`. Compare their fields.

The class route:

9. **`checkProblem` can't prove a problem sound.** It catches an
   evaluation that isn't deterministic, or a `minScore` that is too
   high, only on the rows it tries.
10. **The entry points use `BinnedGrid`, `ArchiveOptions` and `N` 0.** A
    grid or an options type of your own, or an archive held by value,
    means calling the templates over `DynamicProblem` yourself (see
    [The entry points](#the-entry-points)).

Python and the C interface:

11. **The library owns the loop.** You give an `evaluate` and the search
    calls it; there is no way to ask for candidates and tell their
    results, as pycma and pyribs offer. It could be added without
    rewriting the searches and isn't planned. To run your own code
    between generations, use the report.
12. **The Python module is a folder, not a wheel.** Nothing installs it
    with pip: build the library and put `bindings/python` on Python's
    path (see [Setting up](#setting-up)).
13. **What the interface leaves out of the D route**: the archive's
    `evaluated` callback and its per-emitter reports (a report is per
    iteration), [focused runs](#focused-runs), a grid other than
    `BinnedGrid`, a learning spec of the problem's own (a model learns
    the class's defaults), `screenPenalty`, the `noFeasibleSeed` hint,
    and a `Best` carried from one stage to the next (start the next
    stage from the last one's best genes).
14. **The library can't be started again** in one process after
    `solver2_term`. The Python module never ends it.
15. **Only 64-bit, and only the Pythons CI runs**: the newest, and 3.10
    on Linux. The module uses nothing newer than 3.7, but older versions
    aren't tested.
16. **A saved run's bytes are the library's own format**, not the
    checkpoint template's: a run saved from Python doesn't load in a D
    program that uses the template, nor the reverse.
17. **A model's state is saved after every generation or iteration** of
    a run from Python or C, so that a failed evaluation can leave the
    run whole. That is most of what the Python example costs over the D one
    (see [What Python costs](#what-python-costs)).

---

## Module reference

| Module | What it holds |
|---|---|
| `solver2.problem` | `Tier`, `Breach`, `ConstraintSpec`, `breachOf`, `value`, `Ranking`, `Scoring`, `clamp01`, `boundExcess`, `Scale`, `GeneSpec`, `Predicted`, `isProblem`, `hasDescriptors`, `hasSurrogate`, `isSaveable`, `valueOf`, `fitnessOf`, `isGrid`, `Axis`, `BinnedGrid`, `Outcome`, `CauseCounts`, `DynamicProblem`, `checkProblem`, `checkDescription`, `fingerprintMismatch` |
| `solver2.dynamic` | `runStage`, `runArchive`, and the searches' types over `DynamicProblem`: `DynamicBest`, `DynamicStageGeneration`, `DynamicArchiving`, `DynamicRun`, `DynamicArchiveGeneration`, `DynamicBatch`, `DynamicModel` |
| `solver2.stage` | `StageModel`, `StageOptions`, `Best`, `hasCauses`, `StageGeneration`, `StageEnd`, `StageRun`, `cmaStage`, `pickRendered`, `stalled` |
| `solver2.archive` | `Elite`, `Offer`, `Archive`, `copyOf`, `selectDiverse`, `ArchiveOptions`, `ArchiveSurrogate`, `GuideStats`, `ScreenStats`, `RestartReason`, `reasonText`, `NoCauses`, `CausesOf`, `ArchiveGeneration`, `Emitter`, `ArchiveRun`, `RunState`, `runState`, `restoreRun`, `ArchiveBatch`, `NoFeasibleSeed`, `hasMinScore`, `Archiving` (`runArchive`), `ArchiveRules` (and its defaults' former names, `RESTART_*`, `EXPLORE_*`, `AIM_SAMPLES`), the restart helpers |
| `solver2.model` | `isLearning`, `LearnedModel`, `ColumnStats`, `columnStats` |
| `solver2.cmaes` | `CmaEs` (with `State`), `jacobiEigen` |
| `solver2.ensemble` | `EnsembleOptions`, `Ensemble`, `Stats`, `targetStats` |
| `solver2.bytes` | `Bytes` |
| `solver2.rng` | `SplitMix64` |
| `solver2.detmath` | `log2d`, `exp2d`, `log2v`, `exp2v`, `select`, `V`, `VL` |

Outside the core, in the same repository:

| Where | What it holds |
|---|---|
| `bindings/c` | The shared library (`solver2.dll`, `libsolver2.so`): `include/solver2.h`, the interface, and `source/solver2c`, the exported functions over `solver2.dynamic` |
| `bindings/python/solver2` | The Python module: `Problem`, `Gene`, `Limit`, `Axis`, `Ranking`, `Stage`, `Archive`, `Model`, `Error`, and the records `Candidate`, `Elites`, `StageGeneration`, `StageResult`, `ArchiveIteration`, `ArchiveResult` |
| `docs/example`, `docs/example-class`, `bindings/python/examples/beam.py` | The worked example on each route |

---

## References

1. M. C. Fontaine and S. Nikolaidis, "Covariance Matrix Adaptation
   MAP-Annealing," *GECCO 2023*. [arXiv:2205.10752](https://arxiv.org/abs/2205.10752)
2. M. C. Fontaine, J. Togelius, S. Nikolaidis and A. K. Hoover,
   "Covariance Matrix Adaptation for the Rapid Illumination of Behavior
   Space," *GECCO 2020*. [arXiv:1912.02400](https://arxiv.org/abs/1912.02400)
3. B. Tjanaka, M. C. Fontaine, D. H. Lee, Y. Zhang, N. R. Balam,
   N. Dennler, S. S. Garlanka, N. D. Klapsis and S. Nikolaidis, "pyribs: A
   Bare-Bones Python Library for Quality Diversity Optimization,"
   *GECCO 2023*. [arXiv:2303.00191](https://arxiv.org/abs/2303.00191);
   the restart rules in its
   [`EvolutionStrategyEmitter`](https://docs.pyribs.org/en/stable/api/ribs.emitters.EvolutionStrategyEmitter.html)
   and its CMA-ES's `check_stop`.
4. N. Hansen and A. Ostermeier, "Completely Derandomized Self-Adaptation
   in Evolution Strategies," *Evolutionary Computation* 9(2), 2001.
5. N. Hansen, "The CMA Evolution Strategy: A Tutorial," 2016.
   [arXiv:1604.00772](https://arxiv.org/abs/1604.00772)
6. J.-B. Mouret and J. Clune, "Illuminating search spaces by mapping
   elites," 2015. [arXiv:1504.04909](https://arxiv.org/abs/1504.04909)
