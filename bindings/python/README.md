# solver2 from Python

solver2's searches with the problem written in Python: a CMA-ES stage and
the archive (CMA-MAE with explorers), a learned model that pre-screens
and guides them, and runs saved and resumed as bytes. The searches run in
D, in solver2's shared library; this module is `ctypes` and numpy over
it, with no compiled part of its own, so it isn't tied to a Python
version.

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

[`examples/beam.py`](examples/beam.py) is the whole worked example: a
pre-screened stage, the archive with a model, saved, loaded and resumed.
The docstrings of [`solver2/__init__.py`](solver2/__init__.py) are the
reference.

## Setting up

1. Build the library, with LDC 1.42.0, in `bindings/c`:

   ```
   dub build --compiler=ldc2 --build=release
   ```

2. Put this folder on Python's path (`sys.path`, or `PYTHONPATH`) and
   `import solver2`. It needs numpy.

The module looks for `solver2.dll` (`libsolver2.so`) in this order: the
file named by the environment variable `SOLVER2_LIBRARY`; beside the
module, in `solver2/`; in `bindings/c` of the repository the module is
in, which is where dub builds it. So in a checkout nothing needs copying,
and elsewhere copy the library into `solver2/`. The module and the
library must be of the same interface version; the import says so when
they aren't.

## What to know

- **Import it from the main thread** the first time, or from one that
  lives as long as the module is used: the D runtime in the library
  takes that thread for its main one.
- **Any thread may then use it**, several at once. One call at a time per
  object: a run, its model and its problem are one object for this.
  `ctypes` lets other Python threads run while a search does.
- **`evaluate` gets a whole batch**, a stage's population or every
  emitter's in an archive iteration: an array with one row of genes per
  candidate, its own copy. It returns the objectives (NaN: not scored),
  the violations and the descriptors, one row per candidate.
- **Genes come as sampled** and may lie outside [0, 1]. `physical` clamps
  them and maps them over their ranges.
- **An exception raised in `evaluate` or in a report** stops the run and
  comes out of `run`, as itself, a Ctrl-C included. The run and its model
  are then as after their last whole generation or iteration: call `run`
  again and it ends as the run that never failed, bit for bit.
- **Determinism is `evaluate`'s to keep.** The searches give the same
  bits for the same outcomes, on any number of the model's threads. Your
  `evaluate` must give the same bits for the same genes whatever else is
  in its batch: element-wise numpy arithmetic does; a sum over a row is
  the same from run to run but is not added in D's order; BLAS-backed
  operations (`dot`, `@`) may change their last bits with their thread
  count. `problem.check()` catches an evaluation that depends on its
  batch, on the rows it tries.
- **To match a D problem bit for bit**, read the genes through
  `physical`: it is the library's own mapping, and a log scale goes
  through its `exp` and `log`, where numpy's round otherwise.
- **`min_score`** is where the archive's thresholds start: set it at or
  below the lowest value a candidate within the limits can have.
  `check()` refuses a problem whose candidates go below it.
- **A model's cost.** A learned model is trained before every generation
  or iteration, in the library; that, not Python, is where a model-guided
  run spends its time. Give it `threads`: their number doesn't change
  the bits.
- **The call's own cost** is about 25 µs a batch (8 candidates of 3
  genes, measured on Windows with Python 3.15), before whatever
  `evaluate` does.
- **Bounds checks.** The library is built with every array index checked,
  so a bug on the D side comes back as a `solver2.Error`. `solver2.checked`
  tells; `bindings/c/README.md` says how to build without.

## Testing

In this folder, with the library built:

```
python -m unittest discover -s tests
```

- `tests/test_conformance.py` reads `bindings/c/include/solver2.h` and
  holds the module's `ctypes` declarations to it: every struct's fields,
  every function's parameters, every constant.
- `tests/test_solver2.py` runs a problem written in Python through a
  stage, the archive and a model; stops, saves, loads and resumes runs;
  raises in the problem and in a report; and calls from several threads.
- `tests/test_example.py` builds `docs/example-class` with dub, the same
  problem as a D class, and checks that `examples/beam.py` prints the
  same lines, bit patterns included. It needs dub and LDC.

The tests use the standard library's `unittest`; nothing to install but
numpy.

Not run on Linux yet (the work plan's step 6).
