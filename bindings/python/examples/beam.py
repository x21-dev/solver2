"""The worked example from Python: docs/example-class's program, with the
problem's evaluation written in Python. A CMA-ES stage pre-screened by a
learned model, then the archive with a model, saved to bytes part way,
loaded from them, resumed, and checked against the run that never stopped.

It prints what docs/example-class prints, line for line, to the last bit
of the numbers on the lines that start with "bits": the searches are the
same D code, and the evaluation here does the same arithmetic in the same
order.

    python beam.py
"""
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import solver2  # noqa: E402


class Beam(solver2.Problem):
    """Two physical parameters, a width in 1 .. 10 and a log-scaled rate in
    0.1 .. 100. Maximise a smooth objective; width x rate <= 200 is hard,
    width >= 2 terminal, rate <= 40 soft at a cost of 0.05 per unit over."""

    genes = [solver2.Gene("width", 1, 10, neutral=5),
             solver2.Gene("rate", 0.1, 100, neutral=1, scale="log")]
    limits = [solver2.Limit("load", "hard"),
              solver2.Limit("min_width", "terminal"),
              solver2.Limit("rate_cap", "soft", weight=0.05)]
    axes = [solver2.Axis("width", [3, 5, 7], "nmwW", width=2),
            solver2.Axis("rate", [1, 10, 40], "slfF", width=10)]
    # The lowest value a beam within its limits can have is -48.25, at
    # width 2 and rate 100: the archive's thresholds start below it.
    min_score = -50
    prediction_headroom = 0.5
    fingerprint = {"problem": "beam 1", "load_limit": "200"}

    def evaluate(self, genes, feasible_only):
        # physical() clamps, and maps as the library does: the rate's log
        # scale with the library's own exp and log, for the same bits.
        w, r = self.physical(genes).T
        violations = np.column_stack([np.where(w * r > 200, w * r - 200, 0.0),
                                      np.where(w < 2, 2 - w, 0.0),
                                      np.where(r > 40, r - 40, 0.0)])
        objective = 10 - (w - 7) * (w - 7) / 4 - (r - 30) * (r - 30) / 100
        if feasible_only:
            objective[violations[:, 0] != 0] = np.nan       # over the hard limit: not scored
        return objective, violations, np.column_stack([w, r])


def bits(value):
    """A double's 64 bits."""
    return struct.unpack("<Q", struct.pack("<d", value))[0]


def hash_of(values):
    """FNV-1a over the bytes of `values`, least significant first."""
    h = 0xCBF29CE484222325
    for byte in np.ascontiguousarray(values, dtype="<f8").tobytes():
        h = ((h ^ byte) * 0x00000100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


def main(out=print):
    problem = Beam()
    problem.check()
    out("checked: %d genes, %d cells" % (problem.dimension, problem.cells))

    # 1. A CMA-ES stage from the problem's start, pre-screened by a model.
    def every_25th(g):
        if g.generation % 25 == 0:
            out("stage gen %3d  best %.4f  infeasible %d" % (g.generation, g.best, g.infeasible))

    # The models train and predict on four threads of their own. Their
    # number doesn't change the bits: a model is the same on any.
    stage = solver2.Stage(problem, lambda_=8, sigma=0.2, budget=800, seed=42, surrogate_pool=4,
                          model=solver2.Model(problem, 3, threads=4))
    end = stage.run(report=every_25th)
    best = stage.best
    out("stage: %d generations, %d evaluations, stalled %s; best %.4f at width %.3f rate %.3f"
        % (end.generations, end.evaluations, "true" if end.stalled else "false", best.value, best.descriptors[0],
           best.descriptors[1]))
    over = stage.causes
    out("stage: candidates over load %d, min_width %d, rate_cap %d"
        % (over["load"], over["min_width"], over["rate_cap"]))
    out("bits: stage best %016x genes %016x" % (bits(best.value), hash_of(best.genes)))

    # 2. The archive over the problem's grid, seeded with the stage's best:
    #    once straight through, once stopped after 10 iterations, saved,
    #    loaded and resumed.
    options = dict(evals=4000, emitters=4, explorers=2, lambda_=8, seed=9, surrogate_pool=4, surrogate_from=3,
                   surrogate_epochs=5, surrogate_aim=True)
    whole = solver2.Archive(problem, seeds=[best.genes], model=solver2.Model(problem, 5, threads=4), **options)
    whole.run()

    first = solver2.Archive(problem, seeds=[best.genes], model=solver2.Model(problem, 5, threads=4), **options)
    stopped = first.run(report=lambda it: it.iteration == 10)
    saved = first.save()
    out("stopped after %d iterations: %d cells, QD score %.2f" % (stopped.iterations, stopped.coverage,
                                                                 stopped.qd_score))
    first.close()

    run = solver2.Archive.load(problem, saved, model=solver2.Model(problem, 5, threads=4))
    resumed = run.run()
    out("resumed to %d evaluations: %d of %d cells, best %.4f, %d explorer targets reached"
        % (resumed.evaluations, resumed.coverage, resumed.cells, resumed.best_score, resumed.reached))

    a, b = whole.elites(), run.elites()
    same = (np.array_equal(a.occupied, b.occupied)
            and a.value[a.occupied].tobytes() == b.value[b.occupied].tobytes()
            and a.genes[a.occupied].tobytes() == b.genes[b.occupied].tobytes()
            and whole.thresholds.tobytes() == run.thresholds.tobytes()
            and whole.result.evaluations == resumed.evaluations)
    if not same:
        raise SystemExit("the resumed run differs")
    out("the resumed run's archive is the uninterrupted run's, bit for bit")
    for cell in run.diverse(3, 4):
        out("  %s  value %.4f  width %.3f  rate %.3f" % (problem.cell_code(cell), b.value[cell],
                                                       b.descriptors[cell, 0], b.descriptors[cell, 1]))
    out("bits: thresholds %016x values %016x genes %016x QD score %016x"
        % (hash_of(run.thresholds), hash_of(b.value[b.occupied]), hash_of(b.genes[b.occupied]),
           bits(resumed.qd_score)))


if __name__ == "__main__":
    main()
