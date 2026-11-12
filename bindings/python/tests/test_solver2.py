"""The module's tests: a problem written in Python through a stage, the
archive and a learned model; runs stopped, saved, loaded and resumed;
exceptions raised in the problem and in a report; threads.

Two runs are "the same" here when their saved bytes are: those hold a
run's whole state, its model's included, so equal bytes are equal bits of
everything.

    python -m unittest discover -s tests        (in bindings/python)
"""
import os
import threading
import unittest

import numpy as np

import _paths  # noqa: F401
import solver2
from solver2 import Axis, Gene, Limit


class Field(solver2.Problem):
    """Maximise 10 - (x - 1)^2 - (y - 0.8)^2 - z^2 over x and y in 0 .. 2
    and z in -1 .. 1, with x + y <= 1.9 hard, z >= -0.9 terminal and
    x <= 1.5 soft. The descriptors are x and y."""

    genes = [Gene("x", 0, 2, neutral=1), Gene("y", 0, 2, neutral=0.6), Gene("z", -1, 1, neutral=0)]
    limits = [Limit("sum"), Limit("floor", "terminal"), Limit("x_cap", "soft", weight=0.5)]
    axes = [Axis("x", [0.5, 1, 1.5], "abcd", width=0.5), Axis("y", [0.5, 1, 1.5], "abcd", width=0.5)]
    ranking = solver2.Ranking(1000, 500, 20)
    min_score = 5
    prediction_headroom = 0.05
    fingerprint = {"problem": "field 1"}

    def __init__(self, fail_at=0, error=ValueError):
        self.batches = 0
        self.fail_at = fail_at
        self.error = error

    def evaluate(self, genes, feasible_only):
        self.batches += 1
        if self.batches == self.fail_at:
            raise self.error("the problem failed on purpose")
        x, y, z = self.physical(genes).T
        violations = np.column_stack([np.where(x + y > 1.9, x + y - 1.9, 0.0),
                                      np.where(z < -0.9, -0.9 - z, 0.0),
                                      np.where(x > 1.5, x - 1.5, 0.0)])
        objective = 10 - (x - 1) * (x - 1) - (y - 0.8) * (y - 0.8) - z * z
        if feasible_only:
            objective[violations[:, 0] != 0] = np.nan
        genes[:] = 99           # the genes are the evaluation's own copy
        return objective, violations, np.column_stack([x, y])


STAGE = dict(lambda_=8, sigma=0.2, budget=400, seed=3, stall_generations=0)
SCREENED = dict(lambda_=8, sigma=0.2, budget=240, seed=3, stall_generations=0, surrogate_pool=4, surrogate_from=3,
                surrogate_epochs=2)
ARCHIVE = dict(evals=400, emitters=3, explorers=1, lambda_=8, seed=5)
GUIDED = dict(ARCHIVE, surrogate_pool=4, surrogate_from=3, surrogate_epochs=2, surrogate_improvers=True,
              surrogate_aim=True)


def small_model(problem, threads=0, seed=7):
    return solver2.Model(problem, seed, threads=threads, members=2, hidden=16)


def whole_stage(screened=False, threads=0):
    p = Field()
    stage = solver2.Stage(p, model=small_model(p, threads) if screened else None,
                          **(SCREENED if screened else STAGE))
    r = stage.run()
    assert r.finished and not r.stopped
    return stage.save()


def whole_archive(guided=False, threads=0):
    p = Field()
    archive = solver2.Archive(p, model=small_model(p, threads) if guided else None,
                              **(GUIDED if guided else ARCHIVE))
    r = archive.run()
    assert r.finished and not r.stopped
    return archive.save()


class TheLibrary(unittest.TestCase):
    def test_what_it_is(self):
        self.assertEqual(solver2.abi, 2)
        self.assertTrue(os.path.isfile(solver2.library_path))
        self.assertIsInstance(solver2.checked, bool)


class AProblem(unittest.TestCase):
    def test_what_it_says_of_itself(self):
        p = Field()
        self.assertEqual((p.dimension, p.cells), (3, 16))
        self.assertEqual(p.start.tolist(), [0.5, 0.3, 0.5])
        self.assertEqual([p.cell_code(c) for c in (0, 14, 15)], ["aa", "cd", "dd"])
        self.assertEqual(p.batches, 0)
        p.check()
        self.assertGreater(p.batches, 0)

    def test_physical(self):
        p = Field()
        self.assertEqual(p.physical([0.25, 1.5, 0.5]).tolist(), [0.5, 2, 0])            # clamped
        rows = p.physical([[0.25, 1.5, 0.5], [-0.1, 0.3, 1.0]])
        self.assertEqual(rows.tolist(), [[0.5, 2, 0], [0, 0.6, 1]])
        with self.assertRaisesRegex(ValueError, "genes must be 3 wide"):
            p.physical([0.5, 0.5])

        class Scales(solver2.Problem):
            genes = [Gene("rate", 0.1, 100, scale="log"), Gene("count", 1, 8, scale="integer")]
        s = Scales()
        v = s.physical([[0, 0], [1, 1], [0.5, 0.5], [2, -1]])
        self.assertEqual(v[:, 1].tolist(), [1, 8, 5, 1])
        self.assertEqual(v[0, 0], 0.1)
        self.assertAlmostEqual(v[1, 0], 100, places=9)
        self.assertAlmostEqual(v[2, 0], (0.1 * 100) ** 0.5, places=9)
        self.assertEqual(v[3, 0], v[1, 0])
        self.assertEqual(s.start.tolist(), [0.5, 0.5])          # no neutral value: the middle

    def test_a_description_that_does_not_hold(self):
        class Lettered(Field):
            axes = [Axis("x", [0.5, 1, 1.5], "abc"), Field.axes[1]]
        with self.assertRaisesRegex(solver2.Error, "axis x: 3 letters for 4 bins"):
            Lettered().dimension

        class Scaled(Field):
            genes = [Gene("x", 0, 2, scale="lg")] + Field.genes[1:]
        with self.assertRaisesRegex(ValueError, "gene x: its scale 'lg' is none of linear, log, integer"):
            Scaled().dimension

        class Tiered(Field):
            limits = [Limit("sum", "firm")] + Field.limits[1:]
        with self.assertRaisesRegex(ValueError, "limit sum: its tier 'firm'"):
            Tiered().dimension

        class Empty(solver2.Problem):
            pass
        with self.assertRaisesRegex(solver2.Error, "needs genes"):
            Empty().dimension

    def test_check_says_what_is_wrong(self):
        class Drifting(Field):
            def evaluate(self, genes, feasible_only):
                o, v, d = super().evaluate(genes, feasible_only)
                return o + 1e-9 * self.batches, v, d
        with self.assertRaisesRegex(solver2.Error, "other bits when evaluated again"):
            Drifting().check()

        class Sunk(Field):
            min_score = 9.5
        with self.assertRaisesRegex(solver2.Error, "below the problem's minScore"):
            Sunk().check()

        class Unwritten(solver2.Problem):
            genes = Field.genes
        with self.assertRaises(NotImplementedError):
            Unwritten().check()

    def test_what_evaluate_returns_is_checked(self):
        class Short(Field):
            def evaluate(self, genes, feasible_only):
                o, v, d = super().evaluate(genes, feasible_only)
                return o, v[:, :2], d
        with self.assertRaisesRegex(ValueError, r"violations of shape \(\d+, 2\); \d+ candidates by 3 are due"):
            Short().check()

        class Two(Field):
            def evaluate(self, genes, feasible_only):
                return super().evaluate(genes, feasible_only)[:2]
        with self.assertRaisesRegex(TypeError, "must return three arrays"):
            Two().check()

        class Bare(Field):
            def evaluate(self, genes, feasible_only):
                return super().evaluate(genes, feasible_only)[0], None, None
        with self.assertRaisesRegex(ValueError, "returned no violations"):
            Bare().check()

        class Negative(Field):
            def evaluate(self, genes, feasible_only):
                o, v, d = super().evaluate(genes, feasible_only)
                v[0, 1] = -0.25
                return o, v, d
        with self.assertRaisesRegex(solver2.Error, "its violation of the limit floor is negative or NaN"):
            Negative().check()

    def test_without_limits_or_axes(self):
        class Sphere(solver2.Problem):
            genes = [Gene("a", 0, 1), Gene("b", 0, 1)]

            def evaluate(self, genes, feasible_only):
                g = self.physical(genes)
                return -((g - 0.3) * (g - 0.3)).sum(axis=1), None, None
        p = Sphere()
        p.check()
        self.assertEqual(p.cells, 0)
        stage = solver2.Stage(p, **STAGE)
        r = stage.run()
        best = stage.best
        self.assertTrue(r.finished and r.feasible and r.best_value > -1e-3)
        self.assertTrue(abs(best.genes[0] - 0.3) < 0.02 and best.violations.shape == (0,) and stage.causes == {})
        with self.assertRaisesRegex(solver2.Error, "needs a problem with axes"):
            solver2.Archive(p, **ARCHIVE)


class AStage(unittest.TestCase):
    def test_a_stage(self):
        p = Field()
        stage = solver2.Stage(p, **STAGE)
        self.assertEqual(stage.result[:6], (0, 0, False, False, False, False))
        self.assertTrue(np.isnan(stage.result.best_value))
        self.assertIsNone(stage.best)
        reports = []
        r = stage.run(report=reports.append)
        self.assertEqual((r.generations, r.evaluations, r.finished, r.stalled, r.stopped, r.feasible),
                         (50, 400, True, False, False, True))
        self.assertEqual(p.batches, 50)
        self.assertEqual([g.generation for g in reports], list(range(1, 51)))
        self.assertEqual(reports[-1].evaluations, 400)
        self.assertEqual(reports[-1].best, r.best_value)
        self.assertEqual(stage.result, r)

        best = stage.best
        self.assertEqual(best.value, r.best_value)
        self.assertTrue(best.scored and best.hard == 0 and best.terminal == 0 and best.value > 9.9)
        self.assertEqual((best.genes.shape, best.violations.shape, best.descriptors.shape), ((3,), (3,), (2,)))
        self.assertEqual(best.descriptors.tolist(), p.physical(best.genes)[:2].tolist())
        self.assertEqual(best.value, best.objective - 0.5 * best.violations[2])
        self.assertEqual(sorted(stage.causes), ["floor", "sum", "x_cap"])
        self.assertGreater(stage.causes["sum"], 0)

        # The same again; and its budget spent, another call does nothing.
        self.assertEqual(stage.save(), whole_stage())
        stage.run()
        self.assertEqual(p.batches, 50)
        # A larger budget goes on.
        self.assertEqual(stage.run(budget=480).generations, 60)

    def test_its_own_start_and_what_it_refuses(self):
        p = Field()
        stage = solver2.Stage(p, start=[0.2, 0.2, 0.6], **STAGE)
        stage.run()
        self.assertNotEqual(stage.save(), whole_stage())
        with self.assertRaisesRegex(ValueError, r"a stage's start must be of shape \(3,\)"):
            solver2.Stage(p, start=[0.2, 0.2], **STAGE)
        with self.assertRaisesRegex(TypeError, "a stage has no option 'population'. It has: lambda_, sigma"):
            solver2.Stage(p, population=8)
        with self.assertRaisesRegex(solver2.Error, "needs a population"):
            solver2.Stage(p, sigma=0.2, budget=80)
        with self.assertRaisesRegex(solver2.Error, "needs surrogate_pool above 0"):
            solver2.Stage(p, model=small_model(p), **STAGE)

    def test_stopped_saved_loaded_and_gone_on_with(self):
        for screened in (False, True):
            p = Field()
            options = SCREENED if screened else STAGE
            stage = solver2.Stage(p, model=small_model(p) if screened else None, **options)
            r = stage.run(report=lambda g: g.generation == 10)
            self.assertEqual((r.generations, r.stopped, r.finished), (10, True, False))
            saved = stage.save()
            stage.close()
            with self.assertRaisesRegex(solver2.Error, "this stage was closed"):
                stage.result

            back = solver2.Stage.load(p, saved, model=small_model(p, threads=2) if screened else None)
            self.assertEqual(back.save(), saved)
            self.assertEqual(back.result.generations, 10)
            r = back.run()
            self.assertTrue(r.finished and not r.stopped)
            self.assertEqual(back.save(), whole_stage(screened))

    def test_a_model(self):
        p = Field()
        model = small_model(p)
        stage = solver2.Stage(p, model=model, **SCREENED)
        reports = []
        stage.run(report=reports.append)
        self.assertEqual(model.rows, 240)
        self.assertEqual((reports[0].sampled, reports[-1].sampled), (0, 32))
        # The same on any number of the model's threads.
        self.assertEqual(stage.save(), whole_stage(True, threads=3))
        with self.assertRaisesRegex(TypeError, "a model has no option 'layers'"):
            solver2.Model(p, 7, layers=3)
        model.close()
        with self.assertRaisesRegex(solver2.Error, "this model was closed"):
            model.rows


class WhenSomethingIsRaised(unittest.TestCase):
    def test_in_evaluate_the_stage_goes_on(self):
        for screened in (False, True):
            whole = whole_stage(screened)
            for fail_at in (1, 4, 17):
                p = Field(fail_at=fail_at)
                stage = solver2.Stage(p, model=small_model(p) if screened else None,
                                      **(SCREENED if screened else STAGE))
                with self.assertRaisesRegex(ValueError, "the problem failed on purpose"):
                    stage.run()
                # As after its last whole generation.
                self.assertEqual(stage.result.generations, fail_at - 1)
                self.assertTrue(stage.run().finished)
                self.assertEqual(stage.save(), whole)

    def test_in_evaluate_the_archive_goes_on(self):
        for guided in (False, True):
            whole = whole_archive(guided)
            for fail_at in (1, 2, 6):           # the seeds, the first iteration, a later one
                p = Field(fail_at=fail_at)
                archive = solver2.Archive(p, model=small_model(p) if guided else None,
                                          **(GUIDED if guided else ARCHIVE))
                with self.assertRaises(ValueError):
                    archive.run()
                r = archive.result
                self.assertEqual((r.started, r.iterations), (fail_at > 2, max(fail_at - 2, 0)))
                self.assertTrue(archive.run().finished)
                self.assertEqual(archive.save(), whole)

    def test_it_comes_out_as_itself(self):
        class Mine(Exception):
            pass
        p = Field(fail_at=3, error=Mine)
        stage = solver2.Stage(p, **STAGE)
        with self.assertRaises(Mine) as caught:
            stage.run()
        self.assertEqual(str(caught.exception), "the problem failed on purpose")

        # A Ctrl-C in the evaluation is a KeyboardInterrupt out of run, and
        # the run is whole.
        p = Field(fail_at=5, error=KeyboardInterrupt)
        stage = solver2.Stage(p, **STAGE)
        with self.assertRaises(KeyboardInterrupt):
            stage.run()
        self.assertEqual(stage.result.generations, 4)
        stage.run()
        self.assertEqual(stage.save(), whole_stage())

    def test_in_a_report_the_run_stops_there(self):
        def failing(g):
            if g.generation == 4:
                raise RuntimeError("the report failed on purpose")

        p = Field()
        stage = solver2.Stage(p, **STAGE)
        with self.assertRaisesRegex(RuntimeError, "the report failed on purpose"):
            stage.run(report=failing)
        self.assertEqual((stage.result.generations, stage.result.stopped), (4, True))
        stage.run()
        self.assertEqual(stage.save(), whole_stage())

        archive = solver2.Archive(p, **ARCHIVE)
        with self.assertRaisesRegex(RuntimeError, "the report failed on purpose"):
            archive.run(report=lambda it: failing(solver2.StageGeneration(it.iteration, 0, 0, 0, 0, 0, 0, 0)))
        self.assertEqual(archive.result.iterations, 4)
        archive.run()
        self.assertEqual(archive.save(), whole_archive())


class AnArchive(unittest.TestCase):
    def test_an_archive(self):
        p = Field()
        archive = solver2.Archive(p, **ARCHIVE)
        r = archive.result
        self.assertEqual((r.cells, r.started, r.coverage, r.evaluations, r.best_cell), (16, False, 0, 0, None))
        self.assertFalse(archive.elites().occupied.any())
        self.assertEqual(archive.thresholds.tolist(), [5.0] * 16)
        self.assertEqual(archive.diverse(3, 4), [])

        reports = []
        r = archive.run(report=reports.append)
        self.assertEqual((r.started, r.finished, r.stopped, r.iterations, r.evaluations, r.seed_evaluations,
                          r.lambda_), (True, True, False, 17, 1 + 17 * 24, 1, 8))
        self.assertEqual([it.iteration for it in reports], list(range(1, 18)))
        self.assertEqual((reports[-1].evaluations, reports[-1].coverage, reports[-1].qd_score),
                         (r.evaluations, r.coverage, r.qd_score))

        e = archive.elites()
        self.assertEqual((e.occupied.sum(), e.occupied.dtype), (r.coverage, bool))
        self.assertGreaterEqual(r.coverage, 6)
        self.assertEqual((e.genes.shape, e.violations.shape, e.descriptors.shape), ((16, 3), (16, 3), (16, 2)))
        held, empty = e.occupied, ~e.occupied
        self.assertTrue(np.isnan(e.value[empty]).all() and np.isnan(e.genes[empty]).all()
                        and (e.found[empty] == 0).all() and not e.scored[empty].any())
        self.assertTrue((e.value[held] >= 5).all() and e.scored[held].all() and (e.hard[held] == 0).all()
                        and (e.found[held] > 0).all())
        self.assertEqual(e.value[held].tolist(), (e.objective[held] - e.soft[held]).tolist())
        self.assertEqual(e.descriptors[held].tolist(), p.physical(e.genes[held])[:, :2].tolist())
        self.assertEqual(r.best_cell, int(np.nanargmax(e.value)))
        self.assertEqual(r.best_score, e.value[r.best_cell])
        self.assertEqual(r.qd_score, float(np.sum(e.value[held] - 5)))      # in cell order, as the library sums
        # The cell beyond the limit on the sum is never reached.
        self.assertFalse(e.occupied[15])
        t = archive.thresholds
        self.assertTrue(t[15] == 5 and (t >= 5).all() and (t[held] <= e.value[held]).all())
        picks = archive.diverse(3, 4)
        self.assertTrue(picks[0] == r.best_cell and 1 < len(picks) <= 4 and len(set(picks)) == len(picks))
        self.assertGreater(archive.causes["sum"], 0)

        self.assertEqual(archive.save(), whole_archive())
        self.assertEqual(archive.run(evals=520).iterations, 22)

    def test_its_seeds_and_what_it_refuses(self):
        p = Field()
        archive = solver2.Archive(p, seeds=[[0.5, 0.3, 0.5], [0.2, 0.2, 0.6]], **ARCHIVE)
        self.assertEqual(archive.run().seed_evaluations, 2)
        self.assertEqual(solver2.Archive(p, seeds=[0.5, 0.3, 0.5], **ARCHIVE).run().seed_evaluations, 1)

        far = solver2.Archive(p, seeds=[[1.0, 1.0, 0.5]], **ARCHIVE)
        with self.assertRaisesRegex(solver2.Error, "no seed of the archive is within the limits"):
            far.run()
        self.assertFalse(far.result.started)
        with self.assertRaisesRegex(ValueError, "seeds must be gene vectors 3 wide"):
            solver2.Archive(p, seeds=[[0.5, 0.5]], **ARCHIVE)
        with self.assertRaisesRegex(solver2.Error, "needs a model when surrogate_pool is above 0"):
            solver2.Archive(p, **GUIDED)
        with self.assertRaisesRegex(solver2.Error, "infeasibleShare"):
            solver2.Archive(p, rules=dict(infeasible_share=0), **ARCHIVE)
        with self.assertRaisesRegex(TypeError, "an archive's rules has no option 'patience'"):
            solver2.Archive(p, rules=dict(patience=3), **ARCHIVE)
        # Rules as given change the run.
        quick = solver2.Archive(p, rules=dict(idle_generations=1), alpha=1.0, **ARCHIVE)
        quick.run()
        self.assertNotEqual(quick.save(), whole_archive())

    def test_stopped_saved_loaded_and_gone_on_with(self):
        for guided in (False, True):
            p = Field()
            archive = solver2.Archive(p, model=small_model(p) if guided else None,
                                      **(GUIDED if guided else ARCHIVE))
            r = archive.run(report=lambda it: it.iteration == 6)
            self.assertEqual((r.iterations, r.stopped, r.finished), (6, True, False))
            saved = archive.save()
            archive.close()

            back = solver2.Archive.load(p, saved, model=small_model(p, threads=2) if guided else None)
            self.assertEqual(back.save(), saved)
            r = back.run()
            self.assertTrue(r.finished and not r.stopped and r.iterations == 17)
            self.assertEqual(back.save(), whole_archive(guided))

    def test_a_model(self):
        p = Field()
        model = small_model(p)
        archive = solver2.Archive(p, model=model, **GUIDED)
        r = archive.run()
        self.assertEqual(model.rows, r.evaluations)
        self.assertEqual(archive.save(), whole_archive(True, threads=3))

    def test_what_a_load_refuses(self):
        p = Field()
        plain = whole_archive()
        guided = whole_archive(True)
        with self.assertRaisesRegex(solver2.Error, "not a solver2-c stage file"):
            solver2.Stage.load(p, plain)
        with self.assertRaisesRegex(solver2.Error, "cut short"):
            solver2.Archive.load(p, plain[:len(plain) // 2])
        with self.assertRaisesRegex(solver2.Error, "the saved run had a model"):
            solver2.Archive.load(p, guided)
        with self.assertRaisesRegex(solver2.Error, "the saved run had no model"):
            solver2.Archive.load(p, plain, model=small_model(p))
        with self.assertRaisesRegex(solver2.Error, "another seed or other options"):
            solver2.Archive.load(p, guided, model=small_model(p, seed=8))

        class Revised(Field):
            fingerprint = {"problem": "field 2"}
        with self.assertRaisesRegex(solver2.Error, "of another problem: problem: saved field 1, now field 2"):
            solver2.Archive.load(Revised(), plain)

        class Wider(Field):
            genes = Field.genes[:2] + [Gene("z", -1, 2, neutral=0)]
        with self.assertRaisesRegex(solver2.Error, "solver2 genes"):
            solver2.Archive.load(Wider(), plain)
        self.assertEqual(solver2.Archive.load(p, plain).save(), plain)


class Threads(unittest.TestCase):
    def test_the_same_from_any_thread(self):
        stage, archive = whole_stage(), whole_archive(True)
        out = []

        def errand():
            try:
                out.append((whole_stage(), whole_archive(True, threads=2)))
            except BaseException as e:
                out.append(e)

        for count in (4, 1, 1, 1):
            threads = [threading.Thread(target=errand) for _ in range(count)]
            for t in threads:
                t.start()
            for t in threads:
                t.join()
        self.assertEqual(out, [(stage, archive)] * 7)

    def test_an_error_on_another_thread(self):
        out = []

        def errand():
            try:
                solver2.Stage(Field(), sigma=0.2, budget=80)
            except solver2.Error as e:
                out.append(str(e))

        t = threading.Thread(target=errand)
        t.start()
        t.join()
        self.assertEqual(len(out), 1)
        self.assertIn("needs a population", out[0])


if __name__ == "__main__":
    unittest.main()
