"""solver2 from Python: a deterministic black-box optimiser whose searches
run in D and whose problem is yours, written in Python.

A problem is a class with its description and an ``evaluate``:

    import numpy as np
    import solver2

    class Beam(solver2.Problem):
        genes = [solver2.Gene("width", 1, 10, neutral=5),
                 solver2.Gene("rate", 0.1, 100, neutral=1, scale="log")]
        limits = [solver2.Limit("load", "hard"), solver2.Limit("rate_cap", "soft", weight=0.05)]
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
    stage = solver2.Stage(problem, lambda_=8, sigma=0.2, budget=800, seed=42)
    stage.run()
    archive = solver2.Archive(problem, seeds=[stage.best.genes], evals=4000, seed=9)
    archive.run()

The searches are solver2's own, through its C interface (bindings/c): this
module is ctypes and numpy over the shared library, with no compiled part
of its own. What it adds is small: the description as Python objects, the
arrays as numpy's, and exceptions. An exception raised in ``evaluate`` or
in a report stops the run and comes out of the call that started it, as
itself; the run is left whole, and can be gone on with.

Runs are deterministic: the same problem, options and seed give the same
bits, as long as ``evaluate`` gives the same bits for the same genes
whatever else is in its batch.
"""
import ctypes as C
import math
from collections import namedtuple

import numpy as np

from . import _c

__all__ = ["Archive", "ArchiveIteration", "ArchiveResult", "Axis", "Candidate", "Elites", "Error", "Gene", "Limit",
           "Model", "Problem", "Ranking", "Stage", "StageGeneration", "StageResult", "abi", "checked",
           "library_path"]

_lib, library_path = _c.load()
if _lib.solver2_init() != _c.OK:
    raise ImportError("solver2's library didn't start: " + _lib.solver2_last_error().decode("utf-8", "replace"))

#: The interface version of the library and of this module.
abi = _c.ABI
#: Whether the library checks every array index (the usual build).
checked = bool(_lib.solver2_checked())

_double_p = C.POINTER(C.c_double)
_NONE = 2 ** 64 - 1             # the library's "no cell"


class Error(RuntimeError):
    """The library refused a call. Its message says why."""


def _check(code, *sources):
    """Raises unless `code` is the library's OK: the exception one of
    `sources` kept from its own callback when there is one, else an Error
    with the library's message."""
    if code == _c.OK:
        return
    message = _lib.solver2_last_error().decode("utf-8", "replace")
    if code == _c.CALLBACK:
        for source in sources:
            raised = source.__dict__.pop("_raised", None)
            if raised is not None:
                raise raised
    raise Error(message)


def _text(s):
    return None if s is None else str(s).encode("utf-8")


def _doubles(values, shape, what):
    """`values` as a contiguous array of doubles of `shape`."""
    a = np.ascontiguousarray(values, dtype=np.float64)
    if a.shape != shape:
        raise ValueError("%s must be of shape %s, not %s" % (what, shape, a.shape))
    return a


def _pointer(array):
    return array.ctypes.data_as(_double_p)


def _fill(struct, options, what):
    """Sets `options` on a struct of the library's defaults. A keyword's
    last underscore is dropped (lambda_ is lambda)."""
    names = [name for name, _ in struct._fields_]
    for key, value in options.items():
        name = key[:-1] if key.endswith("_") else key
        if name not in names:
            raise TypeError("%s has no option %r. It has: %s" % (what, key, ", ".join(
                n + "_" if n == "lambda" else n for n in names)))
        if isinstance(getattr(struct, name), C.Structure):
            _fill(getattr(struct, name), dict(value), "%s's %s" % (what, name))
        else:
            setattr(struct, name, value)


def _record(name, struct, flags=()):
    """A named tuple with a ctypes struct's fields, and what makes one."""
    cls = namedtuple(name, [n + "_" if n == "lambda" else n for n, _ in struct._fields_])

    def make(s):
        return cls(*[bool(getattr(s, n)) if n in flags else getattr(s, n) for n, _ in struct._fields_])
    return cls, make


# ─────────────────────────────────────────────────────────────────────────────
# A problem's description
# ─────────────────────────────────────────────────────────────────────────────

Gene = namedtuple("Gene", "name lo hi neutral scale", defaults=(math.nan, "linear"))
Gene.__doc__ = """One searched parameter: gene 0 is `lo`, gene 1 is `hi`. `neutral` is where a
search starts (NaN: the middle). `scale`: "linear", "log" (needs lo > 0) or
"integer" (rounds to the nearest whole value)."""

Limit = namedtuple("Limit", "name tier weight group", defaults=("hard", 1.0, None))
Limit.__doc__ = """One limit. `tier`: "hard" (any breach ranks below every candidate within),
"terminal" (ranked after hard, before the objective) or "soft" (costs
`weight` per unit of breach)."""

Axis = namedtuple("Axis", "name edges letters nan_bin width", defaults=(0, 1.0))
Axis.__doc__ = """One axis of the archive's grid, over one descriptor: values below edges[0]
in bin 0, from the last edge up in the last bin, NaN in `nan_bin`.
`letters`: one character per bin (one more than the edges), for a cell's
name. `width`: the unit of distance on this axis, a typical bin's width."""

Ranking = namedtuple("Ranking", "hard_offset terminal_offset bound_penalty")
Ranking.__doc__ = """How candidates rank: over a hard limit at `hard_offset` plus the breach,
over a terminal one at `terminal_offset` plus the breach, else at minus the
value; plus `bound_penalty` per squared gene unit outside [0, 1]. The
offsets must be above every value the problem can have."""

_SCALES = {"linear": _c.LINEAR, "log": _c.LOG, "integer": _c.INTEGER}
_TIERS = {"hard": _c.HARD, "terminal": _c.TERMINAL, "soft": _c.SOFT}


def _choice(table, key, what):
    try:
        return table[key]
    except (KeyError, TypeError):
        raise ValueError("%s %r is none of %s" % (what, key, ", ".join(table))) from None


class Problem:
    """A problem to search. Subclass it: give `genes`, and what applies of
    `limits`, `axes` (for the archive), `ranking`, `min_score`,
    `prediction_headroom` and `fingerprint`; and `evaluate`.

    The description is read once, when the problem is first used.
    """

    #: The searched parameters: `Gene`s.
    genes = ()
    #: The limits: `Limit`s, in the order of `evaluate`'s violations.
    limits = ()
    #: The archive's grid: one `Axis` per descriptor. None: no archive.
    axes = ()
    #: A `Ranking`; None: the library's default (100, 50, 10).
    ranking = None
    #: The value's lower bound: no candidate within the limits has a value
    #: (its objective less its soft breach) below it. The archive's
    #: thresholds start there.
    min_score = 0.0
    #: How far over a limit a model's prediction may be before it counts.
    prediction_headroom = 0.0
    #: What a saved run must match to be loaded on this problem: names and
    #: values of your own (a version, a setting). The description is
    #: checked as well.
    fingerprint = {}

    def evaluate(self, genes, feasible_only):
        """Evaluates a batch: `genes` is an array of one row per candidate,
        as sampled, so they may lie outside [0, 1] (`physical` clamps).

        Returns three arrays, one row per candidate:

        - the objectives, higher is better. NaN: not scored, which a
          candidate over a hard limit may be when `feasible_only`;
        - the violations, one column per limit: how far over, 0 within,
          never negative (None for a problem without limits);
        - the descriptors, one column per axis (None without axes).

        The same genes must give the same bits whatever the batch.
        """
        raise NotImplementedError("a problem needs its evaluate")

    # ── What a caller can ask of it ──────────────────────────────────────

    @property
    def dimension(self):
        """The number of genes."""
        self._handle
        return self._n

    @property
    def cells(self):
        """The number of cells of the grid; 0 without axes."""
        self._handle
        return self._cells

    @property
    def start(self):
        """The genes a search starts from when given none: each gene's
        neutral value."""
        x = np.empty(self.dimension)
        _check(_lib.solver2_problem_start(self._handle, _pointer(x)))
        return x

    def physical(self, genes):
        """What `genes` stand for: each clamped to [0, 1] and mapped over
        its range by its scale. One gene vector, or an array of one per
        row; the result has the same shape.

        These are the library's own numbers (a log scale goes through its
        exp and log), so an evaluation that reads its genes through here
        gets the bits a problem written in D gets. numpy's exp and log
        round otherwise.
        """
        handle = self._handle
        g = np.ascontiguousarray(genes, dtype=np.float64)
        if g.ndim not in (1, 2) or g.shape[-1] != self._n:
            raise ValueError("genes must be %d wide, in one row or several; these are of shape %s"
                             % (self._n, g.shape))
        out = np.empty_like(g)
        _check(_lib.solver2_problem_physical(handle, g.size // self._n, _pointer(g), _pointer(out)))
        return out

    def cell_code(self, cell):
        """A cell's name: one letter per axis."""
        handle = self._handle
        buffer = C.create_string_buffer(self._nd + 1)
        _check(_lib.solver2_problem_cell_code(handle, cell, buffer, len(buffer)))
        return buffer.value.decode("utf-8")

    def check(self, rows=8, seed=1):
        """Checks the evaluation on the start genes, both corners and
        `rows` random rows: the same bits evaluated twice, in another order
        and alone, and no candidate within the limits below `min_score`.
        Raises an Error that says what is wrong. About 5 * (rows + 3)
        evaluations."""
        _check(_lib.solver2_problem_check(self._handle, rows, seed), self)

    def close(self):
        """Lets go of the problem in the library, for good: asked for
        anything afterwards it raises an Error, as a closed run does.

        Its runs hold it in the library and can still be read and saved.
        One gone on with stops, whole, at the Error when its `evaluate`
        asks the problem for something (`physical`): close the runs first.
        """
        handle = self.__dict__.pop("_made", None)
        self._closed = True
        if handle is not None:
            _lib.solver2_problem_free(handle)

    def __del__(self):
        try:
            self.close()
        except Exception:       # the interpreter is going away
            pass

    # ── The library's side ───────────────────────────────────────────────

    @property
    def _handle(self):
        handle = self.__dict__.get("_made")
        if handle is None:
            # Never made again: its runs call the callback of the first
            # making for as long as they live, and read arrays by its sizes.
            if self.__dict__.get("_closed"):
                raise Error("this problem was closed")
            handle = self._make()
        return handle

    def _make(self):
        genes = [g if isinstance(g, Gene) else Gene(*g) for g in self.genes]
        limits = [c if isinstance(c, Limit) else Limit(*c) for c in self.limits]
        axes = [a if isinstance(a, Axis) else Axis(*a) for a in (self.axes or ())]
        spec = _c.ProblemSpec()
        _check(_lib.solver2_problem_spec_default(C.byref(spec)))

        c_genes = (_c.Gene * len(genes))(*[
            _c.Gene(_text(g.name), g.lo, g.hi, g.neutral, _choice(_SCALES, g.scale, "gene %s: its scale" % g.name))
            for g in genes])
        c_limits = (_c.Limit * len(limits))(*[
            _c.Limit(_text(c.name), _choice(_TIERS, c.tier, "limit %s: its tier" % c.name), c.weight,
                     _text(c.group))
            for c in limits])
        edges = [(C.c_double * len(a.edges))(*a.edges) for a in axes]
        c_axes = (_c.Axis * len(axes))(*[
            _c.Axis(_text(a.name), e, len(e), _text(a.letters), a.nan_bin, a.width) for a, e in zip(axes, edges)])
        pairs = []
        for name, value in dict(self.fingerprint).items():
            pairs += [_text(name), _text(value)]
        c_pairs = (C.c_char_p * len(pairs))(*pairs)

        spec.genes, spec.gene_count = c_genes, len(genes)
        spec.limits, spec.limit_count = c_limits, len(limits)
        spec.axes, spec.axis_count = c_axes, len(axes)
        if self.ranking is not None:
            spec.ranking = _c.Ranking(*Ranking(*self.ranking))
        spec.min_score = self.min_score
        spec.prediction_headroom = self.prediction_headroom
        spec.fingerprint, spec.fingerprint_count = c_pairs, len(pairs) // 2
        # The library calls this for as long as the problem is used: kept
        # here, or Python would free it under the library.
        self._callback = _c.EvaluateFn(self._batch)
        spec.evaluate = self._callback

        self._n, self._limits, self._nd = len(genes), len(limits), len(axes)
        self._limit_names = [c.name for c in limits]
        handle, cells = C.c_void_p(), C.c_size_t()
        _check(_lib.solver2_problem_new(C.byref(spec), C.byref(handle)))
        self._made = handle
        # Kept, with the sizes above: a run is read by them after the
        # problem was closed.
        _check(_lib.solver2_problem_cells(handle, C.byref(cells)))
        self._cells = cells.value
        return handle

    def _batch(self, context, rows, genes, feasible_only, scored, objective, violations, descriptors):
        """The library's callback: one batch through `evaluate`. Nothing
        may be raised into the library: an exception is kept, the library
        is told the evaluation failed, and the call that started the run
        raises it."""
        try:
            # A copy: the library's memory is good only until this returns.
            g = np.ctypeslib.as_array(genes, shape=(rows, self._n)).copy()
            result = self.evaluate(g, bool(feasible_only))
            try:
                o, v, d = result
            except (TypeError, ValueError):
                raise TypeError("evaluate must return three arrays: the objectives, the violations and the "
                                "descriptors") from None
            o = np.asarray(o, dtype=np.float64)
            if o.shape != (rows,):
                raise ValueError("evaluate returned objectives of shape %s for %d candidates" % (o.shape, rows))
            np.ctypeslib.as_array(objective, shape=(rows,))[:] = o
            np.ctypeslib.as_array(scored, shape=(rows,))[:] = ~np.isnan(o)
            for given, into, columns, what in ((v, violations, self._limits, "violations"),
                                               (d, descriptors, self._nd, "descriptors")):
                if columns == 0:
                    continue
                if given is None:
                    raise ValueError("evaluate returned no %s, and the problem has %d columns of them"
                                     % (what, columns))
                given = np.asarray(given, dtype=np.float64)
                if given.shape != (rows, columns):
                    raise ValueError("evaluate returned %s of shape %s; %d candidates by %d are due"
                                     % (what, given.shape, rows, columns))
                np.ctypeslib.as_array(into, shape=(rows, columns))[:] = given
            return 0
        except BaseException as e:      # a Ctrl-C too
            self._raised = e
            return 1


# ─────────────────────────────────────────────────────────────────────────────
# A learned model
# ─────────────────────────────────────────────────────────────────────────────

class _Held:
    """Something the library holds a handle for, until it is closed."""

    def __getattr__(self, name):
        if name == "_made":
            raise Error("this %s was closed" % type(self).__name__.lower())
        raise AttributeError(name)


class Model(_Held):
    """An ensemble of small networks that learns every evaluation of a run
    and predicts candidates before they are evaluated: it pre-screens a
    stage (`surrogate_pool`) and guides the archive's emitters.

    `threads`: worker threads of its own for training and predicting. They
    never call `evaluate`, and their number doesn't change the bits.
    Options: members, hidden, batch, learning_rate, l2.
    """

    def __init__(self, problem, seed, *, threads=0, **options):
        o = _c.ModelOptions()
        _check(_lib.solver2_model_options_default(C.byref(o)))
        _fill(o, dict(options, threads=threads), "a model")
        self.problem = problem
        handle = C.c_void_p()
        _check(_lib.solver2_model_new(problem._handle, seed, C.byref(o), C.byref(handle)))
        self._made = handle

    @property
    def rows(self):
        """The evaluations it has learned from."""
        n = C.c_size_t()
        _check(_lib.solver2_model_rows(self._made, C.byref(n)))
        return n.value

    def close(self):
        """Stops its threads and lets go of it. The runs made with it
        hold it still, and go on without the threads, to the same bits."""
        handle = self.__dict__.pop("_made", None)
        if handle is not None:
            _lib.solver2_model_free(handle)

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


# ─────────────────────────────────────────────────────────────────────────────
# Runs
# ─────────────────────────────────────────────────────────────────────────────

StageGeneration, _stage_generation = _record("StageGeneration", _c.StageGeneration)
StageResult, _stage_result = _record("StageResult", _c.StageResult, ("finished", "stalled", "stopped", "feasible"))
ArchiveIteration, _archive_iteration = _record("ArchiveIteration", _c.ArchiveIteration)
ArchiveResult, _archive_result_raw = _record("ArchiveResult", _c.ArchiveResult, ("started", "finished", "stopped"))

Candidate = namedtuple("Candidate", "genes value objective scored hard terminal soft violations descriptors")
Candidate.__doc__ = """One candidate: its genes in [0, 1], its value (the objective less the soft
breach), its objective, each tier's breach, its violations per limit and
its descriptors."""

Elites = namedtuple("Elites", "occupied value found genes objective scored hard terminal soft violations "
                              "descriptors")
Elites.__doc__ = """Every cell of an archive at once, cell c at index c: whether it holds a
candidate, and that candidate's value, the evaluations spent when it was
found, its genes, its objective and breaches, its violations and its
descriptors. An empty cell's numbers are NaN, its `found` 0."""


def _archive_result(s):
    r = _archive_result_raw(s)
    return r._replace(best_cell=None) if r.best_cell == _NONE else r


class _Run(_Held):
    """What a stage and an archive share: a handle, a report, saved bytes."""

    def _report(self, report, kind, make):
        """`report` as the library calls it; a null pointer for none. A
        report that raises stops the run at its generation: the exception
        is kept, and `run` raises it."""
        if report is None:
            return C.cast(None, kind)

        def call(context, pointer):
            try:
                return 1 if report(make(pointer.contents)) else 0
            except BaseException as e:
                self._raised = e
                return -1
        return kind(call)

    def _counts(self, function):
        over = (C.c_size_t * max(len(self.problem._limit_names), 1))()
        _check(function(self._made, over))
        return dict(zip(self.problem._limit_names, over))

    def _bytes(self, function):
        size = C.c_size_t()
        _check(function(self._made, None, 0, C.byref(size)))
        buffer = (C.c_ubyte * size.value)()
        _check(function(self._made, buffer, size.value, C.byref(size)))
        return bytes(buffer)

    @classmethod
    def _loaded(cls, function, problem, data, model):
        self = cls.__new__(cls)
        self.problem, self.model = problem, model
        data = bytes(data)
        buffer = (C.c_ubyte * len(data)).from_buffer_copy(data)
        handle = C.c_void_p()
        _check(function(problem._handle, None if model is None else model._made, buffer, len(data),
                        C.byref(handle)))
        self._made = handle
        return self

    def close(self):
        """Lets go of the run in the library."""
        handle = self.__dict__.pop("_made", None)
        if handle is not None:
            self._free(handle)

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


class Stage(_Run):
    """A CMA-ES stage over `problem`: the best candidate from `start` (the
    problem's start when None), pre-screened by `model` when given.

    Options, as keywords: lambda_ (the population, at least 2), sigma (the
    initial step size, in gene units) and budget (evaluations) must be
    given; seed; stall_generations and stall_tolerance (stop when the best
    rose by less than the tolerance over that many generations; 0
    generations: never); with a model, surrogate_pool (at least 2),
    surrogate_share (from 0 to 1), surrogate_from and surrogate_epochs.
    """

    _free = staticmethod(_lib.solver2_stage_free)

    def __init__(self, problem, *, start=None, model=None, **options):
        o = _c.StageOptions()
        _check(_lib.solver2_stage_options_default(C.byref(o)))
        _fill(o, options, "a stage")
        handle = problem._handle
        x0 = None if start is None else _doubles(start, (problem.dimension,), "a stage's start")
        self.problem, self.model = problem, model
        made = C.c_void_p()
        _check(_lib.solver2_stage_new(handle, C.byref(o), None if x0 is None else _pointer(x0),
                                      None if model is None else model._made, C.byref(made)))
        self._made = made

    @classmethod
    def load(cls, problem, data, model=None):
        """The stage `save` gave the bytes of, on a problem with the same
        fingerprint and description. `model`: one made as the saved
        stage's was (the same seed and options), when it had one."""
        return cls._loaded(_lib.solver2_stage_load, problem, data, model)

    def run(self, budget=None, report=None):
        """Runs until `budget` evaluations are spent in all (None: the
        options'), the stage stalls, or `report` stops it; returns a
        `StageResult`. Called again it goes on, and a stage stopped and
        gone on with ends as the one that never stopped, bit for bit. A
        larger budget extends a finished stage.

        `report` is called after each generation with a `StageGeneration`;
        a true return stops the stage there. It may read the stage and
        `save` it: the stage is then as after that generation.

        An exception raised by `evaluate` or by `report` comes out of this
        call. The stage and its model are then as they were after the last
        whole generation, and `run` can be called again: it goes on to the
        budget last asked for, the failed call's included.
        """
        callback = self._report(report, _c.StageReportFn, _stage_generation)
        r = _c.StageResult()
        _check(_lib.solver2_stage_run(self._made, budget or 0, callback, None, C.byref(r)), self, self.problem)
        return _stage_result(r)

    @property
    def result(self):
        """Where the stage stands: a `StageResult`."""
        r = _c.StageResult()
        _check(_lib.solver2_stage_get(self._made, C.byref(r)))
        return _stage_result(r)

    @property
    def best(self):
        """Its best candidate, a `Candidate`; None while it has none within
        the limits."""
        if not self.result.feasible:
            return None
        p = self.problem
        genes, violations, descriptors = np.empty(p._n), np.empty(p._limits), np.empty(p._nd)
        e = _c.Outcome()
        _check(_lib.solver2_stage_best(self._made, _pointer(genes), C.byref(e), _pointer(violations),
                                       _pointer(descriptors)))
        return Candidate(genes, e.value, e.objective, bool(e.scored), e.hard, e.terminal, e.soft, violations,
                         descriptors)

    @property
    def causes(self):
        """Per limit, by its name, the candidates over it so far."""
        return self._counts(_lib.solver2_stage_causes)

    def save(self):
        """The stage as bytes, its model's state included: for `load`."""
        return self._bytes(_lib.solver2_stage_save)


class Archive(_Run):
    """The archive over `problem`'s grid: the best candidate in each cell,
    found by CMA-MAE with explorers that head for empty cells. From `seeds`
    (gene vectors; the problem's start when None), guided by `model` when
    given.

    Options, as keywords: evals (the budget beyond the seeds), emitters,
    lambda_ (per emitter, at least 2; 0: CMA-ES's default), sigma, alpha
    (the archive learning rate), seed, explorers, explorer_sigma; with a
    model, surrogate_pool (at least 2), surrogate_share (from 0 to 1),
    surrogate_from, surrogate_epochs, surrogate_improvers and
    surrogate_aim; rules (a dict
    of idle_generations, infeasible_generations, infeasible_share,
    min_sigma, stale_generations, min_progress, block_iterations,
    aim_samples).
    """

    _free = staticmethod(_lib.solver2_archive_free)

    def __init__(self, problem, *, seeds=None, model=None, **options):
        o = _c.ArchiveOptions()
        _check(_lib.solver2_archive_options_default(C.byref(o)))
        _fill(o, options, "an archive")
        handle = problem._handle
        n = problem.dimension
        rows = None
        if seeds is not None:
            rows = np.ascontiguousarray(seeds, dtype=np.float64)
            if rows.ndim == 1:
                rows = rows.reshape(1, -1)
            if rows.ndim != 2 or rows.shape[1] != n:
                raise ValueError("an archive's seeds must be gene vectors %d wide; these are of shape %s"
                                 % (n, rows.shape))
        self.problem, self.model = problem, model
        made = C.c_void_p()
        _check(_lib.solver2_archive_new(handle, C.byref(o), None if rows is None else _pointer(rows),
                                        0 if rows is None else len(rows), None if model is None else model._made,
                                        C.byref(made)))
        self._made = made

    @classmethod
    def load(cls, problem, data, model=None):
        """The run `save` gave the bytes of, on a problem with the same
        fingerprint and description. `model`: one made as the saved run's
        was (the same seed and options), when it had one."""
        return cls._loaded(_lib.solver2_archive_load, problem, data, model)

    def run(self, evals=None, report=None):
        """Runs until `evals` evaluations beyond the seeds are spent in all
        (None: the options'), or `report` stops it; returns an
        `ArchiveResult`. The first call evaluates the seeds: one must be
        within the limits. Called again it goes on, and a run stopped and
        gone on with ends as the one that never stopped, bit for bit. A
        larger `evals` extends a finished run.

        `report` is called after each iteration with an `ArchiveIteration`;
        a true return stops the run there. It may read the run and `save`
        it: the run is then as after that iteration.

        An exception raised by `evaluate` or by `report` comes out of this
        call. The run and its model are then as they were after the last
        whole iteration, and `run` can be called again: it goes on to the
        `evals` last asked for, the failed call's included.
        """
        callback = self._report(report, _c.ArchiveReportFn, _archive_iteration)
        r = _c.ArchiveResult()
        _check(_lib.solver2_archive_run(self._made, evals or 0, callback, None, C.byref(r)), self, self.problem)
        return _archive_result(r)

    @property
    def result(self):
        """Where the run stands: an `ArchiveResult`."""
        r = _c.ArchiveResult()
        _check(_lib.solver2_archive_get(self._made, C.byref(r)))
        return _archive_result(r)

    def elites(self):
        """Every cell at once: an `Elites` of arrays, cell c at index c."""
        p = self.problem
        cells = p._cells
        occupied = np.zeros(cells, dtype=np.uint8)
        value, found = np.empty(cells), np.zeros(cells, dtype=np.uintp)
        genes, violations, descriptors = (np.empty((cells, k)) for k in (p._n, p._limits, p._nd))
        outcomes = (_c.Outcome * cells)()
        _check(_lib.solver2_archive_elites(
            self._made, occupied.ctypes.data_as(C.POINTER(C.c_ubyte)), _pointer(value),
            found.ctypes.data_as(C.POINTER(C.c_size_t)), _pointer(genes), outcomes, _pointer(violations),
            _pointer(descriptors)))
        # The structs' layout from their ctypes declaration, padding and
        # all. (numpy's own reading of a ctypes array of structs goes by a
        # format string that Python before 3.11 gets wrong for a padded
        # struct, and numpy then guesses, with a warning.)
        e = np.frombuffer(outcomes, dtype=np.dtype(_c.Outcome))
        return Elites(occupied.astype(bool), value, found, genes, e["objective"].copy(), e["scored"] != 0,
                      e["hard"].copy(), e["terminal"].copy(), e["soft"].copy(), violations, descriptors)

    @property
    def thresholds(self):
        """Each cell's threshold."""
        t = np.empty(self.problem._cells)
        _check(_lib.solver2_archive_thresholds(self._made, _pointer(t)))
        return t

    @property
    def causes(self):
        """Per limit, by its name, the candidates over it so far."""
        return self._counts(_lib.solver2_archive_causes)

    def diverse(self, within, count):
        """The cells to show: the best, then again and again the one
        farthest from those picked, among the cells at most `within` below
        the best. At most `count`: none for 0."""
        cells = (C.c_size_t * max(count, 1))()
        picked = C.c_size_t()
        _check(_lib.solver2_archive_diverse(self._made, within, count, cells, C.byref(picked)))
        return list(cells[:picked.value])

    def save(self):
        """The run as bytes, its model's state included: for `load`."""
        return self._bytes(_lib.solver2_archive_save)
