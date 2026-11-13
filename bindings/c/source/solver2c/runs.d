/// What the C interface's handles stand for besides a problem: a learned
/// model, a CMA-ES stage and a run of the archive, each over a
/// `ForeignProblem`. They are solver2.dynamic's entry points with what a
/// caller in another language can't hold: the state between calls, a report
/// per generation as a plain struct, and the whole of it as bytes.
///
/// Two things here are the binding's own.
///
/// **A failed evaluation leaves a run whole.** The caller's evaluation can
/// fail in the middle of a run (in Python, any exception, a Ctrl-C). The
/// run and its model are then put back as they were after the last whole
/// generation or iteration, so the same call made again goes on, and ends
/// as the run that never failed, bit for bit. A stage keeps its state
/// after every generation as it is; an archive is copied after every
/// iteration (`RunData`), and a model saved (its bytes), since the searches
/// train it and sample from the emitters before they evaluate.
///
/// **The bytes** of a saved run are written by checkpoint.d, a copy of the
/// manual's template: the problem's fingerprint, the run's state as the
/// searches give it, and the model's.
module solver2c.runs;

import std.conv : text;
import std.exception : enforce;
import std.math : isFinite;
import std.parallelism : TaskPool;

import solver2.archive : Archive, ArchiveOptions, ArchiveRules, RestartReason, restoreRun, runState, RunState,
                         selectDiverse;
import solver2.bytes : Bytes;
import solver2.cmaes : CmaEs;
import solver2.dynamic;
import solver2.ensemble : EnsembleOptions;
import solver2.problem;
import solver2.stage : stalled, StageEnd, StageOptions, StageRun;

static import solver2c.checkpoint;
import solver2c.c;
import solver2c.foreign;

/// The saved bytes' magic strings and format version. Bump the version
/// when a saved field changes meaning or type, and teach the loaders the
/// older one; adding a field needs no bump (see checkpoint.d).
enum string STAGE_MAGIC = "solver2-c stage";
enum string ARCHIVE_MAGIC = "solver2-c archive";    /// ditto
enum uint FORMAT = 1;                               /// ditto

// ─────────────────────────────────────────────────────────────────────────────
// The model
// ─────────────────────────────────────────────────────────────────────────────

/// A model as a saved run holds it: how it was made, to check the model it
/// is restored into, and its state.
struct SavedModel
{
    bool present;
    ulong seed;
    EnsembleOptions options;
    size_t rows;
    ubyte[] state;
}

/// A learned model of a foreign problem, with worker threads of its own.
final class ForeignModel
{
    ForeignProblem problem;
    DynamicModel learned;
    immutable ulong seed;
    immutable EnsembleOptions options;      // its seed is `seed`

    this(ForeignProblem problem, ulong seed, const ref solver2_model_options o)
    {
        enforce(o.members > 0 && o.hidden > 0 && o.batch > 0,
                "a model needs members, hidden units and a batch size above 0");
        enforce(isFinite(o.learning_rate) && o.learning_rate > 0 && isFinite(o.l2) && o.l2 >= 0,
                "a model needs a learning rate above 0 and an l2 not below 0");
        EnsembleOptions e = {members: o.members, hidden: o.hidden, batch: o.batch, learningRate: o.learning_rate,
                             l2: o.l2, seed: seed};
        this.problem = problem;
        this.seed = seed;
        options = e;
        if (o.threads > 0)
        {
            pool = new TaskPool(o.threads);
            pool.isDaemon = true;           // they never keep the process from ending
        }
        learned = new DynamicModel(problem, seed, pool, e);
    }

    /// Stops its worker threads. The model works on without them, in the
    /// runs that hold it too: with the same bits, on the calling thread.
    void close()
    {
        learned.usePool(null);      // first: the model holds the pool as well
        if (pool !is null)
            pool.finish(true);
        pool = null;
    }

    /// Its state, to put back with `restore`.
    SavedModel save() const
    {
        return SavedModel(true, seed, options, learned.rows, learned.save(learned.rows));
    }

    /// Puts `s` back: all of it, or nothing and a throw.
    void restore(const ref SavedModel s)
    {
        enforce(s.present, "there is no model's state to restore");
        enforce(s.seed == seed && s.options == options,
                "the model given was made with another seed or other options than the saved run's");
        learned.restore(s.state, s.rows);
    }

private:
    TaskPool pool;
}

/// `model`'s state, or the note that there is none.
private SavedModel savedOf(const ForeignModel model)
{
    return model is null ? SavedModel.init : model.save();
}

/// Whether `lambda` is a population CMA-ES can rank: with one candidate
/// it has no better half to move to, and samples NaN. `orDefault`: 0 too,
/// for the archive's "CMA-ES's own".
private bool isPopulation(size_t lambda, bool orDefault = false) pure nothrow @nogc @safe
{
    return lambda >= 2 || (orDefault && lambda == 0);
}

/// Whether a run with a model can sample `pool` times what it evaluates:
/// twice at least. With 1 nothing is ranked on the model, and a stage
/// would make its CMA-ES anew every generation.
private bool isPool(size_t pool) pure nothrow @nogc @safe
{
    return pool >= 2;
}

/// Whether `share` is one: from 0 to 1 (a NaN is not). The searches take
/// that many of a population by it, without asking.
private bool isShare(double share) pure nothrow @nogc @safe
{
    return share >= 0 && share <= 1;
}

/// Checks that a run is given the model its problem and its saved state
/// call for.
private void checkModel(const ForeignProblem problem, const ForeignModel model, bool saved, bool hasSaved)
{
    enforce(model is null || model.problem is problem, "the model given is of another problem");
    if (hasSaved)
        enforce(saved == (model !is null),
                saved ? "the saved run had a model: give one made as it was"
                      : "the saved run had no model: give none");
}

// ─────────────────────────────────────────────────────────────────────────────
// A stage
// ─────────────────────────────────────────────────────────────────────────────

/// A stage as bytes hold it.
struct SavedStage
{
    string[2][] fingerprint;
    StageOptions options;
    double[] start;
    StageRun run;
    DynamicBest best;
    CauseCounts causes;
    bool stopped;
    SavedModel model;
}

/// A CMA-ES stage between calls.
final class ForeignStage
{
    ForeignProblem problem;
    ForeignModel model;
    StageOptions options;           // its budget: the last one asked for
    double[] start;
    StageRun run;                   // after the last whole generation
    DynamicBest best;
    CauseCounts causes;             // over every candidate so far
    bool stopped;                   // the last call was stopped by its report

    this(ForeignProblem problem, const ref solver2_stage_options o, const(double)[] start, ForeignModel model)
    {
        checkModel(problem, model, false, false);
        enforce(isPopulation(o.lambda), "a stage needs a population (lambda) of at least 2");
        enforce(isFinite(o.sigma) && o.sigma > 0, "a stage needs a step size (sigma) above 0");
        enforce(o.budget > 0, "a stage needs a budget");
        enforce(model is null || isPool(o.surrogate_pool), "a stage with a model needs surrogate_pool of at least 2");
        enforce(isShare(o.surrogate_share), "a stage's surrogate_share must be from 0 to 1");
        StageOptions so = {lambda: o.lambda, sigma: o.sigma, budget: o.budget, seed: o.seed,
                           stallGenerations: o.stall_generations, stallTolerance: o.stall_tolerance,
                           surrogatePool: o.surrogate_pool, surrogateShare: o.surrogate_share,
                           surrogateFrom: o.surrogate_from, surrogateEpochs: o.surrogate_epochs,
                           ranking: problem.ranking};
        this.problem = problem;
        this.model = model;
        options = so;
        this.start = start.length ? start.dup : problem.startGenes;
        enforce(this.start.length == problem.dimension,
                text("the stage's start has ", this.start.length, " genes, the problem ", problem.dimension));
        causes = problem.newCauses();
    }

    private this(ForeignProblem problem, ForeignModel model)
    {
        this.problem = problem;
        this.model = model;
    }

    /// Runs until `budget` evaluations are spent in all (0: as before), the
    /// stage stalls or `report` stops it. Throws what the evaluation
    /// throws, the stage and its model then as after their last whole
    /// generation; throws CallbackFailed when the report says so, the
    /// stage then stopped at that generation. Either way the budget asked
    /// for stays the stage's: it is `finished` against that one, and the
    /// next call goes on to it.
    void advance(size_t budget, solver2_stage_report_fn report, void* context)
    {
        if (budget > 0)
            options.budget = budget;
        const o = options;
        stopped = false;            // a report reads and saves a stage that is running
        SavedModel before = savedOf(model);
        solver2_stage_generation last;
        bool pending;
        int failure;

        // Tells the caller of the generation just done; whether to stop.
        bool tell()
        {
            pending = false;
            if (report is null)
                return false;
            immutable int answer = report(context, &last);
            if (answer < 0)
                failure = answer;
            return answer != 0;
        }
        void done(const ref DynamicStageGeneration g)
        {
            causes += g.causes;
            last = solver2_stage_generation(g.generation, g.before + g.genes.length, g.best, g.populationBest, g.sigma,
                                            g.infeasible, g.sampled, g.predictedFeasible);
            pending = true;
        }
        // Asked after every generation but the last, with the stage's state
        // kept: what an evaluation's failure goes back to.
        bool after(const ref StageRun)
        {
            before = savedOf(model);
            return tell();
        }

        StageEnd end;
        try
            end = runStage(problem, best, start, o, model is null ? null : model.learned, &done, &run, &after);
        catch (Exception e)
        {
            // `run` and `best` are those of the last whole generation: a
            // stage changes them only once a generation is evaluated. The
            // model may have trained for the one that failed.
            if (model !is null)
                model.restore(before);
            stopped = false;
            throw e;
        }
        stopped = end.stopped;
        if (pending)
            tell();
        if (failure != 0)
            throw new CallbackFailed(text("the stage's report reported a failure (it returned ", failure, ")"), failure);
    }

    solver2_stage_result result() const
    {
        solver2_stage_result r;
        r.generations = run.generation;
        r.evaluations = run.spent;
        r.stalled = stalled(run.history, options.stallGenerations, options.stallTolerance);
        r.finished = run.generation > 0 && (run.spent >= options.budget || r.stalled);
        r.stopped = stopped;
        r.feasible = best.feasible;
        r.best_value = best.feasible ? value(best.eval.objective, best.eval.breach) : double.nan;
        return r;
    }

    ubyte[] save() const
    {
        const SavedStage s = {fingerprint: problem.fingerprint, options: options, start: start, run: run, best: best,
                              causes: causes, stopped: stopped, model: savedOf(model)};
        Bytes o;
        solver2c.checkpoint.writeHeader(o, STAGE_MAGIC, FORMAT);
        solver2c.checkpoint.write(o, s);
        return o.data;
    }

    /// The stage `save` wrote, on `problem`, its model's state put back
    /// into `model`. Throws when the bytes aren't a stage's, are of another
    /// problem, or don't fit this one.
    static ForeignStage load(ForeignProblem problem, ForeignModel model, const(ubyte)[] bytes)
    {
        auto r = Bytes(bytes.dup);
        solver2c.checkpoint.readHeader(r, STAGE_MAGIC, FORMAT);
        SavedStage s;
        solver2c.checkpoint.read(r, s);
        enforce(r.done, "the saved stage has bytes left over");
        if (const why = fingerprintMismatch(s.fingerprint, problem.fingerprint))
            throw new Exception("the saved stage is of another problem: " ~ why);
        checkModel(problem, model, s.model.present, true);

        immutable size_t n = problem.dimension;
        enforce(s.start.length == n && s.options.ranking == problem.ranking && isPopulation(s.options.lambda)
                && isShare(s.options.surrogateShare) && (model is null || isPool(s.options.surrogatePool)),
                "the saved stage doesn't fit the problem");
        enforce(s.run.history.length == s.run.generation, "the saved stage lacks its history");
        if (s.run.generation > 0)
            checkState(s.run.es, n);
        // solver2_stage_best copies these into the caller's buffers at their own lengths.
        enforce(!s.best.feasible || (s.best.genes.length == n && s.best.eval.descriptors.length == problem.descriptorCount
                                     && s.best.eval.violations.length == problem.constraints.length),
                "the saved stage's best doesn't fit the problem");
        enforce(s.causes.over.length == problem.constraints.length, "the saved stage's cause counts don't fit the problem");

        auto stage = new ForeignStage(problem, model);
        if (model !is null)
            model.restore(s.model);         // last: all or nothing
        stage.options = s.options;
        stage.start = s.start;
        stage.run = s.run;
        stage.best = s.best;
        stage.causes = s.causes;
        stage.stopped = s.stopped;
        return stage;
    }
}

/// Throws unless a saved CMA-ES state is one of `n` genes: the searches
/// index it without asking.
private void checkState(const ref CmaEs.State s, size_t n)
{
    bool square(const(double[])[] m)
    {
        if (m.length != n)
            return false;
        foreach (row; m)
            if (row.length != n)
                return false;
        return true;
    }
    enforce(s.lambda > 0 && s.mean.length == n && s.pc.length == n && s.ps.length == n && s.D.length == n
            && (s.bestX.length == 0 || s.bestX.length == n) && square(s.C) && square(s.B),
            text("the saved state holds a CMA-ES of another shape than the problem's ", n, " genes"));
}

// ─────────────────────────────────────────────────────────────────────────────
// The archive
// ─────────────────────────────────────────────────────────────────────────────

/// A run of the archive as plain data: what the searches give to save, and
/// what a failed iteration goes back to.
struct RunData
{
    RunState!ArchiveOptions state;
    Archive!Outcome archive;
    CauseCounts causes;
}

/// `run` as plain data that shares nothing with it that the run changes: a
/// cell's elite is replaced whole, never changed in place, so the cells are
/// copied and what they hold is shared.
private RunData dataOf(const ref DynamicRun run)
{
    auto r = cast(DynamicRun*) &run;        // read only
    RunData d;
    d.state = runState(*r);
    d.archive = r.archive;
    d.archive.elites = r.archive.elites.dup;
    d.archive.threshold = r.archive.threshold.dup;
    d.causes = CauseCounts(r.causes.over.dup);
    return d;
}

/// Puts `d` back into `run`.
private void restore(ref DynamicRun run, ref RunData d)
{
    restoreRun(run, d.state);
    run.archive = d.archive;
    run.archive.elites = d.archive.elites.dup;
    run.archive.threshold = d.archive.threshold.dup;
    run.causes = CauseCounts(d.causes.over.dup);
    run.stopped = false;
}

/// An archive's run as bytes hold it.
struct SavedArchive
{
    string[2][] fingerprint;
    ArchiveOptions options;
    double[][] seeds;
    bool started;
    bool stopped;
    RunData run;
    SavedModel model;
}

/// A run of the archive between calls.
final class ForeignArchive
{
    ForeignProblem problem;
    ForeignModel model;
    ArchiveOptions options;         // its budget: the last one asked for
    double[][] seeds;
    DynamicRun run;                 // after the last whole iteration, once started
    bool started;                   // the seeds are evaluated
    bool stopped;                   // the last call was stopped by its report

    this(ForeignProblem problem, const ref solver2_archive_options o, const(double[])[] seeds, ForeignModel model)
    {
        checkModel(problem, model, false, false);
        enforce(problem.axes.length > 0, "the archive needs a problem with axes");
        enforce(o.emitters > 0, "the archive needs emitters");
        enforce(isPopulation(o.lambda, true),
                "the archive needs a population (lambda) of at least 2 per emitter, or 0 for CMA-ES's own");
        enforce(isFinite(o.sigma) && o.sigma > 0 && isFinite(o.explorer_sigma) && o.explorer_sigma > 0,
                "the archive needs step sizes (sigma, explorer_sigma) above 0");
        enforce(o.alpha >= 0 && o.alpha <= 1, "the archive's alpha must be from 0 to 1");
        enforce(o.surrogate_pool == 0 || isPool(o.surrogate_pool),
                "the archive's surrogate_pool must be 0, or at least 2 with a model");
        enforce((model !is null) == (o.surrogate_pool > 0),
                "an archive needs a model when surrogate_pool is above 0, and only then");
        enforce(isShare(o.surrogate_share), "the archive's surrogate_share must be from 0 to 1");
        const r = &o.rules;
        ArchiveRules rules = {idleGenerations: r.idle_generations, infeasibleGenerations: r.infeasible_generations,
                              infeasibleShare: r.infeasible_share, minSigma: r.min_sigma,
                              staleGenerations: r.stale_generations, minProgress: r.min_progress,
                              blockIterations: r.block_iterations, aimSamples: r.aim_samples};
        rules.check();
        ArchiveOptions ao = {evals: o.evals, emitters: o.emitters, lambda: o.lambda, sigma: o.sigma, alpha: o.alpha,
                             seed: o.seed, explorers: o.explorers, explorerSigma: o.explorer_sigma,
                             surrogatePool: o.surrogate_pool, surrogateShare: o.surrogate_share,
                             surrogateFrom: o.surrogate_from, surrogateEpochs: o.surrogate_epochs,
                             surrogateImprovers: o.surrogate_improvers != 0, surrogateAim: o.surrogate_aim != 0,
                             ranking: problem.ranking, rules: rules};
        this.problem = problem;
        this.model = model;
        options = ao;
        foreach (i, s; seeds)
        {
            enforce(s.length == problem.dimension,
                    text("seed ", i, " of the archive has ", s.length, " genes, the problem ", problem.dimension));
            this.seeds ~= s.dup;
        }
        if (this.seeds.length == 0)
            this.seeds = [problem.startGenes];
    }

    private this(ForeignProblem problem, ForeignModel model)
    {
        this.problem = problem;
        this.model = model;
    }

    /// Runs until `evals` evaluations beyond the seeds are spent in all (0:
    /// as before) or `report` stops it; the first call evaluates the seeds.
    /// Throws what the evaluation throws, the run and its model then as
    /// after their last whole iteration (not started, when it was the seeds
    /// or the first iteration that failed); throws CallbackFailed when the
    /// report says so, the run then stopped at that iteration. Either way
    /// the budget asked for stays the run's: it is `finished` against that
    /// one, and the next call goes on to it.
    ///
    /// The report finds the run as after the iteration it is told of, in
    /// the first call too.
    void advance(size_t evals, solver2_archive_report_fn report, void* context)
    {
        if (evals > 0)
            options.evals = evals;
        const o = options;
        stopped = false;            // a report reads and saves a run that is running
        // The first call's run is the searches' own until they return.
        immutable bool first = !started;
        SavedModel modelBefore = savedOf(model);
        RunData before;
        bool haveBefore = started;
        if (started)
            before = dataOf(run);
        solver2_archive_iteration it;
        bool pending;
        int failure;

        bool tell()
        {
            pending = false;
            if (report is null)
                return false;
            immutable int answer = report(context, &it);
            if (answer < 0)
                failure = answer;
            return answer != 0;
        }
        // One emitter's generation: the iteration's report adds them up,
        // and the last emitter's has the archive as the iteration left it.
        void progress(const ref DynamicArchiveGeneration g)
        {
            if (g.emitter == 0)
                it = solver2_archive_iteration.init;
            it.iteration = g.iteration;
            it.evaluations = g.evals;
            it.coverage = g.cells;
            it.qd_score = g.qdScore;
            it.best_score = g.bestScore;
            it.improved += g.improved;
            it.new_cells += g.newCells;
            it.restarts += g.restartReason != RestartReason.none;
            pending = true;
        }
        // Asked after every iteration but the last, with the run as it
        // would be saved: what an evaluation's failure goes back to.
        bool after(const ref DynamicRun r)
        {
            before = dataOf(r);
            haveBefore = true;
            modelBefore = savedOf(model);
            if (first)
            {
                // A copy for the report to read and save, as a later
                // call's finds the run itself.
                restore(run, before);
                started = true;
            }
            return tell();
        }

        try
        {
            auto learned = model is null ? null : model.learned;
            if (started)
                runArchive(problem, null, o, &progress, &run, &after, null, learned);
            else
                run = runArchive(problem, seeds, o, &progress, null, &after, null, learned);
        }
        catch (Exception e)
        {
            // The searches train the model and sample from the emitters
            // before they evaluate: put both back.
            if (haveBefore)
                restore(run, before);
            else
                run = DynamicRun.init;
            started = haveBefore;
            if (model !is null)
                model.restore(modelBefore);
            stopped = false;
            throw e;
        }
        started = true;
        stopped = run.stopped;
        if (pending)
            tell();
        if (failure != 0)
            throw new CallbackFailed(text("the archive's report reported a failure (it returned ", failure, ")"),
                                     failure);
    }

    solver2_archive_result result() const
    {
        solver2_archive_result r;
        r.cells = problem.grid.cells;
        r.started = started;
        r.stopped = stopped;
        r.best_score = double.nan;
        r.best_cell = size_t.max;
        r.qd_score = 0;
        if (!started)
            return r;
        r.coverage = run.archive.coverage;
        r.evaluations = run.evals;
        r.seed_evaluations = run.seedEvals;
        r.iterations = run.iterations;
        r.lambda = run.lambda;
        r.finished = run.evals >= run.seedEvals + options.evals;
        r.qd_score = run.archive.qdScore;
        r.best_cell = run.archive.best;
        if (r.best_cell != size_t.max)
            r.best_score = run.archive.elites[r.best_cell].value;
        r.restarts_idle = run.restartsIdle;
        r.restarts_infeasible = run.restartsInfeasible;
        r.restarts_converged = run.restartsConverged;
        r.reached = run.reached.length;
        r.abandoned = run.abandoned.length;
        return r;
    }

    /// The cells to show (`selectDiverse`); none before the run has started.
    size_t[] diverse(double within, size_t count) const
    {
        if (!started)
            return null;
        const grid = problem.grid;
        return selectDiverse(run.archive, grid, within, count);
    }

    ubyte[] save() const
    {
        SavedArchive s = {fingerprint: problem.fingerprint.dup, options: options, started: started, stopped: stopped,
                          model: savedOf(model)};
        foreach (seed; seeds)
            s.seeds ~= seed.dup;
        if (started)
            s.run = dataOf(run);
        Bytes o;
        solver2c.checkpoint.writeHeader(o, ARCHIVE_MAGIC, FORMAT);
        solver2c.checkpoint.write(o, s);
        return o.data;
    }

    /// The run `save` wrote, on `problem`, its model's state put back into
    /// `model`. Throws when the bytes aren't an archive's, are of another
    /// problem, or don't fit this one.
    static ForeignArchive load(ForeignProblem problem, ForeignModel model, const(ubyte)[] bytes)
    {
        auto r = Bytes(bytes.dup);
        solver2c.checkpoint.readHeader(r, ARCHIVE_MAGIC, FORMAT);
        SavedArchive s;
        solver2c.checkpoint.read(r, s);
        enforce(r.done, "the saved archive has bytes left over");
        if (const why = fingerprintMismatch(s.fingerprint, problem.fingerprint))
            throw new Exception("the saved archive is of another problem: " ~ why);
        checkModel(problem, model, s.model.present, true);

        immutable size_t n = problem.dimension, cells = problem.grid.cells, nd = problem.descriptorCount;
        enforce(cells > 0 && s.options.ranking == problem.ranking && isPopulation(s.options.lambda, true)
                && isShare(s.options.surrogateShare) && (model is null || isPool(s.options.surrogatePool)),
                "the saved archive doesn't fit the problem");
        foreach (seed; s.seeds)
            enforce(seed.length == n, "the saved archive's seeds don't fit the problem");
        if (s.started)
        {
            const st = &s.run.state;
            const a = &s.run.archive;
            enforce(a.elites.length == cells && a.threshold.length == cells && st.blockedUntil.length == cells,
                    text("the saved archive has ", a.elites.length, " cells, the problem's grid ", cells));
            foreach (ref el; a.elites)
                enforce(!el.occupied || (el.genes.length == n && el.eval.descriptors.length == nd),
                        "the saved archive holds a candidate that doesn't fit the problem");
            enforce(st.lambda > 0 && st.emitters.length == st.options.emitters,
                    "the saved archive lacks its emitters");
            foreach (ref e; st.emitters)
            {
                checkState(e.es, n);
                enforce(e.target == size_t.max || e.target < cells, "the saved archive has an explorer aimed off the grid");
            }
            foreach (list; [st.reached, st.abandoned, st.focus])
                foreach (c; list)
                    enforce(c < cells, "the saved archive names a cell off the grid");
            enforce(s.run.causes.over.length == problem.constraints.length,
                    "the saved archive's cause counts don't fit the problem");
        }

        auto archive = new ForeignArchive(problem, model);
        if (model !is null)
            model.restore(s.model);         // last: all or nothing
        archive.options = s.options;
        archive.seeds = s.seeds;
        archive.started = s.started;
        archive.stopped = s.stopped;
        if (s.started)
            restore(archive.run, s.run);
        return archive;
    }
}
