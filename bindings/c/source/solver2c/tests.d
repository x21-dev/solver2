/// The C interface's tests, through its exported functions, in one process
/// with the D searches: so that what comes out of the interface can be held
/// to what the searches give natively, bit for bit. (The library as a
/// caller in another language meets it, loaded as a shared library with no
/// D runtime around it, is test/source/driver.d's part.)
///
/// The problem is one function of the genes, given twice: as a C callback
/// behind the interface, and as a `DynamicProblem` subclass for the native
/// runs.
module solver2c.tests;

version (unittest):

import std.string : fromStringz, indexOf;

import solver2.archive : ArchiveOptions, selectDiverse;
import solver2.dynamic;
import solver2.ensemble : EnsembleOptions;
import solver2.problem;
import solver2.stage : StageOptions, StageRun;

import solver2c.api;
import solver2c.c;

shared static this()
{
    assert(solver2_init() == SOLVER2_OK);
}

// ─────────────────────────────────────────────────────────────────────────────
// The problem, twice
// ─────────────────────────────────────────────────────────────────────────────

private static immutable GeneSpec[3] fieldGenes = [GeneSpec("x", 0, 2, 1), GeneSpec("y", 0, 2, 0.6),
                                                   GeneSpec("z", -1, 1, 0)];
private static immutable ConstraintSpec[3] fieldLimits = [ConstraintSpec("sum", Tier.hard),
                                                          ConstraintSpec("floor", Tier.terminal),
                                                          ConstraintSpec("x_cap", Tier.soft, 0.5)];
private static immutable Axis[2] fieldAxes = [Axis("x", [0.5, 1, 1.5], "abcd", 0, 0.5),
                                              Axis("y", [0.5, 1, 1.5], "abcd", 0, 0.5)];
private enum Ranking fieldRanking = Ranking(1000, 500, 20);
private enum double fieldMinScore = 5, fieldHeadroom = 0.05;

/// One candidate: maximise 10 - (x - 1)² - (y - 0.8)² - z², with
/// x + y ≤ 1.9 hard, z ≥ -0.9 terminal and x ≤ 1.5 soft. The descriptors
/// are x and y.
private struct Row
{
    double x, y;
    double[3] violations;
    double objective;
}

private Row fieldRow(const(double)[] g) pure nothrow @nogc @safe
{
    immutable double x = fieldGenes[0].toPhysical(g[0]), y = fieldGenes[1].toPhysical(g[1]),
                     z = fieldGenes[2].toPhysical(g[2]);
    return Row(x, y, [x + y > 1.9 ? x + y - 1.9 : 0, z < -0.9 ? -0.9 - z : 0, x > 1.5 ? x - 1.5 : 0],
               10 - (x - 1) * (x - 1) - (y - 0.8) * (y - 0.8) - z * z);
}

/// The callback's context: its batches counted, one of them made to fail,
/// and a drift that makes the evaluation depend on what came before.
private struct Calls
{
    size_t batches, rows;
    size_t failAt;
    double drift = 0;
    double sink = 0;        // lowers every objective: below the problem's min_score
}

private extern (C) int fieldEvaluate(void* context, size_t rows, const(double)* genes, int feasibleOnly, ubyte* scored,
                                     double* objective, double* violations, double* descriptors)
{
    auto c = cast(Calls*) context;
    ++c.batches;
    if (c.batches == c.failAt)
        return 7;
    c.rows += rows;
    foreach (k; 0 .. rows)
    {
        const r = fieldRow(genes[3 * k .. 3 * k + 3]);
        violations[3 * k .. 3 * k + 3] = r.violations[];
        if (!feasibleOnly || r.violations[0] == 0)
        {
            scored[k] = 1;
            objective[k] = r.objective + c.drift * c.batches - c.sink;
        }
        descriptors[2 * k] = r.x;
        descriptors[2 * k + 1] = r.y;
    }
    return 0;
}

/// The same problem for the native runs.
private final class NativeField : DynamicProblem
{
    override const(GeneSpec)[] genes() const { return fieldGenes[]; }
    override const(ConstraintSpec)[] constraints() const { return fieldLimits[]; }
    override const(string[2])[] fingerprint() const { return [["problem", "field 1"]]; }
    override const(Axis)[] axes() const { return fieldAxes[]; }
    override Ranking ranking() const { return fieldRanking; }
    override double minScore() const { return fieldMinScore; }
    override double predictionHeadroom() const { return fieldHeadroom; }

    override void evaluateInto(const(double[])[] g, Scoring scoring, Outcome[] o)
    {
        foreach (k, row; g)
        {
            const r = fieldRow(row);
            o[k].violations = r.violations.dup;
            o[k].breach = breachOf(fieldLimits[], o[k].violations);
            o[k].scored = scoring == Scoring.all || r.violations[0] == 0;
            if (o[k].scored)
                o[k].objective = r.objective;
            o[k].descriptors = [r.x, r.y];
        }
    }
}

/// The problem's description for the interface. `edit` changes it before
/// the problem is made.
private solver2_problem_spec fieldSpec(Calls* calls, const(char)* revision = "field 1")
{
    static immutable solver2_gene[3] genes = [solver2_gene("x", 0, 2, 1, SOLVER2_LINEAR),
                                              solver2_gene("y", 0, 2, 0.6, SOLVER2_LINEAR),
                                              solver2_gene("z", -1, 1, 0, SOLVER2_LINEAR)];
    static immutable solver2_limit[3] limits = [solver2_limit("sum", SOLVER2_HARD, 1, null),
                                                solver2_limit("floor", SOLVER2_TERMINAL, 1, null),
                                                solver2_limit("x_cap", SOLVER2_SOFT, 0.5, "caps")];
    static immutable double[3] edges = [0.5, 1, 1.5];
    static immutable solver2_axis[2] axes = [immutable solver2_axis("x", edges.ptr, 3, "abcd", 0, 0.5),
                                             immutable solver2_axis("y", edges.ptr, 3, "abcd", 0, 0.5)];
    solver2_problem_spec s;
    ok(solver2_problem_spec_default(&s));
    assert(s.ranking == solver2_ranking(100, 50, 10) && s.min_score == 0 && s.genes is null && s.evaluate is null);
    s.genes = genes.ptr;
    s.gene_count = 3;
    s.limits = limits.ptr;
    s.limit_count = 3;
    s.axes = axes.ptr;
    s.axis_count = 2;
    s.ranking = solver2_ranking(fieldRanking.hardOffset, fieldRanking.terminalOffset, fieldRanking.boundPenalty);
    s.min_score = fieldMinScore;
    s.prediction_headroom = fieldHeadroom;
    auto pairs = new const(char)*[2];
    pairs[0] = "problem";
    pairs[1] = revision;
    s.fingerprint = pairs.ptr;
    s.fingerprint_count = 1;
    s.evaluate = &fieldEvaluate;
    s.context = calls;
    return s;
}

private solver2_problem* field(Calls* calls, const(char)* revision = "field 1")
{
    const s = fieldSpec(calls, revision);
    solver2_problem* p;
    ok(solver2_problem_new(&s, &p));
    return p;
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

private string lastError()
{
    return fromStringz(solver2_last_error()).idup;
}

/// The call succeeded.
private void ok(int code, string file = __FILE__, size_t line = __LINE__)
{
    if (code != SOLVER2_OK)
        throw new Exception("the call failed: " ~ lastError(), file, line);
}

/// The call failed with `code` and a message that holds `part`.
private void refused(int got, string part, int code = SOLVER2_ERROR, string file = __FILE__, size_t line = __LINE__)
{
    import std.conv : text;
    if (got != code || lastError().indexOf(part) < 0)
        throw new Exception(text("expected code ", code, " and \"", part, "\"; got ", got, " and \"", lastError(), "\""),
                            file, line);
}

private bool sameBits(const(double)[] a, const(double)[] b)
{
    if (a.length != b.length)
        return false;
    foreach (i; 0 .. a.length)
        if (a[i] !is b[i])
            return false;
    return true;
}

/// A report that keeps every generation, and may stop or fail at one.
private struct StageTrace
{
    solver2_stage_generation[] generations;
    size_t stopAt, failAt;
}

private extern (C) int stageReport(void* context, const(solver2_stage_generation)* g)
{
    auto t = cast(StageTrace*) context;
    t.generations ~= *g;
    return g.generation == t.failAt ? -5 : g.generation == t.stopAt;
}

private struct ArchiveTrace
{
    solver2_archive_iteration[] iterations;
    size_t stopAt, failAt;
}

private extern (C) int archiveReport(void* context, const(solver2_archive_iteration)* it)
{
    auto t = cast(ArchiveTrace*) context;
    t.iterations ~= *it;
    return it.iteration == t.failAt ? -5 : it.iteration == t.stopAt;
}

private solver2_stage_options stageOptions(bool screened = false)
{
    solver2_stage_options o;
    ok(solver2_stage_options_default(&o));
    assert(o.lambda == 0 && o.budget == 0 && o.stall_generations == 100 && o.stall_tolerance == 0.01
           && o.surrogate_pool == 0 && o.surrogate_share == 0.5 && o.surrogate_from == 5 && o.surrogate_epochs == 10);
    o.lambda = 8;
    o.sigma = 0.2;
    o.budget = 400;
    o.seed = 3;
    o.stall_generations = 0;
    if (screened)
    {
        o.budget = 240;
        o.surrogate_pool = 4;
        o.surrogate_from = 3;
        o.surrogate_epochs = 2;
    }
    return o;
}

private StageOptions nativeOf(const solver2_stage_options o)
{
    StageOptions s = {lambda: o.lambda, sigma: o.sigma, budget: o.budget, seed: o.seed,
                      stallGenerations: o.stall_generations, stallTolerance: o.stall_tolerance,
                      surrogatePool: o.surrogate_pool, surrogateShare: o.surrogate_share,
                      surrogateFrom: o.surrogate_from, surrogateEpochs: o.surrogate_epochs};
    return s;
}

private solver2_archive_options archiveOptions(bool guided = false)
{
    solver2_archive_options o;
    ok(solver2_archive_options_default(&o));
    immutable ArchiveOptions d;
    assert(o.evals == d.evals && o.emitters == d.emitters && o.lambda == 0 && o.sigma == d.sigma && o.alpha == d.alpha
           && o.seed == 1 && o.explorers == 2 && o.explorer_sigma == d.explorerSigma && o.surrogate_pool == 0
           && !o.surrogate_improvers && !o.surrogate_aim && o.rules.idle_generations == 5
           && o.rules.infeasible_share == 0.9 && o.rules.min_sigma == 1e-4 && o.rules.aim_samples == 200);
    o.evals = 400;
    o.emitters = 3;
    o.explorers = 1;
    o.lambda = 8;
    o.seed = 5;
    if (guided)
    {
        o.surrogate_pool = 4;
        o.surrogate_from = 3;
        o.surrogate_epochs = 2;
        o.surrogate_improvers = 1;
        o.surrogate_aim = 1;
    }
    return o;
}

private ArchiveOptions nativeOf(const solver2_archive_options o)
{
    ArchiveOptions a = {evals: o.evals, emitters: o.emitters, lambda: o.lambda, sigma: o.sigma, alpha: o.alpha,
                        seed: o.seed, explorers: o.explorers, explorerSigma: o.explorer_sigma,
                        surrogatePool: o.surrogate_pool, surrogateShare: o.surrogate_share,
                        surrogateFrom: o.surrogate_from, surrogateEpochs: o.surrogate_epochs,
                        surrogateImprovers: o.surrogate_improvers != 0, surrogateAim: o.surrogate_aim != 0};
    return a;
}

/// A small model: two members of sixteen units, on `threads` threads.
private solver2_model* smallModel(solver2_problem* p, size_t threads = 0, ulong seed = 7)
{
    solver2_model_options o;
    ok(solver2_model_options_default(&o));
    assert(o.members == 5 && o.hidden == 64 && o.batch == 128 && o.learning_rate == 1e-3 && o.l2 == 1e-3
           && o.threads == 0);
    o.members = 2;
    o.hidden = 16;
    o.threads = threads;
    solver2_model* m;
    ok(solver2_model_new(p, seed, &o, &m));
    return m;
}

private ubyte[] bytesOf(solver2_stage* s)
{
    size_t size;
    ok(solver2_stage_save(s, null, 0, &size));
    auto b = new ubyte[size];
    ok(solver2_stage_save(s, b.ptr, b.length, &size));
    assert(size == b.length);
    return b;
}

private ubyte[] bytesOf(solver2_archive* a)
{
    size_t size;
    ok(solver2_archive_save(a, null, 0, &size));
    auto b = new ubyte[size];
    ok(solver2_archive_save(a, b.ptr, b.length, &size));
    assert(size == b.length);
    return b;
}

/// A stage run to its end in one call: its bytes.
private ubyte[] wholeStage(bool screened, size_t threads = 0)
{
    Calls calls;
    auto p = field(&calls);
    const o = stageOptions(screened);
    auto m = screened ? smallModel(p, threads) : null;
    solver2_stage* s;
    ok(solver2_stage_new(p, &o, null, m, &s));
    solver2_stage_result r;
    ok(solver2_stage_run(s, 0, null, null, &r));
    assert(r.finished && !r.stopped && r.evaluations == o.budget && calls.rows == o.budget);
    auto b = bytesOf(s);
    ok(solver2_stage_free(s));
    if (m !is null)
        ok(solver2_model_free(m));
    ok(solver2_problem_free(p));
    return b;
}

private ubyte[] wholeArchive(bool guided, size_t threads = 0)
{
    Calls calls;
    auto p = field(&calls);
    const o = archiveOptions(guided);
    auto m = guided ? smallModel(p, threads) : null;
    solver2_archive* a;
    ok(solver2_archive_new(p, &o, null, 0, m, &a));
    solver2_archive_result r;
    ok(solver2_archive_run(a, 0, null, null, &r));
    assert(r.started && r.finished && !r.stopped && r.evaluations == calls.rows);
    auto b = bytesOf(a);
    ok(solver2_archive_free(a));
    if (m !is null)
        ok(solver2_model_free(m));
    ok(solver2_problem_free(p));
    return b;
}

// ─────────────────────────────────────────────────────────────────────────────
// The library and its handles
// ─────────────────────────────────────────────────────────────────────────────

unittest
{
    // The library says what it is; init and term are counted.
    assert(solver2_abi() == SOLVER2_ABI && solver2_checked() == 0);
    assert(solver2_init() == SOLVER2_OK && solver2_term() == SOLVER2_OK);

    // A problem from its description: what a caller can ask of it.
    Calls calls;
    auto p = field(&calls);
    size_t n, cells;
    ok(solver2_problem_dimension(p, &n));
    ok(solver2_problem_cells(p, &cells));
    double[3] start;
    ok(solver2_problem_start(p, start.ptr));
    assert(n == 3 && cells == 16 && start == [0.5, 0.3, 0.5] && calls.batches == 0);
    // Genes as the values they stand for: GeneSpec's own, clamped.
    double[6] sampled = [0.25, 1.5, 0.5, -0.1, 0.3, 1.0], physical;
    ok(solver2_problem_physical(p, 2, sampled.ptr, physical.ptr));
    assert(physical == [0.5, 2, 0, 0, 0.6, 1]);
    foreach (k; 0 .. 2)
        foreach (i; 0 .. 3)
            assert(physical[3 * k + i] is fieldGenes[i].toPhysical(sampled[3 * k + i]));
    ok(solver2_problem_physical(p, 0, null, null));
    refused(solver2_problem_physical(p, 2, sampled.ptr, null), "values is NULL");
    char[3] code;
    ok(solver2_problem_cell_code(p, 14, code.ptr, code.length));
    assert(code == "cd\0");
    refused(solver2_problem_cell_code(p, 14, code.ptr, 2), "the buffer holds 2 bytes");
    refused(solver2_problem_cell_code(p, 16, code.ptr, code.length), "cell 16 is not one of the grid's 16");
    ok(solver2_problem_check(p, 8, 1));
    assert(calls.batches > 0);

    // A handle is good until it is freed, and only for its kind.
    refused(solver2_problem_dimension(null, &n), "the handle is NULL");
    refused(solver2_problem_dimension(cast(solver2_problem*) 12_345, &n), "not a handle");
    refused(solver2_model_free(cast(solver2_model*) p), "a handle of another kind");
    refused(solver2_problem_dimension(p, null), "dimension is NULL");
    ok(solver2_problem_free(p));
    refused(solver2_problem_dimension(p, &n), "not a handle");
    refused(solver2_problem_free(p), "not a handle");
    // Its slot is given out again, under another handle.
    auto q = field(&calls);
    assert(q !is p);
    refused(solver2_problem_dimension(p, &n), "not a handle");
    ok(solver2_problem_free(q));
}

unittest
{
    // A description that doesn't hold is refused when the problem is made.
    Calls calls;
    solver2_problem* p;
    auto s = fieldSpec(&calls);
    refused(solver2_problem_new(null, &p), "spec is NULL");
    refused(solver2_problem_new(&s, null), "problem is NULL");

    s.evaluate = null;
    refused(solver2_problem_new(&s, &p), "needs its evaluate function");
    s = fieldSpec(&calls);
    s.gene_count = 0;
    refused(solver2_problem_new(&s, &p), "needs genes");

    auto genes = s.genes[0 .. 3].dup;
    s = fieldSpec(&calls);
    genes[1].scale = 9;
    s.genes = genes.ptr;
    refused(solver2_problem_new(&s, &p), "gene y: its scale 9");
    genes[1] = solver2_gene("y", 0, 2, 0.6, SOLVER2_LOG);
    refused(solver2_problem_new(&s, &p), "gene y: a log scale needs lo above 0");
    genes[1] = solver2_gene("y", 0, 2, 7, SOLVER2_LINEAR);
    refused(solver2_problem_new(&s, &p), "gene y: its neutral value is outside its range");
    genes[1] = solver2_gene(null, 2, 0, 1, SOLVER2_LINEAR);
    refused(solver2_problem_new(&s, &p), "gene gene 1: its range");
    // A gene without a neutral value starts in the middle.
    genes[1] = solver2_gene("y", 0, 2, double.nan, SOLVER2_LINEAR);
    ok(solver2_problem_new(&s, &p));
    double[3] start;
    ok(solver2_problem_start(p, start.ptr));
    assert(start == [0.5, 0.5, 0.5]);
    ok(solver2_problem_free(p));

    s = fieldSpec(&calls);
    auto limits = s.limits[0 .. 3].dup;
    limits[2].tier = 3;
    s.limits = limits.ptr;
    refused(solver2_problem_new(&s, &p), "limit x_cap: its tier 3");
    limits[2] = solver2_limit("x_cap", SOLVER2_SOFT, -1, null);
    refused(solver2_problem_new(&s, &p), "limit x_cap: its weight");

    s = fieldSpec(&calls);
    auto axes = s.axes[0 .. 2].dup;
    axes[0].letters = "abc";
    s.axes = axes.ptr;
    refused(solver2_problem_new(&s, &p), "axis x: 3 letters for 4 bins");
    s = fieldSpec(&calls);
    s.min_score = double.infinity;
    refused(solver2_problem_new(&s, &p), "minScore must be finite");
    s = fieldSpec(&calls);
    s.ranking.hard_offset = 1;
    refused(solver2_problem_new(&s, &p), "its ranking");
    assert(calls.batches == 0);

    // A problem without limits or axes is one too: a stage runs on it (its
    // evaluation is handed NULL for what it hasn't), an archive can't.
    static extern (C) int bare(void*, size_t rows, const(double)* genes, int, ubyte* scored, double* objective,
                               double* violations, double* descriptors)
    {
        if (violations !is null || descriptors !is null)
            return 1;
        foreach (k; 0 .. rows)
        {
            scored[k] = 1;
            objective[k] = 0;
            foreach (v; genes[3 * k .. 3 * k + 3])
                objective[k] -= (clamp01(v) - 0.3) * (clamp01(v) - 0.3);
        }
        return 0;
    }
    s = fieldSpec(&calls);
    s.limit_count = 0;
    s.axis_count = 0;
    s.evaluate = &bare;
    ok(solver2_problem_new(&s, &p));
    size_t cells = 9;
    ok(solver2_problem_cells(p, &cells));
    assert(cells == 0);
    ok(solver2_problem_check(p, 8, 1));
    const o = stageOptions();
    solver2_stage* stage;
    ok(solver2_stage_new(p, &o, null, null, &stage));
    solver2_stage_result r;
    ok(solver2_stage_run(stage, 0, null, null, &r));
    double[3] best;
    ok(solver2_stage_best(stage, best.ptr, null, null, null));
    assert(r.finished && r.feasible && r.best_value > -1e-3 && best[0] > 0.28 && best[0] < 0.32);
    ok(solver2_stage_free(stage));
    solver2_archive* a;
    const ao = archiveOptions();
    refused(solver2_archive_new(p, &ao, null, 0, null, &a), "needs a problem with axes");
    ok(solver2_problem_free(p));
}

unittest
{
    // solver2_problem_check: an evaluation that depends on what came
    // before, and values below min_score, are told so.
    Calls drifting = {drift: 1e-9};
    auto p = field(&drifting);
    refused(solver2_problem_check(p, 8, 1), "other bits when evaluated again");
    ok(solver2_problem_free(p));

    Calls sunk = {sink: 9.5};
    p = field(&sunk);
    refused(solver2_problem_check(p, 8, 1), "below the problem's minScore");
    ok(solver2_problem_free(p));

    // What the evaluation gives is checked in every batch: a violation
    // below 0, and a candidate within the limits without its objective.
    static extern (C) int negative(void*, size_t rows, const(double)*, int, ubyte* scored, double* objective,
                                   double* violations, double* descriptors)
    {
        violations[1] = -0.25;
        return 0;
    }
    static extern (C) int unscored(void*, size_t rows, const(double)*, int, ubyte*, double*, double*, double*)
    {
        return 0;
    }
    Calls calls;
    auto s = fieldSpec(&calls);
    s.evaluate = &negative;
    ok(solver2_problem_new(&s, &p));
    refused(solver2_problem_check(p, 8, 1), "its violation of the limit floor is negative or NaN");
    ok(solver2_problem_free(p));
    s.evaluate = &unscored;
    ok(solver2_problem_new(&s, &p));
    refused(solver2_problem_check(p, 8, 1), "within the limits, but without an objective");
    ok(solver2_problem_free(p));
}

// ─────────────────────────────────────────────────────────────────────────────
// A stage
// ─────────────────────────────────────────────────────────────────────────────

unittest
{
    // A stage through the interface is the stage the searches run
    // natively: every generation's report, the best and its outcome.
    Calls calls;
    auto p = field(&calls);
    const o = stageOptions();
    solver2_stage* s;
    ok(solver2_stage_new(p, &o, null, null, &s));
    solver2_stage_result r;
    ok(solver2_stage_get(s, &r));
    assert(r.generations == 0 && !r.finished && !r.feasible && r.best_value != r.best_value && calls.batches == 0);
    refused(solver2_stage_best(s, null, null, null, null), "no candidate within the limits yet");

    StageTrace trace;
    ok(solver2_stage_run(s, 0, &stageReport, &trace, &r));
    assert(r.generations == 50 && r.evaluations == 400 && r.finished && !r.stalled && !r.stopped && r.feasible);
    assert(calls.batches == 50 && calls.rows == 400 && trace.generations.length == 50);

    DynamicProblem native = new NativeField;
    DynamicBest best;
    solver2_stage_generation[] nativeTrace;
    CauseCounts nativeCauses = native.newCauses();
    runStage(native, best, null, nativeOf(o), null, (const ref DynamicStageGeneration g) {
        nativeTrace ~= solver2_stage_generation(g.generation, g.before + g.genes.length, g.best, g.populationBest,
                                                g.sigma, g.infeasible, g.sampled, g.predictedFeasible);
        nativeCauses += g.causes;
    });
    foreach (i, ref g; trace.generations)
    {
        const n = &nativeTrace[i];
        assert(g.generation == n.generation && g.evaluations == n.evaluations && g.best is n.best
               && g.population_best is n.population_best && g.sigma is n.sigma && g.infeasible == n.infeasible);
    }
    double[3] genes, violations;
    double[2] descriptors;
    solver2_outcome e;
    ok(solver2_stage_best(s, genes.ptr, &e, violations.ptr, descriptors.ptr));
    assert(sameBits(genes[], best.genes) && e.scored && e.objective is best.eval.objective
           && e.value is value(best.eval.objective, best.eval.breach) && r.best_value is e.value);
    assert(e.hard == 0 && e.terminal == 0 && e.soft is best.eval.breach.soft
           && sameBits(violations[], best.eval.violations) && sameBits(descriptors[], best.eval.descriptors));
    size_t[3] over;
    ok(solver2_stage_causes(s, over.ptr));
    assert(over[] == nativeCauses.over && over[0] > 0);

    // Its budget spent, another call does nothing; a larger budget goes on,
    // as the stage given it at the start.
    ok(solver2_stage_run(s, 0, &stageReport, &trace, &r));
    assert(calls.batches == 50 && trace.generations.length == 50 && r.finished);
    ok(solver2_stage_run(s, 480, &stageReport, &trace, &r));
    assert(r.generations == 60 && r.evaluations == 480 && r.finished && trace.generations.length == 60);
    DynamicBest longer;
    auto more = nativeOf(o);
    more.budget = 480;
    runStage(native, longer, null, more);
    ok(solver2_stage_best(s, genes.ptr, null, null, null));
    assert(sameBits(genes[], longer.genes));

    ok(solver2_stage_free(s));
    refused(solver2_stage_get(s, &r), "a stage: not a handle");
    ok(solver2_problem_free(p));
}

unittest
{
    // A stage's own start, and what a stage refuses.
    Calls calls;
    auto p = field(&calls);
    auto o = stageOptions();
    solver2_stage* s;
    double[3] start = [0.2, 0.2, 0.6];
    ok(solver2_stage_new(p, &o, start.ptr, null, &s));
    ok(solver2_stage_run(s, 0, null, null, null));
    DynamicProblem native = new NativeField;
    DynamicBest best;
    runStage(native, best, start[], nativeOf(o));
    double[3] genes;
    ok(solver2_stage_best(s, genes.ptr, null, null, null));
    assert(sameBits(genes[], best.genes));
    ok(solver2_stage_free(s));

    refused(solver2_stage_new(p, null, null, null, &s), "options is NULL");
    refused(solver2_stage_new(p, &o, null, null, null), "stage is NULL");
    // A population of one has no better half for CMA-ES to move to.
    foreach (lambda; [0, 1])
    {
        o.lambda = lambda;
        refused(solver2_stage_new(p, &o, null, null, &s), "needs a population (lambda) of at least 2");
    }
    o = stageOptions();
    o.sigma = double.nan;
    refused(solver2_stage_new(p, &o, null, null, &s), "needs a step size");
    o = stageOptions();
    o.budget = 0;
    refused(solver2_stage_new(p, &o, null, null, &s), "needs a budget");
    start[1] = double.nan;
    o = stageOptions();
    ok(solver2_stage_new(p, &o, start.ptr, null, &s));
    refused(solver2_stage_run(s, 0, null, null, null), "gene 1 not finite");
    ok(solver2_stage_free(s));
    // A model needs surrogate_pool, of 2 at least (with 1 nothing is
    // pre-screened), and is of its problem.
    auto m = smallModel(p);
    foreach (pool; [0, 1])
    {
        o.surrogate_pool = pool;
        refused(solver2_stage_new(p, &o, null, m, &s), "a stage with a model needs surrogate_pool of at least 2");
    }
    o.surrogate_pool = 2;
    ok(solver2_stage_new(p, &o, null, m, &s));
    ok(solver2_stage_free(s));
    Calls others;
    auto q = field(&others);
    o = stageOptions(true);
    refused(solver2_stage_new(q, &o, null, m, &s), "of another problem");
    // The share the model picks is one: from 0 to 1, both included.
    foreach (share; [-0.1, 1.5, double.nan, double.infinity])
    {
        o.surrogate_share = share;
        refused(solver2_stage_new(p, &o, null, m, &s), "surrogate_share must be from 0 to 1");
        refused(solver2_stage_new(p, &o, null, null, &s), "surrogate_share must be from 0 to 1");
    }
    foreach (share; [0.0, 1.0])
    {
        o.surrogate_share = share;
        ok(solver2_stage_new(p, &o, null, m, &s));
        solver2_stage_result r;
        ok(solver2_stage_run(s, 0, null, null, &r));
        assert(r.finished && r.evaluations == o.budget);
        ok(solver2_stage_free(s));
    }
    ok(solver2_model_free(m));
    ok(solver2_problem_free(q));
    ok(solver2_problem_free(p));
}

unittest
{
    // Pre-screened by a learned model: the native stage with the same
    // model, and the same bytes on any number of the model's threads.
    Calls calls;
    auto p = field(&calls);
    const o = stageOptions(true);
    auto m = smallModel(p);
    solver2_stage* s;
    ok(solver2_stage_new(p, &o, null, m, &s));
    StageTrace trace;
    solver2_stage_result r;
    ok(solver2_stage_run(s, 0, &stageReport, &trace, &r));
    size_t rows;
    ok(solver2_model_rows(m, &rows));
    assert(r.finished && rows == 240 && trace.generations[$ - 1].sampled == 32 && trace.generations[0].sampled == 0);

    DynamicProblem native = new NativeField;
    auto model = new DynamicModel(native, 7, null, EnsembleOptions(2, 16));
    DynamicBest best;
    size_t at;
    runStage(native, best, null, nativeOf(o), model, (const ref DynamicStageGeneration g) {
        const t = &trace.generations[at++];
        assert(t.best is g.best && t.sigma is g.sigma && t.sampled == g.sampled
               && t.predicted_feasible == g.predictedFeasible);
    });
    double[3] genes;
    ok(solver2_stage_best(s, genes.ptr, null, null, null));
    assert(at == 30 && sameBits(genes[], best.genes));
    const one = bytesOf(s);
    ok(solver2_stage_free(s));
    ok(solver2_model_free(m));
    refused(solver2_model_rows(m, &rows), "a model: not a handle");
    ok(solver2_problem_free(p));

    assert(one == wholeStage(true) && one == wholeStage(true, 3));
}

unittest
{
    // Stopped by its report, saved, loaded and gone on with: the stage
    // that never stopped, to the last byte of its state. With a model too.
    foreach (screened; [false, true])
    {
        const whole = wholeStage(screened);
        Calls calls;
        auto p = field(&calls);
        const o = stageOptions(screened);
        auto m = screened ? smallModel(p) : null;
        solver2_stage* s;
        ok(solver2_stage_new(p, &o, null, m, &s));
        StageTrace trace = {stopAt: 10};
        solver2_stage_result r;
        ok(solver2_stage_run(s, 0, &stageReport, &trace, &r));
        assert(r.generations == 10 && r.stopped && !r.finished && trace.generations.length == 10);
        const saved = bytesOf(s);
        ok(solver2_stage_free(s));
        if (m !is null)
            ok(solver2_model_free(m));

        // Into a new stage, and a new model made as the first was.
        auto m2 = screened ? smallModel(p, 2) : null;
        solver2_stage* back;
        ok(solver2_stage_load(p, m2, saved.ptr, saved.length, &back));
        ok(solver2_stage_get(back, &r));
        assert(r.generations == 10 && r.evaluations == 80 && r.stopped && !r.finished);
        assert(bytesOf(back) == saved);
        // Stopped again at once, then to the end.
        trace = StageTrace(null, 11);
        ok(solver2_stage_run(back, 0, &stageReport, &trace, &r));
        assert(r.generations == 11 && r.stopped);
        trace.stopAt = 0;
        ok(solver2_stage_run(back, 0, &stageReport, &trace, &r));
        assert(r.finished && !r.stopped && trace.generations.length == (screened ? 20 : 40));
        assert(bytesOf(back) == whole);
        ok(solver2_stage_free(back));
        if (m2 !is null)
            ok(solver2_model_free(m2));
        ok(solver2_problem_free(p));
    }
}

unittest
{
    // The evaluation fails in the middle of a stage: the call says so, the
    // stage and its model are as after their last whole generation, and the
    // call made again ends as the stage that never failed. Whichever batch
    // fails, the first included, and past where the model began.
    foreach (screened; [false, true])
    {
        const whole = wholeStage(screened);
        foreach (failAt; [1, 2, 4, 9, 17])
        {
            Calls calls = {failAt: failAt};
            auto p = field(&calls);
            const o = stageOptions(screened);
            auto m = screened ? smallModel(p) : null;
            solver2_stage* s;
            ok(solver2_stage_new(p, &o, null, m, &s));
            StageTrace trace;
            solver2_stage_result r;
            refused(solver2_stage_run(s, 0, &stageReport, &trace, &r), "evaluation reported a failure (it returned 7)",
                    SOLVER2_CALLBACK);
            assert(r.generations == failAt - 1 && r.evaluations == 8 * (failAt - 1) && !r.finished && !r.stopped);
            assert(trace.generations.length == failAt - 1);
            if (m !is null)
            {
                size_t rows;
                ok(solver2_model_rows(m, &rows));
                assert(rows == 8 * (failAt - 1));
            }
            ok(solver2_stage_run(s, 0, &stageReport, &trace, &r));
            assert(r.finished && trace.generations.length == o.budget / 8);
            assert(bytesOf(s) == whole);
            ok(solver2_stage_free(s));
            if (m !is null)
                ok(solver2_model_free(m));
            ok(solver2_problem_free(p));
        }
    }

    // A report that fails stops the stage at its generation, whole.
    const whole = wholeStage(false);
    Calls calls;
    auto p = field(&calls);
    const o = stageOptions();
    solver2_stage* s;
    ok(solver2_stage_new(p, &o, null, null, &s));
    StageTrace trace = {failAt: 4};
    solver2_stage_result r;
    refused(solver2_stage_run(s, 0, &stageReport, &trace, &r), "the stage's report reported a failure (it returned -5)",
            SOLVER2_CALLBACK);
    assert(r.generations == 4 && r.stopped && calls.batches == 4);
    // ... the last generation's report too, when nothing is left to stop.
    trace.failAt = 50;
    refused(solver2_stage_run(s, 0, &stageReport, &trace, &r), "it returned -5", SOLVER2_CALLBACK);
    assert(r.generations == 50 && r.finished && bytesOf(s) == whole);
    ok(solver2_stage_free(s));
    ok(solver2_problem_free(p));
}

// ─────────────────────────────────────────────────────────────────────────────
// The archive
// ─────────────────────────────────────────────────────────────────────────────

unittest
{
    // The archive through the interface is the one the searches fill
    // natively: every cell, every threshold, the counts and the picks.
    Calls calls;
    auto p = field(&calls);
    const o = archiveOptions();
    solver2_archive* a;
    ok(solver2_archive_new(p, &o, null, 0, null, &a));

    // Before it has run: empty, its thresholds at the problem's min_score.
    solver2_archive_result r;
    ok(solver2_archive_get(a, &r));
    assert(r.cells == 16 && !r.started && r.coverage == 0 && r.evaluations == 0 && r.best_cell == size_t.max
           && r.best_score != r.best_score && calls.batches == 0);
    ubyte[16] occupied = 1;
    double[16] values, thresholds;
    ok(solver2_archive_elites(a, occupied.ptr, values.ptr, null, null, null, null, null));
    ok(solver2_archive_thresholds(a, thresholds.ptr));
    size_t picked = 9;
    size_t[4] cells;
    ok(solver2_archive_diverse(a, 3, 4, cells.ptr, &picked));
    foreach (c; 0 .. 16)
        assert(!occupied[c] && values[c] != values[c] && thresholds[c] == fieldMinScore);
    assert(picked == 0);

    ArchiveTrace trace;
    ok(solver2_archive_run(a, 0, &archiveReport, &trace, &r));
    assert(r.started && r.finished && !r.stopped && r.iterations == 17 && r.evaluations == 1 + 17 * 24
           && r.seed_evaluations == 1 && r.lambda == 8 && calls.batches == 18 && trace.iterations.length == 17);

    DynamicProblem native = new NativeField;
    size_t at, improved, newCells;
    auto run = runArchive(native, null, nativeOf(o), (const ref DynamicArchiveGeneration g) {
        improved += g.improved;
        newCells += g.newCells;
        if (g.emitter < 2)
            return;
        const t = &trace.iterations[at++];
        assert(t.iteration == g.iteration && t.evaluations == g.evals && t.coverage == g.cells
               && t.qd_score is g.qdScore && t.best_score is g.bestScore && t.improved == improved
               && t.new_cells == newCells);
        improved = 0;
        newCells = 0;
    });
    assert(at == 17 && r.coverage == run.archive.coverage && r.coverage >= 6 && r.qd_score is run.archive.qdScore
           && r.best_cell == run.archive.best && r.best_score is run.archive.bestScore);
    assert(r.restarts_idle == run.restartsIdle && r.restarts_infeasible == run.restartsInfeasible
           && r.restarts_converged == run.restartsConverged && r.reached == run.reached.length
           && r.abandoned == run.abandoned.length);

    size_t[16] found;
    auto genes = new double[16 * 3], violations = new double[16 * 3], descriptors = new double[16 * 2];
    auto outcomes = new solver2_outcome[16];
    ok(solver2_archive_elites(a, occupied.ptr, values.ptr, found.ptr, genes.ptr, outcomes.ptr, violations.ptr,
                              descriptors.ptr));
    ok(solver2_archive_thresholds(a, thresholds.ptr));
    assert(sameBits(thresholds[], run.archive.threshold));
    foreach (c, ref el; run.archive.elites)
    {
        assert(occupied[c] == el.occupied);
        if (!el.occupied)
        {
            assert(values[c] != values[c] && found[c] == 0 && genes[3 * c] != genes[3 * c]
                   && !outcomes[c].scored && descriptors[2 * c] != descriptors[2 * c]);
            continue;
        }
        assert(values[c] is el.value && found[c] == el.found && sameBits(genes[3 * c .. 3 * c + 3], el.genes));
        assert(outcomes[c].scored && outcomes[c].objective is el.eval.objective && outcomes[c].value is el.value
               && outcomes[c].hard == 0 && outcomes[c].soft is el.eval.breach.soft);
        assert(sameBits(violations[3 * c .. 3 * c + 3], el.eval.violations)
               && sameBits(descriptors[2 * c .. 2 * c + 2], el.eval.descriptors));
    }
    size_t[3] over;
    ok(solver2_archive_causes(a, over.ptr));
    assert(over[] == run.causes.over && over[0] > 0);
    const grid = native.grid;
    ok(solver2_archive_diverse(a, 3, 4, cells.ptr, &picked));
    assert(cells[0 .. picked] == selectDiverse(run.archive, grid, 3, 4) && picked > 1 && cells[0] == r.best_cell);
    // At most `count` of them: the best alone, and none into no room.
    size_t[2] room = [77, 77];
    ok(solver2_archive_diverse(a, 3, 1, room.ptr, &picked));
    assert(picked == 1 && room == [r.best_cell, 77]);
    room[0] = 77;
    ok(solver2_archive_diverse(a, 3, 0, room.ptr, &picked));
    assert(picked == 0 && room == [77, 77]);
    picked = 9;
    ok(solver2_archive_diverse(a, 3, 0, null, &picked));
    assert(picked == 0);
    refused(solver2_archive_diverse(a, 3, 1, null, &picked), "cells is NULL");
    char[3] code;
    ok(solver2_problem_cell_code(p, cells[0], code.ptr, code.length));
    assert(code[0 .. 2] == grid.code(cells[0]));

    // A larger budget goes on, as the run given it at the start.
    ok(solver2_archive_run(a, 520, null, null, &r));
    auto more = nativeOf(o);
    more.evals = 520;
    auto longer = runArchive(native, null, more);
    ok(solver2_archive_thresholds(a, thresholds.ptr));
    assert(r.iterations == 22 && r.finished && sameBits(thresholds[], longer.archive.threshold));

    ok(solver2_archive_free(a));
    refused(solver2_archive_get(a, &r), "an archive: not a handle");
    ok(solver2_problem_free(p));
}

unittest
{
    // Seeds of the caller's own, and what an archive refuses.
    Calls calls;
    auto p = field(&calls);
    auto o = archiveOptions();
    solver2_archive* a;
    double[6] seeds = [0.5, 0.3, 0.5, 0.2, 0.2, 0.6];
    ok(solver2_archive_new(p, &o, seeds.ptr, 2, null, &a));
    solver2_archive_result r;
    ok(solver2_archive_run(a, 0, null, null, &r));
    DynamicProblem native = new NativeField;
    auto run = runArchive(native, [seeds[0 .. 3], seeds[3 .. 6]], nativeOf(o));
    double[16] thresholds;
    ok(solver2_archive_thresholds(a, thresholds.ptr));
    assert(r.seed_evaluations == 2 && sameBits(thresholds[], run.archive.threshold));
    ok(solver2_archive_free(a));

    // No seed within the limits: nothing has started, and it can't.
    double[3] far = [1.0, 1.0, 0.5];
    ok(solver2_archive_new(p, &o, far.ptr, 1, null, &a));
    refused(solver2_archive_run(a, 0, null, null, &r), "no seed of the archive is within the limits");
    assert(!r.started && r.evaluations == 0);
    ok(solver2_archive_free(a));

    refused(solver2_archive_new(p, &o, null, 2, null, &a), "seeds is NULL");
    o.emitters = 0;
    refused(solver2_archive_new(p, &o, null, 0, null, &a), "needs emitters");
    o = archiveOptions();
    o.lambda = 1;
    refused(solver2_archive_new(p, &o, null, 0, null, &a), "needs a population (lambda) of at least 2 per emitter");
    foreach (share; [-0.1, 1.5, double.nan])
    {
        o = archiveOptions();
        o.surrogate_share = share;
        refused(solver2_archive_new(p, &o, null, 0, null, &a), "surrogate_share must be from 0 to 1");
    }
    o = archiveOptions();
    o.rules.infeasible_share = 0;
    refused(solver2_archive_new(p, &o, null, 0, null, &a), "infeasibleShare");
    o = archiveOptions(true);
    refused(solver2_archive_new(p, &o, null, 0, null, &a), "needs a model when surrogate_pool is above 0");
    o = archiveOptions();
    auto m = smallModel(p);
    refused(solver2_archive_new(p, &o, null, 0, m, &a), "and only then");
    // With 1 no emitter would be guided: with a model or without.
    o.surrogate_pool = 1;
    refused(solver2_archive_new(p, &o, null, 0, m, &a), "surrogate_pool must be 0, or at least 2 with a model");
    refused(solver2_archive_new(p, &o, null, 0, null, &a), "surrogate_pool must be 0, or at least 2 with a model");
    o.surrogate_pool = 2;
    ok(solver2_archive_new(p, &o, null, 0, m, &a));
    ok(solver2_archive_free(a));
    ok(solver2_model_free(m));
    ok(solver2_problem_free(p));
}

unittest
{
    // Guided by a learned model: the native run with the same model, and
    // the same bytes on any number of the model's threads.
    Calls calls;
    auto p = field(&calls);
    const o = archiveOptions(true);
    auto m = smallModel(p);
    solver2_archive* a;
    ok(solver2_archive_new(p, &o, null, 0, m, &a));
    solver2_archive_result r;
    ok(solver2_archive_run(a, 0, null, null, &r));
    size_t rows;
    ok(solver2_model_rows(m, &rows));
    assert(rows == r.evaluations);

    DynamicProblem native = new NativeField;
    auto model = new DynamicModel(native, 7, null, EnsembleOptions(2, 16));
    auto run = runArchive(native, null, nativeOf(o), null, null, null, null, model);
    double[16] thresholds, values;
    ok(solver2_archive_thresholds(a, thresholds.ptr));
    ok(solver2_archive_elites(a, null, values.ptr, null, null, null, null, null));
    assert(sameBits(thresholds[], run.archive.threshold) && run.guide.sampled > 0 && run.screen.sampled > 0);
    foreach (c, ref el; run.archive.elites)
        assert(el.occupied ? values[c] is el.value : values[c] != values[c]);
    const one = bytesOf(a);
    ok(solver2_archive_free(a));
    ok(solver2_model_free(m));
    ok(solver2_problem_free(p));

    assert(one == wholeArchive(true) && one == wholeArchive(true, 3));
}

unittest
{
    // Stopped by its report, saved, loaded and gone on with: the run that
    // never stopped, to the last byte of its state. With a model too, and
    // saved before it ever ran.
    foreach (guided; [false, true])
    {
        const whole = wholeArchive(guided);
        Calls calls;
        auto p = field(&calls);
        const o = archiveOptions(guided);
        auto m = guided ? smallModel(p) : null;
        solver2_archive* a;
        ok(solver2_archive_new(p, &o, null, 0, m, &a));

        // Not started yet: its bytes hold its options and its seeds.
        const unstarted = bytesOf(a);
        solver2_archive* fresh;
        ok(solver2_archive_load(p, m, unstarted.ptr, unstarted.length, &fresh));
        solver2_archive_result r;
        ok(solver2_archive_get(fresh, &r));
        assert(!r.started && bytesOf(fresh) == unstarted && calls.batches == 0);
        ok(solver2_archive_free(a));
        a = fresh;

        ArchiveTrace trace = {stopAt: 6};
        ok(solver2_archive_run(a, 0, &archiveReport, &trace, &r));
        assert(r.iterations == 6 && r.stopped && !r.finished && trace.iterations.length == 6);
        const saved = bytesOf(a);
        ok(solver2_archive_free(a));
        if (m !is null)
            ok(solver2_model_free(m));

        auto m2 = guided ? smallModel(p, 2) : null;
        solver2_archive* back;
        ok(solver2_archive_load(p, m2, saved.ptr, saved.length, &back));
        ok(solver2_archive_get(back, &r));
        assert(r.started && r.iterations == 6 && r.stopped && !r.finished && bytesOf(back) == saved);
        trace = ArchiveTrace(null, 7);
        ok(solver2_archive_run(back, 0, &archiveReport, &trace, &r));
        assert(r.iterations == 7 && r.stopped);
        trace.stopAt = 0;
        ok(solver2_archive_run(back, 0, &archiveReport, &trace, &r));
        assert(r.finished && !r.stopped && trace.iterations[$ - 1].iteration == 17);
        assert(bytesOf(back) == whole);
        ok(solver2_archive_free(back));
        if (m2 !is null)
            ok(solver2_model_free(m2));
        ok(solver2_problem_free(p));
    }
}

unittest
{
    // The evaluation fails in the middle of a run: the call says so, the
    // run and its model are as after their last whole iteration, and the
    // call made again ends as the run that never failed. Whichever batch
    // fails: the seeds' (1), the first iteration's (2), or a later one's.
    foreach (guided; [false, true])
    {
        const whole = wholeArchive(guided);
        foreach (failAt; [1, 2, 3, 6, 12])
        {
            Calls calls = {failAt: failAt};
            auto p = field(&calls);
            const o = archiveOptions(guided);
            auto m = guided ? smallModel(p) : null;
            solver2_archive* a;
            ok(solver2_archive_new(p, &o, null, 0, m, &a));
            ArchiveTrace trace;
            solver2_archive_result r;
            refused(solver2_archive_run(a, 0, &archiveReport, &trace, &r), "evaluation reported a failure",
                    SOLVER2_CALLBACK);
            // A run whose seeds or first iteration failed hasn't started.
            immutable size_t done = failAt <= 2 ? 0 : failAt - 2;
            assert(r.started == (done > 0) && r.iterations == done && !r.finished && !r.stopped);
            assert(r.evaluations == (done ? 1 + 24 * done : 0));
            if (m !is null)
            {
                size_t rows;
                ok(solver2_model_rows(m, &rows));
                assert(rows == r.evaluations);
            }
            ok(solver2_archive_run(a, 0, &archiveReport, &trace, &r));
            assert(r.finished && r.iterations == 17);
            assert(bytesOf(a) == whole);
            ok(solver2_archive_free(a));
            if (m !is null)
                ok(solver2_model_free(m));
            ok(solver2_problem_free(p));
        }
    }

    // A report that fails stops the run at its iteration, whole.
    const whole = wholeArchive(false);
    Calls calls;
    auto p = field(&calls);
    const o = archiveOptions();
    solver2_archive* a;
    ok(solver2_archive_new(p, &o, null, 0, null, &a));
    ArchiveTrace trace = {failAt: 4};
    solver2_archive_result r;
    refused(solver2_archive_run(a, 0, &archiveReport, &trace, &r), "the archive's report reported a failure",
            SOLVER2_CALLBACK);
    assert(r.iterations == 4 && r.stopped && r.started);
    trace.failAt = 0;
    ok(solver2_archive_run(a, 0, &archiveReport, &trace, &r));
    assert(r.finished && bytesOf(a) == whole);
    ok(solver2_archive_free(a));
    ok(solver2_problem_free(p));
}

/// A report that reads the run it is told of, and saves it.
private struct Watch
{
    solver2_stage* stage;
    solver2_archive* archive;
    ubyte[][] saved;
    bool agreed = true;
}

private extern (C) int stageWatch(void* context, const(solver2_stage_generation)* g)
{
    auto w = cast(Watch*) context;
    solver2_stage_result r;
    ok(solver2_stage_get(w.stage, &r));
    w.agreed &= r.generations == g.generation && r.evaluations == g.evaluations;
    w.saved ~= bytesOf(w.stage);
    return 0;
}

private extern (C) int archiveWatch(void* context, const(solver2_archive_iteration)* it)
{
    auto w = cast(Watch*) context;
    solver2_archive_result r;
    ok(solver2_archive_get(w.archive, &r));
    ubyte[16] occupied;
    double[16] thresholds;
    size_t[3] over;
    ok(solver2_archive_elites(w.archive, occupied.ptr, null, null, null, null, null, null));
    ok(solver2_archive_thresholds(w.archive, thresholds.ptr));
    ok(solver2_archive_causes(w.archive, over.ptr));
    size_t held, raised;
    foreach (c; 0 .. 16)
    {
        held += occupied[c];
        raised += thresholds[c] > fieldMinScore;
    }
    w.agreed &= r.started && r.iterations == it.iteration && r.evaluations == it.evaluations
                && r.coverage == it.coverage && r.qd_score is it.qd_score && r.best_score is it.best_score
                && held == it.coverage && raised >= held && over[0] > 0;
    w.saved ~= bytesOf(w.archive);
    return 0;
}

unittest
{
    // A report finds the run as after the generation or iteration it is
    // told of, to read and to save: in a run's first call as in a later
    // one. What it saves, loaded and gone on with, ends as the run that
    // never stopped.
    foreach (withModel; [false, true])
    {
        {
            const whole = wholeStage(withModel);
            Calls calls;
            auto p = field(&calls);
            const o = stageOptions(withModel);
            auto m = withModel ? smallModel(p) : null;
            Watch watch;
            ok(solver2_stage_new(p, &o, null, m, &watch.stage));
            ok(solver2_stage_run(watch.stage, 0, &stageWatch, &watch, null));
            assert(watch.agreed && watch.saved.length == o.budget / 8 && watch.saved[$ - 1] == whole);
            foreach (at; [0, 9])
            {
                solver2_stage* back;
                solver2_stage_result r;
                ok(solver2_stage_load(p, m, watch.saved[at].ptr, watch.saved[at].length, &back));
                ok(solver2_stage_get(back, &r));
                assert(r.generations == at + 1 && !r.finished);
                ok(solver2_stage_run(back, 0, null, null, &r));
                assert(r.finished && bytesOf(back) == whole);
                ok(solver2_stage_free(back));
            }
            ok(solver2_stage_free(watch.stage));
            if (m !is null)
                ok(solver2_model_free(m));
            ok(solver2_problem_free(p));
        }
        {
            const whole = wholeArchive(withModel);
            Calls calls;
            auto p = field(&calls);
            const o = archiveOptions(withModel);
            auto m = withModel ? smallModel(p) : null;
            Watch watch;
            ok(solver2_archive_new(p, &o, null, 0, m, &watch.archive));
            // Stopped and gone on with: the first call's reports, then a later call's.
            ArchiveTrace trace = {stopAt: 6};
            ok(solver2_archive_run(watch.archive, 0, &archiveWatch, &watch, null));
            assert(watch.agreed && watch.saved.length == 17 && watch.saved[$ - 1] == whole);
            foreach (at; [0, 5])
            {
                solver2_archive* back;
                solver2_archive_result r;
                ok(solver2_archive_load(p, m, watch.saved[at].ptr, watch.saved[at].length, &back));
                ok(solver2_archive_get(back, &r));
                assert(r.started && r.iterations == at + 1 && !r.finished && r.evaluations == 1 + 24 * (at + 1));
                ok(solver2_archive_run(back, 0, null, null, &r));
                assert(r.finished && bytesOf(back) == whole);
                ok(solver2_archive_free(back));
            }
            ok(solver2_archive_free(watch.archive));

            // A run stopped in its first call reports from its second as from its first.
            Watch later;
            solver2_archive_result r;
            if (m !is null)
            {
                ok(solver2_model_free(m));
                m = smallModel(p);
            }
            ok(solver2_archive_new(p, &o, null, 0, m, &later.archive));
            ok(solver2_archive_run(later.archive, 0, &archiveReport, &trace, &r));
            assert(r.stopped && r.iterations == 6);
            ok(solver2_archive_run(later.archive, 0, &archiveWatch, &later, &r));
            assert(r.finished && later.agreed && later.saved.length == 11 && later.saved == watch.saved[6 .. $]);
            ok(solver2_archive_free(later.archive));
            if (m !is null)
                ok(solver2_model_free(m));
            ok(solver2_problem_free(p));
        }
    }
}

unittest
{
    // A model freed while a run holds it: its threads stop, and the run
    // goes on with the model, on the caller's thread, to the same bytes.
    foreach (threads; [0, 2])
    {
        Calls calls;
        auto p = field(&calls);
        {
            const o = stageOptions(true);
            auto m = smallModel(p, threads);
            solver2_stage* s;
            ok(solver2_stage_new(p, &o, null, m, &s));
            StageTrace trace = {stopAt: 10};
            solver2_stage_result r;
            ok(solver2_stage_run(s, 0, &stageReport, &trace, &r));
            assert(r.stopped && r.generations == 10);
            ok(solver2_model_free(m));
            size_t rows;
            refused(solver2_model_rows(m, &rows), "a model: not a handle");
            ok(solver2_stage_run(s, 0, &stageReport, &trace, &r));
            assert(r.finished && bytesOf(s) == wholeStage(true));
            ok(solver2_stage_free(s));
        }
        {
            const o = archiveOptions(true);
            auto m = smallModel(p, threads);
            solver2_archive* a;
            ok(solver2_archive_new(p, &o, null, 0, m, &a));
            ArchiveTrace trace = {stopAt: 6};
            solver2_archive_result r;
            ok(solver2_archive_run(a, 0, &archiveReport, &trace, &r));
            assert(r.stopped && r.iterations == 6);
            ok(solver2_model_free(m));
            ok(solver2_archive_run(a, 0, &archiveReport, &trace, &r));
            assert(r.finished && bytesOf(a) == wholeArchive(true));
            ok(solver2_archive_free(a));
        }
        ok(solver2_problem_free(p));
    }
}

unittest
{
    // A larger budget whose call fails is the run's all the same: the run
    // is not finished, and the call made again with no budget goes on to
    // it, and ends as the run given that budget at the start.
    {
        Calls calls;
        auto p = field(&calls);
        auto o = stageOptions();
        solver2_stage* s, longer;
        ok(solver2_stage_new(p, &o, null, null, &s));
        solver2_stage_result r;
        ok(solver2_stage_run(s, 0, null, null, &r));
        assert(r.finished && r.evaluations == 400);
        calls.failAt = calls.batches + 5;
        refused(solver2_stage_run(s, 480, null, null, &r), "evaluation reported a failure", SOLVER2_CALLBACK);
        assert(r.evaluations == 432 && !r.finished && !r.stopped);
        ok(solver2_stage_get(s, &r));
        assert(r.evaluations == 432 && !r.finished);
        ok(solver2_stage_run(s, 0, null, null, &r));
        assert(r.evaluations == 480 && r.finished);
        o.budget = 480;
        ok(solver2_stage_new(p, &o, null, null, &longer));
        ok(solver2_stage_run(longer, 0, null, null, null));
        assert(bytesOf(s) == bytesOf(longer));
        ok(solver2_stage_free(longer));
        ok(solver2_stage_free(s));
        ok(solver2_problem_free(p));
    }
    {
        Calls calls;
        auto p = field(&calls);
        auto o = archiveOptions();
        solver2_archive* a, longer;
        ok(solver2_archive_new(p, &o, null, 0, null, &a));
        solver2_archive_result r;
        ok(solver2_archive_run(a, 0, null, null, &r));
        assert(r.finished && r.iterations == 17);
        calls.failAt = calls.batches + 3;
        refused(solver2_archive_run(a, 520, null, null, &r), "evaluation reported a failure", SOLVER2_CALLBACK);
        assert(r.iterations == 19 && !r.finished && !r.stopped);
        ok(solver2_archive_run(a, 0, null, null, &r));
        assert(r.iterations == 22 && r.finished);
        o.evals = 520;
        ok(solver2_archive_new(p, &o, null, 0, null, &longer));
        ok(solver2_archive_run(longer, 0, null, null, null));
        assert(bytesOf(a) == bytesOf(longer));
        ok(solver2_archive_free(longer));
        ok(solver2_archive_free(a));
        ok(solver2_problem_free(p));
    }
}

unittest
{
    // What a load refuses: bytes that aren't a run's, a run of another
    // problem, and a model that isn't the saved run's.
    Calls calls;
    auto p = field(&calls);
    const so = stageOptions(true);
    const ao = archiveOptions(true);
    auto m = smallModel(p);
    solver2_stage* s;
    solver2_archive* a;
    ok(solver2_stage_new(p, &so, null, m, &s));
    ok(solver2_archive_new(p, &ao, null, 0, m, &a));
    StageTrace st = {stopAt: 5};
    ArchiveTrace at = {stopAt: 5};
    ok(solver2_stage_run(s, 0, &stageReport, &st, null));
    const stage = bytesOf(s);
    ok(solver2_archive_run(a, 0, &archiveReport, &at, null));
    const archive = bytesOf(a);
    size_t size;
    ubyte[8] small;
    refused(solver2_archive_save(a, small.ptr, small.length, &size), "the buffer holds 8 bytes");
    assert(size == archive.length);
    refused(solver2_archive_save(a, null, 0, null), "size is NULL");
    ok(solver2_stage_free(s));
    ok(solver2_archive_free(a));

    refused(solver2_archive_load(p, m, stage.ptr, stage.length, &a), "not a solver2-c archive file");
    refused(solver2_stage_load(p, m, archive.ptr, archive.length, &s), "not a solver2-c stage file");
    refused(solver2_archive_load(p, m, archive.ptr, archive.length / 2, &a), "cut short");
    refused(solver2_archive_load(p, m, null, 0, &a), "bytes is NULL");
    ubyte[64] noise = 0xA7;
    refused(solver2_stage_load(p, m, noise.ptr, noise.length, &s), "");

    // A stage whose best has another number of violations than the problem
    // has limits: solver2_stage_best would copy them past the caller's
    // buffer, or leave it short. Written again as it was, the stage loads.
    {
        import solver2.bytes : Bytes;
        static import solver2c.checkpoint;
        import solver2c.runs : FORMAT, SavedStage, STAGE_MAGIC;
        auto r = Bytes(stage.dup);
        solver2c.checkpoint.readHeader(r, STAGE_MAGIC, FORMAT);
        SavedStage saved;
        solver2c.checkpoint.read(r, saved);
        assert(saved.best.feasible && saved.best.eval.violations.length == 3);
        ubyte[] written(const(double)[] violations)
        {
            auto altered = saved;
            altered.best.eval.violations = violations.dup;
            Bytes o;
            solver2c.checkpoint.writeHeader(o, STAGE_MAGIC, FORMAT);
            solver2c.checkpoint.write(o, altered);
            return o.data;
        }
        const whole = written(saved.best.eval.violations);
        assert(whole == stage);
        foreach (violations; [saved.best.eval.violations ~ 1.0, saved.best.eval.violations[0 .. 2], null])
        {
            const bytes = written(violations);
            refused(solver2_stage_load(p, m, bytes.ptr, bytes.length, &s), "the saved stage's best doesn't fit the problem");
        }
    }

    // The model: the saved run's, made as it was.
    refused(solver2_archive_load(p, null, archive.ptr, archive.length, &a), "the saved run had a model");
    refused(solver2_stage_load(p, null, stage.ptr, stage.length, &s), "the saved run had a model");
    auto other = smallModel(p, 0, 8);
    refused(solver2_archive_load(p, other, archive.ptr, archive.length, &a), "another seed or other options");
    ok(solver2_model_free(other));
    const plain = wholeArchive(false);
    refused(solver2_archive_load(p, m, plain.ptr, plain.length, &a), "the saved run had no model");

    // Another problem: by the caller's fingerprint, and by its description.
    Calls others;
    auto revised = field(&others, "field 2");
    auto m2 = smallModel(revised);
    refused(solver2_archive_load(revised, m2, archive.ptr, archive.length, &a),
            "of another problem: problem: saved field 1, now field 2");
    refused(solver2_stage_load(revised, m2, stage.ptr, stage.length, &s), "of another problem");
    ok(solver2_model_free(m2));
    ok(solver2_problem_free(revised));
    auto spec = fieldSpec(&others);
    spec.min_score = 4;
    solver2_problem* lower;
    ok(solver2_problem_new(&spec, &lower));
    refused(solver2_archive_load(lower, null, plain.ptr, plain.length, &a), "solver2 min score");
    ok(solver2_problem_free(lower));
    auto genes = spec.genes[0 .. 3].dup;
    genes[2].hi = 2;
    spec = fieldSpec(&others);
    spec.genes = genes.ptr;
    solver2_problem* wider;
    ok(solver2_problem_new(&spec, &wider));
    refused(solver2_archive_load(wider, null, plain.ptr, plain.length, &a), "solver2 genes");
    ok(solver2_problem_free(wider));

    // After all that the model is as it was, and the saved runs still load.
    ok(solver2_archive_load(p, m, archive.ptr, archive.length, &a));
    assert(bytesOf(a) == archive);
    ok(solver2_archive_free(a));
    ok(solver2_stage_load(p, m, stage.ptr, stage.length, &s));
    assert(bytesOf(s) == stage);
    ok(solver2_stage_free(s));
    ok(solver2_model_free(m));
    ok(solver2_problem_free(p));
    assert(others.batches == 0);
}
