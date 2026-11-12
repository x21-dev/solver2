/// The worked example on the class route: docs/example's problem as a
/// `DynamicProblem`, run with no template argument at the call site
/// (solver2.dynamic), from what the object says of itself. The same flow:
/// a CMA-ES stage pre-screened by a learned model, then the archive with a
/// model, saved to a file part way, resumed from it, and checked against
/// the run that never stopped.
///
/// It is also the reference of the Python example
/// (bindings/python/examples/beam.py), which is this program through the
/// C interface, its evaluation in Python: the two print the same lines, to
/// the last bit of the numbers on the lines that start with "bits".
///
/// What differs from docs/example, which is a compile-time problem with a
/// learning spec of its own: the learning spec here is the class's default
/// (the value, the hard and the terminal breach, the descriptors), which
/// is all a problem behind the C interface can have, and `minScore` is
/// below the problem's lowest value, so `checkProblem` passes. The two
/// examples' runs therefore differ.
module app;

import std.exception : enforce;
import std.file : read, remove, tempDir, write;
import std.parallelism : TaskPool;
import std.path : buildPath;
import std.stdio : writefln;

static import checkpoint;
import solver2.archive : ArchiveOptions, restoreRun, runState, RunState, selectDiverse;
import solver2.bytes : Bytes;
import solver2.dynamic;
import solver2.problem;
import solver2.stage : StageOptions;

/// Two physical parameters, a width in 1 .. 10 and a log-scaled rate in
/// 0.1 .. 100. Maximise a smooth objective; width × rate ≤ 200 is hard,
/// width ≥ 2 terminal, rate ≤ 40 soft at a cost of 0.05 per unit over.
final class Beam : DynamicProblem
{
    static immutable GeneSpec[2] spec = [GeneSpec("width", 1, 10, 5),
                                         GeneSpec("rate", 0.1, 100, 1, Scale.log)];
    static immutable ConstraintSpec[3] limits = [ConstraintSpec("load", Tier.hard),
                                                 ConstraintSpec("min_width", Tier.terminal),
                                                 ConstraintSpec("rate_cap", Tier.soft, 0.05)];
    static immutable Axis[2] grid2 = [Axis("width", [3, 5, 7], "nmwW", 0, 2),
                                      Axis("rate", [1, 10, 40], "slfF", 0, 10)];

    override const(GeneSpec)[] genes() const { return spec[]; }
    override const(ConstraintSpec)[] constraints() const { return limits[]; }
    override const(string[2])[] fingerprint() const { return [["problem", "beam 1"], ["load_limit", "200"]]; }
    override const(Axis)[] axes() const { return grid2[]; }

    /// The lowest value a beam within its limits can have is -48.25, at
    /// width 2 and rate 100: the archive's thresholds start below it.
    override double minScore() const { return -50; }

    /// 0.5 of headroom on a predicted breach.
    override double predictionHeadroom() const { return 0.5; }

    override void evaluateInto(const(double[])[] g, Scoring scoring, Outcome[] o)
    {
        foreach (k, row; g)
        {
            // toPhysical clamps: the searches pass genes as sampled.
            immutable double w = spec[0].toPhysical(row[0]), r = spec[1].toPhysical(row[1]);
            o[k].violations = [w * r > 200 ? w * r - 200 : 0, w < 2 ? 2 - w : 0, r > 40 ? r - 40 : 0];
            o[k].breach = breachOf(limits[], o[k].violations);
            o[k].scored = scoring == Scoring.all || o[k].breach.hard == 0;
            if (o[k].scored)
                o[k].objective = 10 - (w - 7) * (w - 7) / 4 - (r - 30) * (r - 30) / 100;
            o[k].descriptors = [w, r];
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Saving
// ─────────────────────────────────────────────────────────────────────────────

enum string MAGIC = "solver2-example-class checkpoint";
enum uint FORMAT = 1;

/// What a checkpoint holds: as docs/example's, for the class's types.
struct Saved
{
    string[2][] fingerprint;
    RunState!ArchiveOptions state;
    typeof(DynamicRun.archive) archive;
    CauseCounts causes;
    ubyte[] model;
}

void save(string path, DynamicProblem p, const ref DynamicRun run, DynamicModel model)
{
    const Saved s = {fingerprint: p.fingerprint, state: runState(run), archive: run.archive, causes: run.causes,
                     model: model.save(run.evals)};
    Bytes o;
    checkpoint.writeHeader(o, MAGIC, FORMAT);
    checkpoint.write(o, s);
    write(path, o.data);
}

DynamicRun load(string path, DynamicProblem p, DynamicModel model)
{
    auto r = Bytes(cast(ubyte[]) read(path));
    checkpoint.readHeader(r, MAGIC, FORMAT);
    Saved s;
    checkpoint.read(r, s);
    enforce(r.done, "the checkpoint has bytes left over");
    if (const why = fingerprintMismatch(s.fingerprint, p.fingerprint))
        throw new Exception("the checkpoint is of another problem: " ~ why);
    DynamicRun run;
    restoreRun(run, s.state);
    run.archive = s.archive;
    run.causes = s.causes;
    model.restore(s.model, run.evals);
    return run;
}

// ─────────────────────────────────────────────────────────────────────────────
// Bits
// ─────────────────────────────────────────────────────────────────────────────

ulong bitsOf(double v) pure nothrow @nogc @trusted
{
    return *cast(const(ulong)*) &v;
}

/// FNV-1a over the bytes of `values`, least significant first.
ulong hashOf(const(double)[] values) pure nothrow @nogc @safe
{
    ulong h = 0xCBF2_9CE4_8422_2325UL;
    foreach (v; values)
    {
        immutable ulong w = bitsOf(v);
        foreach (i; 0 .. 8)
        {
            h ^= (w >> (8 * i)) & 0xFF;
            h *= 0x0000_0100_0000_01B3UL;
        }
    }
    return h;
}

// ─────────────────────────────────────────────────────────────────────────────
// The run
// ─────────────────────────────────────────────────────────────────────────────

void main()
{
    // The caller knows the problem only as a DynamicProblem.
    DynamicProblem p = new Beam;
    checkProblem(p);
    writefln("checked: %d genes, %d cells", p.dimension, p.grid.cells);

    // The models train and predict on four threads of their own. Their
    // number doesn't change the bits: a model is the same on any.
    auto pool = new TaskPool(4);
    scope (exit)
        pool.finish(true);

    // 1. A CMA-ES stage from the problem's start, pre-screened by a model.
    DynamicBest best;
    StageOptions so = {lambda: 8, sigma: 0.2, budget: 800, seed: 42, surrogatePool: 4};
    auto over = p.newCauses();
    const end = runStage(p, best, null, so, new DynamicModel(p, 3, pool), (const ref DynamicStageGeneration g) {
        over += g.causes;
        if (g.generation % 25 == 0)
            writefln("stage gen %3d  best %.4f  infeasible %d", g.generation, g.best, g.infeasible);
    });
    writefln("stage: %d generations, %d evaluations, stalled %s; best %.4f at width %.3f rate %.3f",
             end.generations, end.spent, end.stalled, valueOf(p, best.eval), best.eval.descriptors[0],
             best.eval.descriptors[1]);
    writefln("stage: candidates over load %d, min_width %d, rate_cap %d", over.over[0], over.over[1], over.over[2]);
    writefln("bits: stage best %016x genes %016x", bitsOf(valueOf(p, best.eval)), hashOf(best.genes));

    // 2. The archive over the problem's grid, seeded with the stage's best:
    //    once straight through, once stopped after 10 iterations, saved,
    //    loaded and resumed.
    ArchiveOptions ao = {evals: 4000, emitters: 4, explorers: 2, lambda: 8, seed: 9,
                         surrogatePool: 4, surrogateFrom: 3, surrogateEpochs: 5, surrogateAim: true};
    const seeds = [best.genes];
    const whole = runArchive(p, seeds, ao, null, null, null, null, new DynamicModel(p, 5, pool));

    immutable path = buildPath(tempDir, "solver2-example-class.state");
    {
        auto model = new DynamicModel(p, 5, pool);
        const stopped = runArchive(p, seeds, ao, null, null, (const ref DynamicRun r) => r.iterations == 10, null,
                                   model);
        save(path, p, stopped, model);
        writefln("stopped after %d iterations: %d cells, QD score %.2f", stopped.iterations,
                 stopped.archive.coverage, stopped.archive.qdScore);
    }
    auto model = new DynamicModel(p, 5, pool);
    auto run = load(path, p, model);
    remove(path);
    const resumed = runArchive(p, null, ao, null, &run, null, null, model);

    const grid = p.grid;
    writefln("resumed to %d evaluations: %d of %d cells, best %.4f, %d explorer targets reached",
             resumed.evals, resumed.archive.coverage, grid.cells, resumed.archive.bestScore, resumed.reached.length);
    double[] values, genes;
    foreach (c, ref el; whole.archive.elites)
    {
        enforce(el.occupied == resumed.archive.elites[c].occupied && el.genes == resumed.archive.elites[c].genes
                && el.value is resumed.archive.elites[c].value, "the resumed run differs");
        if (el.occupied)
        {
            values ~= el.value;
            genes ~= el.genes;
        }
    }
    enforce(whole.archive.threshold == resumed.archive.threshold && whole.evals == resumed.evals,
            "the resumed run differs");
    writefln("the resumed run's archive is the uninterrupted run's, bit for bit");
    foreach (c; selectDiverse(resumed.archive, grid, 3, 4))
    {
        const el = &resumed.archive.elites[c];
        writefln("  %s  value %.4f  width %.3f  rate %.3f", grid.code(c), el.value, el.eval.descriptors[0],
                 el.eval.descriptors[1]);
    }
    writefln("bits: thresholds %016x values %016x genes %016x QD score %016x", hashOf(resumed.archive.threshold),
             hashOf(values), hashOf(genes), bitsOf(resumed.archive.qdScore));
}
