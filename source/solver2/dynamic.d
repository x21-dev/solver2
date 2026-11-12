/// Problems as classes, through the searches: that a `DynamicProblem`
/// (problem.d) used only through the base class gets what a compile-time
/// problem gets from the stage, the archive and the learned model, and the
/// same bits.
///
/// So far this module holds those checks. The entry points that need no
/// template argument at the call site come here next.
module solver2.dynamic;

import solver2.model : isLearning;
import solver2.problem;
import solver2.stage : hasCauses;

// The optional members the searches look for on the static type are all
// there on the base class.
static assert(hasCauses!DynamicProblem && isLearning!DynamicProblem);
static assert(__traits(hasMember, DynamicProblem, "screenPenalty") && __traits(hasMember, DynamicProblem, "noFeasibleSeed")
              && __traits(hasMember, DynamicProblem, "adoptCauses"));

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

version (unittest)
{
    import std.exception : assertThrown, collectExceptionMsg;
    import std.string : indexOf;

    import solver2.archive;
    import solver2.ensemble : EnsembleOptions;
    import solver2.model : LearnedModel;
    import solver2.stage;

    private static immutable GeneSpec[3] fieldGenes = [GeneSpec("x", 0, 2, 1), GeneSpec("y", 0, 2, 0.6),
                                                       GeneSpec("z", -1, 1, 0)];
    private static immutable ConstraintSpec[3] fieldLimits = [ConstraintSpec("sum", Tier.hard),
                                                              ConstraintSpec("floor", Tier.terminal),
                                                              ConstraintSpec("x_cap", Tier.soft, 0.5)];
    private static immutable Axis[2] fieldAxes = [Axis("x", [0.5, 1, 1.5], "abcd", 0, 0.5),
                                                  Axis("y", [0.5, 1, 1.5], "abcd", 0, 0.5)];

    /// The test problem's evaluation: maximise 10 - (x - 1)² - (y - 0.8)² -
    /// z², with x + y ≤ 1.9 hard, z ≥ -0.9 terminal and x ≤ 1.5 soft. The
    /// descriptors are x and y.
    private void fillField(const(double[])[] genes, Scoring scoring, Outcome[] o)
    {
        foreach (k, g; genes)
        {
            immutable double x = fieldGenes[0].toPhysical(g[0]), y = fieldGenes[1].toPhysical(g[1]),
                             z = fieldGenes[2].toPhysical(g[2]);
            o[k].violations = [x + y > 1.9 ? x + y - 1.9 : 0, z < -0.9 ? -0.9 - z : 0, x > 1.5 ? x - 1.5 : 0];
            o[k].breach = breachOf(fieldLimits[], o[k].violations);
            o[k].scored = scoring == Scoring.all || o[k].breach.hard == 0;
            if (o[k].scored)
                o[k].objective = 10 - (x - 1) * (x - 1) - (y - 0.8) * (y - 0.8) - z * z;
            o[k].descriptors = [x, y];
        }
    }

    /// The problem as a class, with only what a subclass must give.
    private class Field : DynamicProblem
    {
        size_t batches, penalties;

        override const(GeneSpec)[] genes() const { return fieldGenes[]; }
        override const(ConstraintSpec)[] constraints() const { return fieldLimits[]; }
        override const(string[2])[] fingerprint() const { return [["problem", "field 1"]]; }
        override const(Axis)[] axes() const { return fieldAxes[]; }

        override void evaluateInto(const(double[])[] g, Scoring scoring, Outcome[] o)
        {
            ++batches;
            fillField(g, scoring, o);
        }
    }

    /// ... and with what a subclass may add.
    private final class Costly : Field
    {
        override double screenPenalty(const(double)[] x)
        {
            ++penalties;
            return 3 * clamp01(x[2]);
        }

        override string noFeasibleSeed() const { return "start nearer the middle"; }
        override double predictionHeadroom() const { return 0.05; }
    }

    /// The same problem as a compile-time one, with none of the optional
    /// members: what the class's defaults must behave as.
    private struct FieldS
    {
        alias Eval = Outcome;

        size_t dimension() const { return 3; }

        Outcome[] evaluate(const(double[])[] g, Scoring scoring)
        {
            auto o = new Outcome[g.length];
            fillField(g, scoring, o);
            return o;
        }

        bool   scored(const ref Outcome e) const { return e.scored; }
        double objective(const ref Outcome e) const { return e.objective; }
        Breach breach(const ref Outcome e) const { return e.breach; }
        size_t descriptorCount() const { return 2; }
        void   descriptors(const ref Outcome e, double[] d) const { d[] = e.descriptors[]; }
    }

    static assert(hasDescriptors!FieldS && !hasCauses!FieldS && !__traits(hasMember, FieldS, "screenPenalty"));

    /// A model that knows the problem exactly.
    private final class Exact : StageModel!Outcome
    {
        void learn(const(double[])[] genes, const(Outcome)[] evals) {}
        void train(size_t epochs) {}

        void predictScore(const(double[])[] genes, double[] score, double[] violation)
        {
            auto o = new Outcome[genes.length];
            fillField(genes, Scoring.all, o);
            foreach (k, ref e; o)
            {
                score[k] = value(e.objective, e.breach);
                violation[k] = e.breach.hard + e.breach.terminal;
            }
        }
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

    /// A stage over `p`, its every generation's numbers, genes and picks
    /// flat, to compare stages bit for bit.
    private double[] stageTrace(P, M)(ref P p, M model, const ref StageOptions o, out Best!Outcome best,
                                      out size_t overLimits)
    {
        double[] t;
        size_t over;
        cmaStage(p, best, [0.5, 0.3, 0.5], o, model, (const ref StageGeneration!P g) {
            t ~= [cast(double) g.generation, g.best, g.populationBest, g.sigma, cast(double) g.infeasible];
            foreach (row; g.genes)
                t ~= row;
            foreach (i; g.picked)
                t ~= i;
            static if (hasCauses!P)
            {
                // Each candidate over a hard or terminal limit is counted
                // under at least one of them.
                assert(g.causes.over.length == 3 && g.causes.over[0] + g.causes.over[1] >= g.infeasible);
                assert(g.causes.over[0] <= o.lambda && g.causes.over[1] <= o.lambda);
                over += g.causes.over[0];
            }
        });
        overLimits = over;
        return t;
    }

    private bool sameArchive(A)(const ref A a, const ref A b)
    {
        if (!sameBits(a.threshold[], b.threshold[]) || a.elites.length != b.elites.length)
            return false;
        foreach (c, ref el; a.elites)
        {
            const other = &b.elites[c];
            if (el.occupied != other.occupied || !sameBits(el.genes, other.genes) || el.value !is other.value
                || el.found != other.found)
                return false;
        }
        return true;
    }

    private ArchiveOptions fieldOptions()
    {
        ArchiveOptions o;
        o.evals = 400;
        o.emitters = 3;
        o.explorers = 1;
        o.lambda = 8;
        o.seed = 5;
        return o;
    }

    private alias DynamicArchiving = Archiving!(DynamicProblem, BinnedGrid);
}

unittest
{
    // What a caller needs comes from the object: its start genes, its grid
    // and its ranking.
    DynamicProblem p = new Field;
    assert(p.dimension == 3 && p.descriptorCount == 2 && p.startGenes == [0.5, 0.3, 0.5]);
    assert(p.grid.cells == 16 && p.ranking == Ranking.init);
    assert(p.screenPenalty([0.5, 0.5, 0.5]) == 0 && p.noFeasibleSeed is null);
    assert(p.featureCount == 3 && p.targetCount == 5);

    // The default learning spec: the value, the hard and the terminal
    // breach, the descriptors; a prediction read back by its tiers.
    auto e = p.evaluate([[0.9, 0.5, 0.5], [0.5, 0.3, 0.0]], Scoring.all);
    auto t = new double[5];
    p.targets(e[0], t);                    // x 1.8, y 1: over the sum by 0.9, x over its cap by 0.3
    assert(t[0] is value(e[0].objective, e[0].breach) && t[0] == e[0].objective - 0.5 * e[0].violations[2]);
    assert(t[1] == e[0].breach.hard && t[1] > 0 && t[2] == 0 && t[3] == 1.8 && t[4] == 1);
    p.targets(e[1], t);                    // z -1: under the floor
    assert(t[1] == 0 && t[2] == e[1].breach.terminal && t[2] > 0);
    const q = p.predicted([7, 0.01, 0.2, 1.25, 0.5]);
    assert(q.objective == 7 && q.breach == Breach(0.01, 0.2) && q.descriptors == [1.25, 0.5]);
    DynamicProblem c = new Costly;
    assert(c.predicted([7, 0.01, 0.2, 1.25, 0.5]).breach == Breach(0, 0.2));     // its headroom
}

unittest
{
    // checkProblem: a sound problem passes; one that isn't is told why.
    checkProblem(new Field);
    checkProblem(new Costly, 3, 9);

    // The evaluation depends on what else was evaluated.
    static class Drifting : Field
    {
        override void evaluateInto(const(double[])[] g, Scoring scoring, Outcome[] o)
        {
            super.evaluateInto(g, scoring, o);
            foreach (ref e; o)
                e.objective += 1e-9 * batches;
        }
    }
    assert(collectExceptionMsg(checkProblem(new Drifting)).indexOf("other bits when evaluated again") > 0);

    static class Batchy : Field
    {
        override void evaluateInto(const(double[])[] g, Scoring scoring, Outcome[] o)
        {
            super.evaluateInto(g, scoring, o);
            foreach (ref e; o)
                e.objective += 1e-9 * g.length;
        }
    }
    assert(collectExceptionMsg(checkProblem(new Batchy)).indexOf("alone than in a batch") > 0);

    // Scoring.feasibleOnly may skip an objective, not change one.
    static class Moody : Field
    {
        override void evaluateInto(const(double[])[] g, Scoring scoring, Outcome[] o)
        {
            super.evaluateInto(g, scoring, o);
            if (scoring == Scoring.feasibleOnly)
                foreach (ref e; o)
                    e.objective += 1;
        }
    }
    assert(collectExceptionMsg(checkProblem(new Moody)).indexOf("Scoring.feasibleOnly") > 0);

    // A description that doesn't fit.
    static class Lettered : Field
    {
        override const(Axis)[] axes() const { return [Axis("x", [0.5, 1, 1.5], "abc", 0, 0.5), fieldAxes[1]]; }
    }
    assert(collectExceptionMsg(checkProblem(new Lettered)).indexOf("3 letters for 4 bins") > 0);

    // What evaluate checks of every batch, with or without checkProblem.
    static class Short : Field
    {
        override void evaluateInto(const(double[])[] g, Scoring scoring, Outcome[] o)
        {
            super.evaluateInto(g, scoring, o);
            foreach (ref e; o)
                e.descriptors = e.descriptors[0 .. 1];
        }
    }
    DynamicProblem s = new Short;
    assert(collectExceptionMsg(s.evaluate([[0.5, 0.3, 0.5]], Scoring.all)).indexOf("1 descriptors for 2") > 0);
    assertThrown(checkProblem(s));

    static class Negative : Field
    {
        override void evaluateInto(const(double[])[] g, Scoring scoring, Outcome[] o)
        {
            super.evaluateInto(g, scoring, o);
            o[0].breach.soft = -1;
        }
    }
    DynamicProblem n = new Negative;
    assert(collectExceptionMsg(n.evaluate([[0.5, 0.3, 0.5]], Scoring.all)).indexOf("negative or NaN") > 0);

    static class Unscored : Field
    {
        override void evaluateInto(const(double[])[] g, Scoring scoring, Outcome[] o)
        {
            super.evaluateInto(g, scoring, o);
            o[0].scored = false;
        }
    }
    DynamicProblem u = new Unscored;
    assert(collectExceptionMsg(u.evaluate([[0.5, 0.3, 0.5]], Scoring.all)).indexOf("without an objective") > 0);
}

unittest
{
    // A stage over a subclass through the base class: the bits of the same
    // problem known at compile time, and its cause counts on top.
    StageOptions o = {lambda: 8, sigma: 0.2, budget: 400, seed: 3, stallGenerations: 0};
    DynamicProblem p = new Field;
    FieldS s;
    Best!Outcome bp, bs;
    size_t over, none;
    const tp = stageTrace(p, null, o, bp, over);
    const ts = stageTrace(s, null, o, bs, none);
    assert(sameBits(tp, ts) && sameBits(bp.genes, bs.genes) && bp.evals == 400);
    assert(bp.feasible && bp.eval.objective > 9.99 && over > 0);
    assert((cast(Field) p).batches == 50);
}

unittest
{
    // Pre-screened: the default screenPenalty behaves as none (the bits of
    // the compile-time problem without one), and a subclass's is called
    // through the base class and changes what is picked.
    StageOptions o = {lambda: 8, sigma: 0.2, budget: 240, seed: 5, stallGenerations: 0,
                      surrogatePool: 4, surrogateShare: 0.5, surrogateFrom: 3, surrogateEpochs: 2};
    DynamicProblem p = new Field, c = new Costly;
    FieldS s;
    Best!Outcome bp, bs, bc;
    size_t over;
    const tp = stageTrace(p, new Exact, o, bp, over);
    const ts = stageTrace(s, new Exact, o, bs, over);
    const tc = stageTrace(c, new Exact, o, bc, over);
    assert(sameBits(tp, ts) && sameBits(bp.genes, bs.genes));
    assert((cast(Field) p).penalties == 0 && (cast(Field) c).penalties == 28 * 32);
    assert(!sameBits(tc, tp));
}

unittest
{
    // The archive through the base class: the compile-time problem's bits,
    // the cause counts, and the subclass's hint when no seed is feasible.
    DynamicProblem p = new Field;
    FieldS s;
    const grid = p.grid;
    const seeds = [p.startGenes];
    auto opts = fieldOptions();
    opts.ranking = p.ranking;
    auto rp = DynamicArchiving.runArchive(p, grid, seeds, opts);
    auto rs = Archiving!(FieldS, BinnedGrid).runArchive(s, grid, seeds, opts);
    assert(sameArchive(rp.archive, rs.archive) && rp.evals == rs.evals && rp.reached == rs.reached);
    assert(rp.archive.coverage >= 6 && rp.causes.over.length == 3 && rp.causes.over[0] > 0);

    DynamicProblem c = new Costly;
    const double[][] far = [[1.0, 1.0, 0.5]];             // x 2, y 2: over the sum
    assert(collectExceptionMsg(DynamicArchiving.runArchive(c, grid, far, opts)).indexOf("start nearer the middle") > 0);
    assert(collectExceptionMsg(DynamicArchiving.runArchive(p, grid, far, opts)).indexOf(":") < 0);
}

unittest
{
    // Stopped, saved as plain data and resumed through the base class: the
    // run that never stopped, its cause counts included.
    DynamicProblem p = new Field;
    const grid = p.grid;
    const seeds = [p.startGenes];
    const opts = fieldOptions();
    auto whole = DynamicArchiving.runArchive(p, grid, seeds, opts);
    auto stopped = DynamicArchiving.runArchive(p, grid, seeds, opts, null, null,
                                              (const ref DynamicArchiving.Run r) => r.iterations == 6);
    assert(stopped.stopped && stopped.iterations == 6);

    const saved = runState(stopped);
    DynamicArchiving.Run back;
    restoreRun(back, saved);
    back.archive = stopped.archive;
    back.archive.elites = stopped.archive.elites.dup;
    back.archive.threshold = stopped.archive.threshold.dup;
    back.causes = CauseCounts(stopped.causes.over.dup);
    auto resumed = DynamicArchiving.runArchive(p, grid, null, opts, null, &back);
    assert(sameArchive(resumed.archive, whole.archive) && resumed.evals == whole.evals);
    assert(resumed.causes.over == whole.causes.over && resumed.reached == whole.reached);

    // Counts that weren't restored are sized, and count from the resume on.
    DynamicArchiving.Run bare;
    restoreRun(bare, saved);
    bare.archive = stopped.archive;
    bare.archive.elites = stopped.archive.elites.dup;
    bare.archive.threshold = stopped.archive.threshold.dup;
    auto again = DynamicArchiving.runArchive(p, grid, null, opts, null, &bare);
    assert(sameArchive(again.archive, whole.archive) && again.causes.over.length == 3);
    assert(again.causes.over[0] == whole.causes.over[0] - stopped.causes.over[0]);
}

unittest
{
    // The learned model over the base class, on the defaults of its
    // learning spec: a pre-screened stage and model-guided explorers, the
    // same twice.
    alias Model = LearnedModel!DynamicProblem;
    StageOptions so = {lambda: 8, sigma: 0.2, budget: 160, seed: 5, stallGenerations: 0,
                       surrogatePool: 4, surrogateShare: 0.5, surrogateFrom: 3, surrogateEpochs: 2};
    double[][2] traces;
    foreach (ref t; traces)
    {
        DynamicProblem p = new Costly;
        auto model = new Model(p, 7, null, EnsembleOptions(2, 16));
        Best!Outcome best;
        size_t over;
        t = stageTrace(p, model, so, best, over);
        assert(model.rows == 160 && best.feasible && (cast(Field) p).penalties == 18 * 32);
    }
    assert(sameBits(traces[0], traces[1]));

    auto ao = fieldOptions();
    ao.surrogatePool = 4;
    ao.surrogateFrom = 3;
    ao.surrogateEpochs = 2;
    ao.surrogateImprovers = true;
    DynamicArchiving.Run[2] runs;
    foreach (ref r; runs)
    {
        DynamicProblem p = new Costly;
        const grid = p.grid;
        auto model = new Model(p, 7, null, EnsembleOptions(2, 16));
        r = DynamicArchiving.runArchive(p, grid, [p.startGenes], ao, null, null, null, null, model);
        assert(model.rows == r.evals && r.guide.sampled > 0 && r.screen.sampled > 0);
        assert((cast(Field) p).penalties > 0);
    }
    assert(sameArchive(runs[0].archive, runs[1].archive));
}

unittest
{
    // The manual's "What a caller needs comes from the object": a problem
    // the caller knows only as a DynamicProblem, checked, then through a
    // stage and the archive on what the object says of itself.
    DynamicProblem p = new Costly;
    checkProblem(p);

    Best!Outcome best;
    StageOptions so = {lambda: 8, sigma: 0.2, budget: 800, seed: 42, ranking: p.ranking};
    cmaStage(p, best, p.startGenes, so, null);

    alias A = Archiving!(DynamicProblem, BinnedGrid);
    const grid = p.grid;
    ArchiveOptions ao = {evals: 400, seed: 9, ranking: p.ranking};
    auto run = A.runArchive(p, grid, [best.genes], ao);
    assert(best.feasible && run.archive.bestScore >= valueOf(p, best.eval) && run.archive.coverage >= 6);
}
