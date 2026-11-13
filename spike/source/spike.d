/// The spike that came before the C interface: a shared library that runs
/// solver2's CMA-ES stage over a problem whose evaluation is a C callback,
/// so that Python can supply it. Throwaway code: it answers whether the D
/// runtime, its collector and its threads behave inside a Python process,
/// and whether a run through a Python callback has the bits of the same
/// run through a D one.
module spike;

import core.memory : GC;
import core.runtime : rt_init, rt_term;
import core.thread : Thread, thread_attachThis, thread_detachThis;
import std.parallelism : TaskPool;

import solver2.model : LearnedModel;
import solver2.problem;
import solver2.stage;

/// The caller's evaluation of a batch: `n` rows of `dim` genes in `genes`,
/// as sampled (the callback clamps them); it fills `objective` and `over`
/// (how far over the hard limit, 0 within). Returns 0, or anything else
/// to stop the run.
alias EvalFn = extern (C) int function(void* ctx, size_t n, size_t dim, const(double)* genes, double* objective,
                                       double* over);

/// What a stage gives back.
extern (C) struct SpikeResult
{
    ulong hash;               // of every generation's numbers and genes, in order
    ulong generations, spent;
    double best;              // the best value
    ulong learned;            // rows the model learned, 0 without one
}

/// The problem: its evaluation is the callback.
struct Callback
{
    struct Eval
    {
        bool   scored;
        double obj = double.nan;
        double over;
    }

    EvalFn f;
    void*  ctx;
    size_t n;

    size_t dimension() const { return n; }

    Eval[] evaluate(const(double[])[] genes, Scoring)
    {
        auto flat = new double[genes.length * n];
        foreach (k, g; genes)
            flat[k * n .. (k + 1) * n] = g[];
        auto objective = new double[genes.length], over = new double[genes.length];
        if (f(ctx, genes.length, n, flat.ptr, objective.ptr, over.ptr) != 0)
            throw new Exception("the evaluation callback failed");
        auto o = new Eval[genes.length];
        foreach (k, ref e; o)
            e = Eval(true, objective[k], over[k]);
        return o;
    }

    bool   scored(const ref Eval e) const { return e.scored; }
    double objective(const ref Eval e) const { return e.obj; }
    Breach breach(const ref Eval e) const { return Breach(e.over); }

    // The learning spec, for a pre-screened stage: computed here, never
    // through the callback (the model calls these from its own threads).
    size_t featureCount() const { return n; }
    void features(const(double)[] genes, double[] x) const { x[] = genes[]; }
    size_t targetCount() const { return 2; }
    void targets(const ref Eval e, double[] t) const { t[0] = e.obj; t[1] = e.over; }
    Predicted predicted(const(double)[] m) const { return Predicted(m[0], Breach(m[1] > 0.05 ? m[1] : 0), null); }
}

static assert(isProblem!Callback && hasSurrogate!Callback);

private string lastError;

/// Runs `f` and keeps anything it throws on this side of the boundary.
private int guarded(scope void delegate() f) nothrow
{
    try
    {
        f();
        return 0;
    }
    catch (Throwable t)
    {
        try
            lastError = t.msg.idup ~ "\0";
        catch (Throwable)
            lastError = "error\0";
        return 1;
    }
}

private ulong bitsOf(double v) pure nothrow @nogc @trusted
{
    return *cast(const(ulong)*) &v;
}

/// FNV-1a, a word at a time.
private void fold(ref ulong h, ulong w) pure nothrow @nogc @safe
{
    foreach (i; 0 .. 8)
    {
        h ^= (w >> (8 * i)) & 0xFF;
        h *= 0x0000_0100_0000_01B3UL;
    }
}

extern (C):

/// Starts the D runtime. Before any other call, from the thread that makes
/// the others.
export int spike_init() nothrow
{
    try
        return rt_init() ? 0 : 1;
    catch (Throwable)
        return 1;
}

export int spike_term() nothrow
{
    try
        return rt_term() ? 0 : 1;
    catch (Throwable)
        return 1;
}

/// The last error's message.
export const(char)* spike_error() nothrow
{
    return lastError.length ? lastError.ptr : "".ptr;
}

/// Makes the calling thread known to the D runtime (and its collector), and
/// forgets it.
export int spike_attach() nothrow
{
    return guarded({ thread_attachThis(); });
}

export int spike_detach() nothrow
{
    return guarded({ thread_detachThis(); });
}

/// Collections the D collector has run in this process.
export ulong spike_collections() nothrow
{
    ulong n;
    guarded({ n = GC.profileStats().numCollections; });
    return n;
}

/// Throws inside D, to show nothing gets out.
export int spike_throw() nothrow
{
    return guarded({ throw new Exception("thrown inside D"); });
}

/// Whether the D runtime knows the calling thread, and how many it knows.
export int spike_known() nothrow
{
    int known = -1;
    guarded({ known = Thread.getThis() !is null; });
    return known;
}

export int spike_threads() nothrow
{
    int n = -1;
    guarded({ n = cast(int) Thread.getAll().length; });
    return n;
}

/// The stronger probe: megabytes in use after a collection with 64 MB held
/// only by a local (`held`), and after one with nothing holding it
/// (`dropped`), and whether the held bytes were intact after 64 MB more
/// had been allocated over them.
export int spike_probe2(double* held, double* dropped) nothrow
{
    int intact = -1;
    guarded({
        enum size_t n = 64 << 20;
        auto a = new ubyte[n];
        a[] = 0xAB;
        GC.collect();
        *held = GC.stats().usedSize / 1048576.0;
        auto b = new ubyte[n];
        b[] = 0xCD;
        intact = 1;
        foreach (i; 0 .. n / 4096)
            intact &= a[i * 4096] == 0xAB;
        a = null;
        b = null;
        GC.collect();
        *dropped = GC.stats().usedSize / 1048576.0;
    });
    return intact;
}

/// Whether the collector sees the calling thread's stack: an array held
/// only by a local, a collection, and is the array still allocated? 1 yes,
/// 0 it was freed under us, -1 an error.
export int spike_probe() nothrow
{
    int alive = -1;
    guarded({
        auto a = new ubyte[1 << 20];
        a[] = 0xAB;
        GC.collect();
        alive = GC.addrOf(a.ptr) !is null && a[$ - 1] == 0xAB;
    });
    return alive;
}

/// A bug inside D, to see what it does to the process: 0 an index out of
/// range, 1 a failed assert, 2 a call through a null class reference.
export int spike_fault(int kind, size_t index) nothrow
{
    return guarded({
        auto a = new int[3];
        Object o;
        final switch (kind)
        {
            case 0: a[index] = 1; break;
            case 1: assert(index < a.length, "an assert failed inside D"); break;
            case 2: lastError = o.toString() ~ "\0"; break;
        }
    });
}

/// A CMA-ES stage over the callback's problem in `dim` genes, for `budget`
/// evaluations from `seed`. `withModel`: pre-screened by a learned model,
/// which trains and predicts on `threads` worker threads (0: none).
export int spike_stage(EvalFn f, void* ctx, size_t dim, ulong seed, size_t budget, int withModel, size_t threads,
                       SpikeResult* result, double* bestGenes) nothrow
{
    return guarded({
        auto p = Callback(f, ctx, dim);
        Best!(Callback.Eval) best;
        auto x0 = new double[dim];
        foreach (j, ref v; x0)
            v = j % 2 ? 0.2 : 0.8;
        StageOptions o = {lambda: 8, sigma: 0.2, budget: budget, seed: seed, stallGenerations: 0};
        TaskPool pool = threads ? new TaskPool(threads) : null;
        scope (exit)
            if (pool !is null)
                pool.finish(true);
        LearnedModel!Callback learned;
        StageModel!(Callback.Eval) model;
        if (withModel)
        {
            o.surrogatePool = 4;
            o.surrogateFrom = 3;
            o.surrogateEpochs = 2;
            learned = new LearnedModel!Callback(p, seed, pool);
            model = learned;
        }
        ulong h = 0xCBF2_9CE4_8422_2325UL;
        const end = cmaStage(p, best, x0, o, model, (const ref StageGeneration!Callback g) {
            fold(h, g.generation);
            fold(h, bitsOf(g.best));
            fold(h, bitsOf(g.populationBest));
            fold(h, bitsOf(g.sigma));
            fold(h, g.infeasible);
            foreach (row; g.genes)
                foreach (v; row)
                    fold(h, bitsOf(v));
        });
        result.hash = h;
        result.generations = end.generations;
        result.spent = end.spent;
        result.best = best.feasible ? valueOf(p, best.eval) : double.nan;
        result.learned = learned is null ? 0 : learned.rows;
        foreach (j; 0 .. dim)
            bestGenes[j] = best.feasible ? best.genes[j] : double.nan;
    });
}
