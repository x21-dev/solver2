/// The spike's runs with the evaluation in D: what spike.py must print,
/// line for line, with the evaluation in Python.
module native;

import std.stdio : writefln;

import solver2.problem : clamp01;
import spike;

/// Maximise -(x - 0.3)² summed, with x[0] ≤ 0.5 hard. In plain IEEE
/// operations in a fixed order, as spike.py's `sphere` has them.
extern (C) int sphere(void* ctx, size_t n, size_t dim, const(double)* genes, double* objective, double* over)
{
    foreach (k; 0 .. n)
    {
        double s = 0;
        foreach (j; 0 .. dim)
        {
            immutable double d = clamp01(genes[k * dim + j]) - 0.3;
            s -= d * d;
        }
        objective[k] = s;
        immutable double x0 = clamp01(genes[k * dim]);
        over[k] = x0 > 0.5 ? x0 - 0.5 : 0;
    }
    return 0;
}

void report(string name, ulong seed, const ref SpikeResult r, const(double)[] genes)
{
    writefln("%s seed %d: hash %016x generations %d spent %d learned %d best %016x genes %(%016x %)", name, seed,
             r.hash, r.generations, r.spent, r.learned, *cast(const(ulong)*) &r.best, cast(const(ulong)[]) genes);
}

int main()
{
    enum size_t dim = 6;
    if (spike_init() != 0)
        return 1;
    SpikeResult r;
    double[dim] genes;
    // Plain stages, enough of them for the collector to run many times.
    foreach (seed; 1 .. 51)
    {
        if (spike_stage(&sphere, null, dim, seed, 4000, 0, 0, &r, genes.ptr) != 0)
            return 1;
        report("plain", seed, r, genes[]);
    }
    // Pre-screened stages: a learned model on its own threads.
    foreach (seed; 1 .. 6)
    {
        if (spike_stage(&sphere, null, dim, seed, 800, 1, 3, &r, genes.ptr) != 0)
            return 1;
        report("model", seed, r, genes[]);
    }
    return 0;
}
