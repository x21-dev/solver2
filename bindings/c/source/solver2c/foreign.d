/// A problem given as data across the C interface: its genes, limits, axes,
/// ranking and lowest value copied once, and its evaluation the caller's
/// function, called once per batch with flat arrays.
///
/// It is a `DynamicProblem` like any other. The learning spec is the class's
/// own, computed on the D side: the model's worker threads never call the
/// caller's function.
module solver2c.foreign;

import std.conv : text;
import std.exception : enforce;
import std.format : format;
import std.string : fromStringz;

import solver2.problem;
import solver2c.c;

/// The caller's own function reported a failure: its evaluation by a
/// return value other than 0, a report by a negative one. The exported
/// function that was running returns SOLVER2_CALLBACK.
final class CallbackFailed : Exception
{
    int code;

    this(string what, int code, string file = __FILE__, size_t line = __LINE__) pure nothrow @safe
    {
        super(what, file, line);
        this.code = code;
    }
}

/// A C string as a D one, copied; `otherwise` when it is NULL.
string textOf(const(char)* s, lazy string otherwise = null)
{
    return s is null ? otherwise : fromStringz(s).idup;
}

final class ForeignProblem : DynamicProblem
{
    /// Batches and candidates sent to the caller's evaluation so far.
    size_t batches, evaluations;

    /// Copies `s` but for its context, which is passed back as it is.
    /// Throws when a count has no array, or a scale or tier isn't one; the
    /// rest of the description is `checkDescription`'s to check.
    this(const ref solver2_problem_spec s)
    {
        enforce(s.evaluate !is null, "a problem needs its evaluate function");
        enforce(s.gene_count > 0 && s.genes !is null, "a problem needs genes");
        enforce(s.limit_count == 0 || s.limits !is null, text(s.limit_count, " limits, but the array is NULL"));
        enforce(s.axis_count == 0 || s.axes !is null, text(s.axis_count, " axes, but the array is NULL"));
        enforce(s.fingerprint_count == 0 || s.fingerprint !is null,
                text(s.fingerprint_count, " fingerprint pairs, but the array is NULL"));
        fn = s.evaluate;
        context = cast(void*) s.context;
        foreach (i, ref g; s.genes[0 .. s.gene_count])
        {
            immutable string name = textOf(g.name, text("gene ", i));
            enforce(g.scale >= SOLVER2_LINEAR && g.scale <= SOLVER2_INTEGER,
                    text("gene ", name, ": its scale ", g.scale, " is none of SOLVER2_LINEAR, _LOG and _INTEGER"));
            genes_ ~= GeneSpec(name, g.lo, g.hi, g.neutral, cast(Scale) g.scale);
        }
        foreach (i, ref c; s.limits[0 .. s.limit_count])
        {
            immutable string name = textOf(c.name, text("limit ", i));
            enforce(c.tier >= SOLVER2_HARD && c.tier <= SOLVER2_SOFT,
                    text("limit ", name, ": its tier ", c.tier, " is none of SOLVER2_HARD, _TERMINAL and _SOFT"));
            limits_ ~= ConstraintSpec(name, cast(Tier) c.tier, c.weight, textOf(c.group));
        }
        foreach (i, ref a; s.axes[0 .. s.axis_count])
        {
            immutable string name = textOf(a.name, text("axis ", i));
            enforce(a.edge_count == 0 || a.edges !is null, text("axis ", name, ": its edges are NULL"));
            axes_ ~= Axis(name, a.edges[0 .. a.edge_count].dup, textOf(a.letters), a.nan_bin, a.width);
        }
        ranking_ = Ranking(s.ranking.hard_offset, s.ranking.terminal_offset, s.ranking.bound_penalty);
        minScore_ = s.min_score;
        headroom_ = s.prediction_headroom;

        // What a saved run must match: the caller's pairs, and the
        // description, since a run saved on other genes, limits or axes
        // can't go on here.
        foreach (i; 0 .. s.fingerprint_count)
        {
            const name = s.fingerprint[2 * i], value = s.fingerprint[2 * i + 1];
            enforce(name !is null && value !is null, text("fingerprint pair ", i, " has a NULL in it"));
            fingerprint_ ~= [textOf(name), textOf(value)];
        }
        string genesText, limitsText, axesText;
        foreach (ref g; genes_)
            genesText ~= format("%s %a %a %a %d; ", g.name, g.lo, g.hi, g.neutral, cast(int) g.scale);
        foreach (ref c; limits_)
            limitsText ~= format("%s %d %a; ", c.name, cast(int) c.tier, c.weight);
        foreach (ref a; axes_)
            axesText ~= format("%s %(%a %) %s %d %a; ", a.name, a.edges, a.letters, a.nanBin, a.width);
        fingerprint_ ~= ["solver2 genes", genesText];
        fingerprint_ ~= ["solver2 limits", limitsText];
        fingerprint_ ~= ["solver2 axes", axesText];
        fingerprint_ ~= ["solver2 ranking", format("%a %a %a", ranking_.hardOffset, ranking_.terminalOffset,
                                                   ranking_.boundPenalty)];
        fingerprint_ ~= ["solver2 min score", format("%a", minScore_)];
    }

    override const(GeneSpec)[] genes() const { return genes_; }
    override const(ConstraintSpec)[] constraints() const { return limits_; }
    override const(string[2])[] fingerprint() const { return fingerprint_; }
    override const(Axis)[] axes() const { return axes_; }
    override Ranking ranking() const { return ranking_; }
    override double minScore() const { return minScore_; }
    override double predictionHeadroom() const { return headroom_; }

    /// One call of the caller's function for the batch. Each candidate's
    /// breach is summed here, by `breachOf`, from the violations it gave:
    /// so the order of the sum and its bits are the library's.
    override void evaluateInto(const(double[])[] g, Scoring scoring, Outcome[] o)
    {
        immutable size_t rows = g.length, n = genes_.length, limits = limits_.length, nd = axes_.length;
        auto flat = new double[rows * n];
        foreach (k, row; g)
            flat[k * n .. (k + 1) * n] = row[];
        auto scored = new ubyte[rows];
        auto objective = new double[rows], violations = new double[rows * limits],
             descriptors = new double[rows * nd];
        violations[] = 0;
        ++batches;
        evaluations += rows;
        immutable int status = fn(context, rows, flat.ptr, scoring == Scoring.feasibleOnly, scored.ptr, objective.ptr,
                                  limits ? violations.ptr : null, nd ? descriptors.ptr : null);
        if (status != 0)
            throw new CallbackFailed(text("the problem's evaluation reported a failure (it returned ", status, ")"),
                                     status);
        foreach (k, ref e; o)
        {
            e.violations = violations[k * limits .. (k + 1) * limits].dup;
            foreach (i, v; e.violations)
                enforce(v >= 0, text("candidate ", k, " of ", rows, ": its violation of the limit ", limits_[i].name,
                                     " is negative or NaN (", v, ")"));
            e.breach = breachOf(limits_, e.violations);
            e.scored = scored[k] != 0;
            e.objective = e.scored ? objective[k] : double.nan;
            e.descriptors = descriptors[k * nd .. (k + 1) * nd].dup;
        }
    }

private:
    GeneSpec[] genes_;
    ConstraintSpec[] limits_;
    Axis[] axes_;
    Ranking ranking_;
    double minScore_ = 0, headroom_ = 0;
    string[2][] fingerprint_;
    solver2_evaluate_fn fn;
    void* context;
}
