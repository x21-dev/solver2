/// The exported functions of include/solver2.h. Each one is a thin shell:
/// it makes the calling thread known to the D runtime, looks its handles
/// up, calls foreign.d or runs.d, copies the results into the caller's
/// buffers, and turns anything thrown into a return code and a message.
///
/// What the shell sees to:
///
/// - **Nothing thrown crosses the boundary**, a bug's `Error` included
///   (with bounds checks on, an index out of range is one).
/// - **Any thread may call.** A thread the D runtime doesn't know is
///   attached to it on its first call, stays attached, and is detached
///   when it ends, by a hook the operating system calls then (see
///   "Threads" below).
/// - **Handles are numbers**, never addresses: a slot of a table and the
///   count of times the slot was given out, so a handle that was freed, or
///   is of another kind, or never was one, is found out and refused. The
///   table is what keeps the objects alive: the caller holds nothing of
///   the collector's.
module solver2c.api;

import core.runtime : rt_init, rt_term;
import core.sync.mutex : Mutex;
import core.thread : Thread, thread_attachThis, thread_detachThis, thread_setThis;

import std.conv : text;
import std.exception : enforce;

import solver2.archive : ArchiveOptions, ArchiveRules;
import solver2.ensemble : EnsembleOptions;
import solver2.problem;
import solver2.stage : StageOptions;

import solver2c.c;
import solver2c.foreign;
import solver2c.runs;

static assert(size_t.sizeof == 8, "the handles are made for 64-bit targets");

// ─────────────────────────────────────────────────────────────────────────────
// Errors
// ─────────────────────────────────────────────────────────────────────────────

/// The calling thread's last error, as C text. A thread's own, and in no
/// memory of the collector's: a thread detached on its way out may still
/// ask for it.
private char[2048] errorText = '\0';

private void setError(scope const(char)[][] parts...) nothrow @nogc
{
    size_t at = 0;
    foreach (part; parts)
        foreach (ch; part)
            if (at + 1 < errorText.length)
                errorText[at++] = ch;
    errorText[at] = '\0';
}

// ─────────────────────────────────────────────────────────────────────────────
// Threads
// ─────────────────────────────────────────────────────────────────────────────
//
// The collector must know every thread that runs D code: it stops them all
// and reads their stacks. A thread made elsewhere (by C, by Python) is
// unknown to it until it is attached, and must be detached before it ends.
//
// A thread is attached once, on its first call, and stays so. Attaching for
// the length of each call and detaching before it returns, which reads as
// the tidier rule, is not safe with this runtime (druntime 2.112, in LDC
// 1.42): detaching leaves the thread as the runtime's "this thread", so
// its next call finds itself known and runs unseen by the collector; and
// with that cleared, a second attach gives the thread a new Thread object
// while its block cache (core.internal.gc.blkcache) stays recorded on the
// first, so the collector no longer clears the cache's entries of blocks it
// frees. The test driver's second thread showed both: other bits in most
// runs, then a corrupted heap in every run.
//
// So that no thread ends while attached, each attached thread leaves a mark
// in a slot of thread-local storage that has a destructor: the operating
// system calls `threadEnds` on the thread as it ends (fiber-local storage
// on Windows, a pthread key elsewhere), and that detaches it.

private extern (C) void rt_moduleTlsCtor();
private extern (C) void rt_moduleTlsDtor();

version (Windows)
{
    private extern (Windows) nothrow @nogc
    {
        uint FlsAlloc(void function(void*) nothrow callback);
        int FlsSetValue(uint index, void* value);
        int FlsFree(uint index);
    }

    private enum uint noKey = uint.max;         // FLS_OUT_OF_INDEXES
    private __gshared uint exitKey = noKey;

    private extern (Windows) void threadEnds(void* mark) nothrow
    {
        leave();
    }

    /// Makes the slot. With the first solver2_init.
    private void makeExitKey()
    {
        exitKey = FlsAlloc(&threadEnds);
        enforce(exitKey != noKey, "no thread-local slot left for the threads' exit hook");
    }

    /// Frees the slot, and with it the hooks. Windows calls the hook of
    /// every thread's mark from here, on this thread: `ending` tells
    /// `leave` to do nothing.
    private void freeExitKey() nothrow @nogc
    {
        if (exitKey != noKey)
            FlsFree(exitKey);
        exitKey = noKey;
    }

    private bool markThread() nothrow @nogc
    {
        return FlsSetValue(exitKey, cast(void*) 1) != 0;
    }
}
else
{
    import core.sys.posix.pthread : pthread_key_create, pthread_key_delete, pthread_key_t, pthread_setspecific;

    private __gshared pthread_key_t exitKey;
    private __gshared bool haveExitKey;

    private extern (C) void threadEnds(void* mark) nothrow
    {
        leave();
    }

    private void makeExitKey()
    {
        enforce(pthread_key_create(&exitKey, &threadEnds) == 0, "no thread-local key left for the threads' exit hook");
        haveExitKey = true;
    }

    private void freeExitKey() nothrow @nogc
    {
        if (haveExitKey)
            pthread_key_delete(exitKey);
        haveExitKey = false;
    }

    private bool markThread() nothrow @nogc
    {
        return pthread_setspecific(exitKey, cast(void*) 1) == 0;
    }
}

/// Set while solver2_term takes the hooks down.
private __gshared bool ending;

/// Makes the calling thread known to the D runtime, when it isn't: for
/// good, until `leave`.
private void enter()
{
    if (Thread.getThis() !is null)
        return;
    thread_attachThis();
    if (!markThread())
    {
        thread_detachThis();
        thread_setThis(null);
        throw new Exception("this thread can't be given its exit hook, so it can't call solver2");
    }
    rt_moduleTlsCtor();         // as a thread the runtime starts gets them
}

/// An attached thread ends: its thread-local destructors, then off the
/// runtime's list, and no longer its "this thread".
private void leave() nothrow
{
    if (ending || inits == 0)
        return;
    try
    {
        if (Thread.getThis() is null)
            return;
        rt_moduleTlsDtor();
        thread_detachThis();
        thread_setThis(null);
    }
    catch (Throwable)
    {
    }
}

/// Runs `f` with the calling thread known to the D runtime, and keeps
/// whatever it throws on this side.
private int guarded(scope void delegate() f) nothrow
{
    if (inits == 0)
    {
        setError("solver2_init has not been called");
        return SOLVER2_ERROR;
    }
    try
    {
        enter();
        f();
        return SOLVER2_OK;
    }
    catch (CallbackFailed e)
    {
        setError(e.msg);
        return SOLVER2_CALLBACK;
    }
    catch (Exception e)
    {
        setError(e.msg);
        return SOLVER2_ERROR;
    }
    catch (Throwable t)
    {
        // An Error: a bug on this side. Say where.
        char[20] digits = void;
        size_t n = digits.length;
        size_t line = t.line;
        do
        {
            digits[--n] = cast(char)('0' + line % 10);
            line /= 10;
        }
        while (line > 0);
        setError("a bug in solver2: ", t.msg, " (", t.file, ":", digits[n .. $], ")");
        return SOLVER2_ERROR;
    }
}

/// Throws unless the caller gave a pointer where one is needed.
private void need(bool given, lazy string what)
{
    enforce(given, what ~ " is NULL");
}

/// Copies `from` into the caller's buffer, when it gave one.
private void give(T)(T* to, const(T)[] from) nothrow @nogc
{
    if (to !is null)
        to[0 .. from.length] = from[];
}

// ─────────────────────────────────────────────────────────────────────────────
// Handles
// ─────────────────────────────────────────────────────────────────────────────

private __gshared int inits;                // solver2_init calls not yet ended
private __gshared Mutex tableLock;
private __gshared Object[] table;           // a slot's object; null when free
private __gshared uint[] issued;            // how many times each slot was given out

/// A handle for `o`: its slot and the slot's count.
private void* handleOf(Object o)
{
    tableLock.lock();
    scope (exit)
        tableLock.unlock();
    size_t slot = table.length;
    foreach (i, held; table)
        if (held is null)
        {
            slot = i;
            break;
        }
    if (slot == table.length)
    {
        table ~= null;
        issued ~= 0;
    }
    table[slot] = o;
    ++issued[slot];
    return cast(void*)((cast(size_t) issued[slot] << 32) | (slot + 1));
}

/// The object `handle` stands for, or a throw that says what was expected.
private T objectOf(T)(const(void)* handle, string what, bool free = false)
{
    immutable size_t v = cast(size_t) handle, slot = (v & 0xFFFF_FFFF) - 1, count = v >> 32;
    tableLock.lock();
    scope (exit)
        tableLock.unlock();
    enforce(v != 0, what ~ ": the handle is NULL");
    enforce(slot < table.length && issued[slot] == count && table[slot] !is null,
            what ~ ": not a handle (it was freed, or never was one)");
    auto o = cast(T) table[slot];
    enforce(o !is null, what ~ ": a handle of another kind");
    if (free)
        table[slot] = null;
    return o;
}

private ForeignProblem problemOf(const(solver2_problem)* h, bool free = false)
{
    return objectOf!ForeignProblem(h, "a problem", free);
}

private ForeignModel modelOf(const(solver2_model)* h, bool free = false)
{
    return objectOf!ForeignModel(h, "a model", free);
}

/// The model of `h`, or none for NULL.
private ForeignModel modelOrNone(const(solver2_model)* h)
{
    return h is null ? null : modelOf(h);
}

private ForeignStage stageOf(const(solver2_stage)* h, bool free = false)
{
    return objectOf!ForeignStage(h, "a stage", free);
}

private ForeignArchive archiveOf(const(solver2_archive)* h, bool free = false)
{
    return objectOf!ForeignArchive(h, "an archive", free);
}

/// An outcome as the caller reads it.
private solver2_outcome outcomeOf(const ref Outcome e) pure nothrow @nogc
{
    return solver2_outcome(e.scored, e.objective, e.scored ? value(e.objective, e.breach) : double.nan, e.breach.hard,
                           e.breach.terminal, e.breach.soft);
}

/// Bytes into the caller's buffer, the way the save functions give them.
private void giveBytes(const(ubyte)[] bytes, ubyte* buffer, size_t capacity, size_t* size)
{
    need(size !is null, "size");
    *size = bytes.length;
    if (buffer is null)
        return;
    enforce(capacity >= bytes.length, text("the buffer holds ", capacity, " bytes, the saved state has ", bytes.length));
    buffer[0 .. bytes.length] = bytes[];
}

extern (C):

// ─────────────────────────────────────────────────────────────────────────────
// The library
// ─────────────────────────────────────────────────────────────────────────────

export int solver2_init() nothrow
{
    try
    {
        if (!rt_init())
        {
            setError("the D runtime didn't start");
            return SOLVER2_ERROR;
        }
        if (inits == 0)
        {
            tableLock = new Mutex;
            ending = false;
            makeExitKey();
        }
        ++inits;
        return SOLVER2_OK;
    }
    catch (Throwable t)
    {
        setError("the D runtime didn't start: ", t.msg);
        return SOLVER2_ERROR;
    }
}

export int solver2_term() nothrow
{
    if (inits == 0)
    {
        setError("solver2_term without a solver2_init before it");
        return SOLVER2_ERROR;
    }
    try
    {
        if (inits == 1)
        {
            // The last one: stop the models' threads, and let go of everything.
            foreach (o; table)
                if (auto m = cast(ForeignModel) o)
                    m.close();
            table = null;
            issued = null;
            // No thread's end calls back into a runtime that is gone.
            ending = true;
            freeExitKey();
        }
        --inits;
        if (!rt_term())
        {
            setError("the D runtime didn't end");
            return SOLVER2_ERROR;
        }
        return SOLVER2_OK;
    }
    catch (Throwable t)
    {
        setError("the D runtime didn't end: ", t.msg);
        return SOLVER2_ERROR;
    }
}

export const(char)* solver2_last_error() nothrow @nogc
{
    return errorText.ptr;
}

export int solver2_abi() nothrow @nogc
{
    return SOLVER2_ABI;
}

export int solver2_checked() nothrow @nogc
{
    version (Solver2Checked)
        return 1;
    else
        return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// A problem
// ─────────────────────────────────────────────────────────────────────────────

export int solver2_problem_spec_default(solver2_problem_spec* spec) nothrow @nogc
{
    if (spec is null)
    {
        setError("spec is NULL");
        return SOLVER2_ERROR;
    }
    immutable Ranking r;
    *spec = solver2_problem_spec.init;
    spec.ranking = solver2_ranking(r.hardOffset, r.terminalOffset, r.boundPenalty);
    spec.min_score = 0;
    spec.prediction_headroom = 0;
    return SOLVER2_OK;
}

export int solver2_problem_new(const(solver2_problem_spec)* spec, solver2_problem** problem) nothrow
{
    return guarded({
        need(spec !is null, "spec");
        need(problem !is null, "problem");
        auto p = new ForeignProblem(*spec);
        checkDescription(p);
        *problem = cast(solver2_problem*) handleOf(p);
    });
}

export int solver2_problem_free(solver2_problem* problem) nothrow
{
    return guarded({ problemOf(problem, true); });
}

export int solver2_problem_check(solver2_problem* problem, size_t rows, ulong seed) nothrow
{
    return guarded({ checkProblem(problemOf(problem), rows, seed); });
}

export int solver2_problem_dimension(solver2_problem* problem, size_t* dimension) nothrow
{
    return guarded({
        need(dimension !is null, "dimension");
        *dimension = problemOf(problem).dimension;
    });
}

export int solver2_problem_cells(solver2_problem* problem, size_t* cells) nothrow
{
    return guarded({
        need(cells !is null, "cells");
        const p = problemOf(problem);
        *cells = p.axes.length ? p.grid.cells : 0;
    });
}

export int solver2_problem_start(solver2_problem* problem, double* genes) nothrow
{
    return guarded({
        need(genes !is null, "genes");
        give(genes, problemOf(problem).startGenes);
    });
}

export int solver2_problem_cell_code(solver2_problem* problem, size_t cell, char* buffer, size_t capacity) nothrow
{
    return guarded({
        need(buffer !is null, "buffer");
        const p = problemOf(problem);
        const grid = p.grid;
        enforce(p.axes.length > 0 && cell < grid.cells, text("cell ", cell, " is not one of the grid's ",
                                                              p.axes.length ? grid.cells : 0));
        immutable string code = grid.code(cell);
        enforce(capacity > code.length, text("the buffer holds ", capacity, " bytes, the cell's name and its end need ",
                                             code.length + 1));
        buffer[0 .. code.length] = code[];
        buffer[code.length] = '\0';
    });
}

// ─────────────────────────────────────────────────────────────────────────────
// A learned model
// ─────────────────────────────────────────────────────────────────────────────

export int solver2_model_options_default(solver2_model_options* options) nothrow @nogc
{
    if (options is null)
    {
        setError("options is NULL");
        return SOLVER2_ERROR;
    }
    immutable EnsembleOptions e;
    *options = solver2_model_options(e.members, e.hidden, e.batch, e.learningRate, e.l2, 0);
    return SOLVER2_OK;
}

export int solver2_model_new(solver2_problem* problem, ulong seed, const(solver2_model_options)* options,
                             solver2_model** model) nothrow
{
    return guarded({
        need(model !is null, "model");
        solver2_model_options o;
        if (options is null)
            solver2_model_options_default(&o);
        else
            o = *options;
        *model = cast(solver2_model*) handleOf(new ForeignModel(problemOf(problem), seed, o));
    });
}

export int solver2_model_free(solver2_model* model) nothrow
{
    return guarded({ modelOf(model, true).close(); });
}

export int solver2_model_rows(solver2_model* model, size_t* rows) nothrow
{
    return guarded({
        need(rows !is null, "rows");
        *rows = modelOf(model).learned.rows;
    });
}

// ─────────────────────────────────────────────────────────────────────────────
// A CMA-ES stage
// ─────────────────────────────────────────────────────────────────────────────

export int solver2_stage_options_default(solver2_stage_options* options) nothrow @nogc
{
    if (options is null)
    {
        setError("options is NULL");
        return SOLVER2_ERROR;
    }
    immutable StageOptions o;
    *options = solver2_stage_options(o.lambda, o.sigma, o.budget, o.seed, o.stallGenerations, o.stallTolerance,
                                     o.surrogatePool, o.surrogateShare, o.surrogateFrom, o.surrogateEpochs);
    return SOLVER2_OK;
}

export int solver2_stage_new(solver2_problem* problem, const(solver2_stage_options)* options, const(double)* start,
                             solver2_model* model, solver2_stage** stage) nothrow
{
    return guarded({
        need(options !is null, "options");
        need(stage !is null, "stage");
        auto p = problemOf(problem);
        const x0 = start is null ? null : start[0 .. p.dimension];
        *stage = cast(solver2_stage*) handleOf(new ForeignStage(p, *options, x0, modelOrNone(model)));
    });
}

export int solver2_stage_free(solver2_stage* stage) nothrow
{
    return guarded({ stageOf(stage, true); });
}

export int solver2_stage_run(solver2_stage* stage, size_t budget, solver2_stage_report_fn report, void* context,
                             solver2_stage_result* result) nothrow
{
    return guarded({
        auto s = stageOf(stage);
        scope (exit)
            if (result !is null)
                *result = s.result;
        s.advance(budget, report, context);
    });
}

export int solver2_stage_get(solver2_stage* stage, solver2_stage_result* result) nothrow
{
    return guarded({
        need(result !is null, "result");
        *result = stageOf(stage).result;
    });
}

export int solver2_stage_best(solver2_stage* stage, double* genes, solver2_outcome* outcome, double* violations,
                              double* descriptors) nothrow
{
    return guarded({
        const s = stageOf(stage);
        enforce(s.best.feasible, "the stage has no candidate within the limits yet");
        give(genes, s.best.genes);
        if (outcome !is null)
            *outcome = outcomeOf(s.best.eval);
        give(violations, s.best.eval.violations);
        give(descriptors, s.best.eval.descriptors);
    });
}

export int solver2_stage_causes(solver2_stage* stage, size_t* over) nothrow
{
    return guarded({
        need(over !is null, "over");
        give(over, stageOf(stage).causes.over);
    });
}

export int solver2_stage_save(solver2_stage* stage, ubyte* buffer, size_t capacity, size_t* size) nothrow
{
    return guarded({ giveBytes(stageOf(stage).save(), buffer, capacity, size); });
}

export int solver2_stage_load(solver2_problem* problem, solver2_model* model, const(ubyte)* bytes, size_t size,
                              solver2_stage** stage) nothrow
{
    return guarded({
        need(bytes !is null, "bytes");
        need(stage !is null, "stage");
        auto s = ForeignStage.load(problemOf(problem), modelOrNone(model), bytes[0 .. size]);
        *stage = cast(solver2_stage*) handleOf(s);
    });
}

// ─────────────────────────────────────────────────────────────────────────────
// The archive
// ─────────────────────────────────────────────────────────────────────────────

export int solver2_archive_options_default(solver2_archive_options* options) nothrow @nogc
{
    if (options is null)
    {
        setError("options is NULL");
        return SOLVER2_ERROR;
    }
    immutable ArchiveOptions o;
    immutable ArchiveRules r;
    *options = solver2_archive_options(o.evals, o.emitters, o.lambda, o.sigma, o.alpha, o.seed, o.explorers,
                                       o.explorerSigma, o.surrogatePool, o.surrogateShare, o.surrogateFrom,
                                       o.surrogateEpochs, o.surrogateImprovers, o.surrogateAim,
                                       solver2_rules(r.idleGenerations, r.infeasibleGenerations, r.infeasibleShare,
                                                     r.minSigma, r.staleGenerations, r.minProgress, r.blockIterations,
                                                     r.aimSamples));
    return SOLVER2_OK;
}

export int solver2_archive_new(solver2_problem* problem, const(solver2_archive_options)* options, const(double)* seeds,
                               size_t seed_count, solver2_model* model, solver2_archive** archive) nothrow
{
    return guarded({
        need(options !is null, "options");
        need(archive !is null, "archive");
        need(seed_count == 0 || seeds !is null, "seeds");
        auto p = problemOf(problem);
        immutable size_t n = p.dimension;
        const(double)[][] rows;
        foreach (k; 0 .. seed_count)
            rows ~= seeds[k * n .. (k + 1) * n];
        *archive = cast(solver2_archive*) handleOf(new ForeignArchive(p, *options, rows, modelOrNone(model)));
    });
}

export int solver2_archive_free(solver2_archive* archive) nothrow
{
    return guarded({ archiveOf(archive, true); });
}

export int solver2_archive_run(solver2_archive* archive, size_t evals, solver2_archive_report_fn report, void* context,
                               solver2_archive_result* result) nothrow
{
    return guarded({
        auto a = archiveOf(archive);
        scope (exit)
            if (result !is null)
                *result = a.result;
        a.advance(evals, report, context);
    });
}

export int solver2_archive_get(solver2_archive* archive, solver2_archive_result* result) nothrow
{
    return guarded({
        need(result !is null, "result");
        *result = archiveOf(archive).result;
    });
}

export int solver2_archive_elites(solver2_archive* archive, ubyte* occupied, double* value, size_t* found, double* genes,
                                  solver2_outcome* outcome, double* violations, double* descriptors) nothrow
{
    return guarded({
        const a = archiveOf(archive);
        const p = a.problem;
        immutable size_t cells = p.grid.cells, n = p.dimension, limits = p.constraints.length, nd = p.descriptorCount;
        immutable Outcome none;
        foreach (c; 0 .. cells)
        {
            const el = a.started ? &a.run.archive.elites[c] : null;
            immutable bool held = el !is null && el.occupied;
            if (occupied !is null)
                occupied[c] = held;
            if (value !is null)
                value[c] = held ? el.value : double.nan;
            if (found !is null)
                found[c] = held ? el.found : 0;
            if (outcome !is null)
                outcome[c] = outcomeOf(held ? el.eval : none);
            if (genes !is null)
            {
                genes[c * n .. (c + 1) * n] = double.nan;
                if (held)
                    genes[c * n .. (c + 1) * n] = el.genes[];
            }
            // An elite's violations may be none (a problem that gives no
            // cause counts): NaN then, as for an empty cell.
            if (violations !is null)
            {
                violations[c * limits .. (c + 1) * limits] = double.nan;
                if (held && el.eval.violations.length == limits)
                    violations[c * limits .. (c + 1) * limits] = el.eval.violations[];
            }
            if (descriptors !is null)
            {
                descriptors[c * nd .. (c + 1) * nd] = double.nan;
                if (held)
                    descriptors[c * nd .. (c + 1) * nd] = el.eval.descriptors[];
            }
        }
    });
}

export int solver2_archive_thresholds(solver2_archive* archive, double* thresholds) nothrow
{
    return guarded({
        need(thresholds !is null, "thresholds");
        const a = archiveOf(archive);
        if (a.started)
            give(thresholds, a.run.archive.threshold);
        else
            thresholds[0 .. a.problem.grid.cells] = a.problem.minScore;
    });
}

export int solver2_archive_causes(solver2_archive* archive, size_t* over) nothrow
{
    return guarded({
        need(over !is null, "over");
        const a = archiveOf(archive);
        over[0 .. a.problem.constraints.length] = 0;
        if (a.started)
            give(over, a.run.causes.over);
    });
}

export int solver2_archive_diverse(solver2_archive* archive, double within, size_t count, size_t* cells,
                                   size_t* picked) nothrow
{
    return guarded({
        need(picked !is null, "picked");
        need(cells !is null || count == 0, "cells");
        const chosen = archiveOf(archive).diverse(within, count);
        give(cells, chosen);
        *picked = chosen.length;
    });
}

export int solver2_archive_save(solver2_archive* archive, ubyte* buffer, size_t capacity, size_t* size) nothrow
{
    return guarded({ giveBytes(archiveOf(archive).save(), buffer, capacity, size); });
}

export int solver2_archive_load(solver2_problem* problem, solver2_model* model, const(ubyte)* bytes, size_t size,
                                solver2_archive** archive) nothrow
{
    return guarded({
        need(bytes !is null, "bytes");
        need(archive !is null, "archive");
        auto a = ForeignArchive.load(problemOf(problem), modelOrNone(model), bytes[0 .. size]);
        *archive = cast(solver2_archive*) handleOf(a);
    });
}
