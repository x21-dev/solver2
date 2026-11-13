/// The C interface as a caller in another language meets it: the shared
/// library loaded at run time into a process with no D runtime of its own,
/// and driven through nothing but the functions it exports.
///
/// This program is built with -betterC, so it has no D runtime, no
/// collector and no exceptions: it is a C program written in D. What it
/// knows of the library is include/solver2.h, read by ImportC, and it
/// looks every function of the header up in the library by its name.
///
/// It checks what only this arrangement can: that the library loads and
/// exports every function, that it refuses calls before solver2_init, that
/// a run through it is deterministic and resumes bit for bit from its
/// bytes, that a failed evaluation leaves it whole, and that threads the D
/// runtime never heard of can call in, at once and one after another. That the numbers are the searches'
/// own is the binding's unittests' part (source/solver2c/tests.d).
///
///     driver [path to the library]        default: ../solver2.dll, ../libsolver2.so
module driver;

import core.stdc.stdio : _IONBF, printf, setvbuf, stdout;
import core.stdc.stdlib : exit, free, malloc;
import core.stdc.string : memcmp, strstr;

import h = solver2_h;

// ─────────────────────────────────────────────────────────────────────────────
// The library, loaded
// ─────────────────────────────────────────────────────────────────────────────

private bool isFunction(string name)()
{
    return name.length > 8 && name[0 .. 8] == "solver2_" && is(typeof(__traits(getMember, h, name)) == function);
}

// One pointer per function of the header, under the function's name: the
// calls below read as C.
static foreach (name; __traits(allMembers, h))
    static if (isFunction!name)
        mixin("__gshared typeof(&h." ~ name ~ ") " ~ name ~ ";");

version (Windows)
{
    import core.sys.windows.winbase : CloseHandle, CreateThread, GetProcAddress, INFINITE, LoadLibraryA,
                                      WaitForSingleObject;
    import core.sys.windows.windef : HANDLE;

    enum defaultLibrary = "../solver2.dll";

    private void* openLibrary(const(char)* path) { return LoadLibraryA(path); }
    private void* symbol(void* library, const(char)* name) { return GetProcAddress(library, name); }
}
else
{
    import core.sys.posix.dlfcn : dlopen, dlsym, RTLD_LOCAL, RTLD_NOW;
    import core.sys.posix.pthread : pthread_create, pthread_join, pthread_t;

    enum defaultLibrary = "../libsolver2.so";

    private void* openLibrary(const(char)* path) { return dlopen(path, RTLD_NOW | RTLD_LOCAL); }
    private void* symbol(void* library, const(char)* name) { return dlsym(library, name); }
}

/// Loads the library and finds every function the header declares.
private size_t load(const(char)* path)
{
    void* library = openLibrary(path);
    if (library is null)
    {
        printf("FAILED: the library %s doesn't load\n", path);
        exit(1);
    }
    size_t found;
    static foreach (name; __traits(allMembers, h))
        static if (isFunction!name)
        {
            mixin(name) = cast(typeof(mixin(name))) symbol(library, name.ptr);
            if (mixin(name) is null)
            {
                printf("FAILED: the library doesn't export %s\n", name.ptr);
                exit(1);
            }
            ++found;
        }
    return found;
}

// ─────────────────────────────────────────────────────────────────────────────
// Checks
// ─────────────────────────────────────────────────────────────────────────────

private void check(bool held, const(char)* what, int line = __LINE__)
{
    if (held)
        return;
    printf("FAILED at line %d: %s\n", line, what);
    if (solver2_last_error !is null)
        printf("  the library's last error: %s\n", solver2_last_error());
    exit(1);
}

private void ok(int code, int line = __LINE__)
{
    check(code == h.SOLVER2_OK, "a call that should have succeeded", line);
}

/// The call failed with `code`, and its message holds `part`.
private void refused(int got, const(char)* part, int code = h.SOLVER2_ERROR, int line = __LINE__)
{
    check(got == code, "a call that should have been refused, or refused otherwise", line);
    check(strstr(solver2_last_error(), part) !is null, part, line);
}

/// A saved state, in memory of this program's own.
private struct Blob
{
    ubyte* data;
    size_t size;
}

private bool same(const Blob a, const Blob b)
{
    return a.size == b.size && a.size > 0 && memcmp(a.data, b.data, a.size) == 0;
}

private Blob bytesOf(h.solver2_stage* s)
{
    Blob b;
    ok(solver2_stage_save(s, null, 0, &b.size));
    b.data = cast(ubyte*) malloc(b.size);
    ok(solver2_stage_save(s, b.data, b.size, &b.size));
    return b;
}

private Blob bytesOf(h.solver2_archive* a)
{
    Blob b;
    ok(solver2_archive_save(a, null, 0, &b.size));
    b.data = cast(ubyte*) malloc(b.size);
    ok(solver2_archive_save(a, b.data, b.size, &b.size));
    return b;
}

// ─────────────────────────────────────────────────────────────────────────────
// The problem
// ─────────────────────────────────────────────────────────────────────────────

/// The evaluation's context: its batches counted, and one made to fail.
private struct Calls
{
    size_t batches, failAt;
}

private double clamp01(double g) { return g < 0 ? 0 : (g > 1 ? 1 : g); }

/// Maximise 10 - (x - 1)² - (y - 0.8)² - z² over x and y in 0 .. 2 and z in
/// -1 .. 1, with x + y ≤ 1.9 hard, z ≥ -0.9 terminal and x ≤ 1.5 soft. The
/// descriptors are x and y.
private extern (C) int evaluate(void* context, size_t rows, const(double)* genes, int feasibleOnly, ubyte* scored,
                                double* objective, double* violations, double* descriptors)
{
    auto c = cast(Calls*) context;
    ++c.batches;
    if (c.batches == c.failAt)
        return 7;
    foreach (k; 0 .. rows)
    {
        immutable double x = 2 * clamp01(genes[3 * k]), y = 2 * clamp01(genes[3 * k + 1]),
                         z = -1 + 2 * clamp01(genes[3 * k + 2]);
        violations[3 * k] = x + y > 1.9 ? x + y - 1.9 : 0;
        violations[3 * k + 1] = z < -0.9 ? -0.9 - z : 0;
        violations[3 * k + 2] = x > 1.5 ? x - 1.5 : 0;
        if (!feasibleOnly || violations[3 * k] == 0)
        {
            scored[k] = 1;
            objective[k] = 10 - (x - 1) * (x - 1) - (y - 0.8) * (y - 0.8) - z * z;
        }
        descriptors[2 * k] = x;
        descriptors[2 * k + 1] = y;
    }
    return 0;
}

private h.solver2_problem* problem(Calls* calls)
{
    static __gshared h.solver2_gene[3] genes = [h.solver2_gene("x", 0, 2, 1, h.SOLVER2_LINEAR),
                                                h.solver2_gene("y", 0, 2, 0.6, h.SOLVER2_LINEAR),
                                                h.solver2_gene("z", -1, 1, 0, h.SOLVER2_LINEAR)];
    static __gshared h.solver2_limit[3] limits = [h.solver2_limit("sum", h.SOLVER2_HARD, 1, null),
                                                  h.solver2_limit("floor", h.SOLVER2_TERMINAL, 1, null),
                                                  h.solver2_limit("x_cap", h.SOLVER2_SOFT, 0.5, null)];
    static __gshared double[3] edges = [0.5, 1, 1.5];
    static __gshared const(char)*[2] fingerprint = ["problem", "driver 1"];
    h.solver2_axis[2] axes = [h.solver2_axis("x", edges.ptr, 3, "abcd", 0, 0.5),
                              h.solver2_axis("y", edges.ptr, 3, "abcd", 0, 0.5)];
    h.solver2_problem_spec spec;
    ok(solver2_problem_spec_default(&spec));
    spec.genes = genes.ptr;
    spec.gene_count = 3;
    spec.limits = limits.ptr;
    spec.limit_count = 3;
    spec.axes = axes.ptr;
    spec.axis_count = 2;
    spec.ranking = h.solver2_ranking(1000, 500, 20);
    spec.min_score = 5;
    spec.prediction_headroom = 0.05;
    spec.fingerprint = fingerprint.ptr;
    spec.fingerprint_count = 1;
    spec.evaluate = &evaluate;
    spec.context = calls;
    h.solver2_problem* p;
    ok(solver2_problem_new(&spec, &p));
    return p;
}

private h.solver2_model* model(h.solver2_problem* p, size_t threads)
{
    h.solver2_model_options o;
    ok(solver2_model_options_default(&o));
    o.members = 2;
    o.hidden = 16;
    o.threads = threads;
    h.solver2_model* m;
    ok(solver2_model_new(p, 7, &o, &m));
    return m;
}

private h.solver2_stage_options stageOptions()
{
    h.solver2_stage_options o;
    ok(solver2_stage_options_default(&o));
    o.lambda = 8;
    o.sigma = 0.2;
    o.budget = 400;
    o.seed = 3;
    o.stall_generations = 0;
    return o;
}

private h.solver2_archive_options archiveOptions()
{
    h.solver2_archive_options o;
    ok(solver2_archive_options_default(&o));
    o.evals = 400;
    o.emitters = 3;
    o.explorers = 1;
    o.lambda = 8;
    o.seed = 5;
    o.surrogate_pool = 4;
    o.surrogate_from = 3;
    o.surrogate_epochs = 2;
    o.surrogate_improvers = 1;
    return o;
}

/// Reports that stop a run at one generation or iteration.
private extern (C) int stopStage(void* context, const(h.solver2_stage_generation)* g)
{
    return g.generation == *cast(size_t*) context;
}

private extern (C) int stopArchive(void* context, const(h.solver2_archive_iteration)* it)
{
    return it.iteration == *cast(size_t*) context;
}

/// A stage run to its end in one call: its bytes.
// (The options are not const here: ImportC drops const from a pointer to a
// struct, so the header's `const solver2_stage_options*` wants a mutable one.)
private Blob wholeStage(h.solver2_problem* p)
{
    auto o = stageOptions();
    h.solver2_stage* s;
    ok(solver2_stage_new(p, &o, null, null, &s));
    h.solver2_stage_result r;
    ok(solver2_stage_run(s, 0, null, null, &r));
    check(r.finished && r.feasible && r.generations == 50 && r.evaluations == 400, "the stage's end");
    auto b = bytesOf(s);
    ok(solver2_stage_free(s));
    return b;
}

/// A model-guided run of the archive to its end in one call, its model on
/// `threads` threads: its bytes.
private Blob wholeArchive(h.solver2_problem* p, size_t threads)
{
    auto o = archiveOptions();
    auto m = model(p, threads);
    h.solver2_archive* a;
    ok(solver2_archive_new(p, &o, null, 0, m, &a));
    h.solver2_archive_result r;
    ok(solver2_archive_run(a, 0, null, null, &r));
    check(r.finished && r.iterations == 17 && r.coverage >= 6 && r.best_score > 9, "the archive's end");
    auto b = bytesOf(a);
    ok(solver2_archive_free(a));
    ok(solver2_model_free(m));
    return b;
}

// ─────────────────────────────────────────────────────────────────────────────
// Other threads
// ─────────────────────────────────────────────────────────────────────────────

/// What another thread does, and what it found. Each has a problem of its
/// own: one call at a time per object is the interface's rule.
private struct Errand
{
    bool withArchive;
    Blob stage, archive;
    bool refusedThere;
}

private void errand(Errand* e)
{
    Calls calls;
    auto p = problem(&calls);
    e.stage = wholeStage(p);
    if (e.withArchive)
        e.archive = wholeArchive(p, 2);
    // An error on this thread is this thread's.
    size_t n;
    e.refusedThere = solver2_problem_dimension(null, &n) == h.SOLVER2_ERROR
                  && strstr(solver2_last_error(), "a problem: the handle is NULL") !is null;
    ok(solver2_problem_free(p));
}

version (Windows)
{
    private alias ThreadId = HANDLE;

    private extern (Windows) uint threadMain(void* arg)
    {
        errand(cast(Errand*) arg);
        return 0;
    }

    private ThreadId start(Errand* e)
    {
        HANDLE t = CreateThread(null, 0, &threadMain, e, 0, null);
        check(t !is null, "another thread");
        return t;
    }

    private void join(ThreadId t)
    {
        WaitForSingleObject(t, INFINITE);
        CloseHandle(t);
    }
}
else
{
    private alias ThreadId = pthread_t;

    private extern (C) void* threadMain(void* arg)
    {
        errand(cast(Errand*) arg);
        return null;
    }

    private ThreadId start(Errand* e)
    {
        pthread_t t;
        check(pthread_create(&t, null, &threadMain, e) == 0, "another thread");
        return t;
    }

    private void join(ThreadId t)
    {
        pthread_join(t, null);
    }
}

// ─────────────────────────────────────────────────────────────────────────────

extern (C) int main(int argc, char** argv)
{
    // Unbuffered, so that a crash doesn't take the lines before it along:
    // into a pipe (dub's, CI's) they would wait in a buffer and be lost.
    setvbuf(stdout, null, _IONBF, 0);
    const(char)* path = argc > 1 ? argv[1] : defaultLibrary.ptr;
    immutable size_t functions = load(path);
    printf("%s: all %d functions of solver2.h found\n", path, cast(int) functions);

    // What needs no runtime answers at once; the rest says what is missing.
    check(solver2_abi() == h.SOLVER2_ABI, "the library's interface version is the header's");
    printf("interface version %d; bounds checks %s\n", solver2_abi(), solver2_checked() ? "on".ptr : "off".ptr);
    size_t n;
    refused(solver2_problem_dimension(null, &n), "solver2_init has not been called");
    refused(solver2_term(), "without a solver2_init");

    ok(solver2_init());
    Calls calls;
    auto p = problem(&calls);
    size_t cells;
    ok(solver2_problem_dimension(p, &n));
    ok(solver2_problem_cells(p, &cells));
    check(n == 3 && cells == 16, "the problem's size");
    ok(solver2_problem_check(p, 8, 1));
    char[3] code;
    ok(solver2_problem_cell_code(p, 14, code.ptr, code.length));
    check(code[0] == 'c' && code[1] == 'd' && code[2] == 0, "a cell's name");
    printf("a problem made and checked: %d genes, %d cells\n", cast(int) n, cast(int) cells);

    // A stage, the same twice; stopped, saved, freed, loaded and gone on
    // with, the same again.
    const stage = wholeStage(p);
    check(same(stage, wholeStage(p)), "the same stage twice");
    {
        auto o = stageOptions();
        h.solver2_stage* s;
        ok(solver2_stage_new(p, &o, null, null, &s));
        size_t stopAt = 10;
        h.solver2_stage_result r;
        ok(solver2_stage_run(s, 0, &stopStage, &stopAt, &r));
        check(r.generations == 10 && r.stopped && !r.finished, "a stage stopped by its report");
        auto saved = bytesOf(s);
        ok(solver2_stage_free(s));
        h.solver2_stage* back;
        ok(solver2_stage_load(p, null, saved.data, saved.size, &back));
        ok(solver2_stage_run(back, 0, null, null, &r));
        check(r.finished && same(bytesOf(back), stage), "a stage resumed from its bytes");
        double[3] genes;
        h.solver2_outcome e;
        ok(solver2_stage_best(back, genes.ptr, &e, null, null));
        check(e.scored && e.value == r.best_value && e.hard == 0 && e.terminal == 0 && r.best_value > 9.9,
              "the stage's best");
        printf("a stage: %d generations, best %.6f at %.4f %.4f %.4f; resumed from its %d bytes, the same bits\n",
               cast(int) r.generations, r.best_value, genes[0], genes[1], genes[2], cast(int) saved.size);
        ok(solver2_stage_free(back));
        free(saved.data);
    }

    // The archive with a learned model: the same on 0 and 3 of the model's
    // threads; stopped, saved, loaded with a new model and gone on with, the
    // same again.
    const archive = wholeArchive(p, 0);
    check(same(archive, wholeArchive(p, 3)), "the same archive, whatever the model's threads");
    {
        auto o = archiveOptions();
        auto m = model(p, 2);
        h.solver2_archive* a;
        ok(solver2_archive_new(p, &o, null, 0, m, &a));
        size_t stopAt = 6;
        h.solver2_archive_result r;
        ok(solver2_archive_run(a, 0, &stopArchive, &stopAt, &r));
        check(r.iterations == 6 && r.stopped && !r.finished, "an archive stopped by its report");
        auto saved = bytesOf(a);
        ok(solver2_archive_free(a));
        ok(solver2_model_free(m));
        auto m2 = model(p, 1);
        h.solver2_archive* back;
        ok(solver2_archive_load(p, m2, saved.data, saved.size, &back));
        ok(solver2_archive_run(back, 0, null, null, &r));
        check(r.finished && same(bytesOf(back), archive), "an archive resumed from its bytes");
        ubyte[16] occupied;
        double[16] values;
        ok(solver2_archive_elites(back, occupied.ptr, values.ptr, null, null, null, null, null));
        size_t held;
        foreach (c; 0 .. 16)
        {
            held += occupied[c];
            check(occupied[c] ? values[c] >= 5 && values[c] <= r.best_score : values[c] != values[c], "a cell's value");
        }
        check(held == r.coverage && !occupied[15], "the cells that hold a candidate");
        printf("an archive with a model: %d of %d cells, best %.6f; resumed from its %d bytes, the same bits\n",
               cast(int) r.coverage, cast(int) r.cells, r.best_score, cast(int) saved.size);
        ok(solver2_archive_free(back));
        ok(solver2_model_free(m2));
        free(saved.data);
    }

    // The evaluation fails in the middle: the call says so, and made again
    // it ends as the run that never failed.
    {
        Calls failing = {failAt: 9};
        auto q = problem(&failing);
        auto o = archiveOptions();
        auto m = model(q, 2);
        h.solver2_archive* a;
        ok(solver2_archive_new(q, &o, null, 0, m, &a));
        h.solver2_archive_result r;
        refused(solver2_archive_run(a, 0, null, null, &r), "it returned 7", h.SOLVER2_CALLBACK);
        check(r.started && r.iterations == 7 && !r.finished, "a run left after its last whole iteration");
        ok(solver2_archive_run(a, 0, null, null, &r));
        check(r.finished && same(bytesOf(a), archive), "a run gone on with after a failed evaluation");
        printf("an evaluation that fails in iteration 8: the run goes on from 7, to the same bits\n");
        ok(solver2_archive_free(a));
        ok(solver2_model_free(m));
        ok(solver2_problem_free(q));
    }

    // A model freed while its run is stopped: the run holds it, and goes on
    // without its threads to the same bits. And no cell is picked into no
    // room.
    {
        auto o = archiveOptions();
        auto m = model(p, 2);
        h.solver2_archive* a;
        ok(solver2_archive_new(p, &o, null, 0, m, &a));
        size_t stopAt = 6;
        h.solver2_archive_result r;
        ok(solver2_archive_run(a, 0, &stopArchive, &stopAt, &r));
        check(r.iterations == 6 && r.stopped, "an archive stopped by its report");
        ok(solver2_model_free(m));
        ok(solver2_archive_run(a, 0, null, null, &r));
        check(r.finished && same(bytesOf(a), archive), "a run gone on with after its model was freed");
        size_t[2] room = [77, 77];
        size_t picked = 9;
        ok(solver2_archive_diverse(a, 3, 0, room.ptr, &picked));
        check(picked == 0 && room[0] == 77, "no cell picked for a count of 0");
        ok(solver2_archive_diverse(a, 3, 1, room.ptr, &picked));
        check(picked == 1 && room[0] == r.best_cell && room[1] == 77, "the best cell alone for a count of 1");
        printf("a model freed in mid-run: the run goes on without its threads, to the same bits\n");
        ok(solver2_archive_free(a));
    }

    // Threads the D runtime never heard of, each making many calls: the
    // same runs, the same bits, and their errors their own. Two at once,
    // while this one runs a stage too; then a string of short ones, each
    // gone before the next (a thread is attached to the runtime on its
    // first call and detached as it ends).
    {
        refused(solver2_stage_free(null), "a stage: the handle is NULL");
        Errand[2] both = [Errand(true), Errand(true)];
        ThreadId[2] running = [start(&both[0]), start(&both[1])];
        auto meanwhile = wholeStage(p);
        foreach (t; running)
            join(t);
        check(same(meanwhile, stage), "this thread's stage, beside two others");
        free(meanwhile.data);
        foreach (ref e; both)
        {
            check(same(e.stage, stage), "a stage from another thread");
            check(same(e.archive, archive), "an archive with a model from another thread");
            check(e.refusedThere, "another thread's error");
            free(e.stage.data);
            free(e.archive.data);
        }
        check(strstr(solver2_last_error(), "a stage: the handle is NULL") !is null,
              "each thread's last error is its own");
        foreach (i; 0 .. 12)
        {
            Errand e;
            join(start(&e));
            check(same(e.stage, stage) && e.refusedThere, "a stage from a short-lived thread");
            free(e.stage.data);
        }
        printf("two threads unknown to the D runtime at once, then twelve one after another: the same bits\n");
    }

    // Handles outlive whatever the collector does meanwhile: forty more
    // stages' worth of garbage, and the first problem still answers.
    auto kept = model(p, 0);
    foreach (i; 0 .. 40)
    {
        auto again = wholeStage(p);
        check(same(again, stage), "the same stage, again and again");
        free(again.data);
    }
    size_t rows = 1;
    ok(solver2_model_rows(kept, &rows));
    ok(solver2_problem_dimension(p, &n));
    check(rows == 0 && n == 3, "the handles held through the collections");

    // A handle freed is refused, not followed.
    ok(solver2_model_free(kept));
    refused(solver2_model_rows(kept, &rows), "a model: not a handle");
    ok(solver2_problem_free(p));
    refused(solver2_problem_dimension(p, &n), "a problem: not a handle");

    ok(solver2_term());
    refused(solver2_problem_dimension(p, &n), "solver2_init has not been called");
    printf("all passed\n");
    return 0;
}
