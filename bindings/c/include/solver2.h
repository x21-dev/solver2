/* solver2.h: the C interface of solver2, a deterministic black-box
 * optimiser: CMA-ES stages and a quality-diversity archive (CMA-MAE) over a
 * problem whose evaluation is the caller's function.
 *
 * The library is solver2's D code built as a shared library (solver2.dll,
 * libsolver2.so). This header is written by hand; the binding's tests hold
 * it to the library, declaration by declaration.
 *
 * The rules:
 *
 * - Call solver2_init before anything else, and solver2_term at the end.
 * - Every other function returns SOLVER2_OK or an error code; nothing is
 *   thrown across the boundary. After an error, solver2_last_error gives
 *   the message, to the thread that made the call.
 * - Objects are handles: the library keeps what they stand for, and a
 *   handle is good until it is freed. A handle that was freed, or never
 *   was one, is an error, not a crash.
 * - Results are copied into buffers the caller owns. Their sizes follow
 *   from the problem: its genes (the dimension), its limits, its axes (the
 *   descriptors) and its grid's cells. A buffer may be NULL to skip it.
 * - Strings are UTF-8 and end in 0. The library copies what it is given.
 * - Any thread may call, and several at once. One call at a time per
 *   object: a run, its model and its problem are one object for this. The
 *   evaluation is called on the thread that called the library. A thread
 *   that calls is made known to the D runtime in the library on its first
 *   call and forgotten when it ends; the library sees to both.
 * - Deterministic: the same problem, options and seed give the same bits
 *   on any thread count of one platform, as long as the evaluation gives
 *   the same bits for the same genes whatever else is in its batch.
 */
#ifndef SOLVER2_H
#define SOLVER2_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* The interface's version: solver2_abi() of a library built from this
 * header returns it. It changes when a struct or a function here does. */
#define SOLVER2_ABI 2

/* Return codes. */
enum
{
    SOLVER2_OK = 0,
    SOLVER2_ERROR = 1,    /* see solver2_last_error */
    SOLVER2_CALLBACK = 2  /* the caller's own function reported a failure */
};

/* A gene's scale: how [0, 1] maps to lo .. hi. */
enum
{
    SOLVER2_LINEAR = 0,
    SOLVER2_LOG = 1,      /* needs lo > 0 */
    SOLVER2_INTEGER = 2   /* rounds to the nearest whole value */
};

/* A limit's tier. */
enum
{
    SOLVER2_HARD = 0,     /* any breach ranks below every candidate within */
    SOLVER2_TERMINAL = 1, /* ranked after hard, before the objective */
    SOLVER2_SOFT = 2      /* costs weight per unit of breach */
};

typedef struct solver2_problem solver2_problem;
typedef struct solver2_model solver2_model;
typedef struct solver2_stage solver2_stage;
typedef struct solver2_archive solver2_archive;

/* ── The library ──────────────────────────────────────────────────────── */

/* Starts the D runtime inside the library. Calls are counted: as many
 * solver2_term as solver2_init end it. Not to be called from two threads
 * at once. The first call is to come from a thread that lives as long as
 * the library is used (the main thread): the runtime takes it for its
 * own main thread. */
int solver2_init(void);

/* Ends it, with the last of the calls that began it: every handle is then
 * void. No other thread may be inside the library then, and the library
 * is not to be started again in the same process. */
int solver2_term(void);

/* The message of the calling thread's last error; "" when none. Good until
 * that thread's next call that fails. */
const char* solver2_last_error(void);

/* SOLVER2_ABI as the library was built. */
int solver2_abi(void);

/* 1 when the library checks every array index (the usual build), 0 when
 * it was built without (the "unchecked" configuration). */
int solver2_checked(void);

/* ── A problem ────────────────────────────────────────────────────────── */

/* One searched parameter: gene 0 is lo, gene 1 is hi. neutral is where a
 * search starts (NaN: the middle). */
typedef struct solver2_gene
{
    const char* name;
    double lo;
    double hi;
    double neutral;
    int scale;            /* SOLVER2_LINEAR, _LOG or _INTEGER */
} solver2_gene;

/* One limit. weight: the cost per unit of breach of a soft one. group: a
 * label of the caller's own, or NULL. */
typedef struct solver2_limit
{
    const char* name;
    int tier;             /* SOLVER2_HARD, _TERMINAL or _SOFT */
    double weight;
    const char* group;
} solver2_limit;

/* One axis of the archive's grid, over one descriptor: values below
 * edges[0] in bin 0, from the last edge up in the last bin, NaN in
 * nan_bin. letters: one character per bin (edge_count + 1), for a cell's
 * name. width: the unit of distance on this axis, a typical bin's width. */
typedef struct solver2_axis
{
    const char* name;
    const double* edges;  /* ascending */
    size_t edge_count;
    const char* letters;
    size_t nan_bin;
    double width;
} solver2_axis;

/* How candidates rank: over a hard limit at hard_offset plus the breach,
 * over a terminal one at terminal_offset plus the breach, else at minus
 * the value; plus bound_penalty per squared gene unit outside [0, 1]. The
 * offsets must be above every value the problem can have. */
typedef struct solver2_ranking
{
    double hard_offset;
    double terminal_offset;
    double bound_penalty;
} solver2_ranking;

/* The evaluation of a batch: `rows` candidates, candidate k's genes at
 * genes[k * dimension ..]. Genes come as sampled and may lie outside
 * [0, 1]: read them clamped. To fill, for candidate k:
 *
 * - violations[k * limits ..]: how far over each limit, 0 within; never
 *   negative. They start at 0.
 * - when it is within its hard and terminal limits: scored[k] = 1,
 *   objective[k] (higher is better) and descriptors[k * axes ..].
 *   A candidate over a hard limit may be left unscored when feasible_only
 *   is not 0. scored starts at 0, the others at NaN.
 *
 * The same genes must give the same bits whatever the batch. Return 0, or
 * anything else to report a failure: the call that asked for the
 * evaluation then returns SOLVER2_CALLBACK. A pointer is NULL when the
 * problem has none of what it holds (no limits, no axes). */
typedef int (*solver2_evaluate_fn)(void* context, size_t rows, const double* genes, int feasible_only,
                                   unsigned char* scored, double* objective, double* violations,
                                   double* descriptors);

/* A problem as data. Start from solver2_problem_spec_default. */
typedef struct solver2_problem_spec
{
    const solver2_gene* genes;
    size_t gene_count;
    const solver2_limit* limits;      /* may be none */
    size_t limit_count;
    const solver2_axis* axes;         /* none: the problem has no archive */
    size_t axis_count;
    solver2_ranking ranking;
    double min_score;                 /* the value's lower bound: the archive's thresholds start there */
    double prediction_headroom;       /* how far over a limit a model's prediction may be */
    const char* const* fingerprint;   /* name, value, name, value, ...: what a saved run must match */
    size_t fingerprint_count;         /* the number of pairs */
    solver2_evaluate_fn evaluate;
    void* context;                    /* passed to evaluate as it is */
} solver2_problem_spec;

/* Fills `spec` with nothing but the defaults: the ranking 100, 50 and 10,
 * min_score 0, prediction_headroom 0. */
int solver2_problem_spec_default(solver2_problem_spec* spec);

/* Makes a problem from `spec`, which is copied but for `context`. Its
 * description is checked: SOLVER2_ERROR says what is wrong with it. */
int solver2_problem_new(const solver2_problem_spec* spec, solver2_problem** problem);
int solver2_problem_free(solver2_problem* problem);

/* Checks the problem's evaluation on its start genes, both corners and
 * `rows` random rows drawn with `seed` (8 and 1 are the usual): the same
 * bits evaluated twice, in another order and alone; and no candidate
 * within the limits below min_score. About 5 * (rows + 3) evaluations. */
int solver2_problem_check(solver2_problem* problem, size_t rows, uint64_t seed);

/* The number of genes, and of cells of the grid (0 without axes). */
int solver2_problem_dimension(solver2_problem* problem, size_t* dimension);
int solver2_problem_cells(solver2_problem* problem, size_t* cells);

/* The genes a search starts from when given none: `dimension` values. */
int solver2_problem_start(solver2_problem* problem, double* genes);

/* What genes stand for: `rows` gene vectors in `genes`, as many rows of
 * `dimension` physical values into `values`, each gene clamped to [0, 1]
 * and mapped over its range by its scale. These are the library's own
 * numbers (a log scale goes through its exp and log), so an evaluation
 * that reads its genes through here gets the bits a problem written in D
 * gets, whatever the caller's own maths library rounds to. It may be
 * called from inside the evaluation. */
int solver2_problem_physical(solver2_problem* problem, size_t rows, const double* genes, double* values);

/* A cell's name, one letter per axis, into `buffer` of `capacity` bytes
 * (the axes and one more). */
int solver2_problem_cell_code(solver2_problem* problem, size_t cell, char* buffer, size_t capacity);

/* One candidate, as a search reads it. value: the objective less the soft
 * breach, NaN when not scored. */
typedef struct solver2_outcome
{
    int scored;
    double objective;
    double value;
    double hard;          /* the breach of each tier, 0 within */
    double terminal;
    double soft;
} solver2_outcome;

/* ── A learned model ──────────────────────────────────────────────────── */

/* An ensemble of small networks that learns every evaluation of a run and
 * predicts candidates before they are evaluated. threads: worker threads
 * of its own for training and predicting (0: none); they never call the
 * evaluation, and their number doesn't change the bits. */
typedef struct solver2_model_options
{
    size_t members;
    size_t hidden;
    size_t batch;
    double learning_rate;
    double l2;
    size_t threads;
} solver2_model_options;

int solver2_model_options_default(solver2_model_options* options);

/* A model of `problem`, seeded with `seed`. options NULL: the defaults.
 * It serves one stage or archive at a time, and keeps what it learned in
 * one for the next. */
int solver2_model_new(solver2_problem* problem, uint64_t seed, const solver2_model_options* options,
                      solver2_model** model);
int solver2_model_free(solver2_model* model);

/* The evaluations it has learned from. */
int solver2_model_rows(solver2_model* model, size_t* rows);

/* ── A CMA-ES stage ───────────────────────────────────────────────────── */

/* lambda, sigma and budget have no defaults: set them. With a model, set
 * surrogate_pool above 0: from generation surrogate_from on the stage
 * samples surrogate_pool times lambda and evaluates lambda of them,
 * surrogate_share of those the ones the model ranks best. */
typedef struct solver2_stage_options
{
    size_t lambda;            /* the population */
    double sigma;             /* the initial step size, in gene units */
    size_t budget;            /* evaluations, rounded up to whole generations */
    uint64_t seed;
    size_t stall_generations; /* stop when the best rose by less than */
    double stall_tolerance;   /* this over that many generations; 0: never */
    size_t surrogate_pool;
    double surrogate_share;
    size_t surrogate_from;
    size_t surrogate_epochs;
} solver2_stage_options;

int solver2_stage_options_default(solver2_stage_options* options);

/* One generation, as a stage reports it. */
typedef struct solver2_stage_generation
{
    size_t generation;        /* from 1 */
    size_t evaluations;       /* of the stage so far */
    double best;              /* the best value so far; -inf while none is within the limits */
    double population_best;   /* this generation's; NaN when none is within the limits */
    double sigma;
    size_t infeasible;        /* candidates over a hard or terminal limit */
    size_t sampled;           /* pre-screened generations: the candidates the model ranked */
    size_t predicted_feasible;
} solver2_stage_generation;

/* Called after each generation. Return 0 to go on, a positive value to
 * stop the stage here (solver2_stage_run goes on from there when called
 * again), a negative one to report a failure: the stage stops here too,
 * and the call returns SOLVER2_CALLBACK. */
typedef int (*solver2_stage_report_fn)(void* context, const solver2_stage_generation* generation);

/* Where a stage stands. */
typedef struct solver2_stage_result
{
    size_t generations;
    size_t evaluations;
    int finished;             /* its budget is spent, or it stalled */
    int stalled;
    int stopped;              /* the last run was stopped by its report */
    int feasible;             /* it has a best within the limits */
    double best_value;        /* NaN while it has none */
} solver2_stage_result;

/* A stage over `problem` from `start` (`dimension` genes; NULL: the
 * problem's start), pre-screened by `model` (or NULL). Nothing is
 * evaluated yet. */
int solver2_stage_new(solver2_problem* problem, const solver2_stage_options* options, const double* start,
                      solver2_model* model, solver2_stage** stage);
int solver2_stage_free(solver2_stage* stage);

/* Runs the stage until `budget` evaluations are spent in all (0: the
 * options'), it stalls or `report` (or NULL) stops it. Called again it
 * goes on, and a stage stopped and gone on with ends as the one that never
 * stopped, bit for bit. A larger budget extends a finished stage.
 *
 * When the evaluation fails (SOLVER2_CALLBACK) the stage and its model
 * are left as they were after the last whole generation, and the call can
 * be made again. `result` (or NULL) is filled either way. */
int solver2_stage_run(solver2_stage* stage, size_t budget, solver2_stage_report_fn report, void* context,
                      solver2_stage_result* result);

int solver2_stage_get(solver2_stage* stage, solver2_stage_result* result);

/* Its best candidate: genes (`dimension`, in [0, 1]), outcome, violations
 * (one per limit) and descriptors (one per axis). SOLVER2_ERROR while it
 * has none within the limits. */
int solver2_stage_best(solver2_stage* stage, double* genes, solver2_outcome* outcome, double* violations,
                       double* descriptors);

/* Per limit, the candidates over it so far. */
int solver2_stage_causes(solver2_stage* stage, size_t* over);

/* The stage as bytes, its model's state included: `size` gets their
 * number, and `buffer` (or NULL) the bytes when `capacity` holds them
 * (SOLVER2_ERROR when it doesn't). */
int solver2_stage_save(solver2_stage* stage, unsigned char* buffer, size_t capacity, size_t* size);

/* A stage from solver2_stage_save's bytes, on a problem with the same
 * fingerprint and description. `model`: one made as the saved stage's was
 * (the same seed and options, but for threads), when it had one. */
int solver2_stage_load(solver2_problem* problem, solver2_model* model, const unsigned char* bytes, size_t size,
                       solver2_stage** stage);

/* ── The archive ──────────────────────────────────────────────────────── */

/* When an emitter restarts and when an explorer gives its target up. */
typedef struct solver2_rules
{
    size_t idle_generations;
    size_t infeasible_generations;
    double infeasible_share;
    double min_sigma;
    size_t stale_generations;
    double min_progress;
    size_t block_iterations;
    size_t aim_samples;
} solver2_rules;

typedef struct solver2_archive_options
{
    size_t evals;             /* the budget beyond the seeds, in whole iterations */
    size_t emitters;          /* explorers included */
    size_t lambda;            /* per emitter; 0: CMA-ES's default for the dimension */
    double sigma;             /* the improvers' initial step size */
    double alpha;             /* the archive learning rate */
    uint64_t seed;
    size_t explorers;         /* the last emitters head for empty cells */
    double explorer_sigma;
    size_t surrogate_pool;    /* with a model: above 0 */
    double surrogate_share;
    size_t surrogate_from;
    size_t surrogate_epochs;
    int surrogate_improvers;  /* pre-screen the improvers too */
    int surrogate_aim;        /* start an explorer where the model puts it nearest its target */
    solver2_rules rules;
} solver2_archive_options;

int solver2_archive_options_default(solver2_archive_options* options);

/* One iteration (every emitter's generation), as the archive reports it. */
typedef struct solver2_archive_iteration
{
    size_t iteration;         /* from 1 */
    size_t evaluations;       /* of the run so far, the seeds included */
    size_t coverage;          /* cells that hold a candidate */
    double qd_score;          /* the sum of their values above min_score */
    double best_score;
    size_t improved;          /* this iteration's candidates that beat their cell's threshold */
    size_t new_cells;         /* ... that opened a cell */
    size_t restarts;          /* emitters that restart or take a new target after it */
} solver2_archive_iteration;

/* Called after each iteration; its return value as solver2_stage_report_fn's. */
typedef int (*solver2_archive_report_fn)(void* context, const solver2_archive_iteration* iteration);

/* Where a run of the archive stands. */
typedef struct solver2_archive_result
{
    size_t cells;             /* the grid's */
    size_t coverage;
    size_t evaluations;       /* the seeds included */
    size_t seed_evaluations;
    size_t iterations;
    size_t lambda;            /* each emitter's population */
    int started;              /* the seeds are evaluated */
    int finished;             /* its budget is spent */
    int stopped;              /* the last run was stopped by its report */
    double qd_score;
    double best_score;        /* NaN while it holds nothing */
    size_t best_cell;         /* (size_t) -1 while it holds nothing */
    size_t restarts_idle;
    size_t restarts_infeasible;
    size_t restarts_converged;
    size_t reached;           /* targets the explorers reached */
    size_t abandoned;         /* ... gave up on */
} solver2_archive_result;

/* An archive over `problem`'s grid, from `seed_count` gene vectors in
 * `seeds` (seed k at seeds[k * dimension ..]; none: the problem's start),
 * guided by `model` (or NULL). Nothing is evaluated yet. */
int solver2_archive_new(solver2_problem* problem, const solver2_archive_options* options, const double* seeds,
                        size_t seed_count, solver2_model* model, solver2_archive** archive);
int solver2_archive_free(solver2_archive* archive);

/* Runs it until `evals` evaluations beyond the seeds are spent in all (0:
 * the options'), or `report` (or NULL) stops it. The first call evaluates
 * the seeds: one of them must be within the limits. Called again it goes
 * on, and a run stopped and gone on with ends as the one that never
 * stopped, bit for bit. A larger `evals` extends a finished run.
 *
 * When the evaluation fails (SOLVER2_CALLBACK) the run and its model are
 * left as they were after the last whole iteration, and the call can be
 * made again. `result` (or NULL) is filled either way. */
int solver2_archive_run(solver2_archive* archive, size_t evals, solver2_archive_report_fn report, void* context,
                        solver2_archive_result* result);

int solver2_archive_get(solver2_archive* archive, solver2_archive_result* result);

/* Every cell at once, cell c at index c: whether it holds a candidate,
 * and that candidate's value, the evaluations spent when it was found, its
 * genes (cells * dimension), outcome, violations (cells * limits) and
 * descriptors (cells * axes). An empty cell's numbers are NaN, its `found`
 * 0. Any pointer may be NULL. */
int solver2_archive_elites(solver2_archive* archive, unsigned char* occupied, double* value, size_t* found,
                           double* genes, solver2_outcome* outcome, double* violations, double* descriptors);

/* Each cell's threshold. */
int solver2_archive_thresholds(solver2_archive* archive, double* thresholds);

/* Per limit, the candidates over it so far. */
int solver2_archive_causes(solver2_archive* archive, size_t* over);

/* The cells to show: the best, then again and again the one farthest from
 * those picked, among the cells at most `within` below the best. At most
 * `count` into `cells`; `picked` gets how many. */
int solver2_archive_diverse(solver2_archive* archive, double within, size_t count, size_t* cells, size_t* picked);

/* The run as bytes, and back: as solver2_stage_save and solver2_stage_load. */
int solver2_archive_save(solver2_archive* archive, unsigned char* buffer, size_t capacity, size_t* size);
int solver2_archive_load(solver2_problem* problem, solver2_model* model, const unsigned char* bytes, size_t size,
                         solver2_archive** archive);

#ifdef __cplusplus
}
#endif

#endif
