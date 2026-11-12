/// The C interface's types as D sees them: include/solver2.h, declaration
/// by declaration, under the same names. The header is the contract and is
/// written by hand; so is this. The binding's tests read the header with
/// ImportC and hold every struct here to its layout and every exported
/// function to its signature (conformance.d), so the two can't drift.
///
/// The library itself needs no C compiler: only the tests read the header.
module solver2c.c;

extern (C):

enum int SOLVER2_ABI = 2;

enum : int
{
    SOLVER2_OK = 0,
    SOLVER2_ERROR = 1,
    SOLVER2_CALLBACK = 2,
}

enum : int
{
    SOLVER2_LINEAR = 0,
    SOLVER2_LOG = 1,
    SOLVER2_INTEGER = 2,
}

enum : int
{
    SOLVER2_HARD = 0,
    SOLVER2_TERMINAL = 1,
    SOLVER2_SOFT = 2,
}

// The handles: pointers to nothing the caller can look into. Their value
// is a number the library gave out (api.d), never an address.
struct solver2_problem;
struct solver2_model;
struct solver2_stage;
struct solver2_archive;

struct solver2_gene
{
    const(char)* name;
    double lo;
    double hi;
    double neutral;
    int scale;
}

struct solver2_limit
{
    const(char)* name;
    int tier;
    double weight;
    const(char)* group;
}

struct solver2_axis
{
    const(char)* name;
    const(double)* edges;
    size_t edge_count;
    const(char)* letters;
    size_t nan_bin;
    double width;
}

struct solver2_ranking
{
    double hard_offset;
    double terminal_offset;
    double bound_penalty;
}

alias solver2_evaluate_fn = int function(void* context, size_t rows, const(double)* genes, int feasible_only,
                                         ubyte* scored, double* objective, double* violations, double* descriptors);

struct solver2_problem_spec
{
    const(solver2_gene)* genes;
    size_t gene_count;
    const(solver2_limit)* limits;
    size_t limit_count;
    const(solver2_axis)* axes;
    size_t axis_count;
    solver2_ranking ranking;
    double min_score;
    double prediction_headroom;
    const(char*)* fingerprint;
    size_t fingerprint_count;
    solver2_evaluate_fn evaluate;
    void* context;
}

struct solver2_outcome
{
    int scored;
    double objective;
    double value;
    double hard;
    double terminal;
    double soft;
}

struct solver2_model_options
{
    size_t members;
    size_t hidden;
    size_t batch;
    double learning_rate;
    double l2;
    size_t threads;
}

struct solver2_stage_options
{
    size_t lambda;
    double sigma;
    size_t budget;
    ulong seed;
    size_t stall_generations;
    double stall_tolerance;
    size_t surrogate_pool;
    double surrogate_share;
    size_t surrogate_from;
    size_t surrogate_epochs;
}

struct solver2_stage_generation
{
    size_t generation;
    size_t evaluations;
    double best;
    double population_best;
    double sigma;
    size_t infeasible;
    size_t sampled;
    size_t predicted_feasible;
}

alias solver2_stage_report_fn = int function(void* context, const(solver2_stage_generation)* generation);

struct solver2_stage_result
{
    size_t generations;
    size_t evaluations;
    int finished;
    int stalled;
    int stopped;
    int feasible;
    double best_value;
}

struct solver2_rules
{
    size_t idle_generations;
    size_t infeasible_generations;
    double infeasible_share;
    double min_sigma;
    size_t stale_generations;
    double min_progress;
    size_t block_iterations;
    size_t aim_samples;
}

struct solver2_archive_options
{
    size_t evals;
    size_t emitters;
    size_t lambda;
    double sigma;
    double alpha;
    ulong seed;
    size_t explorers;
    double explorer_sigma;
    size_t surrogate_pool;
    double surrogate_share;
    size_t surrogate_from;
    size_t surrogate_epochs;
    int surrogate_improvers;
    int surrogate_aim;
    solver2_rules rules;
}

struct solver2_archive_iteration
{
    size_t iteration;
    size_t evaluations;
    size_t coverage;
    double qd_score;
    double best_score;
    size_t improved;
    size_t new_cells;
    size_t restarts;
}

alias solver2_archive_report_fn = int function(void* context, const(solver2_archive_iteration)* iteration);

struct solver2_archive_result
{
    size_t cells;
    size_t coverage;
    size_t evaluations;
    size_t seed_evaluations;
    size_t iterations;
    size_t lambda;
    int started;
    int finished;
    int stopped;
    double qd_score;
    double best_score;
    size_t best_cell;
    size_t restarts_idle;
    size_t restarts_infeasible;
    size_t restarts_converged;
    size_t reached;
    size_t abandoned;
}
