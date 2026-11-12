"""solver2's C interface as ctypes sees it: bindings/c/include/solver2.h,
declaration by declaration, and the shared library, found and loaded.

Written by hand, as the header is. tests/test_conformance.py reads the
header and holds every struct, function and constant here to it, in both
directions, so the two can't drift: a struct that didn't match would
corrupt memory without a word.
"""
import ctypes as C
import os
import sys

#: The interface version these declarations are of (the header's SOLVER2_ABI).
ABI = 2

OK, ERROR, CALLBACK = 0, 1, 2
LINEAR, LOG, INTEGER = 0, 1, 2
HARD, TERMINAL, SOFT = 0, 1, 2

_double_p = C.POINTER(C.c_double)
_size_p = C.POINTER(C.c_size_t)
_ubyte_p = C.POINTER(C.c_ubyte)
_handle = C.c_void_p            # a handle: a number, never looked into
_handle_p = C.POINTER(C.c_void_p)


class Gene(C.Structure):
    _fields_ = [("name", C.c_char_p), ("lo", C.c_double), ("hi", C.c_double), ("neutral", C.c_double),
                ("scale", C.c_int)]


class Limit(C.Structure):
    _fields_ = [("name", C.c_char_p), ("tier", C.c_int), ("weight", C.c_double), ("group", C.c_char_p)]


class Axis(C.Structure):
    _fields_ = [("name", C.c_char_p), ("edges", _double_p), ("edge_count", C.c_size_t), ("letters", C.c_char_p),
                ("nan_bin", C.c_size_t), ("width", C.c_double)]


class Ranking(C.Structure):
    _fields_ = [("hard_offset", C.c_double), ("terminal_offset", C.c_double), ("bound_penalty", C.c_double)]


EvaluateFn = C.CFUNCTYPE(C.c_int, C.c_void_p, C.c_size_t, _double_p, C.c_int, _ubyte_p, _double_p, _double_p,
                         _double_p)


class ProblemSpec(C.Structure):
    _fields_ = [("genes", C.POINTER(Gene)), ("gene_count", C.c_size_t),
                ("limits", C.POINTER(Limit)), ("limit_count", C.c_size_t),
                ("axes", C.POINTER(Axis)), ("axis_count", C.c_size_t),
                ("ranking", Ranking), ("min_score", C.c_double), ("prediction_headroom", C.c_double),
                ("fingerprint", C.POINTER(C.c_char_p)), ("fingerprint_count", C.c_size_t),
                ("evaluate", EvaluateFn), ("context", C.c_void_p)]


class Outcome(C.Structure):
    _fields_ = [("scored", C.c_int), ("objective", C.c_double), ("value", C.c_double), ("hard", C.c_double),
                ("terminal", C.c_double), ("soft", C.c_double)]


class ModelOptions(C.Structure):
    _fields_ = [("members", C.c_size_t), ("hidden", C.c_size_t), ("batch", C.c_size_t),
                ("learning_rate", C.c_double), ("l2", C.c_double), ("threads", C.c_size_t)]


class StageOptions(C.Structure):
    _fields_ = [("lambda", C.c_size_t), ("sigma", C.c_double), ("budget", C.c_size_t), ("seed", C.c_uint64),
                ("stall_generations", C.c_size_t), ("stall_tolerance", C.c_double),
                ("surrogate_pool", C.c_size_t), ("surrogate_share", C.c_double), ("surrogate_from", C.c_size_t),
                ("surrogate_epochs", C.c_size_t)]


class StageGeneration(C.Structure):
    _fields_ = [("generation", C.c_size_t), ("evaluations", C.c_size_t), ("best", C.c_double),
                ("population_best", C.c_double), ("sigma", C.c_double), ("infeasible", C.c_size_t),
                ("sampled", C.c_size_t), ("predicted_feasible", C.c_size_t)]


StageReportFn = C.CFUNCTYPE(C.c_int, C.c_void_p, C.POINTER(StageGeneration))


class StageResult(C.Structure):
    _fields_ = [("generations", C.c_size_t), ("evaluations", C.c_size_t), ("finished", C.c_int),
                ("stalled", C.c_int), ("stopped", C.c_int), ("feasible", C.c_int), ("best_value", C.c_double)]


class Rules(C.Structure):
    _fields_ = [("idle_generations", C.c_size_t), ("infeasible_generations", C.c_size_t),
                ("infeasible_share", C.c_double), ("min_sigma", C.c_double), ("stale_generations", C.c_size_t),
                ("min_progress", C.c_double), ("block_iterations", C.c_size_t), ("aim_samples", C.c_size_t)]


class ArchiveOptions(C.Structure):
    _fields_ = [("evals", C.c_size_t), ("emitters", C.c_size_t), ("lambda", C.c_size_t), ("sigma", C.c_double),
                ("alpha", C.c_double), ("seed", C.c_uint64), ("explorers", C.c_size_t),
                ("explorer_sigma", C.c_double), ("surrogate_pool", C.c_size_t), ("surrogate_share", C.c_double),
                ("surrogate_from", C.c_size_t), ("surrogate_epochs", C.c_size_t), ("surrogate_improvers", C.c_int),
                ("surrogate_aim", C.c_int), ("rules", Rules)]


class ArchiveIteration(C.Structure):
    _fields_ = [("iteration", C.c_size_t), ("evaluations", C.c_size_t), ("coverage", C.c_size_t),
                ("qd_score", C.c_double), ("best_score", C.c_double), ("improved", C.c_size_t),
                ("new_cells", C.c_size_t), ("restarts", C.c_size_t)]


ArchiveReportFn = C.CFUNCTYPE(C.c_int, C.c_void_p, C.POINTER(ArchiveIteration))


class ArchiveResult(C.Structure):
    _fields_ = [("cells", C.c_size_t), ("coverage", C.c_size_t), ("evaluations", C.c_size_t),
                ("seed_evaluations", C.c_size_t), ("iterations", C.c_size_t), ("lambda", C.c_size_t),
                ("started", C.c_int), ("finished", C.c_int), ("stopped", C.c_int), ("qd_score", C.c_double),
                ("best_score", C.c_double), ("best_cell", C.c_size_t), ("restarts_idle", C.c_size_t),
                ("restarts_infeasible", C.c_size_t), ("restarts_converged", C.c_size_t), ("reached", C.c_size_t),
                ("abandoned", C.c_size_t)]


#: The header's names for what is above, for the conformance test.
STRUCTS = {
    "solver2_gene": Gene, "solver2_limit": Limit, "solver2_axis": Axis, "solver2_ranking": Ranking,
    "solver2_problem_spec": ProblemSpec, "solver2_outcome": Outcome, "solver2_model_options": ModelOptions,
    "solver2_stage_options": StageOptions, "solver2_stage_generation": StageGeneration,
    "solver2_stage_result": StageResult, "solver2_rules": Rules, "solver2_archive_options": ArchiveOptions,
    "solver2_archive_iteration": ArchiveIteration, "solver2_archive_result": ArchiveResult,
}
CALLBACKS = {"solver2_evaluate_fn": EvaluateFn, "solver2_stage_report_fn": StageReportFn,
             "solver2_archive_report_fn": ArchiveReportFn}
HANDLES = ("solver2_problem", "solver2_model", "solver2_stage", "solver2_archive")
CONSTANTS = {"SOLVER2_ABI": ABI, "SOLVER2_OK": OK, "SOLVER2_ERROR": ERROR, "SOLVER2_CALLBACK": CALLBACK,
             "SOLVER2_LINEAR": LINEAR, "SOLVER2_LOG": LOG, "SOLVER2_INTEGER": INTEGER, "SOLVER2_HARD": HARD,
             "SOLVER2_TERMINAL": TERMINAL, "SOLVER2_SOFT": SOFT}

#: Every function of the header: its result, then its parameters.
PROTOTYPES = {
    "solver2_init": (C.c_int, []),
    "solver2_term": (C.c_int, []),
    "solver2_last_error": (C.c_char_p, []),
    "solver2_abi": (C.c_int, []),
    "solver2_checked": (C.c_int, []),

    "solver2_problem_spec_default": (C.c_int, [C.POINTER(ProblemSpec)]),
    "solver2_problem_new": (C.c_int, [C.POINTER(ProblemSpec), _handle_p]),
    "solver2_problem_free": (C.c_int, [_handle]),
    "solver2_problem_check": (C.c_int, [_handle, C.c_size_t, C.c_uint64]),
    "solver2_problem_dimension": (C.c_int, [_handle, _size_p]),
    "solver2_problem_cells": (C.c_int, [_handle, _size_p]),
    "solver2_problem_start": (C.c_int, [_handle, _double_p]),
    "solver2_problem_physical": (C.c_int, [_handle, C.c_size_t, _double_p, _double_p]),
    "solver2_problem_cell_code": (C.c_int, [_handle, C.c_size_t, C.c_char_p, C.c_size_t]),

    "solver2_model_options_default": (C.c_int, [C.POINTER(ModelOptions)]),
    "solver2_model_new": (C.c_int, [_handle, C.c_uint64, C.POINTER(ModelOptions), _handle_p]),
    "solver2_model_free": (C.c_int, [_handle]),
    "solver2_model_rows": (C.c_int, [_handle, _size_p]),

    "solver2_stage_options_default": (C.c_int, [C.POINTER(StageOptions)]),
    "solver2_stage_new": (C.c_int, [_handle, C.POINTER(StageOptions), _double_p, _handle, _handle_p]),
    "solver2_stage_free": (C.c_int, [_handle]),
    "solver2_stage_run": (C.c_int, [_handle, C.c_size_t, StageReportFn, C.c_void_p, C.POINTER(StageResult)]),
    "solver2_stage_get": (C.c_int, [_handle, C.POINTER(StageResult)]),
    "solver2_stage_best": (C.c_int, [_handle, _double_p, C.POINTER(Outcome), _double_p, _double_p]),
    "solver2_stage_causes": (C.c_int, [_handle, _size_p]),
    "solver2_stage_save": (C.c_int, [_handle, _ubyte_p, C.c_size_t, _size_p]),
    "solver2_stage_load": (C.c_int, [_handle, _handle, _ubyte_p, C.c_size_t, _handle_p]),

    "solver2_archive_options_default": (C.c_int, [C.POINTER(ArchiveOptions)]),
    "solver2_archive_new": (C.c_int, [_handle, C.POINTER(ArchiveOptions), _double_p, C.c_size_t, _handle,
                                      _handle_p]),
    "solver2_archive_free": (C.c_int, [_handle]),
    "solver2_archive_run": (C.c_int, [_handle, C.c_size_t, ArchiveReportFn, C.c_void_p,
                                      C.POINTER(ArchiveResult)]),
    "solver2_archive_get": (C.c_int, [_handle, C.POINTER(ArchiveResult)]),
    "solver2_archive_elites": (C.c_int, [_handle, _ubyte_p, _double_p, _size_p, _double_p, C.POINTER(Outcome),
                                         _double_p, _double_p]),
    "solver2_archive_thresholds": (C.c_int, [_handle, _double_p]),
    "solver2_archive_causes": (C.c_int, [_handle, _size_p]),
    "solver2_archive_diverse": (C.c_int, [_handle, C.c_double, C.c_size_t, _size_p, _size_p]),
    "solver2_archive_save": (C.c_int, [_handle, _ubyte_p, C.c_size_t, _size_p]),
    "solver2_archive_load": (C.c_int, [_handle, _handle, _ubyte_p, C.c_size_t, _handle_p]),
}


def _candidates():
    """Where the library may be, in the order tried."""
    name = "solver2.dll" if os.name == "nt" else "libsolver2.so"
    here = os.path.dirname(os.path.abspath(__file__))
    given = os.environ.get("SOLVER2_LIBRARY")
    if given:
        yield given
        return
    yield os.path.join(here, name)                                  # beside this module
    yield os.path.join(here, os.pardir, os.pardir, "c", name)       # in the repository, where dub builds it


def load():
    """Loads the library and declares its functions. Returns the library and
    the path it came from.

    The library is the one named by the environment variable
    SOLVER2_LIBRARY when that is set; else solver2.dll (libsolver2.so)
    beside this module, or in bindings/c of the repository this module is
    in.
    """
    tried = []
    for path in _candidates():
        path = os.path.normpath(path)
        if os.path.isfile(path):
            break
        tried.append(path)
    else:
        raise ImportError(
            "solver2's shared library was not found (looked for %s). Build it in bindings/c with "
            "`dub build --compiler=ldc2 --build=release`, or name it in the environment variable SOLVER2_LIBRARY."
            % ", ".join(tried))
    lib = C.CDLL(path)
    for name, (result, parameters) in PROTOTYPES.items():
        try:
            function = getattr(lib, name)
        except AttributeError:
            raise ImportError("%s has no %s: it is another version of solver2's library than this module's "
                              "(interface version %d). Rebuild it." % (path, name, ABI)) from None
        function.restype = result
        function.argtypes = parameters
    built = lib.solver2_abi()
    if built != ABI:
        raise ImportError("%s is of interface version %d, this module of %d. Rebuild the library, or use the "
                          "module that came with it." % (path, built, ABI))
    return lib, path


if sys.maxsize <= 2 ** 32:
    raise ImportError("solver2 needs a 64-bit Python")
