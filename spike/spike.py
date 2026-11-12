"""Step 1 of docs/work-plan.md, a spike: Python calls D.

Loads the spike's shared library with ctypes, gives it a problem whose
evaluation is a Python function, and prints what spike-native prints with
the evaluation in D, line for line. Then the checks that have no native
counterpart: errors on both sides of the boundary, a second thread, and
how often D's collector ran.

    python spike.py            the runs, to compare with spike-native
    python spike.py --checks   the runs (not printed), then the checks
"""
import ctypes
import os
import struct
import sys
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
NAME = "spike.dll" if os.name == "nt" else "libspike.so"
lib = ctypes.CDLL(os.path.join(HERE, NAME))


class SpikeResult(ctypes.Structure):
    _fields_ = [("hash", ctypes.c_uint64), ("generations", ctypes.c_uint64), ("spent", ctypes.c_uint64),
                ("best", ctypes.c_double), ("learned", ctypes.c_uint64)]


EVAL = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_size_t,
                        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double))

lib.spike_init.restype = ctypes.c_int
lib.spike_term.restype = ctypes.c_int
lib.spike_error.restype = ctypes.c_char_p
lib.spike_attach.restype = ctypes.c_int
lib.spike_detach.restype = ctypes.c_int
lib.spike_collections.restype = ctypes.c_uint64
lib.spike_throw.restype = ctypes.c_int
lib.spike_probe2.restype = ctypes.c_int
lib.spike_known.restype = ctypes.c_int
lib.spike_threads.restype = ctypes.c_int
lib.spike_stage.restype = ctypes.c_int
lib.spike_stage.argtypes = [EVAL, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_uint64, ctypes.c_size_t, ctypes.c_int,
                            ctypes.c_size_t, ctypes.POINTER(SpikeResult), ctypes.POINTER(ctypes.c_double)]

DIM = 6
calls = 0
raised = None          # an exception in the callback, kept for the caller


def clamp01(g):
    return 0.0 if g < 0 else (1.0 if g > 1 else g)


def sphere(ctx, n, dim, genes, objective, over):
    """native.d's sphere, in the same operations in the same order."""
    global calls, raised
    try:
        calls += 1
        if fail_at and calls == fail_at:
            raise ValueError("the problem failed on purpose")
        for k in range(n):
            s = 0.0
            for j in range(dim):
                d = clamp01(genes[k * dim + j]) - 0.3
                s -= d * d
            objective[k] = s
            x0 = clamp01(genes[k * dim])
            over[k] = x0 - 0.5 if x0 > 0.5 else 0.0
        return 0
    except BaseException as e:      # nothing may escape into D
        raised = e
        return 1


fail_at = 0
callback = EVAL(sphere)


def bits(v):
    return struct.unpack("<Q", struct.pack("<d", v))[0]


def stage(seed, budget, with_model, threads):
    """One stage: its report line."""
    global raised
    r = SpikeResult()
    genes = (ctypes.c_double * DIM)()
    raised = None
    if lib.spike_stage(callback, None, DIM, seed, budget, with_model, threads, ctypes.byref(r), genes) != 0:
        if raised is not None:
            raise raised
        raise RuntimeError(lib.spike_error().decode())
    return "%s seed %d: hash %016x generations %d spent %d learned %d best %016x genes %s" % (
        "model" if with_model else "plain", seed, r.hash, r.generations, r.spent, r.learned, bits(r.best),
        " ".join("%016x" % bits(g) for g in genes))


def runs():
    lines = [stage(seed, 4000, 0, 0) for seed in range(1, 51)]
    lines += [stage(seed, 800, 1, 3) for seed in range(1, 6)]
    return lines


def checks(lines):
    global fail_at, calls
    print("collections so far:", lib.spike_collections())

    # An error inside D stays inside D.
    assert lib.spike_throw() == 1
    print("a D error, caught at the boundary:", lib.spike_error().decode())

    # An exception in the callback stops the run and comes back out.
    calls, fail_at = 0, 7
    try:
        stage(1, 4000, 0, 0)
        raise SystemExit("the callback's exception was lost")
    except ValueError as e:
        print("a Python exception, raised again after %d callbacks: %s" % (calls, e))
    fail_at = 0
    # ... and the library goes on working: the same stage, the same line.
    assert stage(1, 4000, 0, 0) == lines[0]
    print("the stage after the failed one: as before")

    # A second thread, made known to the D runtime first: the same lines.
    out = {}

    def worker():
        try:
            assert lib.spike_attach() == 0, lib.spike_error().decode()
            out["plain"] = [stage(seed, 4000, 0, 0) for seed in range(1, 11)]
            out["model"] = stage(1, 800, 1, 3)
            assert lib.spike_detach() == 0, lib.spike_error().decode()
        except BaseException as e:
            out["error"] = e

    t = threading.Thread(target=worker)
    t.start()
    t.join()
    if "error" in out:
        raise out["error"]
    assert out["plain"] == lines[:10] and out["model"] == lines[50]
    print("a second thread, attached: the same 11 lines")

    # The main thread again, after the other came and went.
    assert stage(2, 4000, 0, 0) == lines[1]
    print("collections at the end:", lib.spike_collections())


def sphere_numpy(ctx, n, dim, genes, objective, over):
    """The sphere on numpy views of D's arrays: nothing is copied."""
    import numpy as np
    x = np.clip(np.ctypeslib.as_array(genes, shape=(n, dim)), 0.0, 1.0)
    d = x - 0.3
    np.ctypeslib.as_array(objective, shape=(n,))[:] = -(d * d).sum(axis=1)
    np.ctypeslib.as_array(over, shape=(n,))[:] = np.where(x[:, 0] > 0.5, x[:, 0] - 0.5, 0.0)
    return 0


def unattached():
    """Stages from a thread the D runtime was never told of."""
    out = {}

    def worker():
        out["lines"] = [stage(seed, 4000, 0, 0) for seed in range(1, 51)] + [stage(1, 800, 1, 3)]

    t = threading.Thread(target=worker)
    t.start()
    t.join()
    return out["lines"]


def stress(lines):
    """Six times the runs, while another Python thread makes garbage."""
    import gc
    stop = False
    made = [0]

    def churn():
        while not stop:
            made[0] += len([[i, str(i)] for i in range(2000)])
            gc.collect()

    import time
    t = threading.Thread(target=churn)
    t.start()
    t0 = time.perf_counter()
    try:
        again = [stage(seed, 4000, 0, 0) for seed in range(1, 6)] + [stage(1, 800, 1, 3)]
    finally:
        stop = True
        t.join()
    assert again == lines[:5] + [lines[50]]
    print("6 stages (2600 callbacks) beside a busy Python thread (%d objects made): the same lines, in %.1f s; "
          "%d collections in D" % (made[0], time.perf_counter() - t0, lib.spike_collections()))


def probe():
    """What D's runtime knows of a thread, and whether its collector works
    from it: the main thread, a second one unattached, and attached."""
    def one(tag):
        held, dropped = ctypes.c_double(), ctypes.c_double()
        known, threads = lib.spike_known(), lib.spike_threads()
        intact = lib.spike_probe2(ctypes.byref(held), ctypes.byref(dropped))
        print("%-20s known %d of %d threads; 64 MB held by a local: intact %d, %.0f MB in use after a collection"
              % (tag, known, threads, intact, held.value))

    one("main thread")

    def worker(attach):
        if attach:
            lib.spike_attach()
        one("second, attached" if attach else "second, unattached")
        if attach:
            lib.spike_detach()

    for attach in (False, True, False):
        t = threading.Thread(target=worker, args=(attach,))
        t.start()
        t.join()


if __name__ == "__main__":
    assert lib.spike_init() == 0
    if "--fault" in sys.argv:
        kind, index = int(sys.argv[2]), int(sys.argv[3])
        print("fault", kind, "returned", lib.spike_fault(kind, index), lib.spike_error().decode(), flush=True)
        sys.exit(0)
    if "--numpy" in sys.argv:
        callback = EVAL(sphere_numpy)
        print("\n".join(runs()))
        sys.exit(0)
    if "--probe" in sys.argv:
        probe()
        sys.exit(0)
    if "--unattached" in sys.argv:
        print("\n".join(unattached()[:55]))
        sys.exit(0)
    lines = runs()
    if "--checks" in sys.argv:
        checks(lines)
    elif "--stress" in sys.argv:
        stress(lines)
    else:
        print("\n".join(lines))
    assert lib.spike_term() == 0
