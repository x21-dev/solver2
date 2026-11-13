# solver2 {version}

A deterministic black-box optimiser: CMA-ES stages and a
quality-diversity archive (CMA-MAE, with explorers that head for empty
cells), with a learned surrogate model. The searches are written in D and
compiled into the shared library in this folder; the problem is yours,
written in Python or behind a C function.

This is a built release: nothing here needs compiling, and no D compiler
is needed to use it. The source, and the developer manual, are at
<https://github.com/x21-dev/solver2>; the manual for this release is
<https://github.com/x21-dev/solver2/blob/v{version}/docs/manual.md>.

## What is here

| | |
|---|---|
| `solver2/` | The Python module, with the library inside it: {beside} |
| `include/solver2.h` | The C interface, version {abi}, and its documentation |
| `examples/beam.py` | The worked example in Python |
| `LICENSE.md` | The licence |

## What it needs

- **Linux**: x86-64, glibc 2.35 or newer (Ubuntu 22.04 or newer; the
  library is built on Ubuntu 22.04).
- **Windows**: x86-64.
- **From Python**: a 64-bit Python 3.10 or newer, and numpy.

## From Python

Put this folder (the one that holds `solver2/`) on Python's path, with
`PYTHONPATH` or `sys.path`, and `import solver2`. To see that it works:

```
python examples/beam.py
```

It runs for some seconds and ends with a line that starts `bits:`.

Keep the files of `solver2/` together: the module looks for the library
beside itself, and on Linux the library looks there for the two runtime
libraries that came with it.

The manual's "Interfacing with Python" is the documentation, and the
docstrings of `solver2/__init__.py` the reference.

## From C, or another language

Include `include/solver2.h` and link against `solver2/{library}` (on
Windows through `solver2/solver2.lib`). At run time the library must be
where the loader finds it: beside the program on Windows; on Linux by a
run path, or `LD_LIBRARY_PATH`, that names the folder. On Linux keep the
two runtime libraries beside it.

The header's opening comment has the rules, and the manual's "The C
interface, for other languages" a small program.

## Determinism

The same problem, options and seed give the same bits on any thread count
of one platform. Runs differ between platforms, and may differ between
versions of the C library on Linux: keep reference results per platform.
