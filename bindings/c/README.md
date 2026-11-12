# solver2's C interface

solver2's searches as a shared library (`solver2.dll`, `libsolver2.so`)
that a program in another language loads: a CMA-ES stage and the archive
over a problem whose evaluation is the caller's function, with a learned
model, and a run saved and resumed as bytes. The Python module is built
on it.

[`include/solver2.h`](include/solver2.h) is the interface and its
documentation: read its opening comment for the rules (start the library
first, error codes, handles, buffers, threads), then each function's.

## Building

LDC 1.42.0, in this folder:

```
dub build --compiler=ldc2 --build=release
```

That is the `library` configuration: the D runtime linked into the
library, so it needs no other file beside it, and **bounds checks on**,
so a bug on the D side comes back as an error and a message instead of
corrupting the caller's process. `solver2_checked()` tells which build a
library is.

To build without the checks (the evaluation of a real problem dwarfs
what they cost):

```
dub build --compiler=ldc2 --build=release --config=unchecked
```

Both write the same file, so the last one built is the one in place.

## Testing

```
dub test --compiler=ldc2                        # in this folder
dub build --compiler=ldc2 --build=release       # the library the driver loads
cd test
dub run --compiler=ldc2                         # the driver
```

- **`dub test`** runs the interface's functions in one process with the D
  searches, and holds what comes out of them to what the searches give
  natively, bit for bit: a stage, an archive, a learned model; a run
  stopped, saved, loaded and resumed; a run whose evaluation fails part
  way. It also holds the header to the library: `solver2.h` is read by
  ImportC and compared with the D declarations, struct by struct and
  function by function, and a mismatch fails the build. This needs a C
  compiler for its preprocessor (MSVC on Windows, `cc` elsewhere); the
  library itself doesn't.
- **The driver** (`test/source/driver.d`) is the caller in another
  language: a program built without the D runtime that loads the library
  at run time and uses nothing but the header and the exported
  functions. It checks what only that shows: the library loads and
  exports everything, refuses calls before `solver2_init`, resumes bit
  for bit from its bytes, survives a failed evaluation, and serves
  threads it never heard of, two at once and twelve one after another.
  It takes the library's path as its argument (`../solver2.dll` or
  `../libsolver2.so` without one) and ends with `all passed`.

## What is where

| File | What it holds |
|---|---|
| `include/solver2.h` | The interface, written by hand |
| `source/solver2c/c.d` | The header's types as D sees them, under the same names |
| `source/solver2c/api.d` | The exported functions: threads, handles, errors, copies |
| `source/solver2c/foreign.d` | `ForeignProblem`: a problem from data, its evaluation the caller's function |
| `source/solver2c/runs.d` | A model, a stage and an archive between calls; saving; what a failed evaluation goes back to |
| `source/solver2c/checkpoint.d` | The format of the saved bytes: a copy of the manual's checkpoint template |
| `source/solver2c/conformance.d` | The header held to the library (tests only) |
| `source/solver2c/tests.d` | The tests `dub test` runs |
| `test/` | The driver, a dub package of its own |

Not run on Linux yet: the thread-exit hook's `pthread` half, the driver's
`dlopen` half and the build flags are written as the documentation says
and wait for a Linux machine (the work plan's step 6).
