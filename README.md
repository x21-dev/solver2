# solver2

A deterministic black-box optimiser written in D: CMA-ES stages and a
quality-diversity archive (CMA-MAE, with explorers that head for empty
cells) over a problem interface, with a learned surrogate model. It is
the next version of [solver](https://github.com/x21-dev/solver), the
optimiser extracted from [pleaser](https://github.com/x21-dev/pleaser).

A problem can be written in three ways, and all three run the same
searches:

| Route | The problem is | Worked example |
|---|---|---|
| Compile time | A D struct; the searches are templates over it, so calls inline | [`docs/example`](docs/example/source/app.d) |
| The class | A subclass of `DynamicProblem`, chosen or loaded at run time, run through plain functions (`solver2.dynamic`) | [`docs/example-class`](docs/example-class/source/app.d) |
| Python, or C | A Python class, or data and one function, over the shared library in [`bindings/c`](bindings/c/README.md) | [`bindings/python/examples/beam.py`](bindings/python/examples/beam.py) |

The [developer manual](docs/manual.md) is the documentation: the model,
each route, saving and resuming, determinism, performance, what the
library doesn't do, and how to move a project from solver. solver itself
stays as it is, and the projects that use it need not change.

Runs are bit-identical on any thread count of one platform, not across
platforms (the C maths library differs), so keep golden runs per
platform; see the manual's Determinism section.

## Using it

- **From D**: solver2 is a dub source library with no dependencies, built
  with LDC (developed with 1.42.0). A dependent compiles its modules with
  its own. See the manual's "Adding solver2 to a project".
- **From Python**: unpack a
  [release](https://github.com/x21-dev/solver2/releases) (Linux and Windows, x86-64; no D
  compiler needed), or build the shared library in `bindings/c`; put the
  folder on Python's path and `import solver2`. See
  [`bindings/python`](bindings/python/README.md) and the manual's
  "Interfacing with Python".
- **From C or another language**: the same library, with its header, from
  a release or from `bindings/c`.

solver2 is proprietary: using it needs a licence granted individually.
See [`LICENSE.md`](LICENSE.md).

## Testing

```
dub test --compiler=ldc2
```

runs the core's tests. The examples, the C interface and the Python
module have their own, listed in [`CLAUDE.md`](CLAUDE.md); CI runs all of
them on Linux and Windows for every pull request.
