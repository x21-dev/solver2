# solver2

The next version of [solver](https://github.com/x21-dev/solver), the
deterministic black-box optimiser extracted from
[pleaser](https://github.com/x21-dev/pleaser): CMA-ES stages and a
quality-diversity archive (CMA-MAE) over a problem interface, with a
learned surrogate model.

solver2 is being built to add two things:

- a problem as a class that need not be known at compile time;
- a C interface and a Python module over it, so a problem written in
  Python runs through the D searches.

[`docs/work-plan.md`](docs/work-plan.md) has the design, the steps and
where the work stands. solver2 started as solver's code at commit
`ae68dbb`, with the package and its modules renamed: `solver2.problem`,
`solver2.stage`, `solver2.archive`, `solver2.model`, `solver2.ensemble`,
`solver2.cmaes`, `solver2.bytes`, `solver2.rng`, `solver2.detmath`. So
far it adds the first of the two: `DynamicProblem`, a problem as a class
that gets what a compile-time problem gets, and `solver2.dynamic`, the
stage and the archive over it as plain functions (`runStage`,
`runArchive`). The C interface and the Python module are next.

solver itself stays as it is, and the projects that use it need not
change. The manual will say how to move a project to solver2.

Runs are bit-identical on any thread count of one platform, not across
platforms (the C maths library differs), so keep golden runs per
platform; see the manual's Determinism section.

solver2 is a source library: a dependent compiles its modules with its
own, so calls into it inline as calls within one package do.

```
dub test --compiler=ldc2
```

To use it, see the [developer manual](docs/manual.md) and its worked
example in `docs/example`. The manual is still solver's, with the module
names changed.
