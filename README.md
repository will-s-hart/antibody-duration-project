# Duration of antibody protection project

This repository supports the Duration of Antibody Protection project at the
Nagoya-Oxford Workshop 2026.

A vaccine trial measures an antibody marker repeatedly after vaccination, and
the marker falls quickly at first and then more slowly. At each planned trial
review we want to estimate when it will fall below a level associated with
protection — often beyond the period observed so far — and to say how much
would be gained from further follow-up. This repository holds the simulation
study that answers that question.

[notes/handout.pdf](notes/handout.pdf) is the participant handout and the
methodological source of truth: the question, the estimand, the candidate
models and their starting values.

Supporting material is collected under `notes/`:

- [Dry-run summary](notes/dry-run/summary.pdf): a synthesis of the earlier
  feasibility analyses, using their original calibration.
- [Analysis figure schematic](notes/figures/analysis-schematic-v2.png): the
  latest illustrative figure plan. The
  [first draft](notes/figures/analysis-schematic-v1.png) is retained alongside
  it. These schematics illustrate intended figures, not computed results.

## Getting started

Requires R 4.3 or later. If you do not have R, [rig](https://github.com/r-lib/rig)
is the least painful way to install and switch between versions:

```sh
brew install --cask rig && rig add release   # macOS; see rig's README for Linux and Windows
```

Then, from the repository root:

```r
renv::restore()        # install the pinned package versions (first time only)
pkgload::load_all()    # load the project
testthat::test_local() # run the tests
```

`renv` bootstraps itself from `.Rprofile` when you open the project, so
`renv::restore()` works in a fresh clone with nothing else installed.

The cheap end-to-end check, which is also what to run before opening a pull
request:

```sh
Rscript scripts/00-integration-gate.R
```

## Single-participant prior pilot

The first three steps of handout Section 6 can be illustrated with one
simulated participant, monthly measurements and reviews through 24 months:

```sh
Rscript scripts/01-prior-pilot.R
Rscript scripts/02-prior-pilot-report.R
```

Install the JAGS command-line sampler separately, then use `renv::restore()`
for the R dependencies. The pilot calls JAGS directly; it does not need the
`rjags` connector. Four independent chains fit all four families under both
Section 5 prior sets. Failed convergence checks trigger longer sampling;
fits, attempts and errors are retained in `cache/prior-pilot/` and
`results/prior-pilot/`. Resume reuses fits only when data, configuration and
fitting implementation match. Add `--refit` to repeat all fits.

Open `results/prior-pilot/report.html` for the self-contained report and
figures. Both threshold durations come from each fitted trajectory, and all
reviews use the same noisy history. This is an individual-duration
illustration, not a repeated-simulation coverage or adequacy result. The
hierarchical, least-squares and general sweep workstreams remain separate.

Repeat with the common day-90/day-365 anchor priors using:

```sh
Rscript scripts/01-prior-pilot.R --two-anchor
Rscript scripts/02-prior-pilot-report.R --two-anchor
```

This uses the same simulated history and writes to
`results/prior-pilot-two-anchor/` and `cache/prior-pilot-two-anchor/`.
The report includes a comparison with the original pilot when its results
are available. The exponential retains its original early-decline prior;
the other three families share the joint prior on initial, day-90 and
day-365 latent titres. See the handout for the remaining biphasic parameter
and the explicit relaxation of the initial component-weight restriction.

## Layout

| Path | Contents |
| --- | --- |
| [notes/](notes/) | Handout, bibliography, dry-run summary and figure plans |
| [R/](R/) | Package source: shared schemas, models, one file per workstream |
| [inst/config/default.yaml](inst/config/default.yaml) | Calibration, design grid and adequacy rules |
| [scripts/](scripts/) | Entry points, starting with the integration gate |
| [tests/testthat/](tests/testthat/) | Test suite |
| [INTERFACES.md](INTERFACES.md) | The shared contract — read this before writing code |
| [AGENTS.md](AGENTS.md) | Conventions, for coding agents and for people |

## Workstreams

Work is split four ways so that it can proceed in parallel. Each workstream
owns its own files and builds against the shared schemas; nobody waits for
anybody, because `mock_sim_dataset()` and `mock_fit_result()` stand in for
work that has not landed yet.

| | Workstream | File |
| --- | --- | --- |
| (a) | Coordination | `R/driver.R`, `R/score.R`, `R/plots.R`, `scripts/` |
| (b) | Simulation, inference, MCMC | `R/simulate.R`, `R/fit-ls.R`, `R/fit-mcmc.R` |
| (c) | Mechanistic models | `R/calibrate.R` |
| (d) | Nonlinear mixed-effects models | `R/fit-nlme.R` |

The contract between them — the four schemas, the units convention and the
rule that a trajectory which never crosses the threshold is `Inf` — is in
[INTERFACES.md](INTERFACES.md).

## Contributing

Fork the [main repository](https://github.com/will-s-hart/antibody-duration-project),
create a branch in your fork, and open a pull request from that branch against
`main` in the main repository.

Please run `testthat::test_local()` and the integration gate before opening a
pull request, and see [AGENTS.md](AGENTS.md) for the conventions the project
relies on.
