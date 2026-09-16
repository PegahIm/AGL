# Quantum and Bayesian Models of Artificial Grammar Learning

## Overview

This repository contains all materials required to reproduce the analyses and
model comparisons reported in the manuscript. Across four experiments, we
compare a quantum probability model and a Bayesian model of artificial grammar
learning (AGL), assessing which account provides the better fit to human
grammaticality judgements under different learning conditions.

Statistical analyses of the behavioural data were carried out in R; model
specification, parameter estimation, model comparison, and model recovery were
carried out in Mathematica.

## Repository structure

```
.
├── Data/           Cleaned and classified data for Experiments 1–4 (.xlsx)
├── Figures/        All figures reported in the manuscript
├── Mathematica/    Model fitting, model comparison, and model recovery
└── R/              Experiment-wise statistical analyses (.Rmd)
```

### `Data/`

These files are the common input to both the R analyses and
the Mathematica model-fitting notebooks; no other data source is required.

### `Mathematica/`

One notebook per experiment, each implementing the Bayesian and quantum models,
estimating their free parameters, and comparing fit via {G², BIC}:

| File | Experiment |
| --- | --- |
| `{AGL-Exp1-Models}.nb` | Experiment 1 |
| `{AGL-Exp2-Models}.nb` | Experiment 2 |
| `{AGL-Exp3-Models}.nb` | Experiment 3 |
| `{AGL-Exp4-Models}.nb` | Experiment 4 |
| `{AGL-Exp2-ModelRecoverySimulation}.wl` | Model recovery simulation |

The four notebooks share a common structure and differ only in the number of
participants and, for Experiment 1, in the set of test items. The `.wl` script
implements the model recovery simulation.

### `R/`

One R Markdown file per experiment, containing the behavioural analyses
reported in the manuscript.

### `Figures/`

All figures appearing in the manuscript.

## Contact

Pegah Imannezhad — pegah.imannezhad@city.ac.uk
