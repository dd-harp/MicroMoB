# Micro-MoB (Microsimulation for mosquito-borne pathogens)

<!-- badges: start -->
[![R-CMD-check](https://github.com/dd-harp/MicroMoB/workflows/R-CMD-check/badge.svg)](https://github.com/dd-harp/MicroMoB/actions)
[![codecov](https://codecov.io/gh/dd-harp/MicroMoB/branch/main/graph/badge.svg?)](https://codecov.io/gh/dd-harp/MicroMoB)
<!-- badges: end -->

## What is Micro-MoB?

**Micro-MoB** was made to simplify the task of model building for mosquito-borne pathogen transmission (MBPT) systems. 
It stands for "microsimulation for mosquito-borne pathogens". It is a modular
framework to build discrete time MBPT models. It uses R's [S3 object system](http://adv-r.had.co.nz/S3.html)
to define a set of _components_ and _interfaces_ which can be filled by any specific _model_ that
implements the interface. These parts, along with certain _invariants_ can be put
together to define a full simulation model. Definitions for all these terms can be found
in the documentation.

We hope it proves useful. Please visit the [website](https://dd-harp.github.io/MicroMoB/) to learn more.

## Installation

```
remotes::install_github('dd-harp/MicroMoB')
library(MicroMoB)
```

## Documentation

To start learning more about the software design, the problems it was designed to solve,
and how to build new models in **Micro-MoB**, please read `vignette("MicroMoB")`. 

Next, `vignette("bloodmeal")` describes how the bloodmeal algorithm computes
the distribution of bites using each component's interface, allowing different models
to be linked in a consistent framework.

We also have articles describing some well-known models of specific components
of MBPT models that are implemented in **Micro-MoB**:

  * `vignette("RM_mosquito")`: read about an implementation of a flexible [Ross-Macdonald
  model](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1002588) of adult mosquito dynamics.
  * `vignette("MOI_human")`: read about an implementation of the $(M/M/\infty)$ queuing
  model for superinfection in humans.
  * `vignette("BH_aqua")`: read about a simple non-linear model of aquatic (immature)
  stage mosquito dynamics based on the well known [Beverton-Holt model](https://en.wikipedia.org/wiki/Beverton-Holt_model) from ecology.
  * `vignette("RM_transmission")`: read about how we put together models fulfilling
  each component to run a simple Ross-Macdonald style model of pathogen transmission
  between human hosts and mosquito vectors.

## Contributing

Thank you for your interest in **Micro-MoB**! If you have a bug to report, please
open an [issue on GitHub](https://github.com/dd-harp/MicroMoB/issues). If you would like
to open a pull request or have further questions, please see our guide to
contributing to the project at `vignette("Contributing")`.

## Code of Conduct
  
Please note that the **Micro-MoB** project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
