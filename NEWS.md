# MicroMoB 0.1.2

  * skip plumber API tests on CRAN due to occasional failures out of developer
  control.

# MicroMoB 0.1.1

  * coerce initial state to `integer` type in `setup_humans_SIP()`
  * add getters and setters for parameters for mosquito RM, aquatic BH, and
  aquatic trace models.

# MicroMoB 0.1.0

  * adds ability to model larval breeding pools separately from patches
  * adds human SIP (Susceptible-Infected-Protected) model
  * adds prevalence sampling with test sens/spec `observe_pfpr()`
  * adds simple bloodmeal computation `compute_bloodmeal_simple()`
  * significantly simplifies `step_mosquitoes.RM_stochastic()`
  * update vignettes with better use of `data.table` and discrete time
  equilibrium calculations

# MicroMoB 0.0.12

  * move `plumber` dependency to "Suggests" to avoid NOTE in checks.
  * fix bug in `step_mosquitoes.RM_deterministic` where mortality applied to
  incubating mosquitoes would not correctly be dependent on patch.
  * use `plumber:::findPort()` in test of plumber API 

# MicroMoB 0.0.11

  * add behavioral state model of adult mosquito dynamics.

# MicroMoB 0.0.10

  * add functions to each model allowing parameters to be read in from specially
  formatted JSON files
  * add `inst/extdata` folder with JSON files for testing config
  * add `jsonlite` package to `Imports` for JSON config capability
  * add `plumber` package to `Imports` to expose a web API
  * add `inst/plumber/` subfolders to store Plumber APIs
  * add `callr` package to `Suggests` for testing API
  * add `httr` package to `Suggests` for testing API
  * add `callr` package to `Suggests` for using CSV serialization to return simulation output from API
  * `p` parameter for RM adult mosquito model is now allowed to be patch and time varying
  * add new vignette "Advanced topics" describing how to extend the package and
  the Plumber web API.
  * add `withr` package to `Suggests` to clean up after testing API.

# MicroMoB 0.0.9

  * add information on return value to all function documentation
  * reset `par` options in vignettes
  * improve diagrams in vignettes

# MicroMoB 0.0.8

  * add SIR model of human infection.

# MicroMoB 0.0.7

  * fix moved URL in README.md as per CRAN request.

# MicroMoB 0.0.6

  * improve all vignettes and documentation for clarity.
  * fix DESCRIPTION according to CRAN checks.

# MicroMoB 0.0.5

  * higher quality graphics (https://github.com/dd-harp/MicroMoB/issues/47)
  * bloodmeal directed wiring diagram (https://github.com/dd-harp/MicroMoB/issues/41)

# MicroMoB 0.0.4

  * faster multinomial draws, using the algorithm presented in ["An asymptotically optimal, online algorithm for weighted random sampling with replacement"](https://arxiv.org/abs/1611.00532)

# MicroMoB 0.0.3

  * add Code of Conduct
  * improvements to vignettes (bloodmeal.Rmd and MicroMoB.Rmd)
  * improvements to function reference page

# MicroMoB 0.0.2

  * initial release of minimum viable product: https://github.com/dd-harp/MicroMoB/issues/28
