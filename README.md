# MicroMoB

<!-- badges: start -->
[![R-CMD-check](https://github.com/dd-harp/MicroMoB/workflows/R-CMD-check/badge.svg)](https://github.com/dd-harp/MicroMoB/actions)
[![codecov](https://codecov.io/gh/dd-harp/MicroMoB/branch/main/graph/badge.svg?)](https://codecov.io/gh/dd-harp/MicroMoB)
<!-- badges: end -->

## Software design

The model object will be an environment with a class attribute.
All interior structure will be named lists. Each function will be passed the entire
model object, and will dispatch on classes of objects within the model object.

## Components

The model is broken into components, for humans, immature and adult mosquitoes (and some others).
Each component has an _interface_, which are methods which must be defined for that
component. A component's interface is stored in file, for example, R/humans_interface.R
shows the user what methods must be defined for any human model. Other components (e.g. the bloodmeal)
will call generic methods not knowing what specific code is implementing them, and so
they must return values consistent with their definition.

We call a specific implementation of a component a _model_.
Specific implementations are found in files that replace _\_interface_ with the
model name, for example R/humans_SIR.R. Their accompanying test files are located in
tests/testthat. If you are creating a new model, please remember to test it
adequately.

We list the components which require interfaces below and specific models
to implement them.

### Mosquitoes

The mosquito component is responsible for all dynamics which update adult mosquito
populations. The interface is defined in R/mosquito_interface.R.

### Aquatic

The aquatic component is responsible for all dynamics which update immature (aquatic
stage) mosquito populations. The interface is defined in R/aquatic_interface.R.

### Humans

The human component updates human populations. The interface is defined in R/humans_interface.R.

### Visitor

The human component updates human populations outside of the resident population of the geographic area being simulated. 
The interface is defined in R/visitor_interface.R.

### Alternative blood hosts

The alternative blood host component is responsible for other blood hosts for 
mosquitoes (livestock, dogs, etc). The interface is defined in R/althost_interface.R.

## Update

To update the model, a function is called which gathers information from the various
components to calculate rates which couple components (i.e. infection) together, which
are then passed back to each individual component, and updated using the generic
interface for each component.
