
#' @title Check a model object
#' @description Check that a model object has minimum objects required to execute
#' a time step.
#' @param Mod a model object ([list] or [environment])
#' @export
checkModObject <- function(Mod) {
  stopifnot(inherits(Mod, "environment"))

  stopifnot(is.finite(Mod$nMosquito))
  stopifnot(Mod$nMosquito > 0)

  stopifnot(is.finite(Mod$nHost))
  stopifnot(Mod$nHost > 0)

  stopifnot(is.finite(Mod$nPathogen))
  stopifnot(Mod$nPathogen > 0)

  stopifnot(!is.null(Mod$pathogen))
  stopifnot(inherits(Mod$pathogen, "list"))
  stopifnot(length(Mod$pathogen) == Mod$nPathogen)

  # maybe can remove this for specific checkers that run b4 1st iteration of
  # sim loop
  stopifnot(!is.null(Mod$mosquito))
  stopifnot(inherits(Mod$mosquito, "list"))
  stopifnot(length(Mod$mosquito) == Mod$nMosquito)

  stopifnot(!is.null(Mod$NI))
  stopifnot(inherits(Mod$NI, "matrix"))
  stopifnot(dim(Mod$NI) == c(Mod$nHost, Mod$nPathogen))
}



# checkMosquitoObject <- function(obj) {
#   stopifnot()
# }
