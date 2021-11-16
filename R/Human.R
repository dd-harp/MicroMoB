




#' @title Setup a human object
#' @description Setup a human object
#' @param type a character in `c("simple", "strata")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup.human <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(type %in% c("strata"))
  pop <- structure(list(), class = type)
  UseMethod("setup.human", pop)
}

#' @rdname setup.human
#' @method setup.human strata
#' @param H a vector of human population sizes
#' @param J a matrix whose columns assign human strata to patches (rows); the
#' columns must all sum to one.
#' @export
setup.human.strata <- function(type, model, H, J = NULL, ...) {

  stopifnot(length(H) > 0)
  stopifnot(is.finite(H))
  stopifnot(H >= 0)

  if (is.null(J)) {
    J <- diag(length(H))
  }

  p <- nrow(J)
  n <- length(H)

  stopifnot(n == ncol(J))
  stopifnot(p > 0)

  if (is_binary(J)) {

    pop$H <- H
    pop$J <- J

  } else {

    stopifnot(colSums(J) == 1)

    np <- n * p

    # J does not directly assign H
    assignment_indices <- matrix(data = 1:np, nrow = p, ncol = n, byrow = FALSE)

    pop$J <- matrix(0, nrow = p, ncol = np)
    for (j in 1:p) {
      pop$J[j, assignment_indices[j, ]] <- 1
    }
    pop$H <- as.vector(J %*% diag(H))

  }

  model$human <- pop
}




#' @title Setup a time spent model
#' @description Setup a time spent model. The model object must have already
#' been initialized with a human object (see [MicroMoB::setup.human]).
#' @section Null time spent model:
#' The null model will simply assign each person to spend all time at their
#' home residence patch.
#' @param type a character in `c("null", "matrix")`
#' @param model a model object (an [environment])
#' @param ... other arguments to be passed to type methods
#' @export
setup.timespent <- function(type, model, ...) {
  stopifnot(inherits(model, "environment"))
  stopifnot(!is.null(model$human))
  stopifnot(type %in% c("null", "matrix"))
  tisp <- structure(list(), class = type)
  UseMethod("setup.human", tisp)
}


#' @rdname setup.timespent
#' @method setup.timespent null
#' @export
setup.timespent.null <- function(type, model) {

  tisp$theta <- diag(length(model$human$H))
  model$tisp <- tisp
}

#' #' @rdname setup.timespent
#' #' @method setup.timespent matrix
#' #' @export
#' setup.timespent.matrix <- function(type, model) {
#'
#'   tisp$theta <- diag(length(model$human$H))
#'   model$tisp <- tisp
#' }

