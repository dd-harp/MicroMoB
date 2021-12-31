#' # a Beverton-Holt style model of aquatic stages
#'
#' #' @title Setup aquatic (immature) mosquito model with Beverton-Holt dynamics
#' #' @description A single compartment for all aquatic stages is modeled which
#' #' suffers density dependent mortality like the [https://en.wikipedia.org/wiki/Beverton%E2%80%93Holt_model](Beverton-Holt model).
#' #' @details All parameters can be passed either as a vector of length equal to `p`, a matrix with `p` rows
#' #' and `tmax` columns, or a matrix with `p` rows and `365` columns.
#' #' @param model an object from [MicroMoB::make_MicroMoB]
#' #' @param molt proportion of immature stages which will mature and emerge as adults each day
#' #' @param surv daily survival probability
#' #' @param K carrying capacity
#' #' @param stochastic should the model update deterministically or stochastically?
#' #' @export
#' setup_aqua_trace <- function(model, molt, surv, K, stochastic) {
#'   stopifnot(inherits(model, "MicroMoB"))
#'
#'   tmax <- model$global$tmax
#'   p <- model$global$p
#'
#'   stopifnot(is.finite(molt))
#'   stopifnot(molt >= 0)
#'   stopifnot(is.finite(surv))
#'   stopifnot(surv >= 0)
#'   stopifnot(is.finite(K))
#'   stopifnot(K >= 0)
#'
#'   if (inherits(molt, "matrix")) {
#'     stopifnot(nrow(molt) == p)
#'     if (ncol(molt) == 365L) {
#'       ix <- (1:tmax) %% 365L
#'       ix[which(ix == 0L)] <- 365L
#'       molt_mat <- molt[, ix, drop = FALSE]
#'     } else if (ncol(molt) == tmax) {
#'       molt_mat <- molt
#'     } else {
#'       stop("incorrect dimensions of molt matrix")
#'     }
#'   } else {
#'     stopifnot(length(molt) == p)
#'     if (p > 1) {
#'       molt_mat <- replicate(n = tmax, expr = molt)
#'     } else {
#'       molt_mat <- matrix(data = molt, nrow = 1, ncol = tmax)
#'     }
#'   }
#'
#'   stopifnot(nrow(molt_mat) == p)
#'   stopifnot(ncol(molt_mat) == tmax)
#'
#'   aqua_class <- c("trace")
#'   if (stochastic) {
#'     aqua_class <- c(aqua_class, "trace_stochastic")
#'   } else {
#'     aqua_class <- c(aqua_class, "trace_deterministic")
#'   }
#'
#'   model$aqua <- structure(list(), class = aqua_class)
#'   model$aqua$lambda <- lambda_mat
#'
#' }
