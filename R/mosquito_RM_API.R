# functions to modify parameter values of the RM mosquito model

#' @title Set feeding rate for Ross-Macdonald mosquito model
#' @description Change the feeding rate parameter `f`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param f new blood feeding rate
#' @return no return value
#' @export
set_f_mosquito_RM <- function(model, f) {
  stopifnot(is.numeric(f), is.finite(f), f >= 0)
  stopifnot(length(f) == 1L)

  model$mosquito$f <- f
}

#' @title Get feeding rate for Ross-Macdonald mosquito model
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector
#' @export
get_f_mosquito_RM <- function(model) {
  return(model$mosquito$f)
}

#' @title Set human blood feeding fraction for Ross-Macdonald mosquito model
#' @description Change the human blood feeding fraction parameter `q`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param q new human blood feeding fraction
#' @return no return value
#' @export
set_q_mosquito_RM <- function(model, q) {
  stopifnot(is.numeric(q), is.finite(q), q >= 0)
  stopifnot(length(q) == 1L)

  model$mosquito$q <- q
}

#' @title Get human blood feeding fraction for Ross-Macdonald mosquito model
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector
#' @export
get_q_mosquito_RM <- function(model) {
  return(model$mosquito$q)
}

#' @title Set extrinsic incubation period for Ross-Macdonald mosquito model
#' @description Change the extrinsic incubation period parameter `eip` for some
#' set of times. The new values `eip` should either be a scalar or a vector
#' of length equal to the length of `times`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param eip new extrinsic incubation period values
#' @param times vector of times to set the new values
#' @return no return value
#' @export
set_eip_mosquito_RM <- function(model, eip, times) {
  stopifnot(is.numeric(eip), is.finite(eip), eip >= 0)
  if (length(eip) == 1L) {
    eip <- rep(x = eip, times = length(times))
  }
  stopifnot(length(eip) == length(times))
  stopifnot(all(times <= model$global$tmax), all(times > 0))

  model$mosquito$eip[times] <- eip

  # if maximum allowed EIP has become longer, we need to rebuild data structures
  if (max(model$mosquito$eip) > model$mosquito$maxEIP) {
    maxEIP <- max(model$mosquito$eip)
    model$mosquito$maxEIP <- maxEIP

    ZZ_shift <- matrix(0, nrow = maxEIP, ncol = maxEIP)
    ZZ_shift[1:(maxEIP-1), 2:maxEIP] <- diag(maxEIP-1)
    model$mosquito$ZZ_shift <- ZZ_shift

    model$mosquito$ZZ <- rbind(model$mosquito$ZZ, matrix(data = 0, nrow = maxEIP - nrow(model$mosquito$ZZ), ncol = model$global$p))
  }

}

#' @title Get extrinsic incubation period for Ross-Macdonald mosquito model
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param times vector of times to return
#' @return no return value
#' @export
get_eip_mosquito_RM <- function(model, times) {
  stopifnot(all(times <= model$global$tmax), all(times > 0))

  return(model$mosquito$eip[times])
}

#' @title Set daily survival probability for Ross-Macdonald mosquito model
#' @description Change the daily survival probability parameter `p` for some times
#' and places. The parameter `p` is stored internally as a matrix so that `times`
#' and `places` are used to modify a submatrix, therefore the new value `p` should
#' either be a scalar value to update the entire submatrix or a matrix of `places` rows
#' and `times` columns.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param p new human blood feeding fraction
#' @param times vector of times to set the new values
#' @param places vector of places to set the new values
#' @return no return value
#' @export
set_p_mosquito_RM <- function(model, p, times, places) {
  stopifnot(is.numeric(p), is.finite(p), p >= 0)
  if (length(p) == 1L) {
    p <- matrix(data = p, nrow = length(places), ncol = length(times))
  }
  stopifnot(is.matrix(p), nrow(p) == length(places), ncol(p) == length(times))
  stopifnot(all(times <= model$global$tmax), all(times > 0))
  stopifnot(all(places <= model$global$p), all(places > 0))

  model$mosquito$p[places, times] <- p
}

#' @title Get daily survival probability for Ross-Macdonald mosquito model
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param times vector of times to get values
#' @param places vector of places to get values
#' @return a matrix
#' @export
get_p_mosquito_RM <- function(model, times, places) {
  stopifnot(all(times <= model$global$tmax), all(times > 0))
  stopifnot(all(places <= model$global$p), all(places > 0))

  return(model$mosquito$p[places, times])
}

#' @title Set mosquito dispersal matrix for Ross-Macdonald mosquito model
#' @description Change the mosquito dispersal matrix parameter `psi`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param psi new mosquito dispersal matrix
#' @return no return value
#' @export
set_psi_mosquito_RM <- function(model, psi) {
  stopifnot(is.numeric(psi), is.finite(psi), psi >= 0)
  stopifnot(dim(psi) == c(model$global$p, model$global$p))
  stopifnot(approx_equal(rowSums(psi), 1))

  model$mosquito$psi <- psi
}

#' @title Get mosquito dispersal matrix for Ross-Macdonald mosquito model
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a matrix
#' @export
get_psi_mosquito_RM <- function(model) {
  return(model$mosquito$psi)
}

#' @title Set number of eggs laid per oviposition for Ross-Macdonald mosquito model
#' @description Change the number of eggs laid per oviposition parameter `nu`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param nu new number of eggs laid per oviposition
#' @return no return value
#' @export
set_nu_mosquito_RM <- function(model, nu) {
  stopifnot(is.numeric(nu), is.finite(nu), nu >= 0)
  stopifnot(length(nu) == 1L)

  model$mosquito$nu <- nu
}

#' @title Get number of eggs laid per oviposition for Ross-Macdonald mosquito model
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector
#' @export
get_nu_mosquito_RM <- function(model) {
  return(model$mosquito$nu)
}

#' @title Set kappa for Ross-Macdonald mosquito model
#' @description Change `kappa`.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param kappa new value of `kappa`
#' @return no return value
#' @export
set_kappa_mosquito_RM <- function(model, kappa) {
  stopifnot(is.numeric(kappa), is.finite(kappa), kappa >= 0)
  stopifnot(length(kappa) == model$global$p)

  model$mosquito$kappa <- kappa
}

#' @title Get kappa for Ross-Macdonald mosquito model
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @return a vector
#' @export
get_kappa_mosquito_RM <- function(model) {
  return(model$mosquito$kappa)
}
