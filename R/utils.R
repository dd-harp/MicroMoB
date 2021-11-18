is_binary <- function(x) {
  uniq <- unique(as.vector(x))
  if (all(length(uniq) == 2L)) {
    return(all(uniq %in% c(0, 1)))
  } else {
    return(FALSE)
  }
}


Ff_0 <- function(x, fx = 1, sf = 1){
  fx * (sf * x / (1 + sf * x))
}
