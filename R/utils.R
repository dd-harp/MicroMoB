is_binary <- function(x) {
  UseMethod("is_binary", x)
}

is_binary.matrix <- function(x) {
  uniq <- unique(as.vector(x))
  if (all(length(uniq) == 2L)) {
    return(all(uniq %in% c(0, 1)))
  } else {
    return(FALSE)
  }
}

is_binary.numeric <- function(x) {
  uniq <- unique(x)
  if (all(length(uniq) == 2L)) {
    return(all(uniq %in% c(0, 1)))
  } else {
    return(FALSE)
  }
}

is_binary.integer <- function(x) {
  uniq <- unique(x)
  if (all(length(uniq) == 2L)) {
    return(all(uniq == c(0L, 1L)))
  } else {
    return(FALSE)
  }
}
