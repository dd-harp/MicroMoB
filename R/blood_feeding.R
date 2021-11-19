

# bf.bloodfeeding.Ff0

compute_beta <- function(human) {
  UseMethod("compute_beta", human$tisp)
}
