library(spatstat)


win <- unit.square()

cluster_rad <- 0.2
buffer <- 0.1

n_parents <- 6

# cif, par required
parents_model_pars <- rmhmodel(
  cif = 'hardcore',
  par = list(
    beta = 1,
    hc = cluster_rad + buffer
  ),
  w = win,
  types = 1:n_parents
)

parents_model_rmhcontrol <- rmhcontrol(p = 1 , fixall = TRUE)

parents_model_start <- rmhstart(n.start = rep(1, n_parents))

parents_model <- rmh(model = parents_model_pars, start = parents_model_start, control = parents_model_rmhcontrol)

pairdist(parents_model)
