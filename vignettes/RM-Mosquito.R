## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(ramp.dts)
#devtools::load_all()

## ----eval=F-------------------------------------------------------------------
#  devtools::load_all()

## -----------------------------------------------------------------------------
rm1 <- dts_setup(Xname = "trace")

## -----------------------------------------------------------------------------
dts_solve(rm1, Tmax=200) -> rm1

## -----------------------------------------------------------------------------
dts_steady(rm1) -> rm1

## -----------------------------------------------------------------------------
with(rm1$outputs$orbits$MYZ[[1]], {
  plot(time, M, type = "l") 
  lines(time, U, col = "darkblue") 
  lines(time, Y, col = "purple") 
  lines(time, Z, col = "darkred") 
})

## -----------------------------------------------------------------------------
dts_plot_YZ(rm1)

## -----------------------------------------------------------------------------
compute_MYZ_equil = function(pars, Lambda, kappa, i=1){
  Omega = make_Omega(0, pars$MYZpar[[i]])
  with(pars$MYZpar[[i]],{

    Mbar <-  Lambda/(1-p) 
    Pbar <- f*Mbar/(1 - p + f)
    Ubar <- Lambda/(1-p*exp(-f*q*kappa))
    
    Y1 <- Omega*(1-exp(-f*q*kappa))*Ubar 
    Yi = Y1
    Y = Yi
    for(i in 2:(max_eip-1)){
      Yi <- Omega*Yi*(1-G[i])
      Y = Y+Yi
    }
    Yn <- Omega*Yi
    Y = Y + Yn
    Zbar <- Omega*Yn/(1-p) 
    return(c(M = Mbar, P=Pbar, U=Ubar, Y=Y, Z= Zbar)) 
})}

## -----------------------------------------------------------------------------
compute_MYZ_equil(rm1, 1000, .1) 

## -----------------------------------------------------------------------------
c(
M = tail(rm1$outputs$orbits$MYZ[[1]]$M, 1), 
P = tail(rm1$outputs$orbits$MYZ[[1]]$P, 1), 
U = tail(rm1$outputs$orbits$MYZ[[1]]$U, 1),
Y = tail(rm1$outputs$orbits$MYZ[[1]]$Y, 1),
Z = tail(rm1$outputs$orbits$MYZ[[1]]$Z, 1)
)

## -----------------------------------------------------------------------------
rm10 <- dts_setup(Xname = "trace", nPatches = 10, membership = 1:10)
rm10$Lpar[[1]]$scale = 1.5^c(1:10) 

## -----------------------------------------------------------------------------
dts_solve(rm10, Tmax=200) -> rm10

## -----------------------------------------------------------------------------
with(rm10$outputs$orbits$MYZ[[1]], {
  plot(time, M[,10], type = "l")
  for(i in 2:10)
    lines(time, M[,i]) 
#  lines(time, U, col = "darkblue") 
#  lines(time, Y, col = "purple") 
#  lines(time, Z, col = "darkred") 
})

