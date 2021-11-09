


# setup.timespent.null = function(txt, Mod, dk=1, scl=1, pwr=1, stay=0.9){
#   print("Setting up timespent")
#   spendtimePars = list()
#   class(spendtimePars) = "null"
#   Mod$host$spendtimePars = spendtimePars
#
#   Fxpker = function(d){xpker(d, dk=dk, scl=scl, pwr=pwr)}
#   stay = rep(stay, Mod$landscape$nHaunts)
#   TimeSpent = diag(stay) + makePsiXX(Mod$landscape$haunts$XY, wts=1, Fd=Fxpker) %*% diag(1-stay)
#   Mod$TimeSpent = TimeSpent
#   return(Mod)
# }
