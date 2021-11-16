#
# checkBloodFeeding.basic <- function(Mod) {
#   stopifnot(!is.null(Mod$NI))
#   stopifnot(inherits(Mod$NI, "matrix"))
#   stopifnot(dim(Mod$NI) == c(Mod$nHost, Mod$nPathogen))
#
#   stopifnot(length(Mod$mosquito) == Mod$nMosquito)
#   stopifnot(
#     vapply(X = lapply(X = Mod$mosquito, FUN = "[[", "Z.b"), FUN = function(Zb) {
#       all(dim(Zb) == c(Mod$nHost, Mod$nPathogen))
#     }, FUN.VALUE = logical(1), USE.NAMES = FALSE)
#   )
#
#
#
# }
#
#
# prop_bite <- c("n" = 0.9, "d" = 0.1)
# bite_day <- 3
#
#
# BloodFeeding.basic = function(t, Mod) {
#   # Get NI, the probability a mosquito would become infected after
#   # blood feeding on a host of each type.
#   # + for each human sub-population
#   # + for each pathogen
#   # + NI has dimensions |nHost| X |nPathogen|
#
#   # # Move to setup
#   # Mod$NI = matrix(0, Mod$nHost, Mod$nPathogen)
#   # for(i in 1:Mod$nMosquito) Mod$mosquito[[i]]$Z.b = matrix(0,Mod$nHost, Mod$nPathogen)
#   # for(j in 1:Mod$nPathogen) Mod$NI[,j] = getKappa(t, Mod, j)
#
#   # Get the NI for each human sub-population
#   AltHosts     <- AltHostTime(t, Mod)
#
#   ######################################################
#   # For each mosquito species...
#   ######################################################
#   for(i in 1:Mod$nMosquito){
#     ######################################################
#     # Skin in the Game
#     ######################################################
#     TaR    <- makeTaR(t, TimeSpent, mosquito[[i]])
#     SiG    <- makeSiG(t, TaR, mosquito[[i]])
#     AltSiG <- makeAltSiG(t, AltHosts, mosquito[[i]])
#
#     ######################################################
#     # Host Encounters
#     ######################################################
#     ModSiG <- ChooseHost(t, SiG, AltSiG, Mod, mosquito[[i]])
#     pb     <- SurviveBloodFeeding(t, SiG, ModSiG, AltSiG, Mod, mosquito[[i]])
#     W      <- rowSums(cbind(ModSiG, AltSiG))
#     psi    <- FeedingSuccess(t, W, mosquito[[i]])
#     DiB    <- ModSiG/W
#
#     ######################################################
#     # Save Mosquito Bionomic Parameters
#     ######################################################
#     Mod$mosquito[[i]]$adultPars$psi.b <- psi
#     Mod$mosquito[[i]]$adultPars$p.b <- pb
#     Mod$mosquito[[i]]$Q <- rowSums(DiB)
#
#     ######################################################
#     # Transmission from hosts to mosquitoes: kappa
#     # Dimensions of kappa: |nHaunts| X |nPathogen|
#     # For each mosquito species, the proportion that would
#     # become infected at each place on the landscape by
#     # each pathognen species.
#     ######################################################
#     Mod$mosquito[[i]]$kappa <- psi*(DiB %*% Mod$NI)
#
#     ######################################################
#     # Transmission from mosquitoes to hosts: EIR
#     # Dimensions of EIR are |nHost| X |nMosquito|
#     # For each sub-population, from each mosquito species,
#     # the expected number of infective bites per person.
#     #####################################################
#     for(j in 1:Mod$nPathogen){
#       Zb <- getZ(t, mosquito[[i]], j)
#       Mod$mosquito[[i]]$Z.b[, j] <- Zb
#       Mod$pathogen[[j]]$EIR <- colSums((psi*DiB*as.vector(Zb)) %*% diag(1/host$popSize))
#     }
#   }
# }
