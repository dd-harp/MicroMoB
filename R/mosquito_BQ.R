rm(list=ls());gc()
suppressMessages(require(diagram))
suppressMessages(require(av))
suppressMessages(require(viridis))




frame_bq = function(b,q,mtl=NULL){
  plot(rbind(b,q), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", frame.plot=T,
       main = mtl)
}

frame_bqs = function(b,q,s, mtl=NULL){
  plot(rbind(b,q,s), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", frame.plot=T,
       main = mtl)
}

addP.b = function(b, wts=1, pw=1, adj=2, clr="red"){
  points(b, col = clr, pch = 15, cex=adj*wts^pw/max(wts^pw))
}

addP.q = function(q, wts=1, pw=1, adj=2, clr="blue"){
  points(q, col=clr, pch=19, cex=adj*wts^pw/max(wts^pw))
}

addP.s = function(s, wts=1, pw=1, adj=2, clr="orange"){
  #if(is.matrix(wts)){wts=rowSums(wts)}
  points(s, col=clr, pch=17, cex=adj*wts^pw/max(wts^pw))
}

addP.bb = function(b, M, pw=1, adj=2){
  addP.b(b, as.vector(rowSums(M)), pw, adj, "red")
  diag(M) <- 0
  addP.b(b, as.vector(rowSums(M)), pw, adj, "darkred")
}

addP.qq = function(q, M, pw=1, adj=2){
  addP.q(q, as.vector(rowSums(M)), pw, adj, "blue")
  diag(M)<-0
  addP.q(q, as.vector(rowSums(M)), pw, adj, "darkblue")
}

addP.ss = function(s, M, pw=1, adj=2){
  addP.s(s, as.vector(rowSums(M)), pw, adj, "orange")
  diag(M) <- 0
  addP.s(s, as.vector(rowSums(M)), pw, adj, "brown")
}

plotPoints_bq = function(b, q,
                         bwts=1, qwts=1, pw=1, adj=2,
                         mtl=NULL){
  frame_bq(b,q,mtl)
  addP.b(b, bwts, pw, adj)
  addP.q(q, qwts, pw, adj)
}

plotPoints_bqs = function(b, q, s,
                          bwts=1, qwts=1, swts=1, pw=1, adj=2,
                          mtl=NULL){

  frame_bqs(b,q,s,mtl)
  addP.b(b, bwts, pw, adj)
  addP.q(q, qwts, pw, adj)
  addP.s(s, swts, pw, adj)
}





lattice = function(n, mn, mx){
  points = seq(mn,mx,length.out=n)
  m = matrix(points, n, n)
  cbind(as.vector(m), as.vector(t(m)))
}



clusters = function(xy, nc, vr){
  np = dim(xy)[1]
  #nc = 1+rnbinom(np, mu=nc-1, size=1)
  nc = 1+rpois(np, nc-1)
  for(i in 1:np){
    xi = xy[i,1] + rnorm(nc[i], 0, vr)
    yi = xy[i,2] + rnorm(nc[i], 0, vr)
    xyi = cbind(xi,yi)
    xy = rbind(xy, xyi)
  }
  xy
}







kFmake.exp = function(k=1, s=2, gamma=1,  w=1){
  return(function(dd, w=1){
    wij = w*(exp(-k*(dd/s)^gamma))
    wij/max(wij)
  })
}

kFmake.pwr = function(delta=1, s=1, w=1){
  return(function(dd, w=1){
    wij = w/(dd+s)^delta
    wij/max(wij)
  })
}

kFmake.mix = function(p=0.9, k=1, s=2, gamma=1, delta=1, w=1){
  return(function(dd, w=1){
    wij = p*w*(exp(-k*(dd/s)^gamma)) + (1-p)*w/(dd+s)^delta
    wij/max(wij)
  })
}





arrowsX2Y = function(xy1, xy2, M, mnwd=0.05, bbend=1, endd=0.75, adj=2, clr="red"){
  n1 = dim(xy1)[1]
  n2 = dim(xy2)[1]
  adj = adj/max(M)
  for(i in 1:n1)
    for (j in c(1:n2))
      if (M[j,i] > mnwd*max(M)){
        fac = M[j,i]*adj
        curvedarrow(xy1[i,], xy2[j,],
                    segment=c(0.1, endd), lwd=fac, curve = 0.01*bbend,
                    arr.pos=endd, arr.length=0.15*fac, arr.width=0.1*fac,
                    arr.type="curved", arr.col = clr, lcol = clr)
      }
}



arrowsX2X = function(xy, M, mnwd=0.05, bbend=1, adj=2, endd=0.75, clr = "red" ){
  n = dim(xy)[1]
  diag(M) <- 0
  adj = adj/max(M)
  for(i in 1:n)
    for(j in c(1:n)[-i])
      if (M[j,i]> mnwd*max(M)){
        fac = M[j,i]*adj
        curvedarrow(xy[i,], xy[j,],
                    segment=c(0.1, endd), lwd=fac, curve = 0.01*bbend,
                    arr.pos=endd, arr.length=0.15*fac, arr.width=0.1*fac,
                    arr.type="curved", arr.col = clr, lcol = clr)
      }
}



avgX2X = function(xy, M, mnwd=0.05, bbend=1, adj=2, clr = grey(0.7) ){
  n = dim(xy)[1]
  diag(M)<-0
  M <- (M + t(M))/2
  scl <- max(M)
  adj = adj/max(M)
  for(i in 1:n)
    for(j in c(1:n)[-i])
      if(M[j,i]> mnwd*scl){
        fac = M[i,j]*adj
        segments(xy[i,1], xy[i,2], xy[j,1], xy[j,2], col=clr, lwd=fac)
      }
}



decompM=function(M){
  diagonal = diag(M)
  diag(M) <- 0
  connect = pmin(M, t(M))
  flow = t(M-connect)
  list(connect=connect, flow=flow, diagonal=diagonal)
}

flowX2X = function(xy, M, mnwd=0.05,  bbend=1, adj=2, endd=0.75, clr = "red" ){
  n = dim(xy)[1]
  diag(M)<-0
  scl <- max(M - t(M))
  adj = adj/scl

  for(i in 1:n)
    for(j in c(1:n)[-i])
      if(M[j,i] > M[i,j])
        if(M[j,i]> mnwd*scl){
          fac = (M[j,i]-M[i,j])*adj
          curvedarrow(xy[i,], xy[j,],
                      segment=c(0.1, endd), lwd=fac, curve = 0.01*bbend,
                      arr.pos=endd, arr.length=0.15*fac, arr.width=0.1*fac,
                      arr.type="curved", arr.col = clr, lcol = clr)
        }
}

symX2X = function(xy, M, mnwd=0.05, bbend=1, adj=2, clr = grey(0.7) ){
  n = dim(xy)[1]
  diag(M)<-0
  M <- pmin(M, t(M))
  scl <- max(M)
  adj = adj/max(M)
  for(i in 1:n)
    for(j in c(1:n)[-i])
      if(M[j,i]> mnwd*scl){
        fac = M[i,j]*adj
        segments(xy[i,1], xy[i,2], xy[j,1], xy[j,2], col=clr, lwd=fac)
      }
}



makePsi = function(S, D, kF, w=1){
  lS = length(S[,1])
  lD = length(D[,1])
  K = matrix(0, lD, lS)
  if(length(w)==1) w=rep(w, lD)
  for(i in 1:lS){
    K[,i] = kF(sqrt((S[i,1]-D[,1])^2 + (S[i,2]-D[,2])^2), w)
    K[,i] = K[,i]/sum(K[,i])
  }
  K
}

makePsi_stay = function(S, kF, w=1, stay=0){
  lS = dim(S)[1]
  K = matrix(0, lS, lS)
  if(length(w)==1) w=rep(w, lS)
  for(i in 1:lS){
    K[,i] = kF(sqrt((S[i,1]-S[,1])^2 + (S[i,2]-S[,2])^2), w)
    K[i,i] = 0
    K[,i] = (1-stay)*K[,i] /sum(K[,i])
    K[i,i] = stay
  }
  return(K)
}



PsiProfile_BQ = function(b,q,
                         Psi_bb, Psi_qb,
                         Psi_bq, Psi_qq){
  par (mfrow = c(2,2), mar = c(2,2,2,2))
  # b -> b
  frame_bq(b,q, mtl=expression(Psi[b%->%b]))
  symX2X(b, Psi_bb)
  flowX2X(b, Psi_bb)
  addP.bb(b, Psi_bb, adj=3)

  # q -> b
  frame_bq(b,q, mtl = expression(Psi[q%->%b]))
  arrowsX2Y(q, b, Psi_qb)
  addP.q(q, adj=1)
  addP.b(b, rowSums(Psi_qb), adj=3)

  # b -> q
  frame_bq(b,q, mtl = expression(Psi[b%->%q]))
  arrowsX2Y(b, q, Psi_bq, clr="blue")
  addP.q(q, rowSums(Psi_bq), adj=3)
  addP.b(b, adj=1)

  #q -> q
  frame_bq(b,q, mtl = expression(Psi[q%->%q]))
  symX2X(q,Psi_qq)
  flowX2X(q,Psi_qq, clr="blue")
  addP.qq(q, Psi_qq, adj=3)
}

PsiProfile.BQ = function(simObject){with(simObject,{
  PsiProfile_BQ(b,q,Psi_bb, Psi_qb, Psi_bq, Psi_qq)
})}

PsiProfile_BQS = function(b,q,s,
                          Psi_bb, Psi_qb, Psi_sb,
                          Psi_bq, Psi_qq, Psi_sq,
                          Psi_bs, Psi_qs, Psi_ss){
  par (mfrow = c(3,3), mar = c(2,2,2,2))
  # b -> b
  frame_bqs(b,q,s, mtl = expression(Psi[b%->%b]))
  symX2X(b, Psi_bb)
  flowX2X(b, Psi_bb)
  addP.bb(b, Psi_bb, adj=3)

  # q -> b
  frame_bqs(b,q,s, mtl = expression(Psi[q%->%b]))
  arrowsX2Y(q, b, Psi_qb)
  addP.q(q, adj=1)
  addP.b(b, rowSums(Psi_qb), adj=3)

  # s -> b
  frame_bqs(b,q,s, mtl = expression(Psi[s%->%b]))
  arrowsX2Y(s, b, Psi_sb, clr="red")
  addP.s(s, adj=1)
  addP.b(b, rowSums(Psi_sb), adj=3)

  # b -> q
  frame_bqs(b,q,s, mtl = expression(Psi[b%->%q]))
  arrowsX2Y(b, q, Psi_bq, clr="blue")
  addP.q(q, rowSums(Psi_bq), adj=3)
  addP.b(b, adj=1)

  #q -> q
  frame_bqs(b,q,s, mtl = expression(Psi[q%->%q]))
  symX2X(q,Psi_qq)
  flowX2X(q, Psi_qq, clr = "blue")
  addP.qq(q, Psi_qq, adj=3)

  #s -> q
  frame_bqs(b,q,s, mtl = expression(Psi[s%->%q]))
  arrowsX2Y(s, q, Psi_sq, clr="blue")
  addP.q(q, rowSums(Psi_sq), adj=3)
  addP.s(s, adj=1)

  # b -> s
  frame_bqs(b,q,s, mtl = expression(Psi[b%->%s]))
  arrowsX2Y(b, s, Psi_bs, clr="orange")
  addP.b(b, adj=1)
  addP.s(s, rowSums(Psi_bs), adj=3)

  #q -> s
  plotPoints_bqs(b,q,s, mtl = expression(Psi[q%->%s]))
  arrowsX2Y(q, s, Psi_qs, clr="orange")
  addP.q(q, adj=1)
  addP.s(s, rowSums(Psi_qs), adj=3)

  #s -> s
  plotPoints_bqs(b,q,s, mtl = expression(Psi[s%->%s]))
  symX2X(s, Psi_ss)
  flowX2X(s, Psi_ss, clr="orange")
  addP.ss(s, Psi_ss, adj=3)
}

PsiProfile.BQS = function(simObject){with(simObject,{
  PsiProfile_BQS(b,q,s, Psi_bb, Psi_qb, Psi_sb,
                 Psi_bq, Psi_qq, Psi_sq,
                 Psi_bs, Psi_qs, Psi_ss)
})}





dispersalPMF = function(xy, M){
  dxy = as.matrix(dist(xy, diag=TRUE, upper=TRUE))
  dxy = as.vector(dxy)
  M = as.vector(M)
  ot = order(dxy)
  dxy = dxy[ot]
  PMF = M[ot]/sum(M)
  CMF = cumsum(PMF)
  list(dxy=dxy, pmf=PMF, cmf=CMF, mean = sum(dxy*PMF))
}

plotDDpmf = function(xMFobj, mtl=NULL, clr="black", thresh=.99){with(xMFobj,{
  ix = 1:min(which(cmf>thresh))
  dxy = dxy[ix]
  pmf = pmf[ix]
  plot(dxy, pmf, main = mtl, type ="h", xlab = "Distance", ylab = "PMF", lwd=2)
  segments(mean, 0, mean, max(pmf), col = clr, lwd=2, lty=2)
})}

plotDDcmf = function(xMFobj, mtl=NULL, clr="black", thresh=.99){with(xMFobj,{
  ix = 1:min(which(cmf>thresh))
  dxy = dxy[ix]
  cmf = cmf[ix]
  plot(dxy, cmf, main = mtl, type ="l", xlab = "Distance", ylab = "CMF", ylim = c(0,1), lwd=2)
  segments(mean, 0, mean, 1, col = clr, lwd=2, lty=2)
})}



setupL = function(theta=.9, zeta=.01, xi=0, pL=.9, ova=20){
  p.L = list(theta=theta, zeta=zeta, xi=xi, pL=pL, ova=ova)
  class(p.L) <- "L"
  return(p.L)
}

setupL4P = function(theta=.9, zeta=.01, xi=0, pL=.9, ova=20){
  p.L = list(theta=theta, zeta=zeta, xi=xi, pL=pL, ova=ova)
  class(p.L) <- "L4P"
  return(p.L)
}

setup_BQ = function(pB=.98, pQ = 0.98, psiB=.8, psiQ=0.9){
  par = list(pB=pB,pQ=pQ,psiB=psiB,psiQ=psiQ)
  class(par) <- "BQ"
  return(par)
}

setup_BQS = function(pB=0.98, pQ=0.98, pS=0.98,
                     psiB=0.8, psiQ=0.9, psiS=0.8,
                     sigb=0.1, sigq=0.1, sigf=0.1, sigL=0.1){
  par = list(pB=pB,pQ=pQ,pS=pS,
             psiB=psiB,psiQ=psiQ,psiS=psiS,
             sigb=sigb,sigq=sigq,sigf=sigf,sigL=sigL)
  class(par) <- "BQS"
  return(par)
}

makeM_BQ = function(Pbb, Pqb, Pbq, Pqq, pars){with(pars,{

  # "hardened" adults
  Mbb = Pbb %*% diag(pB*(1-psiB), dim(Pbb)[1])
  Mbq = Pbq %*% diag(pB*psiB, dim(Pbb)[1])
  Mqb = Pqb %*% diag(pQ*psiQ, dim(Pqq)[2])
  Mqq = Pqq %*% diag(pQ*(1-psiQ), dim(Pqq)[2])

  bigM = rbind(
    cbind(Mbb, Mqb),
    cbind(Mbq, Mqq)
  )
  # recently emerged adults
  Mlb = Pqb %*% diag(pQ, dim(Pqq)[2])

  list(Mbb=Mbb, Mbq=Mbq, Mqb=Mqb, Mqq=Mqq, Mlb=Mlb, bigM=bigM)
})}

makeM_BQS = function(Pbb, Pqb, Psb,
                     Pbq, Pqq, Psq,
                     Pbs, Pqs, Pss,
                     pars){with(pars,{
                       nb = dim(Pbb)[1]
                       nq = dim(Pqq)[1]
                       ns = dim(Pss)[1]

                       Mbb = Pbb %*% diag(pB*(1-sigb)*(1-psiB), nb)
                       Mbq = Pbq %*% diag(pB*psiB, nb)
                       Mbs = Pbs %*% diag(pB*sigb*psiB, nb)
                       Mqb = Pqb %*% diag(pQ*(1-sigf)*psiQ, nq)
                       Mqq = Pqq %*% diag(pQ*(1-sigq)*(1-psiQ), nq)
                       Mqs = Pqs %*% diag(pQ*(sigf*psiQ + sigq*(1-psiQ)), nq)
                       Msb = Psb %*% diag(pS*psiS, ns)
                       Msq = 0*t(Mqs)
                       Mss = Pss %*% diag(1-psiS, ns)

                       Mlb = Pqb %*% diag(pQ*(1-sigL), nq)
                       Mls = Pqs %*% diag(pQ*sigL, nq)


                       bigM = rbind(
                         cbind(Mbb, Mqb, Msb),
                         cbind(Mbq, Mqq, Msq),
                         cbind(Mbs, Mqs, Mss)
                       )
                       return(list(Mbb=Mbb, Mbq=Mbq, Mbs=Mbs, Mqb=Mqb, Mqq=Mqq, Mqs=Mqs, Msb=Msb, Msq=Msq, Mss=Mss, Mlb=Mlb, bigM=bigM))
                     })}

makeSimObj_BQ = function(b, q,
                         # Kernel Shapes, Search Weights
                         kFb, kFq,
                         wb=1, wq=1,
                         stayB=0, stayQ=0,
                         # Adult Parameters
                         adultPars=setup_BQ(),
                         # Aquatic Parameters
                         aquaticPars=setupL(),
                         # Parasite Parameters
                         eip=15){
  Psi_bb = makePsi_stay(b,kFb,wb,stayB)
  Psi_bq = makePsi(b,q,kFq,wq)
  Psi_qb = makePsi(q,b,kFb,wb)
  Psi_qq = makePsi_stay(q,kFq,wq,stayQ)

  simObject = makeM_BQ(Psi_bb, Psi_qb, Psi_bq, Psi_qq, adultPars)

  simObject$b = b
  simObject$q = q
  simObject$nb = dim(b)[1]
  simObject$nq = dim(q)[1]

  simObject$Psi_bb = Psi_bb
  simObject$Psi_bq = Psi_bq
  simObject$Psi_qb = Psi_qb
  simObject$Psi_qq = Psi_qq

  class(simObject) <- "BQ"
  simObject$adultPars = adultPars
  simObject$pB = adultPars$pB
  simObject$pQ = adultPars$pQ
  simObject$psiB = adultPars$psiB
  simObject$psiQ = adultPars$psiQ

  simObject$aquaticPars = aquaticPars
  simObject$pL = aquaticPars$pL
  simObject$theta = aquaticPars$theta
  simObject$zeta = aquaticPars$zeta
  simObject$xi = aquaticPars$xi
  simObject$ova = aquaticPars$ova

  simObject$eip = eip

  return(simObject)
}

makeSimObj_BQS = function(b, q, Kfb, Kfq, Kfs, wb=1, wq=1, ws=1,
                          stayB=0.1, stayQ=0.1, stayS=0.1,
                          # Adult Parameters
                          adultPars=setup_BQS(),
                          # Aquatic Parameters
                          aquaticPars = setupL(),
                          # Parasite Parameters
                          eip=15){

  Psi_bb = makePsi_stay(b,kFb,wb,stayB)
  Psi_bq = makePsi(b,q,kFq,wq)
  Psi_bs = makePsi(b,s,kFs,ws)
  Psi_qb = makePsi(q,b,kFb,wb)
  Psi_qq = makePsi_stay(q,kFq,wq,stayQ)
  Psi_qs = makePsi(q,s,kFs,ws)
  Psi_sb = makePsi(s,b,kFs,wb)
  Psi_sq = makePsi(q,q,kFs,wq)
  Psi_ss = makePsi_stay(s,kFs,ws,stayS)

  simObject = makeM_BQS(Psi_bb, Psi_bq, Psi_bs, Psi_qb, Psi_qq, Psi_qs, Psi_sb, Psi_sq, Psi_ss)

  simObject$b=b
  simOqj$q=q
  simOsj$s=s

  nb = dim(b)[2]
  nq = dim(q)[2]
  ns = dim(s)[2]

  simObject$Psi_bb = Psi_bb
  simObject$Psi_bq = Psi_bq
  simObject$Psi_bs = Psi_bs
  simObject$Psi_qb = Psi_qb
  simObject$Psi_qq = Psi_qq
  simObject$Psi_qs = Psi_qs
  simObject$Psi_sb = Psi_sb
  simObject$Psi_sq = Psi_sq
  simObject$Psi_ss = Psi_ss

  class(simObject) <- "BQS"
  simObject$adultPars= adultPars
  simObject$pB = adultPars$pB
  simObject$pQ = adultPars$pQ
  simObject$pS = adultPars$pS
  simObject$psiB = adultPars$psiB
  simObject$psiQ = adultPars$psiQ
  simObject$psiS = adultPars$psiS
  simObject$sigb = adultPars$sigb
  simObject$sigq = adultPars$sigq
  simObject$sigf = adultPars$sigf
  simObject$sigL = adultPars$sigL

  simObject$aquaticPars= aquaticPars
  simObject$pL = aquaticPars$pL
  simObject$theta = aquaticPars$theta
  simObject$zeta = aquaticPars$zeta
  simObject$xi = aquaticPars$xi
  simObject$ova = aquaticPars$ova
  simObject$eip = aquaticPars$eip

  return(simObject)
}

init = function(simObject){
  UseMethod("init", simObject)
}

initL = function(simObject){
  UseMethod("initL", simObject$aquaticPars)
}

init.BQ = function(simObject, B0=10, Q0=10){
  simObject$B = matrix(B0, simObject$nb, 1)
  simObject$Q = matrix(Q0, simObject$nq, 1)
  simObject$eggs = matrix(0, simObject$nq, 1)
  simObject = initL(simObject)
  return(simObject)
}

init.BQS = function(simObject, B0=10, Q0=10, S0=10){
  simObject$B = matrix(B0, simObject$nb, 1)
  simObject$Q = matrix(Q0, simObject$nq, 1)
  simObject$S = matrix(S0, simObject$ns, 1)
  simObject = initL(simObject)
  return(simObject)
}

initL.L = function(simObject, L0=10){
  simObject$L = matrix(L0, simObject$nq, 1)
  simObject$Lambda = matrix(0, simObject$nq, 1)
  return(simObject)
}

initL.L4 = function(simObject, L1=0, L2=0, L3=0, L4=0, P=0){
  simObject$L1 = matrix(L1, simObject$nq, 1)
  simObject$L2 = matrix(L2, simObject$nq, 1)
  simObject$L3 = matrix(L3, simObject$nq, 1)
  simObject$L4 = matrix(L4, simObject$nq, 1)
  simObject$P = matrix(P, simObject$nq, 1)
  simObject$Lambda = matrix(0, simObject$nq, 1)
  simObject$L = simObject$L1 + simObject$L2 + simObject$L3 + simObject$L4
  return(simObject)
}

AquaticDynamics = function(simObject){
  UseMethod("AquaticDynamics", simObject$aquaticPars)
}

AquaticDynamics.L = function(simObject){ with(simObject,{
  survive = pL*exp(-zeta*L)
  mature = theta*exp(-xi*L)
  simObject$Lambda = mature*survive*L
  Lt = (1-mature)*survive*L
  Lt = Lt + eggs
  simObject$L= Lt
  return(simObject)
})}

AquaticDynamics.L4P = function(simObject){ with(simObject,{
  simObject$Lambda = pL*P
  survive = pL*exp(-zeta*L)
  mature = theta*exp(-xi*L)
  L1t = (1-mature)*survive*L1
  L2t = (1-mature)*survive*L2 + survive*mature*L1
  L3t = (1-mature)*survive*L3 + survive*mature*L2
  L4t = (1-mature)*survive*L4 + survive*mature*L3
  Pt = surv*mature*L4
  L1t = L1t + eggs
  simObject$L1 = L1t
  simObject$L2 = L2t
  simObject$L3 = L3t
  simObject$L4 = L4t
  simObject$P = Pt
  simObject$L = L1t+L2t+L3t+L4t
  return(simObject)
})}

AdultDynamics = function(simObject){
  UseMethod("AdultDynamics", simObject)
}

AdultDynamics.BQ = function(simObject){ with(simObject,{
  simObject$eggs = ova*psiQ*Q
  Bt = Mbb %*% B + Mqb %*% Q + Mlb %*% Lambda
  Qt = Mbq %*% B + Mqq %*% Q
  simObject$B = Bt
  simObject$Q = Qt
  return(simObject)
})}

AdultDynamics.BQS = function(simObject){ with(simObject,{
  simObject$eggs = ova*psiQ*Q
  Bt = Mbb %*% B + Mqb %*% Q + Msb%*%S + Mlb %*% Lambda
  Qt = Mbq %*% B + Mqq %*% Q + Msq%*%S
  St = Mbs %*% B + Mqs %*% Q + Mss%*%S + Mls %*% Lambda
  simObject$B = Bt
  simObject$Q = Qt
  simObject$S = St
  return(simObject)
})}

saveStates=function(states, simObject){
  UseMethod("saveStates", simObject)
}

saveLStates=function(states, simObject){
  UseMethod("saveLStates", simObject$aquaticPars)
}

saveStates.BQ = function(states, simObject){
  if(is.null(states))
    return(list(Bt = simObject$B, Qt = simObject$Q))

  states = with(states,{
    Bt = cbind(Bt, simObject$B)
    Qt = cbind(Qt, simObject$Q)
    list(Bt=Bt, Qt=Qt)
  })

  return(states)
}

saveStates.BQS = function(states, simObject){
  if(is.null(states))
    return(list(Bt = simObject$B, Qt=simObject$Q, St=simObject$St))

  states = with(states,{
    Bt = cbind(Bt, simObject$B)
    Qt = cbind(Qt, simObject$Q)
    St = cbind(St, simObject$S)
    list(Bt=Bt, Qt=Qt, St=St)
  })

  return(states)
}

saveLStates.L = function(states, simObject){
  if(is.null(states))
    return(list(Lt = simObject$L))

  states = with(states,{
    Lt = cbind(Lt, simObject$L)
    list(Lt=Lt)
  })

  return(states)
}

saveLStates.L4P = function(states, simObject){
  if(is.null(states)){
    states = with(simObject,list(L1t=L1,L2t=L2,L3t=L3,L4t=L4,Pt=P))
    return(states)
  }
  states = with(states,{
    L1t = cbind(L1t, simObject$L1)
    L2t = cbind(L2t, simObject$L2)
    L3t = cbind(L3t, simObject$L3)
    L4t = cbind(L4t, simObject$L4)
    Pt = cbind(Pt, simObject$P)
    list(L1t=L1,L2t=L2,L3t=L3,L4t=L4,Pt=P)
  })

  return(states)
}

SIM = function(simObject, Tmax=200, TS=FALSE){
  states = saveStates(NULL, simObject)
  Lstates = saveLStates(NULL, simObject)

  for(i in 1:Tmax){
    simObject = AdultDynamics(simObject)
    simObject = AquaticDynamics(simObject)
    if(TS==TRUE) {
      states = save(states, simObject)
      Lstates = saveLStates(Lstates, simObject)
    }
  }

  return(list(simObject = simObject, states=states, Lstates=Lstates))
}

steadyState = function(simObject, Titer=50, tol=.001){
  simObject = init(simObject)
  simObject = SIM(simObject, Titer)$simObject

  err=10*tol
  while(err>tol){
    Bi = simObject$B; Qi = simObject$Q; Li = simObject$L
    simObject = SIM(simObject, Titer)$simObject
    err = sum(abs(simObject$B - Bi)) +
      sum(abs(simObject$Q - Qi)) +
      sum(abs(simObject$L - Li))
  }

  steadyState = list()
  steadyState$B = simObject$B
  steadyState$Q = simObject$Q
  steadyState$L = simObject$L
  simObject$steadyState = steadyState
  return(simObject)
}







makeKqb = function(simObject, Tmax){
  UseMethod("makeKqb", simObject)
}

makeKbq = function(simObject, Tmax){
  UseMethod("makeKbq", simObject)
}



makeKqb.BQ = function(simObject, Tmax=100){with(simObject,{
  Bt = Mqb %*% diag(1, nq)
  Kqb = 0*Bt
  for(i in 1:Tmax){
    Kqb = Kqb + Psi_bb %*% diag(pB*psiB, nb) %*% Bt
    Bt = Mbb%*%Bt
  }
  simObject$Kqb = Kqb
  return(simObject)
})}





makeKqb.BQLS = function(Pbb, Pqb, Psb,
                        Pbq, Pqq, Psq,
                        Pbs, Pqs, Pss,
                        sigf = 0.5, sigq = 0.5, sigb=0.5,
                        pB=.98, pS =0.99, pQ = 0.98,
                        psiB=.8, psiQ=0.9, psiS = 0.98){
  nb = dim(Pbb)[2]
  nq = dim(Pqq)[2]
  ns = dim(Pss)[2]

  Mbb = pB*(1-sigb)*(1-psiB)*Pbb
  Mbq = pB*psiB*Pbq
  Mbs = pB*sigb*psiB*Pbs
  Mqb = pQ*(1-sigf)*psiQ*Pqb
  Mqq = pQ*(1-sigq)*(1-psiQ)*Pqq
  Mqs = pQ*(sigf*psiQ + sigq*(1-psiQ))*Pqs
  Msb = pS*psiS*Psb
  Msq = 0*t(Mqs)
  Mss = (1-psiS)*Pss

  M1 = cbind(Mbb, Mqb, Msb, 0*Mbb)
  M2 = cbind(0*Mbq, Mqq, Msq, 0*Mbq)
  M3 = cbind(0*Mbs, Mqs, Mss, 0*Mbs)
  M4 = cbind(diag(1,nb), 0*Mqb, 0*Msb, diag(1,nb))
  M = rbind(M1, M2, M3, M4)

  Kt = rbind(Mqb %*% diag(1,nq), 0*Mqq, 0*Mqs, 0*Mqb)
  for(i in 1:100) Kt = M%*%Kt
  Kt[-c(1:(nb+nq+ns)),]
}



makeKbq.BQ = function(simObject, Tmax=100){with(simObject,{
  Qt = Psi_bq %*% diag(pB, nb)
  Kbq = 0*Qt
  for(i in 1:Tmax){
    Kbq = Kbq + Psi_qq %*% diag(pQ*psiQ, nq) %*% Qt
    Qt = Mqq%*%Qt
  }
  simObject$Kbq = Kbq
  return(simObject)
})}

makeKbq_BQS = function(Pbb, Pqb, Psb,
                       Pbq, Pqq, Psq,
                       Pbs, Pqs, Pss,
                       sigf = 0.5, sigq = 0.5, sigb = 0.5,
                       pB =.98, pS = 0.99, pQ = 0.98,
                       psiB =.8, psiQ = 0.9, psiS = 0.98){
  nb = dim(b)[1]
  nq = dim(q)[1]
  ns = dim(s)[1]

  Mbb = pB*(1-sigb)*(1-psiB)*Pbb
  Mbq = pB*psiB*Pbq
  Mbs = pB*sigb*psiB*Pbs
  Mqb = pQ*(1-sigf)*psiQ*Pqb
  Mqq = pQ*(1-sigq)*(1-psiQ)*Pqq
  Mqs = pQ*(sigf*psiQ + sigq*(1-psiQ))*Pqs
  Msb = pS*psiS*Psb
  Msq = 0*t(Mqs)
  Mss = (1-psiS)*Pss

  M1 = cbind(Mbb, 0*Mqb, Msb, 0*Mqb)
  M2 = cbind(Mbq, Mqq, Msq, 0*Mqq)
  M3 = cbind(Mbs, 0*Mqs, Mss, 0*Mqs)
  M4 = cbind(0*Mbq, diag(1,nq), 0*Msq, diag(1,nq))
  M = rbind(M1, M2, M3, M4)

  Kt = rbind(0*Mbb, Mbq %*% diag(1,nb), 0*Mbs, 0*Mbq)
  for(i in 1:100) Kt = M%*%Kt
  Kt[-c(1:(nb+nq+ns)),]
}





makeModel_BQ = function(b, q,
                        # Kernel Shapes, Search Weights
                        kFb, kFq,
                        wb=1, wq=1,
                        stayB=0, stayQ=0,
                        # Adult Parameters
                        adultPars = setup_BQ(),
                        # Aquatic Parameters
                        aquaticPars =  setupL(),
                        # Parasite Parameters
                        eip=15){

  model = makeSimObj_BQ(b, q, kFb, kFq, wb, wq, stayB, stayQ, adultPars, aquaticPars, eip)
  model = steadyState(model)
  model = makeKbq(model)
  model = makeKqb(model)
  Kbb = with(model, Kqb %*% Kbq)
  model$Kbb = Kbb
  model$KBB = Kbb %*% diag(as.vector(model$steadyState$B))
  Kqq = with(model, Kbq %*% Kqb)
  model$Kqq = Kqq
  model$KQQ = Kqq %*% diag(as.vector(model$steadyState$Q))
  model = computeG(model)
  model = computeGG(model)
  model = computeV(model)
  model = computeVC(model)
  model
}

makeModel_BQS = function(b, q, s,
                         # Kernel Shapes, Search Weights
                         kFb, kFq, kFs,
                         wb=1, wq=1, ws=1,
                         stayB=0, stayQ=0, stayS=0,
                         # Adult Parameters
                         adultPars = setup_BQS(),
                         # Aquatic Parameters
                         aquaticPars =  setupL(),
                         # Parasite Parameters
                         eip=15){

  model = makeSimObj_BQS(b, q, s, kFb, kFq, kFs, wb, wq, ws, stayB, stayQ, stayS, adultPars, aquaticPars, eip)
  model = steadyState(model)
  model = makeKbq(model)
  model = makeKqb(model)
  Kbb = with(model, Kqb %*% Kbq)
  model$Kbb = Kbb
  model$KBB = Kbb %*% diag(as.vector(model$steadyState$B))
  Kqq = with(model, Kbq %*% Kqb)
  model$Kqq = Kqq
  model$KQQ = Kqq %*% diag(as.vector(model$steadyState$Q))
  model = computeG(model)
  model = computeGG(model)
  model = computeV(model)
  model = computeVC(model)
  model
}

makeMovie.BQ = function(mod, B0, Q0, L0, moviestem="BQLmovie", Tmax=200){
  mod = init(mod, B0, Q0, L0, 0, 0, 0)
  out = SIM(mod, Tmax=200)
  Bt = out$Bt; Qt=out$Qt
  sclB = max(sqrt(Bt))
  sclQ = max(sqrt(Qt))
  framenames = c()
  with(mod,{
    for(i in 1:Tmax){
      fname = paste("./tmp/movie", 100+i, ".png", sep="")
      framenames = c(framenames, fname)
      png(fname, 480, 480)
      plotPoints_bq(b, q, mtl = paste("t =", i), add=FALSE)
      points(b, pch=15, col= "red", cex=2*sqrt(Bt[,i])/sclB)
      points(q, pch=19, col= "blue", cex=2*sqrt(Qt[,i])/sclQ)
      dev.off(dev.cur())
    }
    moviename = paste(moviestem, ".mp4", sep = "")
    av_encode_video(framenames, moviename, verbose=F) -> out
  })
}





xmin=ymin=-20
xmax=ymax=20
nq = 20
q = cbind(x=runif(nq,xmin,xmax), y=runif(nq,ymin,ymax))
nb = 15
b = cbind(x=runif(nb,xmin,xmax), y=runif(nb,ymin,ymax))
ns = 10
s = cbind(x=runif(ns,xmin,xmax), y=runif(ns,ymin,ymax))

# Searching for {b}
kFb = kFmake.exp(k=2, s=4, gamma=1.5)
Psi_bb = makePsi_stay(b,kFb, stay=0.2)
Psi_qb = makePsi(q,b,kFb)
Psi_sb = makePsi(s,b,kFb)

# Searching for {q}
kFq = kFmake.exp(k=2, s=2, gamma=2)
Psi_bq = makePsi(b,q,kFq)
Psi_qq = makePsi_stay(q,kFq, stay=0.2)
Psi_sq = makePsi(s,q,kFq)

# Searching for {s}
kFs = kFmake.exp(k=2, s=1, gamma=1.5)
Psi_qs = makePsi(q,s,kFs)
Psi_bs = makePsi(b,s,kFs)
Psi_ss = makePsi_stay(s, kFs, stay = 0.2)



# Mod1 = makeSimObj_BQ(b, q, kFb, kFq)
# Mod2 = makeSimObj_BQ(bc, qc, kFb, kFq)
# Mod1 = init(Mod1)
# out = SIM(Mod1, 20)
# Mod1 = steadyState(Mod1)
#
# fac = 2.5
# par(mfrow = c(2,2))
# with(Mod1,{
#   plotPoints_bq(b, q, B, 0)
#   plotPoints_bq(b, q, 0, Q)
#   plotPoints_bq(b, q, 0, 0)
#   addP.q(q, L)
#   plot(steadyState$Q, steadyState$L/(steadyState$Q), xlab = "Adults (Q)", ylab = "Larvae : Adults")
# })
#
# par(mfrow = c(3,1))
# with(Mod1,{
#   hist(steadyState$B, 40, col = "red", main = "B")
#   hist(steadyState$Q, 40, col = "blue", main = "Q")
#   hist(steadyState$L, 40, col = "blue", main = "L")
# })




simObject = makeSimObj_BQ(b, q, kFb, kFq)

B0 <- 100
Q0 <- 50
simObject$B = matrix(B0, simObject$nb, 1)
simObject$Q = matrix(Q0, simObject$nq, 1)
simObject$eggs = matrix(0, simObject$nq, 1)
simObject = initL(simObject)

simObject$eggs = simObject$ova*simObject$psiQ*simObject$Q
# Bt = Mbb %*% B + Mqb %*% Q + Mlb %*% Lambda
Bt = simObject$Mbb %*% simObject$B + simObject$Mqb %*% simObject$Q
Qt = simObject$Mbq %*% simObject$B + simObject$Mqq %*% simObject$Q
simObject$B = Bt
simObject$Q = Qt


simObject1 = makeSimObj_BQ(b, q, kFb, kFq)
B0 <- 100
Q0 <- 50
simObject1$B = matrix(B0, simObject1$nb, 1)
simObject1$Q = matrix(Q0, simObject1$nq, 1)
simObject1$eggs = matrix(0, simObject1$nq, 1)
simObject1 = initL(simObject1)

BQ <- rbind(simObject1$B, simObject1$Q)

simObject1$bigM %*% BQ




# code

xmin=ymin=-20
xmax=ymax=20
nq = 5
q = cbind(x=runif(nq,xmin,xmax), y=runif(nq,ymin,ymax))
nb = 4
b = cbind(x=runif(nb,xmin,xmax), y=runif(nb,ymin,ymax))

# Searching for {b}
kFb = kFmake.exp(k=2, s=4, gamma=1.5)
Psi_bb = makePsi_stay(b,kFb, stay=0.2)
Psi_qb = makePsi(q,b,kFb)

# Searching for {q}
kFq = kFmake.exp(k=2, s=2, gamma=2)
Psi_bq = makePsi(b,q,kFq)
Psi_qq = makePsi_stay(q,kFq, stay=0.2)


psiB <- c(0.95, 0.9, 0.875, 0.85)
psiQ <- c(0.99, 0.98, 0.95, 0.90, 0.925)

pB <- c(0.9, 0.95, 0.925, 0.91)
pQ <- c(0.9, 0.95, 0.925, 0.91, 0.075)

# kappa P(infection)
kappa <- c(0.01, 0.05, 0.075, 0.1)

Mbq_inf = Psi_bq %*% diag(pB*psiB*kappa, dim(Psi_bb)[1])
Mbq_noinf = Psi_bq %*% diag(pB*psiB*(1-kappa), dim(Psi_bb)[1])

Mbb = Psi_bb %*% diag(pB*(1-psiB), dim(Psi_bb)[1])
Mbq = Psi_bq %*% diag(pB*psiB, dim(Psi_bb)[1])
Mqb = Psi_qb %*% diag(pQ*psiQ, dim(Psi_qq)[2])
Mqq = Psi_qq %*% diag(pQ*(1-psiQ), dim(Psi_qq)[2])

bigM = rbind(
  cbind(Mbb, Mqb),
  cbind(Mbq, Mqq)
)

bigM[nq:(nq+nb), 1:nb]




maxEIP <- 4
EIP_shift <- matrix(data = 0, nrow = maxEIP + 1, ncol = maxEIP + 1)
EIP_shift[1, 1] <- 1
EIP_shift[2:(maxEIP+1), 1:maxEIP] <- diag(x = 1, nrow = maxEIP, ncol = maxEIP)

Y <- matrix(data = rpois(n = (maxEIP+1)*2, lambda = 10), nrow = 2, ncol = maxEIP+1)

Y %*% EIP_shift






# classes and methods to implement a simple behavioral state model of mosquitoes with infection and incubation

#' @title Setup blood feeding & oviposition (BQ) behavioral state mosquito model
#' @description This is a behavioral state model which allows for time varying EIP and
#' survival probability. Mosquitoes transition between blood feeding (B) and
#' oviposition (Q) depending on teh success (or not) of those biological activities.
#' It complies with the mosquito component interface, and may be simulated deterministically or stochastically.
#' @param model an object from [MicroMoB::make_MicroMoB]
#' @param stochastic should the model update deterministically or stochastically?
#' @param f the blood feeding rate, should be a function which accepts a vector `B` and returns a vector of feeding rates of the same length
#' @param eip the Extrinsic Incubation Period (may be time varying see [MicroMoB::time_varying_parameter])
#' @param pB daily survival probability during blood feeding (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param pQ daily survival probability during oviposition (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param psiQ oviposition success probability (may be time and patch varying see [MicroMoB::time_patch_varying_parameter])
#' @param Psi_bb movement matrix from blood feeding haunts to blood feeding haunts (columns must sum to 1, `p` rows and columns)
#' @param Psi_bq movement matrix from blood feeding haunts to aquatic habitats (columns must sum to 1, `l` rows and `p` columns)
#' @param Psi_qb movement matrix from aquatic habitats to blood feeding haunts (columns must sum to 1, `p` rows and `l` columns)
#' @param Psi_qq movement matrix from aquatic habitats to aquatic habitats (columns must sum to 1, `l` rows and columns)
#' @param nu number of eggs laid per oviposition
#' @param M number of susceptible mosquitoes (vector of length `p + l`)
#' @param Y number of incubating mosquitoes (vector of length `p + l`)
#' @param Z number of infectious mosquitoes (vector of length `p + l`)
#' @return no return value
#' @export
setup_mosquito_BQ <- function(model, stochastic, f, eip, pB, pQ, psiQ, Psi_bb, Psi_bq, Psi_qb, Psi_qq, nu = 25, M, Y) {
  stopifnot(inherits(model, "MicroMoB"))
  stopifnot(is.logical(stochastic))

  p <- model$global$p
  l <- model$global$l
  tmax <- model$global$tmax

  eip_vec <- time_varying_parameter(param = eip, tmax = tmax)

  maxEIP <- max(eip_vec)

  pB_mat <- time_patch_varying_parameter(param = pB, p = p, tmax = tmax)
  pQ_mat <- time_patch_varying_parameter(param = pQ, p = l, tmax = tmax)
  psiQ_mat <- time_patch_varying_parameter(param = psiQ, p = l, tmax = tmax)

  stopifnot(dim(Psi_bb) == p)
  stopifnot(dim(Psi_bq) == c(l, p))
  stopifnot(dim(Psi_qb) == c(p, l))
  stopifnot(dim(Psi_qq) == l)
  stopifnot(approx_equal(colSums(Psi_bb), 1))
  stopifnot(approx_equal(colSums(Psi_bq), 1))
  stopifnot(approx_equal(colSums(Psi_qb), 1))
  stopifnot(approx_equal(colSums(Psi_qq), 1))

  stopifnot(inherits(f, 'function'))

  stopifnot(length(nu) == 1L)

  mosy_class <- c("BQ")
  if (stochastic) {
    mosy_class <- c(mosy_class, "BQ_stochastic")
  } else {
    mosy_class <- c(mosy_class, "BQ_deterministic")
  }

  model$mosquito <- structure(list(), class = mosy_class)
  model$mosquito$f <- f
  model$mosquito$psiQ_mat <- psiQ_mat
  model$mosquito$nu <- nu
  model$mosquito$eip <- eip_vec
  model$mosquito$maxEIP <- maxEIP

  model$mosquito$pB_mat <- pB_mat
  model$mosquito$pQ_mat <- pQ_mat

  model$mosquito$Psi_bb <- Psi_bb
  model$mosquito$Psi_bq <- Psi_bq
  model$mosquito$Psi_qb <- Psi_qb
  model$mosquito$Psi_qq <- Psi_qq

  model$mosquito$kappa <- rep(0, p)

  stopifnot(length(M) == p+l)
  stopifnot(nrow(Y) == p+l)
  stopifnot(ncol(Y) == maxEIP+1)
  stopifnot(is.finite(M))
  stopifnot(is.finite(Y))
  stopifnot(M >= 0)
  stopifnot(Y >= 0)

  model$mosquito$M <- matrix(data = M, ncol = 1)
  model$mosquito$Y <- Y

  # matrix which multiplies Y on the right to shift all by one day
  EIP_shift <- matrix(data = 0, nrow = maxEIP + 1, ncol = maxEIP + 1)
  EIP_shift[1, 1] <- 1
  EIP_shift[2:(maxEIP+1), 1:maxEIP] <- diag(x = 1, nrow = maxEIP, ncol = maxEIP)
  model$mosquito$EIP_shift <- EIP_shift

}



# update mosquitoes over one time step

#' @title Update blood feeding & oviposition (BQ) behavioral state mosquitoes
#' @description This function dispatches on the second argument of `model$mosquito`
#' for stochastic or deterministic behavior.
#' @inheritParams step_mosquitoes
#' @return no return value
#' @details see [MicroMoB::step_mosquitoes.BQ_deterministic] and [MicroMoB::step_mosquitoes.BQ_stochastic]
#' @export
step_mosquitoes.BQ <- function(model) {
  NextMethod()
}



Mbq_inf = Psi_bq %*% diag(pB*psiB*kappa, dim(Psi_bb)[1])
Mbq_noinf = Psi_bq %*% diag(pB*psiB*(1-kappa), dim(Psi_bb)[1])

Mbb = Psi_bb %*% diag(pB*(1-psiB), dim(Psi_bb)[1])
Mbq = Psi_bq %*% diag(pB*psiB, dim(Psi_bb)[1])
Mqb = Psi_qb %*% diag(pQ*psiQ, dim(Psi_qq)[2])
Mqq = Psi_qq %*% diag(pQ*(1-psiQ), dim(Psi_qq)[2])

#' @title Update blood feeding & oviposition (BQ) behavioral state mosquitoes (deterministic)
#' @inheritParams step_mosquitoes
#' @return no return value
#' @importFrom stats pexp
#' @export
step_mosquitoes.BQ_deterministic <- function(model) {

  # parameters
  tnow <- model$global$tnow
  EIP <- model$mosquito$eip[tnow]
  p <- model$global$p
  l <- model$global$l

  psiB <- pexp(q = model$mosquito$f * model$mosquito$q)
  psiQ <- model$mosquito$psiQ_mat[, tnow]

  pB <- model$mosquito$pB_mat[, tnow]
  pQ <- model$mosquito$pQ_mat[, tnow]

  # construct update matrices
  Mbq_inf <- model$mosquito$Psi_bq %*% diag(x = pB*psiB*kappa, ncol = p)
  Mbq_noinf <- model$mosquito$Psi_bq %*% diag(x = pB*psiB*(1-kappa), ncol = p)

  Mbb <- model$mosquito$Psi_bb %*% diag(x = pB*(1-psiB), ncol = p)
  Mbq <- model$mosquito$Psi_bq %*% diag(x = pB*psiB, ncol = p)
  Mqb <- model$mosquito$Psi_qb %*% diag(x = pQ*psiQ, ncol = l)
  Mqq <- model$mosquito$Psi_qq %*% diag(x = pQ*(1-psiQ), ncol = l)

  M <- matrix(data = 0, nrow = l+p, ncol = l+p)
  M[1:p, 1:p] <- Mbb
  M[1:p, l:(l+p)] <- Mqb
  M[l:(l+p), 1:p] <- Mbq
  M[l:(l+p), l:(l+p)] <- Mqq

  M_noinf <- matrix(data = 0, nrow = l+p, ncol = l+p)
  M_noinf[1:p, 1:p] <- Mbb
  M_noinf[1:p, l:(l+p)] <- Mqb
  M_noinf[l:(l+p), 1:p] <- Mbq_noinf
  M_noinf[l:(l+p), l:(l+p)] <- Mqq

  M_inf <- matrix(data = 0, nrow = l+p, ncol = l+p)
  M_inf[l:(l+p), 1:p] <- Mbq_inf

  # update
  model$mosquito$M <- M_noinf %*% model$mosquito$M
  for (i in 1:maxEIP) {
    model$mosquito$Y[, i] <- M %*% model$mosquito$Y[, i]
  }
  model$mosquito$Y <- model$mosquito$Y %*% model$mosquito$EIP_shift
  model$mosquito$Y[, EIP] <- M_inf %*% model$mosquito$M

  # newly emerging adults
  lambda <- compute_emergents(model)

  model$mosquito$M[1:p, 1] <- lambda

}

