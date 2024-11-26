build_submodel <- function(parTable, data) {
  aVs   <- getAVs(parTable)
  iVs   <- getIVs(parTable)
  etas  <- getEtas(parTable, sorted=TRUE)
  xis   <- getXis(parTable)
  mVYs  <- getMVYs(parTable)
  mVXs  <- getMVXs(parTable)
  oVs   <- getOVs(parTable)
  lVs   <- getLVs(parTable)
  p     <- length(oVs)

  Gamma <- buildGamma(etas=etas, xis=xis, mVXs=mVXs, mVYs=mVYs)
  A     <- buildA(etas=etas, mVYs=mVYs)
  Phi   <- buildPhi(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs)
  G     <- buildG(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs)
  Tau   <- buildTau(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs)
  B     <- buildB(A=A) 
  
  IGamma    <- buildIGamma(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs)
  BStar     <- buildBStar(A=A, IGamma=IGamma) 
  BStarInv  <- solve(BStar)
  GammaStar <- buildGammaStar(Gamma=Gamma, IGamma=IGamma)

  Sigma <- G %*% BStar %*% GammaStar %*% Phi %*% t(GammaStar) %*% t(BStar) %*% t(G)
  S     <- stats::cov(data[ , colnames(Sigma)])
  Nu    <- as.matrix(colMeans(data[ , colnames(Sigma)]))
  Mu    <- G %*% BStar %*% GammaStar %*% Tau

  matrices <- list(
    A=A, 
    Gamma=Gamma, 
    Phi=Phi, 
    Tau=Tau,
    G=G,
    B=B,
    BStar=BStar,
    BStarInv=BStarInv,
    GammaStar=GammaStar,
    Sigma=Sigma,
    S=S,
    Nu=Nu,
    Mu=Mu,
    p=p
  )

  info <- list(
    etas=etas,
    xis=xis,
    mVYs=mVYs,
    mVXs=mVXs,
    aVs=aVs,
    iVs=iVs,
    lVs=lVs,
    oVs=oVs
  )

  list(matrices=matrices, info=info, parTable=parTable)
}


build_model <- function(parTable, data) {
  oVs    <- getOVs(parTable)
  groups <- unique(parTable$group)
  
  if (all(groups == 1)) {
    models <- list(build_submodel(parTable, data=data[, oVs]))

  } else {
    warning("GROUPS ARE NOT TESTED YET:\nAtleast these must be fixed!\n  1. Starting params are wrong! \n  2. Labels are wrong!")
    models <- namedList(n=length(groups), names=groups)
    
    for (group in groups) {
      models[[group]] <- build_submodel(parTable[parTable$group == group, ],
                                        data=data[data$group == group, oVs])
    }
  }
  mVYs <- models[[1]]$info$mVYs
  mVXs <- models[[1]]$info$mVXs
  xis  <- models[[1]]$info$xis
  etas <- models[[1]]$info$etas
  lVs  <- models[[1]]$info$lVs

  parTable.d    <- getDetailedParTable(models, parTable=parTable)
  plsEstimates  <- pls(parTable=parTable[parTable$group == 1, ], data=data)
  start         <- getStartingParams(parTable.d, pls=plsEstimates, xis=xis,
                                     etas=etas, mVYs=mVYs, mVXs=mVXs, lVs=lVs)

  list(models=models, parTable.b=parTable, parTable.d=parTable.d,
       data=data, groups=groups, start=start)
}
