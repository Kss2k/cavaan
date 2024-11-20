build_submodel <- function(parTable, data) {
  aVs   <- getAVs(parTable)
  iVs   <- getIVs(parTable)
  etas  <- getEtas(parTable, sorted=TRUE)
  xis   <- getXis(parTable)
  mVYs  <- getMVYs(parTable)
  mVXs  <- getMVXs(parTable)
  oVs   <- getOVs(parTable)
  p     <- length(oVs)

  Gamma <- buildGamma(etas=etas, xis=xis, mVXs=mVXs, mVYs=mVYs)
  A     <- buildA(etas=etas, mVYs=mVYs)
  Phi   <- buildPhi(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs)
  G     <- buildG(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs)
  B     <- buildB(A=A) 
  
  IGamma    <- buildIGamma(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs)
  BStar     <- buildBStar(A=A, IGamma=IGamma) 
  BStarInv  <- solve(BStar)
  GammaStar <- buildGammaStar(Gamma=Gamma, IGamma=IGamma)

  Sigma <- G %*% BStar %*% GammaStar %*% Phi %*% t(GammaStar) %*% t(BStar) %*% t(G)
  S     <- cov(data[ , colnames(Sigma)])

  matrices <- list(
    A=A, 
    Gamma=Gamma, 
    Phi=Phi, 
    G=G,
    B=B,
    BStar=BStar,
    BStarInv=BStarInv,
    GammaStar=GammaStar,
    Sigma=Sigma,
    S=S,
    p=p
  )

  info <- list(
    etas=etas,
    xis=xis,
    mVYs=mVYs,
    mVXs=mVXs,
    aVs=aVs,
    iVs=iVs
  )

  list(matrices=matrices, info=info, parTable=parTable)
}


build_model <- function(parTable, data) {
  oVs    <- getOVs(parTable)
  groups <- unique(parTable$group)
  
  if (all(groups == 1)) {
    models <- list(build_submodel(parTable, data=data[, oVs]))

  } else {
    models <- namedList(n=length(groups), names=groups)
    
    for (group in groups) {
      models[[group]] <- build_submodel(parTable[parTable$group == group, ],
                                        data=data[data$group == group, oVs])
    }
  }

  parTable.d    <- getDetailedParTable(models, parTable=parTable)
  baseGradients <- getBaseGradients(models, parTable.d=parTable.d)
  start         <- getStartingParams(parTable.d)

  list(models=models, parTable.b=parTable, parTable.d=parTable.d,
       baseGradients=baseGradients, data=data, groups=groups, start=start)
}
