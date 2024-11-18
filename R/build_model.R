build_submodel <- function(parTable) {
  aVs   <- getAVs(parTable)
  iVs   <- getIVs(parTable)
  etas  <- getEtas(parTable, sorted=TRUE)
  xis   <- getXis(parTable)
  mVYs  <- getMVYs(parTable)
  mVXs  <- getMVXs(parTable)

  Gamma <- buildGamma(etas=etas, xis=xis, mVXs=mVXs, mVYs=mVYs)
  A     <- buildA(etas=etas, mVYs=mVYs)
  Phi   <- buildPhi(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs)
  G     <- buildG(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs)
  B     <- buildB(A=A) 
  
  IGamma    <- buildIGamma(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs)
  BStar     <- buildBStar(A=A, IGamma=IGamma) 
  GammaStar <- buildGammaStar(Gamma=Gamma, IGamma=IGamma)

  matrices <- list(
    A=A, 
    Gamma=Gamma, 
    Phi=Phi, 
    G=G,
    BStar=BStar,
    GammaStar=GammaStar
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


build_model <- function(parTable) {
  groups <- unique(parTable$group)

  if (all(groups == "")) {
    models <- list(build_submodel(parTable))

  } else {
    models <- namedList(n=length(groups), names=groups)
    
    for (group in groups) {
      models[[group]] <- build_submodel(parTable[parTable$group == group, ])
    }
  }

  parTable.d <- getDetailedParTable(models, parTable=parTable)
  list(models=models, parTable.b=parTable, parTable.d=parTable.d)
}
