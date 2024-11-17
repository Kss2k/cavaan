build_model <- function(parTable) {
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


  list(matrices=matrices)
}
