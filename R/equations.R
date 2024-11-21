logLik <- function(theta, model) {
  model.f <- fillModel(model, params=theta, calc.sigma=TRUE)
  logLik  <- 0
  
  for (group in model$groups) {
    matrices <- model.f$models[[group]]$matrices
    logLik  <- logLik + logLikMatrices(matrices)
  }

  logLik 
}


logLikMatrices <- function(matrices) {
  Sigma <- matrices$Sigma
  S     <- matrices$S
  p     <- matrices$p

  ln(det(Sigma)) + tr(solve(Sigma) %*% S) - ln(det(S)) - p
}


ln <- function(x) {
  log(x, base=exp(1))
}


tr <- function(X) {
  sum(diag(X))
}


gradientLogLik <- function(model, theta) {
  model.f     <- fillModel(model, params=theta, calc.sigma=TRUE)
 
  gradients.b <- model$baseGradients
  groups      <- model$groups
  gradient    <- rep(0, length(theta))

  Qs <- vector("list", length(groups))
  for (group in groups) {
    S        <- model.f$models[[group]]$matrices$S
    Sigma    <- model.f$models[[group]]$matrices$Sigma
    SigmaInv <- solve(Sigma)

    Qs[[group]] <- SigmaInv - SigmaInv %*% S %*% SigmaInv
  }

  for (i in seq_along(theta)) {
    t <- theta[i]
    gradients.bt <- gradients.b[[i]]
    
    for (group in groups) {
      matrices <- model.f$models[[group]]$matrices
      Sigma.g  <- getGradientSigma(matrices, gradients.bt[[group]])
      gradient[i] <- gradient[i] + tr(Qs[[group]] %*% Sigma.g)
    }
  }
  
  gradient
}


getGradientSigma <- function(matrices, gradients.bt) {
  # Normal Matrices
  G         <- matrices$G
  GammaStar <- matrices$GammaStar
  Phi       <- matrices$Phi
  BStarInv  <- matrices$BStarInv

  # Gradients
  BStar.g     <- gradients.bt$BStar
  GammaStar.g <- gradients.bt$GammaStar
  Phi.g       <- gradients.bt$Phi
  BStarInv.g  <- - BStarInv %*% BStar.g %*% BStarInv

  G %*% (
      BStarInv.g %*% GammaStar   %*% Phi   %*% t(GammaStar)   %*% t(BStarInv) +
      BStarInv   %*% GammaStar.g %*% Phi   %*% t(GammaStar)   %*% t(BStarInv) +
      BStarInv   %*% GammaStar   %*% Phi.g %*% t(GammaStar)   %*% t(BStarInv) +
      BStarInv   %*% GammaStar   %*% Phi   %*% t(GammaStar.g) %*% t(BStarInv) +
      BStarInv   %*% GammaStar   %*% Phi   %*% t(GammaStar)   %*% t(BStarInv.g)
  ) %*% t(G)
}
