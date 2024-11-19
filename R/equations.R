logLik <- function(model, theta) {
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
