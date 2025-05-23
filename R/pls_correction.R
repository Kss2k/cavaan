getConsistenCorrMat <- function(model, P) {
  C <- model$matrices$C
  for (i in seq_len(nrow(C))) {
    for (j in seq_len(i - 1)) {
      x <- rownames(C)[[i]]
      y <- colnames(C)[[j]]
      C[x, y] <- C[y, x] <- C[x, y] / sqrt(P[[x]] * P[[y]])
    }
  }
  C
}


getConsistentLoadings <- function(model, P) {
  lVs <- model$info$lVs
  indsLvs <- model$info$indsLvs
  lambda <- model$matrices$lambda 
  for (lV in lVs) {
    wq <- lambda[indsLvs[[lV]], lV]
    lambda[indsLvs[[lV]], lV] <- wq %*% (sqrt(P[[lV]]) / t(wq) %*% wq)
  }
  lambda
}


getReliabilityCoefs <- function(model) {
  gamma    <- model$matrices$gamma
  lVs      <- model$info$lVs
  indsLvs  <- model$info$indsLvs
  lambda   <- model$matrices$lambda
  S        <- model$matrices$S
  P        <- vector("numeric", length(lVs))
  names(P) <- lVs

  for (lV in lVs) {
    if (length(indsLvs[[lV]]) <= 1) {
      P[[lV]] <- 1
      next
    }
    wq   <- lambda[indsLvs[[lV]], lV]
    subS <- S[indsLvs[[lV]], indsLvs[[lV]]]
    A    <- subS
    diag(A) <- 0

    B       <- wq %*% t(wq)
    diag(B) <- 0

    D <- t(wq) %*% A %*% wq 
    N <- t(wq) %*% B %*% wq 
    s <- (t(wq) %*% wq) ^ 2
    P[[lV]] <- s * D / N
  }

  P
}

rescaleVCOV <- function(vcov, vars, sigmas) {
  for (i in seq_along(sigmas)) {
    var <- vars[i]
    sigma1 <- sigmas[i]
    
    # Original variance of the variable
    sigma0 <- vcov[var, var]
    scalingFactor <- sqrt(sigma1 / sigma0)
    
    vcov[var, ] <- vcov[var, ] * scalingFactor
    vcov[, var] <- vcov[, var] * scalingFactor
    
    vcov[var, var] <- sigma1
  }
  
  vcov
}
