step0 <- function(model) {
  lambda <- model$matrices$lambda
  for (i in 1:ncol(lambda)) {
    lambda[, i] <- lambda[, i] / sum(lambda[, i])
  }

  partLambda <- cbind(model$matrices$Ip, lambda)
  S  <- model$matrices$S 
  C  <- model$matrices$C 
  SC <- model$matrices$SC
  
  # get new expected matrices
  model$matrices$C <- t(lambda) %*% S %*% lambda
  model$matrices$SC <- t(partLambda) %*% S %*% partLambda
  model 
}


# step 1: using path scheme
step1 <- function(model) {
  lVs   <- model$info$lVs
  succs <- model$matrices$succs
  preds <- model$matrices$preds
  gamma <- model$matrices$gamma
  S     <- model$matrices$S 
  C     <- model$matrices$C
  SC    <- model$matrices$SC 

  for (lV in lVs) {
    predsLv <- lVs[preds[ , lV, drop = TRUE]]
    succsLv <- lVs[succs[ , lV, drop = TRUE]]
    for (succ in succsLv) {
      gamma[succ, lV] <- C[lV, succ]
    }
    if (length(predsLv) > 0) {
      gamma[predsLv, lV] <- solve(SC[predsLv, predsLv]) %*% SC[predsLv, lV]
    }
    # standardize 
    gamma[, lV] <- gamma[, lV, drop = TRUE] / c(sqrt(t(gamma[, lV]) %*% C %*% gamma[, lV]))
  }
  model$matrices$gamma <- gamma
  model
}


step2 <- function(model) {
  lVs <- model$info$lVs
  Ip <- model$matrices$Ip
  lambda <- model$matrices$lambda
  partLambda <- cbind(Ip, lambda)
  gamma <- model$matrices$gamma
  partGamma <- rbind(cbind(Ip, matrix(0, nrow = nrow(Ip), ncol = ncol(gamma))),
                     cbind(matrix(0, nrow = nrow(gamma), ncol = ncol(Ip)), gamma))
  S <- model$matrices$S 
  C <- model$matrices$C 
  SC <- model$matrices$SC 
  model$matrices$C <- t(gamma) %*% C %*% gamma
  model$matrices$SC <- t(partGamma) %*% t(partLambda) %*% S %*% partLambda %*% partGamma
  dimnames(model$matrices$SC) <- dimnames(SC)
  model
}


# assuming that all lVs are reflective (i.e., using mode A)
step3 <- function(model) {
  lVs <- model$info$lVs
  indsLvs <- model$info$indsLvs
  lambda <- model$matrices$lambda
  S <- model$matrices$S
  SC <- model$matrices$SC 
  for (lV in lVs) {
    indsLv <- indsLvs[[lV]]
    for (indLv in indsLv) {
      lambda[indLv, lV] <- SC[indLv, lV]
    }
    lambda[, lV] <- lambda[, lV, drop = TRUE] / c(sqrt(t(lambda[, lV]) %*% S %*% lambda[, lV]))
  }
  model$matrices$lambda <- lambda
  model
}


# basically just step 0
step4 <- function(model) {
  lambda <- model$matrices$lambda
  partLambda <- cbind(model$matrices$Ip, lambda)
  S <- model$matrices$S 
  C <- model$matrices$C 
  SC <- model$matrices$SC
  
  # get new expected matrices
  model$matrices$C <- t(lambda) %*% S %*% lambda
  model$matrices$SC <- t(partLambda) %*% S %*% partLambda
  model 
}


step5 <- function(model, convergence = 1e-5)  {
  oldOuterWeights <- model$matrices$outerWeights
  newOuterWeights <- getNonZeroElems(model$matrices$lambda)
  if (all(abs((oldOuterWeights - newOuterWeights) / oldOuterWeights) < convergence)) {
    model$convergence <- TRUE 
  } else {
    model$convergence <- FALSE
  }
  model$matrices$outerWeights <- newOuterWeights
  model
}
