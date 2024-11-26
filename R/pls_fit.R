getFitPLS <- function(model, consistent = TRUE, unstandardized = TRUE) {
  model$matrices$C <- model$matrices$C
   
  lambda   <- model$matrices$lambda
  gamma    <- model$matrices$gamma
  preds    <- model$matrices$preds
  lVs      <- model$info$lVs
  indsLvs  <- model$info$indsLvs
  parTable <- model$parTable
  SC       <- model$matrices$SC
  etas     <- parTable[parTable$op == "~", "lhs"] |> unique()
  xis      <- lVs[!lVs %in% etas]
  
  # measurement model 
  Lambda <- lambda
  Lambda[TRUE] <- 0
  for (lV in model$info$lVs) {
    for (indsLv in indsLvs[[lV]]) {
      Lambda[indsLv, lV] <- SC[indsLv, lV]
    }
  }

  # Caluculate consistent weights, based on measurement model
  if (consistent) {
    P                 <- getReliabilityCoefs(model)
    Lambda            <- getConsistentLoadings(model, P = P)
    Sff               <- getConsistenCorrMat(model, P = P)

    # Sxx <- Lambda %*% Sff %*% t(Lambda)
    Sxx <- model$matrices$S
    Sfx <- Sff %*% t(Lambda)
    Sxf <- Lambda %*% Sff
  
    model$matrices$C  <- Sff
    model$matrices$SC <- rbind(cbind(Sxx, Sxf),
                               cbind(Sfx, Sff))
  }

  if (unstandardized) {
    firstLoadings     <- getFirstLoadingsLVs(Lambda=Lambda, lVs=lVs)
    newVariances      <- firstLoadings ^ 2
    model$matrices$SC <- rescaleVCOV(model$matrices$SC, vars=lVs, sigmas=newVariances)
    model$matrices$C  <- model$matrices$SC[lVs, lVs]
  }
  
  # structural model
  Gamma <- gamma
  Gamma[TRUE] <- 0
  for (lV in lVs) {
    predsLv <- model$info$lVs[preds[ , lV, drop = TRUE]]
    if (length(predsLv) > 0) {
      Gamma[predsLv, lV] <- getPathCoefs(lV, predsLv,  model$matrices$C)
    }
  }

  Sigma <- model$matrices$SC[!grepl("__tmp", rownames(model$matrices$SC)), 
                             !grepl("__tmp", colnames(model$matrices$SC))]
  Lambda <- Lambda[!grepl("__tmp", rownames(Lambda)), 
                   !grepl("__tmp", colnames(Lambda))]
  Gamma <- Gamma[!grepl("__tmp", rownames(Gamma)),
                 !grepl("__tmp", colnames(Gamma))]
  
  list(Lambda = Lambda, Gamma = Gamma, Sigma = Sigma, info=model$info)
}


getParamVecNames <- function(model) {
  selectLambda <- model$matrices$selectLambda
  lambda <- model$matrices$lambda 
  for (j in colnames(lambda)) {
    for (i in rownames(lambda)) {
      lambda[i, j] <- paste0(j, " =~ ", i)
    }
  }

  selectGamma <- model$matrices$selectGamma
  gamma <- model$matrices$gamma
  for (j in colnames(gamma)) {
    for (i in rownames(gamma)) {
      gamma[i, j] <- paste0(j, " ~ ", i)
    }
  }
  
  selectCov <- model$matrices$selectCovXis
  covXis <- selectCov
  for (j in colnames(covXis)) {
    for (i in rownames(covXis)) {
      covXis[i, j] <- paste0(j, " ~~ ", i)
    }
  }

  c(lambda[selectLambda], gamma[selectGamma], covXis[selectCov]) 
}


extractCoefs <- function(model) {
  fit          <- model$fit 
  lambda       <- fit$Lambda
  selectLambda <- model$matrices$selectLambda
  gamma        <- fit$Gamma 
  selectGamma  <- model$matrices$selectGamma
  Sigma        <- fit$Sigma
  selectCov    <- model$matrices$selectCovXis

  c(lambda[selectLambda], gamma[selectGamma], Sigma[selectCov])
}
