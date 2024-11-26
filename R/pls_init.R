initPLS <- function(parTable, data, S=NULL) {
  matricesAndInfo      <- getMatricesAndInfoPLS(parTable)
  matrices             <- matricesAndInfo$matrices
  sortedData           <- sortData(data, matricesAndInfo$info$allInds) 
  matrices$S           <- if (is.null(S)) cov(as.data.frame(sortedData)) else S
  matrices$C           <- diag(nrow(matrices$gamma))
  dimnames(matrices$C) <- dimnames(matrices$gamma)

  matrices$SC <- rbind(
    cbind(matrices$S, matrix(0, nrow=nrow(matrices$S), ncol=nrow(matrices$C))),
    cbind(matrix(0, nrow=nrow(matrices$C), ncol=nrow(matrices$S)), matrices$C)
  )

  colnames(matrices$SC) <- rownames(matrices$SC) <- 
    c(colnames(matrices$S), colnames(matrices$C))
                
  model <- list(parTable=parTable, matrices=matrices, data=sortedData, 
                factorScores=NULL, info=matricesAndInfo$info, 
                params=NULL, fit=NULL)

  model$params <- list(names=getParamVecNames(model), 
                       values=rep(NA, length(getParamVecNames(model))))
  model
}


getMatricesAndInfoPLS <- function(parTable) {
  lVs <- parTable[parTable$op == "=~", "lhs"] |> unique()

  allInds <- vector("character", 0)
  indsLvs <- vector("list", length(lVs))
  names(indsLvs) <- lVs
  for (lV in lVs) {
    indsLv <- parTable[parTable$lhs == lV & parTable$op == "=~", "rhs"]
    allInds <- c(allInds, indsLv)
    indsLvs[[lV]] <- indsLv
  }
  
  lambda <- matrix(0, nrow=length(allInds), ncol=length(lVs),
                   dimnames=list(allInds, lVs))
  selectLambda <- matrix(FALSE, nrow=length(allInds), ncol=length(lVs),
                         dimnames=list(allInds, lVs))
  for (lV in lVs) {
    lambda[indsLvs[[lV]], lV] <- 1
    selectLambda[indsLvs[[lV]], lV] <- TRUE
  }

  gamma <- matrix(0, nrow=length(lVs), ncol=length(lVs),
                  dimnames=list(lVs, lVs))
  selectGamma <- matrix(FALSE, nrow=length(lVs), ncol=length(lVs),
                  dimnames=list(lVs, lVs))

  # Predecessors and successors
  preds <- succs <- matrix(FALSE, nrow=length(lVs), ncol=length(lVs),
                  dimnames=list(lVs, lVs))
  for (lV in lVs) {
    predsLv <- parTable[parTable$lhs == lV & parTable$op == "~", "rhs"]
    succsLv <- parTable[parTable$rhs == lV & parTable$op == "~", "lhs"]
    preds[predsLv, lV] <- TRUE
    succs[succsLv, lV] <- TRUE
    # selectionmatrix 
    selectGamma[predsLv, lV] <- TRUE
  }

  xis <- lVs[!lVs %in% parTable[parTable$op == "~", "lhs"]]
  selectCovXis <- matrix(FALSE, nrow=length(xis), ncol=length(xis),
                   dimnames=list(xis, xis))
  selectCovXis[lower.tri(selectCovXis)] <- TRUE
  
  Ip <- diag(nrow=nrow(lambda))
  colnames(Ip) <- rownames(Ip) <- rownames(lambda)

  matrices <- list(lambda=lambda, gamma=gamma, 
                   preds=preds, succs=succs, 
                   outerWeights=getNonZeroElems(lambda), 
                   Ip=Ip, 
                   selectLambda=selectLambda, selectGamma=selectGamma,
                   selectCovXis=selectCovXis)
  info <- list(indsLvs=indsLvs, allInds=allInds, lVs=lVs)

  list(matrices=matrices, info=info)
}
