sortData <- function(data, indicators) {
  if (!is.matrix(data)) data <- as.matrix(data)
  data[, indicators]
}


getNonZeroElems <- function(x) {
  as.vector(x[!is.na(x) & x != 0])
}


getPathCoefs <- function(y, x, S) {
  # y: dependent varibales
  # x: independent variables
  # S: correlation matrix
  Sxx <- S[x, x]
  Syx <- S[y, x]
  Syx %*% solve(Sxx)
}


diagPartitioned <- function(X, Y) {
  out <- rbind(cbind(X, matrix(0, nrow = nrow(X), ncol = ncol(Y))),
               cbind(matrix(0, nrow = nrow(Y), ncol = ncol(X)), Y))
  colnames(out) <- c(colnames(X), colnames(Y))
  rownames(out) <- c(rownames(X), rownames(Y)) 
  out
}


getFirstLoadingsLVs <- function(Lambda, lVs, tol = 1e-8) {
  out <- c()
  for (lV in lVs) {
    loadings <- Lambda[ , lV, drop = TRUE]
    out <- c(out, loadings[min(which(abs(loadings) > tol))])
  }

  out
}


redefOVsParTable <- function(parTable, oVs) {
  lVs  <- paste0("__tmp__", oVs) 
  rows <- data.frame(lhs=lVs, op="=~", rhs=oVs)
  for (i in seq_along(oVs)) {
    oV <- oVs[i]
    lV <- lVs[i]
    parTable[parTable$op == "~" & parTable$rhs == oV, "rhs"] <- lV
    parTable[parTable$op == "~" & parTable$lhs == oV, "rhs"] <- lV
  }

  rbindNA(parTable, rows)
}


getResidualsVCOV <- function(S) {
  p       <- ncol(S)
  R       <- S
  R[TRUE] <- 0
  
  for (i in seq_len(p)) {
    others <- setdiff(seq_len(p), i)
    
    S_xx <- S[others, others]      
    S_yx <- S[i, others, drop = FALSE]
    
    beta <-  S_yx %*% solve(S_xx)

    e <- S[i, i] - S_yx %*% t(beta)
    R[i, i] <- e
  }

  R
}


getResidualsLVs <- function(fit) {
  lVs   <- fit$info$lVs
  Sigma <- fit$Sigma
  Sigma <- Sigma[colnames(Sigma) %in% lVs, colnames(Sigma) %in% lVs]
  getResidualsVCOV(Sigma) 
}


getResidualsOVs <- function(fit) {
  oVs   <- fit$info$oVs
  Sigma <- fit$Sigma
  R <- getResidualsVCOV(Sigma)
  R[colnames(R) %in% oVs, colnames(R) %in% oVs]
}
