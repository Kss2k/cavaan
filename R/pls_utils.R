sortData <- function(data, indicators) {
  if (!is.matrix(data)) data <- as.matrix(data)
  data[, indicators]
}


simpleLm <- function(x, y) {
  lm(y ~ x)
}


getCoefs <- function(lmObject) {
  lmObject$coefficients[-1]
}


getNonZeroElems <- function(x) {
  as.vector(x[!is.na(x) & x != 0])
}


standardizeMatrix <- function(x) {
  apply(x, MARGIN = 2, FUN = standardizeAtomic)
}


scaleMatrix <- function(x, weights) {
  for (i in seq_len(ncol(x))) x[, i] <- scaleAtomic(x[, i], sd = weights[[i]])
  x
}


standardizeAtomic <- function(x) {
  (x - mean(x)) / sd(x)
}


scaleAtomic <- function(x, sd = 1, center = TRUE) {
  if (center) x <- x - mean(x)
  (x / sd(x)) * sd
}


getScalingSDs <- function(lVs, indsLvs, data) {
  sds <- vector("numeric", length(lVs))
  names(sds) <- lVs
  for (lV in lVs) {
    sds[[lV]] <- sd(data[, indsLvs[[lV]][[1]]])
  }
  sds
}


getPathCoefs <- function(y, X, C) {
  # y: dependent factor 
  # X: independent factors
  # C: correlation matrix
  solve(C[X, X]) %*% C[X, y]
}


normalizeVec <- function(x) {
  # return x with length = 1 
  x / sqrt(t(x) %*% x)[[1]]
}


covProdInds <- function(x, y, data) {
  combos <- as.data.frame(expand.grid(x, y))
  colnames(combos) <- c("x", "y")
  prodInds <- apply(combos, MARGIN = 1, FUN = function(row) 
                    data[, row[[1]]] * data[, row[[2]]])
  colnames(prodInds) <- paste0(combos$x, ":", combos$y)
  cov(prodInds)
}


weightsProdInds <- function(wx, wy) {
  combos <- as.data.frame(expand.grid(wx, wy))
  colnames(combos) <- c("wx", "wy")
  w <- apply(combos, MARGIN = 1, FUN = function(row) 
        row[[1]] * row[[2]])

  if (!is.null(names(wx)) && !is.null(names(wy))) {
    comboNames <- as.data.frame(expand.grid(names(wx), names(wy)))
    colnames(comboNames) <- c("wx", "wy")
    names(w) <- apply(comboNames, MARGIN = 1, FUN = function(row) 
                      paste0(row[[1]], ":", row[[2]]))
  }

  w
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
