warning2 <- function(...) {
  warning(..., call. = FALSE, immediate. = TRUE)
}


stop2 <- function(...) {
  stop(..., call. = FALSE)
}


stopif <- function(cond, ...) {
  if (cond) stop2(...)
}


warnif <- function(cond, ...) {
  if (cond) warning2(...)
}


dmatrix <- function(n, dimnames) {
  out <- diag(n)
  dimnames(out) <- dimnames
  out
}


superDiag <- function(...) {
  matrices <- list(...)

  Y <- NULL
  for (X in matrices) {
    if (is.null(Y)) {
      Y <- X
      next
    }
    
    n <- nrow(X)
    k <- ncol(X)
    m <- nrow(Y)
    p <- ncol(Y)

    rnames <- c(rownames(Y), rownames(X))
    cnames <- c(colnames(Y), colnames(X))
  
    YX <- matrix(0, nrow=m, ncol=k)
    XY <- matrix(0, nrow=n, ncol=p)
    Y <- rbind(cbind(Y, YX),
               cbind(XY, X))
    if (length(rnames) == n + m && length(cnames) == p + k) {
      dimnames(Y) <- list(rnames, cnames)
    }
  }

  Y
}


getDiagLike <- function(X) {
  I <- diag(nrow(X))
  dimnames(I) <- dimnames(X)
  I
}


lapplyNamed <- function(X, FUN, ..., names=names(X)) {
  structure(lapply(X, FUN, ...), names=names)
}


namedList <- function(n=0L, names=NULL) {
  structure(vector("list", n), names=names)
}
