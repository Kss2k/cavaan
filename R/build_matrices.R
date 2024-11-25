buildA <- function(etas, mVYs) {
  m <- length(etas)
  p <- length(mVYs)
  Y <- c(etas, mVYs)
  matrix(0, nrow=m+p, ncol=m+p, dimnames=list(Y, Y))
}


buildB <- function(A) {
  diag(nrow(A)) - A
}


buildBStar <- function(A, IGamma) {
  B <- diag(nrow(A)) - A
  superDiag(B, IGamma)
}


buildGammaEpsilon <- function(etas, mVYs) {
  m <- length(etas)
  p <- length(mVYs)
  Y <- c(etas, mVYs)
  out <- matrix(0, nrow=m+p, ncol=m+p, dimnames=list(Y, Y))
  getDiagLike(out)
}


buildGammaXi <- function(etas, xis, mVYs) {
  n <- length(xis)
  m <- length(etas)
  p <- length(mVYs)
  Y <- c(etas, mVYs)
  matrix(0, nrow=m+p, ncol=n, dimnames=list(Y, xis))
}


buildGammaX <- function(etas, mVXs, mVYs) {
  q <- length(mVXs)
  m <- length(etas)
  p <- length(mVYs)
  Y <- c(etas, mVYs)
  matrix(0, nrow=m+p, ncol=q, dimnames=list(Y, mVXs))
}


buildIGamma <- function(xis, mVXs, etas, mVYs) {
  getDiagLike(buildPhi(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs))
}


buildGamma <- function(etas, xis, mVXs, mVYs) {
  Gamma <- cbind(buildGammaXi(etas=etas, xis=xis, mVYs=mVYs),
                 buildGammaX(etas=etas, mVXs=mVXs, mVYs=mVYs),
                 buildGammaEpsilon(etas=etas, mVYs=mVYs))
}


buildGammaStar <- function(Gamma, IGamma) {
  rbind(Gamma, IGamma)
}


buildPhi <- function(xis, mVXs, etas, mVYs) {
  n <- length(xis) 
  q <- length(mVXs)
  m <- length(etas)
  p <- length(mVYs)
  Y <- c(xis, mVXs, etas, mVYs)
  matrix(0, nrow=n+m+p+q, ncol=n+m+p+q, dimnames=list(Y, Y))
}


buildGy <- function(etas, mVYs) {
  p <- length(mVYs)
  m <- length(etas)
  Z <- matrix(0, nrow=p, ncol=m, dimnames=list(mVYs, etas))
  I <- dmatrix(n=p, dimnames=list(mVYs, mVYs))
  cbind(Z, I)
}


buildGx <- function(xis, mVXs, etas, mVYs) {
  n <- length(xis)
  q <- length(mVXs)
  m <- length(etas)
  p <- length(mVYs)
  Y <- c(etas, mVYs)

  Z1 <- matrix(0, nrow=q, ncol=n, dimnames=list(mVXs, xis))
  I  <- dmatrix(n=q, dimnames=list(mVXs, mVXs))
  Z2 <- matrix(0, nrow=q, ncol=m+p, dimnames=list(mVXs, Y))
  cbind(Z1, I, Z2)
}


buildG <- function(xis, etas, mVXs, mVYs) {
  Gy <- buildGy(etas=etas, mVYs=mVYs)
  Gx <- buildGx(xis=xis, mVXs=mVXs, etas=etas, mVYs=mVYs)
  superDiag(Gy, Gx)
}


buildTau <- function(xis, mVXs, etas, mVYs) {
  n <- length(xis)
  q <- length(mVXs)
  m <- length(etas)
  p <- length(mVYs)
  Y <- c(xis, mVXs, etas, mVYs)
  matrix(0, nrow=n+q+m+p, ncol=1, dimnames=list(Y, "~1"))
}
