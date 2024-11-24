devtools::load_all()
library(lavaan)

tpb <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Causal Relationsships
  INT ~ a * ATT + b * SN + c * PBC
  BEH ~ INT + PBC 
'

msr <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2
'

# using a cfa
lsem <- function(y, x, S) {
  Sxx <- S[x, x]
  Syx <- S[y, x]

  Syx %*% solve(Sxx)
}

matrices <- lavInspect(cfa(msr, TPB), "estimates")
est <- sem(msr, TPB, num.grad = FALSE)
est
psi <- matrices$psi
lambda <- matrices$lambda
theta <- matrices$theta

Sff <- psi
Sxf <- lambda %*% psi
Sfx <- psi %*% t(lambda)
Sxx <- lambda %*% psi %*% t(lambda) + theta
S <- rbind(cbind(Sff, Sfx),
           cbind(Sxf, Sxx))

lsem(y="BEH", x=c("PBC", "INT"), S=S)
lsem(y="INT", x=c("PBC", "SN", "ATT"), S=S)
lsem(y=c("att1", "att2", "att3", "att4", "att5"), x="ATT", S=S)

# using efa
X <- cov(TPB)
E <- eigen(X)
EVal <- E$values
EVec <- E$vectors
k <- 5
LambdaUnrot <- EVec[, 1:5] %*% diag(sqrt(EVal[1:5]))
rotation <- GPArotation::oblimin(LambdaUnrot)
phi <- rotation$Phi[c(1, 4, 2, 3, 5)]
lambda <- rotation$loadings[, c(1, 4, 2, 3, 5)]
names(phi) <- c("ATT", "SN", "PBC", "INT", "BEH")
names(lambda) <- c("ATT", "SN", "PBC", "INT", "BEH")
Sff <- psi
Sxf <- lambda %*% psi
Sfx <- psi %*% t(lambda)
Sxx <- lambda %*% psi %*% t(lambda) + theta
S <- rbind(cbind(Sff, Sfx),
           cbind(Sxf, Sxx))

lsem(y="BEH", x=c("PBC", "INT"), S=S)
lsem(y="INT", x=c("PBC", "SN", "ATT"), S=S)
lsem(y=c("att1", "att2", "att3", "att4", "att5"), x="ATT", S=S)


replace_with_residuals <- function(S) {
  p <- ncol(S)
  residual_cov <- S
  
  for (i in seq_len(p)) {
    others <- setdiff(seq_len(p), i)
    
    S_xx <- S[others, others]      
    S_yx <- S[i, others, drop = FALSE]
    
    beta <-  S_yx %*% solve(S_xx)

    residual_variance <- S[i, i] - S_yx %*% t(beta)
    residual_cov[i, i] <- residual_cov[i, i] - residual_variance
  }
  
  residual_cov
}


fscores <- function(inds, data) {
  X <- cov(data[, inds])
  # S <- replace_with_residuals(X)
  E <- eigen(X)
  Eval <- E$values[1]
  Evec <- E$vectors[, 1]
  if (all(Evec < 0)) Evec <- -Evec
  
  lambda <- Evec * sqrt(Eval)
  print(lambda)
  Sff    <- as.matrix(X[1, 1])
  Sxf    <- lambda %*% Sff
  Sxx    <- X
  B      <- t(Sxf) %*% solve(Sxx)

  apply(data[, inds], MARGIN=1, FUN=\(X_i) B %*% X_i)
}

PBC <- fscores(c("pbc1", "pbc2", "pbc3"), data=TPB )
INT <- fscores(c("int1", "int2", "int3"), data=TPB )
SN  <- fscores(c("sn1", "sn2"), data=TPB )
BEH <- fscores(c("b1", "b2"), data=TPB )
ATT <- fscores(c("att1", "att2", "att3", "att4", "att5"), data=TPB )
lm(INT ~ PBC + ATT + SN)
lm(BEH ~ PBC + INT)
