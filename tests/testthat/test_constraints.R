devtools::load_all()
library(modsem)

m1 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  
  # Inner model
  Y ~ X + Z + X:Z 
'

time <- function(expr) {
  start <- Sys.time()
  result <- expr
  end <- Sys.time()
  cat("Elapsed = ", end - start, "\n")
  result
}

data   <- get_pi_data(m1, method="ca", data=oneInt)
syntax <- get_pi_syntax(m1, method="ca")

# for now it is sligthly faster
est_lav <- time(lavaan::sem(syntax, data))
est_cav <- time(sem(syntax, data, num.grad=TRUE))
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

