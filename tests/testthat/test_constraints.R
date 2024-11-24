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
syntax <- "
X =~ lambda_x1_X*x1
X =~ lambda_x2_X*x2
X =~ lambda_x3_X*x3
Y =~ lambda_y1_Y*y1
Y =~ lambda_y2_Y*y2
Y =~ lambda_y3_Y*y3
Z =~ lambda_z1_Z*z1
Z =~ lambda_z2_Z*z2
Z =~ lambda_z3_Z*z3
Y ~ Gamma_X_Y*X
Y ~ Gamma_Z_Y*Z
Y ~ Gamma_XZ_Y*XZ
XZ =~ lambda_x1z1_XZ*x1z1
XZ =~ lambda_x2z1_XZ*x2z1
XZ =~ lambda_x3z1_XZ*x3z1
XZ =~ lambda_x1z2_XZ*x1z2
XZ =~ lambda_x2z2_XZ*x2z2
XZ =~ lambda_x3z2_XZ*x3z2
XZ =~ lambda_x1z3_XZ*x1z3
XZ =~ lambda_x2z3_XZ*x2z3
XZ =~ lambda_x3z3_XZ*x3z3
X ~~ Var_X*X
Y ~~ Zeta_Y*Y
Z ~~ Var_Z*Z
XZ ~~ Var_XZ*XZ
x1 ~~ Var_x1*x1
x2 ~~ Var_x2*x2
x3 ~~ Var_x3*x3
y1 ~~ Var_y1*y1
y2 ~~ Var_y2*y2
y3 ~~ Var_y3*y3
z1 ~~ Var_z1*z1
z2 ~~ Var_z2*z2
z3 ~~ Var_z3*z3
x1z1 ~~ Var_x1z1*x1z1
x2z1 ~~ Var_x2z1*x2z1
x3z1 ~~ Var_x3z1*x3z1
x1z2 ~~ Var_x1z2*x1z2
x2z2 ~~ Var_x2z2*x2z2
x3z2 ~~ Var_x3z2*x3z2
x1z3 ~~ Var_x1z3*x1z3
x2z3 ~~ Var_x2z3*x2z3
x3z3 ~~ Var_x3z3*x3z3
X ~~ Cov_X_Z*Z
X ~~ Cov_X_XZ*XZ
Z ~~ Cov_Z_XZ*XZ
x2z1 ~~ Cov_x2z1_x1z1*x1z1
x3z1 ~~ Cov_x3z1_x1z1*x1z1
x3z1 ~~ Cov_x3z1_x2z1*x2z1
x1z2 ~~ Cov_x1z2_x1z1*x1z1
x1z2 ~~ 0*x2z1
x1z2 ~~ 0*x3z1
x2z2 ~~ 0*x1z1
x2z2 ~~ Cov_x2z2_x2z1*x2z1
x2z2 ~~ 0*x3z1
x2z2 ~~ Cov_x2z2_x1z2*x1z2
x3z2 ~~ 0*x1z1
x3z2 ~~ 0*x2z1
x3z2 ~~ Cov_x3z2_x3z1*x3z1
x3z2 ~~ Cov_x3z2_x1z2*x1z2
x3z2 ~~ Cov_x3z2_x2z2*x2z2
x1z3 ~~ Cov_x1z3_x1z1*x1z1
x1z3 ~~ 0*x2z1
x1z3 ~~ 0*x3z1
x1z3 ~~ Cov_x1z3_x1z2*x1z2
x1z3 ~~ 0*x2z2
x1z3 ~~ 0*x3z2
x2z3 ~~ 0*x1z1
x2z3 ~~ Cov_x2z3_x2z1*x2z1
x2z3 ~~ 0*x3z1
x2z3 ~~ 0*x1z2
x2z3 ~~ Cov_x2z3_x2z2*x2z2
x2z3 ~~ 0*x3z2
x2z3 ~~ Cov_x2z3_x1z3*x1z3
x3z3 ~~ 0*x1z1
x3z3 ~~ 0*x2z1
x3z3 ~~ Cov_x3z3_x3z1*x3z1
x3z3 ~~ 0*x1z2
x3z3 ~~ 0*x2z2
x3z3 ~~ Cov_x3z3_x3z2*x3z2
x3z3 ~~ Cov_x3z3_x1z3*x1z3
x3z3 ~~ Cov_x3z3_x2z3*x2z3
Cov_x2z1_x1z1 == lambda_x2_X * lambda_x1_X * (Var_X) * Var_z1
Cov_x3z1_x1z1 == lambda_x3_X * lambda_x1_X * (Var_X) * Var_z1
Cov_x3z1_x2z1 == lambda_x3_X * lambda_x2_X * (Var_X) * Var_z1
Cov_x1z2_x1z1 == lambda_z2_Z * lambda_z1_Z * (Var_Z) * Var_x1
Cov_x2z2_x2z1 == lambda_z2_Z * lambda_z1_Z * (Var_Z) * Var_x2
Cov_x2z2_x1z2 == lambda_x2_X * lambda_x1_X * (Var_X) * Var_z2
Cov_x3z2_x3z1 == lambda_z2_Z * lambda_z1_Z * (Var_Z) * Var_x3
Cov_x3z2_x1z2 == lambda_x3_X * lambda_x1_X * (Var_X) * Var_z2
Cov_x3z2_x2z2 == lambda_x3_X * lambda_x2_X * (Var_X) * Var_z2
Cov_x1z3_x1z1 == lambda_z3_Z * lambda_z1_Z * (Var_Z) * Var_x1
Cov_x1z3_x1z2 == lambda_z3_Z * lambda_z2_Z * (Var_Z) * Var_x1
Cov_x2z3_x2z1 == lambda_z3_Z * lambda_z1_Z * (Var_Z) * Var_x2
Cov_x2z3_x2z2 == lambda_z3_Z * lambda_z2_Z * (Var_Z) * Var_x2
Cov_x2z3_x1z3 == lambda_x2_X * lambda_x1_X * (Var_X) * Var_z3
Cov_x3z3_x3z1 == lambda_z3_Z * lambda_z1_Z * (Var_Z) * Var_x3
Cov_x3z3_x3z2 == lambda_z3_Z * lambda_z2_Z * (Var_Z) * Var_x3
Cov_x3z3_x1z3 == lambda_x3_X * lambda_x1_X * (Var_X) * Var_z3
Cov_x3z3_x2z3 == lambda_x3_X * lambda_x2_X * (Var_X) * Var_z3
Var_XZ == (Var_X) * (Var_Z) + (Cov_X_Z) ^ 2
Cov_X_XZ == 0
Cov_Z_XZ == 0
Var_x1z1 == lambda_x1_X ^ 2 * (Var_X) * Var_z1 + lambda_z1_Z ^ 2 * (Var_Z) * Var_x1 + Var_x1 * Var_z1
Var_x2z1 == lambda_x2_X ^ 2 * (Var_X) * Var_z1 + lambda_z1_Z ^ 2 * (Var_Z) * Var_x2 + Var_x2 * Var_z1
Var_x3z1 == lambda_x3_X ^ 2 * (Var_X) * Var_z1 + lambda_z1_Z ^ 2 * (Var_Z) * Var_x3 + Var_x3 * Var_z1
Var_x1z2 == lambda_x1_X ^ 2 * (Var_X) * Var_z2 + lambda_z2_Z ^ 2 * (Var_Z) * Var_x1 + Var_x1 * Var_z2
Var_x2z2 == lambda_x2_X ^ 2 * (Var_X) * Var_z2 + lambda_z2_Z ^ 2 * (Var_Z) * Var_x2 + Var_x2 * Var_z2
Var_x3z2 == lambda_x3_X ^ 2 * (Var_X) * Var_z2 + lambda_z2_Z ^ 2 * (Var_Z) * Var_x3 + Var_x3 * Var_z2
Var_x1z3 == lambda_x1_X ^ 2 * (Var_X) * Var_z3 + lambda_z3_Z ^ 2 * (Var_Z) * Var_x1 + Var_x1 * Var_z3
Var_x2z3 == lambda_x2_X ^ 2 * (Var_X) * Var_z3 + lambda_z3_Z ^ 2 * (Var_Z) * Var_x2 + Var_x2 * Var_z3
Var_x3z3 == lambda_x3_X ^ 2 * (Var_X) * Var_z3 + lambda_z3_Z ^ 2 * (Var_Z) * Var_x3 + Var_x3 * Var_z3
lambda_x1z1_XZ == lambda_x1_X * lambda_z1_Z
lambda_x2z1_XZ == lambda_x2_X * lambda_z1_Z
lambda_x3z1_XZ == lambda_x3_X * lambda_z1_Z
lambda_x1z2_XZ == lambda_x1_X * lambda_z2_Z
lambda_x2z2_XZ == lambda_x2_X * lambda_z2_Z
lambda_x3z2_XZ == lambda_x3_X * lambda_z2_Z
lambda_x1z3_XZ == lambda_x1_X * lambda_z3_Z
lambda_x2z3_XZ == lambda_x2_X * lambda_z3_Z
lambda_x3z3_XZ == lambda_x3_X * lambda_z3_Z
# XZ ~ Mean_XZ*1
# Mean_XZ == (Cov_X_Z)
"

# for now it is sligthly faster
est_lav <- time(lavaan::sem(syntax, data))
est_cav <- time(sem(syntax, data))
