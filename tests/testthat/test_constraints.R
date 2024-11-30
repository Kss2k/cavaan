devtools::load_all()
library(modsem)
library(rbenchmark)

m1 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  
  # Inner model
  Y ~ X + Z + X:Z 
'

data   <- get_pi_data(m1, method="ca", data=oneInt)
syntax <- get_pi_syntax(m1, method="ca")

# benchmark(lavaan::sem(syntax, data), replications=10)
#>        test replications elapsed relative user.self sys.self ...
#> lavaan::sem           10 230.962        1   230.923     0.06 ...
# benchmark(cavaan::sem(syntax, data, num.grad=true), replications=10)
#>        test replications elapsed relative user.self sys.self ...
#> cavaan::sem           10  68.511        1    68.516        0 ...
time(est <- cavaan::sem(syntax, data, num.grad=TRUE))
time(est_lav <- lavaan::sem(syntax, data))
