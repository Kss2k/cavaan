# `cavaan`
`cavaan` is an `R`-package for estimating structural equation models (SEMs)
`cavaan` is modelled after the `lavaan` package, but is designed to be faster and more efficient.
This is not a drop-in replacement for `lavaan`, but it is designed to be as similar as possible to `lavaan` in terms of syntax and output.
The main difference is that `cavaan` is designed to be faster and more efficient than `lavaan`, for
models with complicated non-linear equatlity constraints.

# To Install 
To install the latest version from GitHub, use the following code:
```
# Latest version from GitHub
install.packages("devtools")
devtools::install_github("kss2k/cavaan")
```
Here is an example of how to use `cavaan` to estimate a model, and showing an 
example where `cavaan` is faster than `lavaan`. It uses the `modsem` package to generate a model and data,
for estimating an interaction effect between `X` and `Z` on `Y`, using the constrained approach,
which has a lot of non-linear equality constraints.

```{r}
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

benchmark(lavaan::sem(syntax, data), replications=10)
#>        test replications elapsed relative user.self sys.self ...
#> lavaan::sem           10 230.962        1   230.923     0.06 ...
benchmark(cavaan::sem(syntax, data, num.grad=TRUE), replications=10)
#>        test replications elapsed relative user.self sys.self ...
#> cavaan::sem           10  68.511        1    68.516        0 ...
```

# To Do

- [ ] Improve calculation of derivatives of latent intercepts (broken)
- [ ] Improve calculation of derivatives of equatlity constraints: 
      let $`\theta`$ ($`(J \times 1)`$ be the vector of free parameters, and $`f(\theta)`$ ($`(I \times 1)`$ be the vector of equality constraints,
      then the derivative of the equality constraints is $`\frac{\partial f(\theta)}{\partial \theta}`$. 
      the derivative of the equality constraints is a (Jacobian) ($`I \times J`$) matrix, where the element of the $`i`$-th row is 
      the derivative of the $`i`$-th element of $`f(\theta)`$ with respect to the $`j`$-th element of $`\theta`$.
      The derivative of the log-likelihood with respect to $`f(\theta)`$ is $`\frac{\partial F(f(\theta))}{\partial f(\theta})`$.
      The derivative of the log-likelihood with respect to $`\theta`$ is can then be writeen as $`{\frac{\partial f(\theta)}{\partial \theta}}^T\frac{\partial F(f(\theta))}{\partial f(\theta})`$.
