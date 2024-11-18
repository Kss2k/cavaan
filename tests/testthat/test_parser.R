devtools::load_all()
m1 <- '
  # Outer Model
  X2 + X =~ x1 + x2 +x3 +
    x4
  Y =~ y1 + y2 + y3
  X =~ x5 + x6
  Z =~ 1*z1 + c(l1, l2) * z2 + z3
  
  # Inner model
  Y ~ a*X + c*Z + u1
  Y ~ 1
  X ~1

  Z ~ 0* 1
  h := 1 * a / 8 ^ 2
  e == 1 * a / 8 ^ 2
  f >= 1 * a / 8 ^ 2
  g <= 1 * a / 8 ^ 2
'

parTable <- cavaanify(m1, groups=c("g1", "g2"))


m1 <- '
  # Outer Model
  X2 + X =~ x1 + x2 +x3 +
    x4
  Y =~ y1 + y2 + y3
  X =~ x5 + x6
  Z =~ 1*z1 + z2 + z3
  
  # Inner model
  Y ~ a*X + c*Z + u1
  Y ~ 1
  X ~1

  Z ~ 0* 1
  h := 1 * a / 8 ^ 2
  e == 1 * a / 8 ^ 2
  f >= 1 * a / 8 ^ 2
  g <= 1 * a / 8 ^ 2
'

parTable <- cavaanify(m1)
