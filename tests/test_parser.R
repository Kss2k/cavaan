devtools::load_all()
m1 <- '
  # Outer Model
  X2 + X =~ x1 + x2 +x3 +
    x4
  Y =~ y1 + y2 + y3
  Z =~ 1*z1 + c(l1, l2) * z2 + z3
  
  # Inner model
  Y ~ a*X + c*Z
'

tokens <- tokenizer(m1)
parseTokens(tokens)
