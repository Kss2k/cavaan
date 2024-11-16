devtools::load_all()
m1 <- '
  # Outer Model
  X =~ x1 + x2 +x3
  Y =~ y1 + y2 + y3
  Z =~ z1 + z2 + z3
  
  # Inner model
  Y ~ a*X + Z
'

tokens <- tokenizer(m1)
parseTokens(tokens)
