devtools::load_all()

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

fit2 <- lavaan::sem(tpb, TPB)
fit <- sem(tpb, data=TPB, num.grad=TRUE)
fit

tpb_mean <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

  BEH ~1
  INT ~1
  PBC ~1
  ATT ~1
  SN ~1
  
  b1   ~ 0*1
  int1 ~ 0*1
  pbc1 ~ 0*1
  att1 ~ 0*1
  sn1  ~ 0*1

# Inner Model (Based on Steinmetz et al., 2011)
  # Causal Relationsships
  INT ~ a * ATT + b * SN + c * PBC
  BEH ~ INT + PBC 
'

fit2 <- lavaan::sem(tpb_mean, TPB)
matrices <- lavaan::lavInspect(fit2, 'estimates')
fit <- sem(tpb_mean, data=TPB, num.grad=FALSE)
