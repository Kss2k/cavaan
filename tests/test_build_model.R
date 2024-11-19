devtools::load_all()

tpb <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Covariances
  ATT ~~ SN + PBC
  PBC ~~ SN 
  # Causal Relationsships
  INT ~ a * ATT + b * SN + c * PBC
  BEH ~ INT + PBC 
'

parTable <- cavaanify(tpb)
model <- sem(tpb, data=TPB)
model$start
start <- parTable.d[parTable.d$free, 'est']
logLik(model, start)
est <- nlminb(start, logLik, model=model)
fit <- fillModel(model, est$par)

filled <- fillModel(model, model$start)
matrices <- filled$models[[1]]$matrices
matrices$Sigma
tpb <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Covariances
  ATT ~~ SN + PBC
  PBC ~~ SN 
  # Causal Relationsships
  INT ~ a * ATT + b * SN + c * PBC
  BEH ~ INT + PBC 
  BEH ~ c(d1, d2) * income

  a.b := a * b * d
'

parTable <- cavaanify(tpb, groups=c(1, 2))
