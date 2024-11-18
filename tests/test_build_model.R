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
  BEH ~ d * income

  a.b := a * b * d
'

parTable <- cavaanify(tpb)
model <- build_model(parTable)
parTable.d <- model$parTable.d
fillModel(model, params=parTable.d[parTable.d$free, "est"])
matrices <- model$models[[1]]$matrices
G     <- matrices$G
Gamma <- matrices$Gamma
Phi   <- matrices$Phi
GammaStar <- matrices$GammaStar
BS    <- matrices$BStar
BSi   <- solve(BS)

Sigma <- G %*% BSi %*% GammaStar %*% Phi %*% t(GammaStar) %*% t(BSi) %*% t(G)


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

parTable <- cavaanify(tpb, groups=c("g1", "g2"))
model <- build_model(parTable)

matrices <- model$models[["g2"]]$matrices
G     <- matrices$G
Gamma <- matrices$Gamma
Phi   <- matrices$Phi
GammaStar <- matrices$GammaStar
BS    <- matrices$BStar
BSi   <- solve(BS)

Sigma <- G %*% BSi %*% GammaStar %*% Phi %*% t(GammaStar) %*% t(BSi) %*% t(G)
