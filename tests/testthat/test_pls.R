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


parTable <- cavaanify(tpb, groups=1)
est_pls <- pls(parTable, 10 * TPB)

# var(y) = l^2 * var(x)

lambda <- est_pls$fit$Lambda
Sigma <- est_pls$fit$Sigma
l1a <- lambda["att1", "ATT"]
va1 <- Sigma["att1", "att1"]
resvar <- va1 - l1a^2 # latent-variable var = 1
varATT <- l1a^2

rescaleVCOV(Sigma, vars="ATT", sigmas=l1a^2)
getParamVecNames(est_pls)
