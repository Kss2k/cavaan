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


tpb_strct_ind <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Causal Relationsships
  INT ~ a * ATT + b * SN + c * pbc1 + PBC
  BEH ~ INT + pbc1 + PBC
'
est_pls_strct_ind <- pls(cavaanify(tpb_strct_ind, groups=1), 10 * TPB)


tpb_ov <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Causal Relationsships
  INT ~ a * ATT + b * SN + c * pbc1
  BEH ~ INT + pbc1
'
est_pls_ov <- pls(cavaanify(tpb_ov, groups=1), 10 * TPB)
est_pls_ov

getPathCoefs("INT", c("ATT", "pbc1", "SN"), est_pls_ov$fit$Sigma)
getPathCoefs(c("int1", "int2", "int3"), c("INT"), est_pls_ov$fit$Sigma)
getPathCoefs("INT", c("ATT", "pbc1", "SN"), est_pls_ov$fit$Sigma)
getResidualsVCOV(est_pls_ov$fit$Sigma)


tpb_ov <- ' 
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  # Causal Relationsships
  INT ~ a * ATT + b * SN + c * pbc1
  BEH ~ INT + pbc1
'
est_pls_ov <- pls(cavaanify(tpb_ov, groups=1), TPB)
est_pls_ov

getPathCoefs("INT", c("ATT", "pbc1", "SN"), est_pls_ov$fit$Sigma)
getPathCoefs(c("int1", "int2", "int3"), c("INT"), est_pls_ov$fit$Sigma)
getPathCoefs("INT", c("ATT", "pbc1", "SN"), est_pls_ov$fit$Sigma)
getResidualsVCOV(est_pls_ov$fit$Sigma)
getResidualsLVs(model=est_pls_ov)
getResidualsOVs(model=est_pls_ov)
