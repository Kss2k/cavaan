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

fit <- sem(tpb, data=TPB)
fit
ViewModelCreation(fit)



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
#     lhs op  rhs          est      label row col matrix matrix.label  free fill continue group
# 22  INT  ~  ATT  0.178816343          a   1   1      1    GammaStar  TRUE TRUE    FALSE     1
# 1   ATT =~ att1  1.000000000  ATT=~att1   3   1      1    GammaStar FALSE TRUE    FALSE     1
# 2   ATT =~ att2  0.689265295  ATT=~att2   4   1      1    GammaStar  TRUE TRUE    FALSE     1
# 3   ATT =~ att3  0.683198526  ATT=~att3   5   1      1    GammaStar  TRUE TRUE    FALSE     1
# 4   ATT =~ att4  0.676722150  ATT=~att4   6   1      1    GammaStar  TRUE TRUE    FALSE     1
# 5   ATT =~ att5  0.689891807  ATT=~att5   7   1      1    GammaStar  TRUE TRUE    FALSE     1
# 27  ATT ~~  ATT  0.708782907   ATT~~ATT   1   1      2          Phi  TRUE TRUE    FALSE     1
# 18  ATT ~~  PBC -0.230473549   ATT~~PBC   1   3      2          Phi  TRUE TRUE    FALSE     1
# 19  ATT ~~  PBC  0.000000000   ATT~~PBC   3   1      2          Phi FALSE TRUE     TRUE     1
# 16  ATT ~~   SN -0.245606878    ATT~~SN   1   2      2          Phi  TRUE TRUE    FALSE     1
# 17  ATT ~~   SN  0.000000000    ATT~~SN   2   1      2          Phi FALSE TRUE     TRUE     1
# 32 att1 ~~ att1  0.656731349 att1~~att1   6   6      2          Phi  TRUE TRUE    FALSE     1
# 33 att2 ~~ att2  0.653888308 att2~~att2   7   7      2          Phi  TRUE TRUE    FALSE     1
# 34 att3 ~~ att3  0.651774164 att3~~att3   8   8      2          Phi  TRUE TRUE    FALSE     1
# 35 att4 ~~ att4  0.650196054 att4~~att4   9   9      2          Phi  TRUE TRUE    FALSE     1
# 36 att5 ~~ att5  0.654732067 att5~~att5  10  10      2          Phi  TRUE TRUE    FALSE     1
# 23  INT  ~   SN -0.462588002          b   1   2      1    GammaStar  TRUE TRUE    FALSE     1
# 45   b1 ~~   b1  0.663668412     b1~~b1  19  19      2          Phi  TRUE TRUE    FALSE     1
# 46   b2 ~~   b2  0.661010584     b2~~b2  20  20      2          Phi  TRUE TRUE    FALSE     1
# 14  BEH =~   b1  1.000000000    BEH=~b1  16   2      0        BStar FALSE TRUE    FALSE     1
# 15  BEH =~   b2  0.678707963    BEH=~b2  17   2      0        BStar  TRUE TRUE    FALSE     1
# 31  BEH ~~  BEH  0.678337617   BEH~~BEH   5   5      2          Phi  TRUE TRUE    FALSE     1
# 25  BEH  ~  INT -0.002109853    BEH~INT   2   1      0        BStar  TRUE TRUE    FALSE     1
# 26  BEH  ~  PBC -0.261666188    BEH~PBC   2   3      1    GammaStar  TRUE TRUE    FALSE     1
# 24  INT  ~  PBC -0.043084712          c   1   3      1    GammaStar  TRUE TRUE    FALSE     1
# 11  INT =~ int1  1.000000000  INT=~int1  13   1      0        BStar FALSE TRUE    FALSE     1
# 12  INT =~ int2  0.681992062  INT=~int2  14   1      0        BStar  TRUE TRUE    FALSE     1
# 13  INT =~ int3  0.679478728  INT=~int3  15   1      0        BStar  TRUE TRUE    FALSE     1
# 30  INT ~~  INT  0.686252622   INT~~INT   4   4      2          Phi  TRUE TRUE    FALSE     1
# 42 int1 ~~ int1  0.660740200 int1~~int1  16  16      2          Phi  TRUE TRUE    FALSE     1
# 43 int2 ~~ int2  0.658589974 int2~~int2  17  17      2          Phi  TRUE TRUE    FALSE     1
# 44 int3 ~~ int3  0.655606678 int3~~int3  18  18      2          Phi  TRUE TRUE    FALSE     1
# 8   PBC =~ pbc1  1.000000000  PBC=~pbc1  10   3      1    GammaStar FALSE TRUE    FALSE     1
# 9   PBC =~ pbc2  0.687602569  PBC=~pbc2  11   3      1    GammaStar  TRUE TRUE    FALSE     1
# 10  PBC =~ pbc3  0.682290916  PBC=~pbc3  12   3      1    GammaStar  TRUE TRUE    FALSE     1
# 29  PBC ~~  PBC  0.695788025   PBC~~PBC   3   3      2          Phi  TRUE TRUE    FALSE     1
# 20  PBC ~~   SN -0.254938720    PBC~~SN   3   2      2          Phi  TRUE TRUE    FALSE     1
# 21  PBC ~~   SN  0.000000000    PBC~~SN   2   3      2          Phi FALSE TRUE     TRUE     1
# 39 pbc1 ~~ pbc1  0.666751197 pbc1~~pbc1  13  13      2          Phi  TRUE TRUE    FALSE     1
# 40 pbc2 ~~ pbc2  0.664674859 pbc2~~pbc2  14  14      2          Phi  TRUE TRUE    FALSE     1
# 41 pbc3 ~~ pbc3  0.658804581 pbc3~~pbc3  15  15      2          Phi  TRUE TRUE    FALSE     1
# 6    SN =~  sn1  1.000000000    SN=~sn1   8   2      1    GammaStar FALSE TRUE    FALSE     1
# 7    SN =~  sn2  0.680500054    SN=~sn2   9   2      1    GammaStar  TRUE TRUE    FALSE     1
# 28   SN ~~   SN  0.687420483     SN~~SN   2   2      2          Phi  TRUE TRUE    FALSE     1
# 37  sn1 ~~  sn1  0.680155968   sn1~~sn1  11  11      2          Phi  TRUE TRUE    FALSE     1
# 38  sn2 ~~  sn2  0.672631627   sn2~~sn2  12  12      2          Phi  TRUE TRUE    FALSE     1
