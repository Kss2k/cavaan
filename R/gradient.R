getGradientFunction <- function(num.grad, parTable.d, lVs = getLVs(parTable.d)) {
  if (num.grad) return(gradLogLikNumericCpp)

  intercepts   <- parTable.d$op == "~1"
  lVIntercepts <- intercepts & parTable.d$lhs %in% lVs

  if (!any(intercepts)) gradLogLikCppNoMeans
  else if (!any(lVIntercepts)) gradLogLikCppOVMeans
  else gradLogLikCppLVMeans
}
