#' Fit a Structural Equation Model (SEM)
#'
#' This function fits a structural equation model (SEM) using a specified model syntax and dataset.
#'
#' @param syntax A character string specifying the SEM model syntax.
#' @param data A data frame containing the observed variables used in the model.
#' @param group (Optional) A character string specifying the grouping variable in the data.
#'              If provided, the function performs multi-group SEM.
#' @param start (Optional) A numeric vector of starting values for the optimization. If `NULL`, default values are used.
#' @param num.grad Logical, whether to use numerical gradients for the optimization (`TRUE`) or analytical gradients (`FALSE`).
#'
#' @return A list of class `"cavaan"` containing the fitted model, parameter estimates, and additional details:
#' \describe{
#'   \item{`coef`}{Estimated parameters of the SEM model.}
#'   \item{`nlminb`}{Optimization results from the `nlminb` function.}
#'   \item{`parTable.d`}{Cleaned parameter table with estimates.}
#' }
#'
#' @details
#' - The function supports single-group and multi-group SEM. 
#' - Model syntax is parsed using `cavaanify()` and a model is built with `build_model()`.
#' - The optimization is performed using `nlminb` with numerical or analytical gradients.
#'
#' @examples
#
#' 
#' tpb <- "
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN  =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#' 
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC
#' "
#' 
#' fit2 <- lavaan::sem(tpb, TPB)
#' fit <- sem(tpb, data = TPB, num.grad = TRUE)
#' fit
#' fit <- sem(tpb, data = TPB, num.grad = FALSE)
#' fit
#' 
#' tpb_mean <- "
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#' 
#'   BEH ~1
#'   INT ~1
#'   PBC ~1
#'   ATT ~1
#'   SN ~1
#' 
#'   b1   ~ 0*1
#'   int1 ~ 0*1
#'   pbc1 ~ 0*1
#'   att1 ~ 0*1
#'   sn1  ~ 0*1
#' 
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   # Causal Relationsships
#'   INT ~ a * ATT + b * SN + c * PBC
#'   BEH ~ d * INT + e * PBC
#' "
#' 
#' fit2 <- lavaan::sem(tpb_mean, TPB)
#' matrices <- lavaan::lavInspect(fit2, "estimates")
#' fit <- sem(tpb_mean, data = TPB, num.grad = TRUE)
#' fit
#' fit <- sem(tpb_mean, data = TPB, num.grad = FALSE)
#' fit
#' 
#' 
#' tpb_mean_ov <- "
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#' 
#'   b1   ~ 1
#'   int1 ~ 1
#'   pbc1 ~ 1
#'   att1 ~ 1
#'   sn1  ~ 1
#' 
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   # Causal Relationsships
#'   INT ~ a * ATT + b * SN + c * PBC
#'   BEH ~ INT + PBC
#' "
#' 
#' fit <- sem(tpb_mean_ov, data = TPB, num.grad = TRUE)
#' fit
#' fit <- sem(tpb_mean_ov, data = TPB, num.grad = FALSE)
#' fit
#' @export
sem <- function(syntax, data, group=NULL, start=NULL, num.grad=TRUE) {
  data   <- as.data.frame(data)
  
  if (!is.null(group)) {
    checkGroupData(group=group, data=data)
    groups        <- unique(data[[group]])
    valuesGroup   <- data[[group]]
    data$group    <- as.numeric(as.factor(valuesGroup))
    data[[group]] <- NULL

  } else {
    groups     <- 1
    data$group <- 1
  }

  parTable <- cavaanify(syntax, groups=groups)
  oVs      <- getOVs(parTable)
  checkOVsData(oVs, data)

  model      <- build_model(parTable, data=data)
  parTable.d <- model$parTable.d 

  # PUT IN A SEPERATE FUNCTION -------------------------
  start  <- if (is.null(start)) model$start else start
  isfree <- parTable.d$free
  isVar  <- parTable.d$op[isfree] == "~~" & 
    parTable.d$lhs[isfree] == parTable.d$rhs[isfree]
  upper  <- rep(Inf, length(start)) 
  lower  <- ifelse(isVar, yes=0, no=-Inf)

  # ------------------------------------------------------
  RcppModel <- createRcppModel(model)
  # debugCppModel(RcppModel, start)
  gradient  <- getGradientFunction(num.grad, parTable.d=parTable.d)

  analytic.lv.g <- structure(as.vector(gradLogLikCppLVMeans(start, RcppModel)), names = names(start))
  analytic.ov.g <- structure(as.vector(gradLogLikCppOVMeans(start, RcppModel)), names = names(start))
  numeric.g <- structure(gradLogLikNumericCpp(start, RcppModel), names = names(start))
  browser()
  est       <- #suppressWarnings(
    stats::nlminb(start, objective=logLikCpp, xptr=RcppModel, lower=lower, upper=upper,
                  gradient=gradient)
  #)
 
  par            <- est$par
  model.f        <- fillModel(model, par)
  model.f$coef   <- par
  model.f$nlminb <- est
  model.f$parTable.d <- cleanParTable.d(model.f$parTable.d, theta=par, model=model.f)
  
  class(model.f) <- "cavaan"
  model.f
}
