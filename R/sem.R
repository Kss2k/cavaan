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
  gradient  <- if (num.grad) NULL else gradLogLikNumericCpp # gradLogLikCppLVMeans
  analytic.g <- structure(as.vector(gradLogLikCppLVMeans(start, RcppModel)), names = names(start))
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
