sem <- function(syntax, data, group=NULL, start=NULL) {
  data   <- as.data.frame(data)
  
  if (!is.null(group)) {
    checkGroupData(group=group, data=data)
    groups        <- unique(data[[group]])
    valuesGroup   <- data[[group]]
    data$group    <- as.numeric(as.factor(valuesGroup))
    data[[group]] <- NULL

  } else {
    group      <- "group"
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

  RcppModel <- createRcppModel(model)
  est <- suppressWarnings(nlminb(start, objective=logLikCpp, gradient=gradLogLikCpp, xptr=RcppModel))
  
  # R and C++/R
  par <- est$par
  model.f <- fillModel(model, par)
  model.f$coef   <- par
  model.f$nlminb <- est

  # ----------------------------------------------------

  # CLEAN THIS UP! -------------------------------------
  parTable.d <- model.f$parTable.d  
  parTable.d[model.f$parTable.d$free, 'est'] <- par
  parTable.d[parTable.d$op == "~", "est"] <- 
      - parTable.d[parTable.d$op == "~", "est"]
  model.f$parTable.d <- parTable.d
  # ----------------------------------------------------
  
  class(model.f) <- "cavaan"
  model.f
}


viewModelCreation <- function(model) {
  ViewModelCreation(fit)
}
