sem <- function(syntax, data, group=NULL) {
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
  
  model <- build_model(parTable, data=data)

  est <- nlminb(model$start, logLik, model=model)
  
  model.f <- fillModel(model, est$par)
  parTable.d <- model.f$parTable.d  
  parTable.d[model.f$parTable.d$free, 'est'] <- est$par
  parTable.d[parTable.d$matrix.label == "BStar" &
             parTable.d$op == "~", "est"] <- 
      - parTable.d[parTable.d$matrix.label == "BStar" &
                   parTable.d$op == "~", "est"]
  model.f$parTable.d <- parTable.d
  
  class(model.f) <- "cavaan"
  model.f
}
