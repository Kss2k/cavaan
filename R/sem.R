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

  model
}
