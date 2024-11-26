cavaanify <- function(syntax, groups=1) {
  tokenizeSyntax(syntax) |> 
    parseTokens() |> 
    fixIntercepts() |> 
    pivotGroups(groups) |> 
    getIsExpression()
}


fixIntercepts <- function(parTable) {
  isIntercept <- parTable$op == "~" & parTable$rhs == "1"
  parTable[isIntercept, "rhs"] <- ""
  parTable[isIntercept, "op"]  <- "~1"

  parTable
}


checkGroupLabels <- function(parTable, groups=NULL) {
  stopif(is.null(groups) && any(parTable$func == "c"), 
         "c() function is not allowed in the absence of groups")

  ngroups <- length(groups)
  subParTable <- parTable[parTable$func == "c", ]

  stopif(!all(countElemsClosure(subParTable$closure) == ngroups),
         "Number of groups in c() function does not match the number of groups")
}


pivotGroups <- function(parTable, groups=NULL) {
  checkGroupLabels(parTable, groups) # throws error if group-labels are not valid

  if (is.null(groups) || length(groups) == 1)
    return(cbind(parTable, data.frame(group=1)))

  hasGroupLabel <- parTable$func == "c"
  elemsClosure  <- parseClosures(parTable$closure)
  parTableLong  <- NULL
  for (i in seq_along(groups)) {
    group <- groups[i]
    subParTable <- cbind(parTable, data.frame(group=group))
    subParTable[hasGroupLabel, "label"] <- elemsClosure[hasGroupLabel, i]
    parTableLong <- rbind(parTableLong, subParTable)
  }

  parTableLong
}


getAVs <- function(parTable) { # AVs = All Variables
  subParTable <- parTable[!parTable$op %in% LARGE_MATH_OPS, ]
  aVs <- unique(c(subParTable$lhs, subParTable$rhs))
  aVs[aVs != ""]
}


getLVs <- function(parTable) { # LVs = Latent Variables
  unique(parTable[parTable$op == "=~", "lhs", drop=TRUE])
}


getOVs <- function(parTable) { # OVs = Observed Variables
  aVs <- getAVs(parTable)
  lVs <- getLVs(parTable)
  aVs[!aVs %in% lVs]
}


getIndsLVs <- function(parTable, lVs = getLVs(parTable)) { # Inds = indicators (observed)
  sortedPT <- NULL
  lmabdaLV <- NULL
  for (lV in lVs) {
    lambdaLV <- parTable[parTable$lhs == lV & parTable$op == "=~", ]
    sortedPT <- rbind(sortedPT, lambdaLV)
  }

  unique(sortedPT$rhs)
}


getIndsX <- function(parTable) {
  xis <- getXis(parTable)
  getIndsLVs(parTable, lVs=xis)
}


getIndsY <- function(parTable) {
  etas <- getEtas(parTable)
  getIndsLVs(parTable, lVs=etas)
}


getMVYs <- function(parTable) { # MVs = Manifest Variables, Y = endogenous
  unique(c(getIndsX(parTable), getIndsY(parTable), getODVs(parTable)))
}


getMVXs <- function(parTable) { # MVs = Manifest Variables, X = exogenous
  # a wrapper of getOIVs 
  oVs <- getOVs(parTable)
  iVs <- getIVs(parTable)
  iVs[iVs %in% oVs]
}


getDVs <- function(parTable, sorted=TRUE) { # DVs = Dependent Variables
  if (!sorted) {
    unique(c(parTable[parTable$op == "~", "lhs"],
             parTable[parTable$op == "=~", "rhs"]))
  } else getSortedDVs(parTable)
}


getODVs <- function(parTable) { # Obseved Dependent Variables
  dVs <- getDVs(parTable, sorted=TRUE)
  oVs <- getOVs(parTable)
  dVs[dVs %in% oVs]
}


getIDVs <- function(parTable) { # Independent Dependent Variables
  dVs <- getDVs(parTable, sorted=TRUE)
  oVs <- getOVs(parTable)
  dVs[dVs %in% oVs]
}


getIVs <- function(parTable) { # IVs = Independed Variables
  aVs <- getAVs(parTable)
  dVs <- getDVs(parTable, sorted=FALSE)
  aVs[!aVs %in% dVs]
}


getEtas <- function(parTable, sorted=TRUE) {
  if (!sorted) {
    lVs <- getLVs(parTable)
    dVs <- getDVs(parTable, sorted=FALSE)
    lVs[lVs %in% dVs]
  } else getSortedEtas(parTable)
}


getXis <- function(parTable) {
  lVs  <- getLVs(parTable)
  etas <- getEtas(parTable, sorted=FALSE)
  lVs[!lVs %in% etas]
}


getSortedEtas <- function(parTable) {
  getSortedEtasOrDVs(parTable, isLV=TRUE)
}


getSortedDVs <- function(parTable) {
  getSortedEtasOrDVs(parTable, isLV=FALSE)
}


getSortedEtasOrDVs <- function(parTable, isLV) {
  rows     <- parTable[parTable$op == "~", ]
  sorted   <- character(0L)
  if (isLV) unsorted <- getEtas(parTable, sorted=FALSE)
  else      unsorted <- getDVs(parTable, sorted=FALSE)

  while (length(sorted) < length(unsorted) && NROW(rows) > 0) {
    stopif(all(unique(rows$lhs) %in% rows$rhs), "Model is non-recursive")

    for (i in seq_len(nrow(rows))) {
      if ((x <- rows[i, "lhs"]) %in% rows$rhs) next

      sorted  <- c(x, sorted)
      rows <- rows[!grepl(x, rows$lhs), ]
      break
    }
  }

  if (!all(sorted %in% unsorted) && length(sorted) < length(unsorted)) {
    sorted <- c(sorted, unsorted[!unsorted %in% sorted])
  }

  sorted
}


getIsExpression <- function(parTable) {
  parTable$isEquation <- parTable$op %in% LARGE_MATH_OPS
  parTable
}


cleanParTable.d <- function(parTable.d, theta, model) {
  etas <- model$models[[1]]$info$etas
  
  isNegative <- parTable.d$op == "~" | 
    (parTable.d$op == "~1" & parTable.d$lhs %in% etas)
  isFree     <- parTable.d$free

  parTable.d[isFree, "est"] <- ifelse(isNegative[isFree], yes=-theta, no=theta)
 
  parTable.d <- sortDfBy(parTable.d, x="continue", decreasing=TRUE) |>
    sortDfBy(x="free", decreasing=TRUE) |> sortDfBy(x="lhs") |>
    sortDfBy(x="rhs") |> sortDfBy(x="op")

  parTable.d[!duplicated(parTable.d$label), ]
}
