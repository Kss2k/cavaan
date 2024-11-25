getDetailedParTable <- function(models, parTable, seed=pi) {
  set.seed(seed)

  groups <- unique(parTable$group)
  groups <- ifelse(all(groups == ""), yes=1, no=groups)
  for (group in groups) {
    aVs <- models[[group]]$info$aVs
    iVs <- models[[group]]$info$iVs
    oVs <- models[[group]]$info$oVs
    lVs <- models[[group]]$info$lVs

    parTable <- getMissingVariances(parTable, aVs=aVs)
    parTable <- getMissingCovariances(parTable, iVs=iVs)
    parTable <- getMissingIntercepts(parTable, oVs=oVs, lVs=lVs)
  }

  parTableY <- NULL
  for (i in seq_len(NROW(parTable))) {
    row <- parTable[i, , drop=FALSE]

    tbl_r <- switch(row$op,
      `~`  = buildParamGammaB(models, row),
      `=~` = buildParamGammaB(models, row),
      `~~` = buildParamPhi(models, row),
      `:=` = buildCustomParam(models, row),
      `==` = buildCustomParam(models, row),
      `~1` = buildParamTau(models, row),
      stop2("unrecoginized operator: ", row$op)
    )

    parTableY <- rbind(parTableY, tbl_r)
  }

  parTableY <- constrainParams(parTableY)

  sortParTable(parTableY)
}


buildParamGammaB <- function(models, row) {
  group <- ifelse(row$group == "", yes=1, no=row$group)
  model <- models[[group]]

  G <- model$matrices$Gamma
  B <- model$matrices$B

  lhs   <- row$lhs
  op    <- row$op
  rhs   <- row$rhs
  label <- row$label
  etas  <- model$info$etas
  mVYs  <- model$info$mVYs
  Y     <- c(etas, mVYs)

  if (op == "~") {
    r     <- lhs
    c     <- rhs
    est   <- ifelse(row$const == "", yes=stats::runif(1), no=as.numeric(row$const))
    isEta <- rhs %in% Y
  } else if (op == "=~") {
    r     <- rhs
    c     <- lhs
    est   <- ifelse(row$const == "", yes=stats::runif(1), no=as.numeric(row$const))
    isEta <- lhs %in% etas
  } else {
    stop2("unrecoginized operator: ", row$op)
  }

  M        <- if (isEta) B else G
  i        <- which(rownames(M) == r)
  j        <- which(colnames(M) == c)
  free     <- ifelse(row$const == "", yes=TRUE, no=FALSE)
  continue <- FALSE
  fill     <- TRUE
  label    <- ifelse(label != "", label, sprintf("%s%s%s", lhs, op, rhs))

  matrix       <- ifelse(isEta, yes=0, no=1)
  matrix.label <- ifelse(isEta, yes="BStar", no="GammaStar")
  isEquation <- FALSE

  stopif(!validRowColMatch(i=i, j=j), "row and col must be unique")
  
  data.frame(lhs=lhs, op=op, rhs=rhs, est=est,
             label=label, row=i, col=j, matrix=matrix,
             matrix.label=matrix.label, free=free, fill=fill,
             continue=continue, group=group, isEquation=isEquation)
}


buildParamPhi <- function(models, row) {
  group <- ifelse(row$group == "", yes=1, no=row$group)
  model <- models[[group]]

  G <- model$matrices$Phi

  lhs   <- row$lhs
  op    <- row$op
  rhs   <- row$rhs
  label <- row$label

  r <- lhs
  c <- rhs

  i     <- which(rownames(G) == r)
  j     <- which(colnames(G) == c)
  label <- ifelse(label != "", label, sprintf("%s%s%s", lhs, op, rhs))
  fill  <- TRUE

  stopif(!validRowColMatch(i=i, j=j), "row and col must be unique")

  if (j != i) {
    rows     <- c(i, j)
    cols     <- c(j, i)
    free     <- ifelse(row$const == "", yes=TRUE, no=FALSE)
    free     <- c(free, FALSE)
    continue <- c(FALSE, TRUE)
    est      <- ifelse(row$const == "", yes=0, no=as.numeric(row$const))
    est      <- c(est, est)
  } else {
    rows     <- i
    cols     <- j
    free     <- ifelse(row$const == "", yes=TRUE, no=FALSE)
    continue <- FALSE
    est      <- ifelse(row$const == "", yes=stats::runif(1), no=as.numeric(row$const))
  }
  
  isEquation <- FALSE

  data.frame(lhs=lhs, op=op, rhs=rhs, est=est,
             label=label, row=rows, col=cols, matrix=2,
             matrix.label="Phi", free=free, fill=fill,
             continue=continue, group=group, isEquation=isEquation)
}


buildParamTau <- function(models, row) {
  group <- ifelse(row$group == "", yes=1, no=row$group)
  model <- models[[group]]

  T <- model$matrices$Tau

  lhs   <- row$lhs
  op    <- row$op
  label <- row$label

  r <- lhs

  i     <- which(rownames(T) == r)
  label <- ifelse(label != "", label, sprintf("%s%s", lhs, op))
  fill  <- TRUE

  stopif(length(i) > 1, "Found multiple matches for intercept row!")

  free       <- ifelse(row$const == "", yes=TRUE, no=FALSE)
  continue   <- FALSE
  isEquation <- FALSE
  isOv       <- r %in% model$info$oVs
  est        <- ifelse(row$const != "", yes=as.numeric(row$const),
                       no=ifelse(isOv, yes=models[[group]]$matrices$Nu[r, 1],
                                 no=stats::runif(1)))
  
  data.frame(lhs=lhs, op=op, rhs="", est=est, label=label, 
             row=i, col=1, # there is only one column
             matrix=3, matrix.label="Tau", free=free, fill=fill,
             continue=continue, group=group, isEquation=isEquation)
}


buildCustomParam <- function(models, row) {
  lhs      <- row$lhs
  op       <- row$op
  rhs      <- row$rhs
  i        <- NA
  j        <- NA
  free     <- FALSE
  continue <- FALSE
  label    <- row$lhs
  group    <- row$group
  fill     <- FALSE

  isEquation <- TRUE

  data.frame(lhs=lhs, op=op, rhs=rhs, est=NA, label=label, row=i, col=j, 
             matrix=NA, matrix.label=NA, free=free, fill=fill,
             continue=continue, group=group, isEquation=isEquation)
}


getMissingVariances <- function(parTable, aVs) {
  group <- unique(parTable$group)
  for (x in aVs) {
    if (any(parTable$lhs == x & parTable$op == "~~" & parTable$rhs == x)) next

    row <- data.frame(lhs=x, op="~~", rhs=x, const="", label="", func="",
                      closure="", group=group, isEquation=FALSE)
    parTable <- rbind(parTable, row)
  }

  parTable
}


getMissingCovariances <- function(parTable, iVs) {
  group <- unique(parTable$group)

  for (i in seq_along(iVs)) {
    for (j in seq_along(iVs)[1:i]) {
      if (i == j) next
      x <- iVs[i]
      y <- iVs[j]
      match <- parTable$lhs == x & parTable$op == "~~" & parTable$rhs == y
      match <- match | (parTable$lhs == y & parTable$op == "~~" & parTable$rhs == x)

      if (any(match)) next

      row <- data.frame(lhs=y, op="~~", rhs=x, const="", label="", func="",
                        closure="", group=group, isEquation=FALSE)
      parTable <- rbind(parTable, row)
    }
  }

  parTable
}


needConstrainedIntercepts <- function(parTable, oV) {
  lVRow <- parTable$rhs == oV & parTable$op == "=~"
  
  if (!any(lVRow)) return(FALSE)

  lVs <- unique(parTable[lVRow, "lhs"])

  lVIntercepts <- parTable[parTable$lhs %in% lVs & parTable$op == "~1", ]
  if (!NROW(lVIntercepts) || any(lVIntercepts$const != "" | 
                                 lVIntercepts$label != "")) return(FALSE)

  inds <- parTable[parTable$op == "=~" & parTable$lhs %in% lVs, "rhs"]

  indIntercepts <- parTable[parTable$lhs %in% inds & parTable$op == "~1", ]
  if (!NROW(indIntercepts)) return(TRUE)

  all(indIntercepts$const == "" & indIntercepts$label == "")
}


getMissingIntercepts <- function(parTable, lVs, oVs) {
  group <- unique(parTable$group)

  hasInterceptsOVs <- any(parTable$lhs %in% oVs & parTable$op == "~1")
  hasInterceptsLVs <- any(parTable$lhs %in% lVs & parTable$op == "~1")

  if (!hasInterceptsOVs & !hasInterceptsLVs) return(parTable)
  
  if (hasInterceptsLVs) for (lV in lVs) {
    if (any(parTable$lhs == lV & parTable$op == "~1")) next

    row <- data.frame(lhs=lV, op="~1", rhs=NA, const="", label="", func="",
                      closure="", group=group, isEquation=FALSE)
    parTable <- rbind(parTable, row)
  }

  for (oV in oVs) {
    if (any(parTable$lhs == oV & parTable$op == "~1")) next
    needConstrained <- needConstrainedIntercepts(parTable, oV)
    const <- ifelse(needConstrained, yes="0", no="")
    row <- data.frame(lhs=oV, op="~1", rhs=NA, const=const, label="", func="",
                      closure="", group=group, isEquation=FALSE)
    parTable <- rbind(parTable, row)
  }

  parTable
}


constrainParams <- function(parTable, fix.first=TRUE) {
  lVs <- getLVs(parTable)

  if (fix.first) {
    for (lV in lVs) {
      parTable[parTable$lhs == lV, "free"][1] <- FALSE
      parTable[parTable$lhs == lV, "est"][1]  <- 1
    }
  }

  labels <- unique(parTable$label)
  for (label in labels) {
    rows <- parTable$label == label
    if (!any(parTable[rows, "free"])) next

    if (NROW(rows) > 1) {
      parTable[rows, "free"]    <- FALSE
      parTable[rows, "free"][1] <- TRUE
    }
  }

  parTable
}


validRowColMatch <- function(i, j) {
  nunique(c(length(i), length(j))) == 1 &&
    unique(c(length(i), length(j))) == 1
}


sortParTable <- function(parTable) {
  partiallySorted <- sortDfBy(parTable, x=free, decreasing=TRUE) |>
    sortDfBy(x=continue, decreasing=FALSE) |>
    sortDfBy(x=isEquation, decreasing=TRUE)

  equations <- parTable[parTable$op %in% LARGE_MATH_OPS, , drop=FALSE]
  unsortedLabels <- unique(c(equations$lhs, parTable$label))

  parTableEqLabels <- NULL
  for (i in seq_len(NROW(equations))) {
    lhs <- equations[i, "lhs"]
    rhs <- getVariablesEquation(equations[i, "rhs"])

    if (!length(rhs)) next
    row <- data.frame(lhs=lhs, op="~", rhs=rhs, label=NA)
    parTableEqLabels <- rbind(parTableEqLabels, row)
  }

  depEqLabels   <- getSortedDVs(parTableEqLabels)
  indepEqLabels <- getIVs(parTableEqLabels)

  sortedEqLabels <- c(indepEqLabels, depEqLabels) 
  sortedLabels <- unique(c(sortedEqLabels, unsortedLabels))

  sortedParTable <- NULL
  for (label in sortedLabels) {
    row <-  partiallySorted[partiallySorted$label == label, , drop=FALSE]
    sortedParTable <- rbind(sortedParTable, row)
  }

  sortedParTable
}


getStartingParams <- function(parTable.d) {
  theta <- parTable.d[parTable.d$free, "est"]
  names(theta) <- parTable.d[parTable.d$free, "label"]
  theta
}
