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

  lhs     <- row$lhs
  op      <- row$op
  rhs     <- row$rhs
  label   <- row$label
  etas    <- model$info$etas
  mVYs    <- model$info$mVYs
  Y       <- c(etas, mVYs)
  free    <- row$const == "" & row$label == ""
  isconst <- row$const != ""
  if (op == "~") {
    r     <- lhs
    c     <- rhs
    est   <- ifelse(row$const == "", yes=stats::runif(1), no=const2num(row$const))
    isEta <- rhs %in% Y
  } else if (op == "=~") {
    r     <- rhs
    c     <- lhs
    est   <- ifelse(row$const == "", yes=stats::runif(1), no=const2num(row$const))
    isEta <- lhs %in% etas
  } else {
    stop2("unrecoginized operator: ", row$op)
  }

  M        <- if (isEta) B else G
  i        <- which(rownames(M) == r)
  j        <- which(colnames(M) == c)
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
             continue=continue, group=group, isEquation=isEquation, 
             isconst=isconst)
}


buildParamPhi <- function(models, row) {
  group <- ifelse(row$group == "", yes=1, no=row$group)
  model <- models[[group]]
  xis   <- models[[group]]$info$xis

  G <- model$matrices$Phi

  lhs   <- row$lhs
  op    <- row$op
  rhs   <- row$rhs
  label <- row$label

  r <- lhs
  c <- rhs

  i       <- which(rownames(G) == r)
  j       <- which(colnames(G) == c)
  label   <- ifelse(label != "", label, sprintf("%s%s%s", lhs, op, rhs))
  fill    <- TRUE
  free    <- row$const == "" & row$label == ""
  isconst <- row$const != ""
  stopif(!validRowColMatch(i=i, j=j), "row and col must be unique")

  if (j != i) {
    rows     <- c(i, j)
    cols     <- c(j, i)
    est      <- ifelse(row$const == "", yes=0, no=const2num(row$const))
    free     <- c(free, FALSE)
    continue <- c(FALSE, TRUE)
    est      <- c(est, est)
  } else {
    rows     <- i
    cols     <- j
    est      <- ifelse(row$const != "", yes=const2num(row$const),
                       no=ifelse(lhs %in% xis, yes=1, no=stats::runif(1) + 0.2))
    continue <- FALSE
  }
  
  isEquation <- FALSE
  data.frame(lhs=lhs, op=op, rhs=rhs, est=est,
             label=label, row=rows, col=cols, matrix=2,
             matrix.label="Phi", free=free, fill=fill,
             continue=continue, group=group, isEquation=isEquation,
             isconst=isconst)
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

  free       <- row$const == "" & row$label == ""
  isconst    <- row$const != ""
  continue   <- FALSE
  isEquation <- FALSE
  isOv       <- r %in% model$info$oVs
  est        <- ifelse(row$const != "", yes=const2num(row$const),
                       no=ifelse(isOv, yes=models[[group]]$matrices$Nu[r, 1],
                                 no=stats::runif(1)))

  data.frame(lhs=lhs, op=op, rhs="", est=est, label=label, 
             row=i, col=1, # there is only one column
             matrix=3, matrix.label="Tau", free=free, fill=fill,
             continue=continue, group=group, isEquation=isEquation,
             isconst=isconst)
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
             continue=continue, group=group, isEquation=isEquation,
             isconst=FALSE)
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
    row <- data.frame(lhs=lV, op="~1", rhs=NA, const="0", label="", func="",
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

  labels <- unique(parTable$label)
  for (label in labels) {
    rows <- parTable$label == label
    if (any(rows & parTable$lhs == label & parTable$op %in% LARGE_MATH_OPS)) {
      parTable[rows, "free"]     <- FALSE
      parTable[rows, "continue"] <- TRUE
      parTable[rows & parTable$lhs == label, "continue"] <- FALSE

    } else if (NROW(rows) > 1) {
      parTable[rows, "free"]        <- FALSE
      parTable[rows, "free"][1]     <- TRUE
      parTable[rows, "continue"]    <- TRUE
      parTable[rows, "continue"][1] <- FALSE
    } else if (NROW(rows) && !any(rows & parTable$isconst)) {
      parTable[rows, "free"] <- TRUE
    }
  }

  if (fix.first) {
    for (lV in lVs) {
      if (!all(parTable[parTable$lhs == lV & parTable$op == "=~", "free"])) next
      parTable[parTable$lhs == lV & parTable$op == "=~", "free"][1]    <- FALSE
      parTable[parTable$lhs == lV & parTable$op == "=~", "isconst"][1] <- TRUE
      parTable[parTable$lhs == lV & parTable$op == "=~", "est"][1]  <- 1
    }
  }

  parTable[parTable$isconst, "free"] <- FALSE
  # parTable[parTable$free | parTable$isconst, "continue"] <- FALSE
  # garuantee failure if trying to fill without initializing first!
  parTable[!parTable$free & !parTable$isconst, "est"] <- NA 

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
    row <- partiallySorted[partiallySorted$label == label, , drop=FALSE]
    sortedParTable <- rbind(sortedParTable, row)
  }

  sortedParTable
}


getStartingParams <- function(parTable.d) {
  theta <- parTable.d[parTable.d$free, "est"]
  names(theta) <- parTable.d[parTable.d$free, "label"]
  theta
}
