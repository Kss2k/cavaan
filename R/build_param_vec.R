getDetailedParTable <- function(models, parTable, seed=pi) {
  set.seed(seed)

  warning("B is not correctly filled, FIX ME!")

  groups <- unique(parTable$group)
  groups <- ifelse(all(groups == ""), yes=1, no=groups)
  for (group in groups) {
    aVs <- models[[group]]$info$aVs
    iVs <- models[[group]]$info$iVs
    parTable <- getMissingVariances(parTable, aVs)
    parTable <- getMissingCovariances(parTable, iVs)
  }

  parTableY <- NULL
  for (i in seq_len(NROW(parTable))) {
    row <- parTable[i, , drop=FALSE]

    tbl_r <- switch(row$op,
      `~`  = buildParamGamma(models, row),
      `=~` = buildParamGamma(models, row),
      `~~` = buildParamPhi(models, row),
      `:=` = buildCustomParam(models, row),
      stop2("unrecoginized operator: ", row$op)
    )

    parTableY <- rbind(parTableY, tbl_r)
  }

  parTableY <- constrainParams(parTableY)

  parTableY
}


buildParamGamma <- function(models, row) {
  group <- ifelse(row$group == "", yes=1, no=row$group)
  model <- models[[group]]

  G <- model$matrices$Gamma

  lhs   <- row$lhs
  op    <- row$op
  rhs   <- row$rhs
  label <- row$label

  if (op == "~") {
    r <- lhs
    c <- rhs
  } else if (op == "=~") {
    r <- rhs
    c <- lhs
  } else {
    stop2("unrecoginized operator: ", row$op)
  }

  i     <- which(rownames(G) == r)
  j     <- which(colnames(G) == c)
  free  <- TRUE
  label <- ifelse(label != "", label, sprintf("%s%s%s", lhs, op, rhs))

  stopif(nunique(c(length(i), length(j))) != 1, "row and col must be unique")
  
  data.frame(lhs=lhs, op=op, rhs=rhs, est=runif(1), label=label, row=i, col=j, 
             matrix=1, matrix.label="Gamma", free=free, group=group)
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

  stopif(nunique(c(length(i), length(j))) != 1, "row and col must be unique")

  if (j != i) {
    rows <- c(i, j)
    cols <- c(j, i)
    free <- c(TRUE, FALSE)
  } else {
    rows <- i
    cols <- j
    free <- TRUE
  }
  
  data.frame(lhs=lhs, op=op, rhs=rhs, est=runif(1), label=label, row=rows, col=cols, 
             matrix=2, matrix.label="Phi", free=free, group=group)
}


buildCustomParam <- function(models, row) {
  lhs   <- row$lhs
  op    <- row$op
  rhs   <- row$rhs
  i     <- NA
  j     <- NA
  free  <- FALSE
  label <- row$lhs
  group <- row$group

  data.frame(lhs=lhs, op=op, rhs=rhs, est=NA, label=label, row=i, col=j, 
             matrix=NA, matrix.label=NA, free=free, group=group)
}


getMissingVariances <- function(parTable, aVs) {
  group <- unique(parTable$group)
  for (x in aVs) {
    if (!any(parTable$lhs == x & parTable$op == "~~" & parTable$rhs == x)) {
      row <- data.frame(lhs=x, op="~~", rhs=x, const="", label="", func="",
                        closure="", group=group)
      parTable <- rbind(parTable, row)
    }
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

      if (!any(match)) {
        row <- data.frame(lhs=y, op="~~", rhs=x, const="", label="", func="",
                          closure="", group=group)
        parTable <- rbind(parTable, row)
      }
    }
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

  parTable
}
