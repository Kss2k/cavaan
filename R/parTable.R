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
