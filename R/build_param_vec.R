buildParamVec <- function(models, parTable) {
  for (i in seq_len(NROW(parTable))) {
    row <- parTable[i, , drop=FALSE]

        
  }
}


buildParamGamma <- function(models, row) {
  model  <- models[[row$group]]
  G <- model$matrices$GammaStar

  if (row$op == "~") {
    r <- row$lhs
    c <- row$rhs
  } else if (row$op == "=~") {
    r <- row$lhs
    c <- row$lhs
  } else {
    stop2("unrecoginized operator: ", row$op)
  }

  row  <- which(rownames(G) == r)
  col  <- which(colnames(G) == c)
  free <- TRUE

  
}


