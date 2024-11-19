fillModel <- function(model, params, parTable.d=NULL, calc.sigma=TRUE) {
  parTable.d <- if (is.null(parTable.d)) model$parTable.d else parTable.d
  parTable.d[parTable.d$free, "est"] <- params
  parTable.d <- parTable.d[parTable.d$fill, ] # quick and dirty for now
  groups     <- unique(parTable.d$group)

  est <- NA
  for (i in seq_len(NROW(parTable.d))) {
    row    <- parTable.d[i, , drop=FALSE]
    matrix <- row$matrix.label
    rowi   <- row$row
    colj   <- row$col
    group  <- row$group
    
    if (!row$continue) est <- parTable.d[i, "est"]
    if (!row$fill || is.na(rowi) || is.na(colj)) next
  
    model$models[[group]]$matrices[[matrix]][rowi, colj] <- est
  }

  if (calc.sigma) {
    for (group in groups) {
      matrices <- model$models[[group]]$matrices
      G         <- matrices$G
      Gamma     <- matrices$Gamma
      Phi       <- matrices$Phi
      GammaStar <- matrices$GammaStar
      B         <- matrices$B
      BStar     <- matrices$BStar
      BStarInv  <- solve(BStar)

      Sigma <- G %*% BStarInv %*% GammaStar %*% Phi %*% 
        t(GammaStar) %*% t(BStarInv) %*% t(G)

      model$models[[group]]$matrices$Sigma <- Sigma
    }
  }
  
  model
}


getBaseGradients <- function(models, parTable.d) {
  model      <- list(models=models)
  params     <- unique(parTable.d[parTable.d$free, "label"])
  groups     <- unique(parTable.d[parTable.d$fill, "group"])

  parTable.f      <- parTable.d
  parTable.f$free <- FALSE
  parTable.f$fill <- FALSE

  gradients <- namedList(n=length(params), names=params)

  for (param in params) {
    parTable.fp <- parTable.f

    rows <- parTable.fp$label == param
    parTable.fp[!rows, "est"] <- 0
    parTable.fp[rows,  "est"] <- 1
    parTable.fp[rows, "fill"] <- TRUE
    parTable.fp[rows, "free"] <- TRUE

    theta   <- parTable.fp[parTable.fp$free, "est"]
    
    model.f <- fillModel(model, params=theta, parTable.d=parTable.fp, 
                         calc.sigma=FALSE)
                                
    gradients_param <- namedList(n=length(groups), names=groups)
    for (group in groups) {
      gradients_param[[group]] <- model.f$models[[group]]$matrices
    }

    gradients[[param]] <- gradients_param
  }

  gradients
}
