fillModel <- function(model, params) {
  parTable.d <- model$parTable.d
  parTable.d[parTable.d$free, "est"] <- params
  for (i in seq_len(NROW(parTable.d[!parTable.d$free, ]))) {
    row <- parTable.d[i, , drop=FALSE]
    replacement <- parTable.d[parTable.d$label == row$label, "est"]
    if (NROW(replacement)) row$est <- replacement
  }

  for (i in seq_len(NROW(parTable.d))) {
    row <- parTable.d[i, , drop=FALSE]
    matrix <- row$matrix.label
    rowi   <- row$row
    colj   <- row$col
    if (is.na(rowi) || is.na(colj)) next
    group  <- ifelse(row$group=="1", yes=1, no=row$group)
    model$models[[group]]$matrices[[matrix]][rowi, colj] <- parTable.d[i, "est"]
  }
  
  model
}
