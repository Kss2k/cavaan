#' @export
print.cavaan <- function(x, ...) {
  cat("NLMINB:", x$nlminb$message, "\n")
  print(x$parTable.d[c("lhs", "op", "rhs", "est", "label")])
}

