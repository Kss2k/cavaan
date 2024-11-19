checkOVsData <- function(oVs, data) {
  stopif(!is.data.frame(data) && !is.matrix(data), "data must be a data.frame or a matrix")
  stopif(!all(oVs %in% colnames(data)), "missing observed variables in data: ",
          oVs[!oVs %in% colnames(data)])
}


checkGroupData <- function(group, data) {
  stopif(!is.data.frame(data) && !is.matrix(data), "data must be a data.frame or a matrix")
  stopif(!group %in% colnames(data), "group must be a column in data")
}
