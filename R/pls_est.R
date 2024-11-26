pls <- function(parTable, data, maxIter = 100, consistent = TRUE) {
  if (!is.matrix(data)) data <- as.matrix(data)
  
  model <- initPLS(parTable, data) 
  model <- estimatePLS(model, maxIter = maxIter)
  getFitPLS(model, consistent = consistent)
}


estimatePLS <- function(model, maxIter = 100, standardize = TRUE) {
  model <- step0(model)
  for (i in seq_len(maxIter)) {
    model <- model |> 
      step1() |>
      step2() |>
      step3() |>
      step4() |>
      step5() 
  } 
  model 
}
