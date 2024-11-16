tokenizer <- function(string) {
  chars <- stringr::str_split_1(string, "")

  tokens  <- character(0L)
  token   <- NULL
  tokenT <- NULL
  for (char in chars) {
    if (fitsToken(char=char, token=token, tokenT=tokenT)) 
      token <- paste0(token, char)
    else if (char == " ") {
      if (!is.null(tokenT)) 
        tokens  <- c(tokens, token)
      token   <- NULL
      tokenT <- NULL
    } else {
      tokens  <- c(tokens, token)
      token   <- char
      tokenT <- getTokenType(char=char)
    }
  }

  tokens
}


getTokenType <- function(char) {
  if      (grepl("[[:alpha:]]", char)) "name"
  else if (grepl("[[:digit:]]", char)) "numeric"
  else switch(char,
              "+" = "simple_op", "*" = "simple_op", 
              "~" = "large_op", "=" = "large_op",
              "#" = "comment", "\n" = "newline")
}


fitsToken <- function(char, token, tokenT) {
  if (is.null(tokenT)) return(FALSE)
  switch(tokenT, 
         name      = grepl("[[:alpha:]]|[[:digit:]]", char),
         numeric   = grepl("[[:alpha:]]", char),
         simple_op = FALSE,
         large_op  = paste0(token, char) %in% c("=~", "~~"),
         comment   = char != "\n",
         newline   = FALSE,
         stop("unrecoginized token type: ", tokenT)
  ) 
}


getNextToken <- function(tokens, i, N) {
  if (i + 1 > N) return(NULL)
  tokens[i + 1]
}


parseTokens <- function(tokens) {
  parTable <- NULL
  # parTable <- data.frame(lhs=NULL, op=NULL, rhs=NULL, label=NULL, numeric=NULL)
  emptyRow   <- list(lhs=NA, op=NA, rhs=NA, label=NA, numeric=NA)
  currentRow <- emptyRow

  position <- "lhs"
  for (i in seq_along(tokens)) {
    token      <- tokens[i]
    tokenT     <- getTokenType(token[1])
    nextToken  <- getNextToken(tokens, i=i, N=length(tokens))
    NextTokenT <- getTokenType(token[1])
    if (nextToken == "*" && position == "rhs") {
      if (NextTokenT == "name") 
        currentRow$label <- c(currentRow$label, nextToken)
      else if (NextTokenT == "numeric") 
        currentRow$numeric <- c(currentRow$numeric, nextToken)
    } else if (nextToken == "*") {
      stop("* not allowed for lhs")
    } else if (tokenT %in% c("name", "numeric")) {
      currentRow[[position]] <- c(currentRow[[position]], token)
    } else if (tokenT == "large_op") {
  
    }
    else if (tokenT == "newline") {
      parTable <- c(parTable, currentRow)
      currentRow <- emptyRow
    }
  }
}
