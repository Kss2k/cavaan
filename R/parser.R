LARGE_MATH_OPS <- c(":=", ">=", "<=", "==")
LARGE_OPS <- c("~", "~~", "=~", ":=", ">=", "<=", "==")


tokenizeSyntax <- function(string) {
  chars <- stringr::str_split_1(string, "")

  tokens <- character(0L)
  token  <- NULL
  tokenT <- NULL
  for (char in chars) {
    if (fitsToken(char=char, token=token, tokenT=tokenT)) {
      token <- paste0(token, char)

    } else if (char == " ") {
      if (!is.null(tokenT)) tokens  <- c(tokens, token)
      token  <- NULL
      tokenT <- NULL

    } else {
      tokens <- c(tokens, token)
      token  <- char
      tokenT <- getTokenType(char=char)
    }
  }
  if (last(tokens) != "\n") c(tokens, "\n")
  else tokens
}


getTokenType <- function(char) { 
  if (is.null(char)) return(NULL)
  char <- stringr::str_split_1(char, "")[1]
  if      (grepl("[[:alpha:]]", char)) "name"
  else if (grepl("[[:digit:]]", char)) "numeric"
  else switch(char,
              "+" = "simple.op", "*" = "simple.op", 
              "/" = "simple.op", "^" = "simple.op", 
              "-" = "simple.op", 
              "~" = "large.op", "=" = "large.op",
              "<" = "large.op", ">" = "large.op",
              ":" = "large.op", "(" = "closure", 
              "#" = "comment", "\n" = "newline")
}


last <- function(x) {
  x[length(x)]
}


lastChar <- function(string) {
  if (is.null(string)) return(NULL)
  last(stringr::str_split_1(string, ""))
}


fitsToken <- function(char, token, tokenT) {
  if (is.null(tokenT)) return(FALSE)
  switch(tokenT, 
         name      = grepl("[[:alpha:]]|[[:digit:]]", char),
         numeric   = grepl("[[:alpha:]]", char),
         simple.op = FALSE,
         large.op  = paste0(token, char) %in% LARGE_OPS,
         comment   = char != "\n",
         newline   = FALSE,
         closure   = lastChar(token) != ")",
         stop("unrecoginized token type: ", tokenT)
  ) 
}


getNextToken <- function(tokens, i, N) {
  if (i + 1 > N) return(NULL)
  tokens[i + 1]
}


tokenExprsToParTable <- function(exprs) {
  parTable <- NULL

  for (lhs in exprs$lhs) {
    row <- data.frame(lhs=lhs, op=exprs$op, rhs=exprs$rhs, 
                      label=exprs$label, const=exprs$const, 
                      func=exprs$func, closure=exprs$closure)
    parTable <- rbind(parTable, row)
  }

  parTable
}


parseTokens <- function(tokens) {
  parTable <- NULL
  # parTable <- data.frame(lhs=NULL, op=NULL, rhs=NULL, label=NULL, numeric=NULL)
  emptyTokenExprs   <- list(lhs=NULL, op=NULL, rhs=NULL, label=NULL, const=NULL,
                            func=NULL, closure=NULL)
  currentTokenExprs <- NULL 

  position <- "lhs"
  state    <- "open"
  exprType <- "standard"
  labelMod <- constMod <- closureMod <- funcMod <- NA
  skipNext <- FALSE

  for (i in seq_along(tokens)) {
    if (skipNext) {
      skipNext <- FALSE
      next
    }

    token       <- tokens[i]
    tokenT      <- getTokenType(token)
    nextToken   <- getNextToken(tokens, i=i, N=length(tokens))
    nextTokenT  <- getTokenType(nextToken)
    nextToken2  <- getNextToken(tokens, i=i+1, N=length(tokens))
    nextToken2T <- getTokenType(nextToken2)

    if (exprType == "math") {
      if (token == "\n" && state == "closed") {
        currentTokenExprs$label   <- c(currentTokenExprs$label, NA)
        currentTokenExprs$const   <- c(currentTokenExprs$const, NA)
        currentTokenExprs$func    <- c(currentTokenExprs$func,  NA)
        currentTokenExprs$closure <- c(currentTokenExprs$closure, NA)

        parTable <- rbind(parTable, tokenExprsToParTable(currentTokenExprs))
        exprType <- "standard" 
        position <- "lhs"
        state    <- "open"
        currentTokenExprs <- NULL 
          
      } else {
        mathExpr <- paste0(currentTokenExprs[[position]], 
                           ifelse(token == "\n", yes="", no=token))
        currentTokenExprs[[position]] <- mathExpr
        state <- ifelse(tokenT == "simple.op", yes="open", no="closed")
      }

    } else if (!is.null(nextToken) && nextToken == "*" && position == "rhs") {
      if (tokenT == "name") 
        labelMod <- token
      else if (tokenT == "numeric") 
        constMod <- token
      skipNext <- TRUE

    } else if (!is.null(nextToken) && !is.null(nextToken2T) && 
               nextTokenT == "closure" && nextToken2 == "*") {
      funcMod    <- token
      closureMod <- nextToken
      skipNext   <- TRUE
  
    } else if (!is.null(nextToken) && nextToken == "*" && position == "lhs") {
      stop("* not allowed for lhs") 
      
    } else if (tokenT %in% c("name", "numeric")) {
      if (state == "closed") stop("Unexpected token: ", token)
      if (is.null(currentTokenExprs)) currentTokenExprs <- emptyTokenExprs
      state <- "closed"

      currentTokenExprs[[position]] <- c(currentTokenExprs[[position]], token)
      if (position == "rhs") {
        currentTokenExprs$label    <- c(currentTokenExprs$label, labelMod)
        currentTokenExprs$const    <- c(currentTokenExprs$const, constMod)
        currentTokenExprs$func     <- c(currentTokenExprs$func, funcMod)
        currentTokenExprs$closure  <- c(currentTokenExprs$closure, closureMod)
        labelMod <- constMod <- closureMod <- funcMod <- NA
      }

    } else if (tokenT == "large.op") {
      currentTokenExprs$op <- token
      position <- "rhs"
      state    <- "open"
      exprType <- ifelse(token %in% LARGE_MATH_OPS, yes="math", no="standard")

    } else if (token == "+") {
      if (state == "open") stop("Unexpected `+`")
      state <- "open"

    } else if (tokenT == "newline" && state == "closed") {
      parTable          <- rbind(parTable, tokenExprsToParTable(currentTokenExprs))
      currentTokenExprs <- NULL 
      position   <- "lhs"
      state      <- "open"

    } else if (i >= length(tokens) && is.null(currentTokenExprs)) {
      stop("Unexpected end of input")
    }
  }

  parTable
}


cavaanify <- function(syntax) {
  parseTokens(tokenizeSyntax(syntax))
}
