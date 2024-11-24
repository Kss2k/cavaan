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
         name      = grepl("[[:alpha:]]|[[:digit:]]|_|\\.", char),
         numeric   = grepl("[[:alpha:]]", char),
         simple.op = FALSE,
         large.op  = paste0(token, char) %in% LARGE_OPS,
         comment   = char != "\n",
         newline   = FALSE,
         closure   = lastChar(token) != ")",
         stop2("unrecoginized token type: ", tokenT)
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
  empty    <- list(lhs=NULL, op=NULL, rhs=NULL, label=NULL, 
                   const=NULL, func=NULL, closure=NULL)
  parTable <- NULL
  current  <- NULL 
  pos      <- "lhs"
  state    <- "open"
  exprType <- "standard"
  labelMod <- constMod <- closureMod <- funcMod <- ""
  skipNext <- FALSE

  for (i in seq_along(tokens)) {
    if (skipNext) {skipNext <- FALSE; next}

    token       <- tokens[i]
    tokenT      <- getTokenType(token)
    nextToken   <- getNextToken(tokens, i=i, N=length(tokens))
    nextTokenT  <- getTokenType(nextToken)
    nextToken2  <- getNextToken(tokens, i=i+1, N=length(tokens))
    nextToken2T <- getTokenType(nextToken2)

    if (exprType == "math") {
      if (token == "\n" && state == "closed") {
        current$label   <- c(current$label, "")
        current$const   <- c(current$const, "")
        current$func    <- c(current$func,  "")
        current$closure <- c(current$closure, "")

        parTable <- rbind(parTable, tokenExprsToParTable(current))
        exprType <- "standard" 
        pos      <- "lhs"
        state    <- "open"
        current  <- NULL 
         
      } else {
        mathExpr <- paste0(current[[pos]], ifelse(token == "\n", yes="", no=token))
        current[[pos]] <- mathExpr
        state <- ifelse(tokenT == "simple.op", yes="open", no="closed")
      }

    } else if (!is.null(nextToken) && nextToken == "*" && pos == "rhs") {
      if      (tokenT == "name")    labelMod <- token
      else if (tokenT == "numeric") constMod <- token
      skipNext <- TRUE

    } else if (tokenT == "name" && !is.null(nextToken) && !is.null(nextToken2T) && 
               nextTokenT == "closure" && nextToken2 == "*") {
      funcMod    <- token
      closureMod <- nextToken
      skipNext   <- TRUE
  
    } else if (!is.null(nextToken) && nextToken == "*" && pos == "lhs") {
      stop2("* not allowed for lhs") 
      
    } else if (tokenT %in% c("name", "numeric")) {
      if (state == "closed") stop2("Unexpected token: ", token)
      if (is.null(current)) current <- empty
      
      state          <- "closed"
      current[[pos]] <- c(current[[pos]], token)

      if (pos == "rhs") {
        current$label   <- c(current$label, labelMod)
        current$const   <- c(current$const, constMod)
        current$func    <- c(current$func, funcMod)
        current$closure <- c(current$closure, closureMod)
        labelMod <- constMod <- closureMod <- funcMod <- ""
      }

    } else if (tokenT == "large.op") {
      current$op <- token
      pos        <- "rhs"
      state      <- "open"
      exprType   <- ifelse(token %in% LARGE_MATH_OPS, yes="math", no="standard")

    } else if (token == "+") {
      if (state == "open") stop2("Unexpected `+`")
      state <- "open"

    } else if (tokenT == "newline" && state == "closed") {
      parTable <- rbind(parTable, tokenExprsToParTable(current))
      current  <- NULL 
      pos      <- "lhs"
      state    <- "open"

    } else stopif(i >= length(tokens) && !is.null(current), "Unexpected end of input")
  }

  parTable
}


parseClosure <- function(string) {
  if (is.na(string) || string == "") return(NULL)
  out <- stringr::str_split_1(string, pattern="\\(|,|\\)") |>
    stringr::str_remove_all(pattern=" ") |> (\(x) x[x != ""])()
  if (!length(out)) return(NULL)
  out
}


parseClosures <- function(strings) {
  nelems <- unique(countElemsClosure(strings[strings != ""]))
  stopif(length(nelems) > 1, "Number of elements in c() functions do not match")

  labels <- matrix(NA, nrow=length(strings), ncol=nelems)
  
  for (i in seq_along(strings)) {
    string <- strings[i]
    elems  <- parseClosure(string)
    if (is.null(elems)) next

    labels[i, seq_along(elems)] <- elems
  
  }

  labels
}


countElemsClosure <- function(string) {
  vapply(string, FUN.VALUE=numeric(1L), FUN=\(s) length(parseClosure(s)))
}
