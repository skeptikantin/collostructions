#' @export
#' @importFrom stats reshape
reshape.cca <- function(x, cond = "shorter", str.dir = TRUE, value = "COLL.STR",
                        abs.dev = FALSE, max.assoc = FALSE,
                        sorton = "abs.dev", decimals = 5) {

  if (str.dir==TRUE & value != "COLL.STR") {
    stop("\nThe default 'str.dir = TRUE' only makes sense if measure = 'COLL.STR'.\nSet 'str.dir = FALSE' or use a different measure.")
  }

  # checks if sensible argument for measure has been entered
  na.value <- match(value, c("COLL.STR", "OBS"))
  if (is.na(na.value)) {
    stop("Invalid 'measure' argument. See ?reshape.cca().")
  }

  # determine positions
  s1 <- which(names(x)=="SLOT1")
  s2 <- which(names(x)=="SLOT2")
  assoc <- which(names(x)=="ASSOC")
  value <- grep(value, names(x))

  # determine condition and
  if (cond == "shorter") {
    cond = ifelse(length(unique(x[,s1])) > length(unique(x[,s2])), 1, 2)
  } else if (cond %in% c(1:2)) {
    cond = cond
  } else {
    stop("Please enter either '1', '2' or 'shorter' as a condition.")
  }

  # for the column names later, fetch the column names
  cols <- as.character(unique(x[,ifelse(cond == 1, 2, 1)]))

  # create directed association if str.dir == TRUE
  if (str.dir) {
    x[,value] <- ifelse(x[,assoc] == "rep", -x[,value], x[,value])
  }

  # now subset to have only the wanted columns
  y <- subset(x, select = c(s1, s2, value))

  if (cond == 2) {
    y <- reshape(y, idvar = names(x)[s2], timevar = names(x)[s1], direction = "wide")
  } else if (cond == 1) {
    y <- reshape(y, idvar = names(x)[s1], timevar = names(x)[s2], direction = "wide")
  } else {
    stop("Please enter either '1', '2' or 'shorter' as a condition.")
  }

  # calculate the absoulte deviations for every item:
  y$abs.dev <- apply(y[2:ncol(y)], 1, function(x) sum(abs(x), na.rm = TRUE))

  # sorting by user-specified value
  if (sorton == "abs.dev") {
    y <- y[with(y, order(abs.dev, decreasing = T)),]
  } else if (sorton == "alphabetical") {
    y <- y[with(y, order(y[,1], decreasing = F)),]
  } else {
    stop("Please use 'abs.dev' or 'alphabetical' as sorton argument.")
  }

  # set the condition column as rownames
  rownames(y) <- y[,1]
  y <- y[,-1]
  # set the condition 2 as column names
  names(y) <- c(cols, "abs.dev")

  # round table
  y[,2:ncol(y)] <- round(y[,2:ncol(y)], decimals)


  # Should the type of maximum association be displayed?
  if (max.assoc) {
    y$max.assoc <- as.factor(apply(y[,-ncol(y)], 1, function(x) {
      x.p <- which.max(abs(x))
      x.p <- names(x)[x.p]
      x.p}))
  }

  # Should the abs.dev column be displayed in the output
  if (abs.dev == FALSE) {
    y <- subset(y, select = -abs.dev)
  }
  y

}
