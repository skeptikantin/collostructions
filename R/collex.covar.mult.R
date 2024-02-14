#' @export
collex.covar.mult <- function(x, am = "t", raw = TRUE, all = FALSE, reverse = FALSE, threshold = 1, decimals = 5) {

  # checks if sensible argument for association measure has been entered
  na.am <- match(am, c("pmi", "poiss", "t", "tmi", "z", "random"))

  # checks validity of threshold argument
  if (threshold <= 0 || !is.numeric(threshold)) {
    stop("Invalid 'threshold' argument. Must be a positive number.")
  }

  # Check if data frame and convert to one
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }

  # Checks if input df is the result of some tidyverse something-something
  if (any(class(x) %in% c("tbl_df", "tbl", "grouped_df"))) {
    x <- as.data.frame(x)
  }

  # If raw obs are inserted, sum the counts
  if (raw == TRUE) {
    d <- as.data.frame(table(x))
  } else {
    d <- x
  }

  # If only attested combinations are to be calculated, reduce
  if (all == FALSE) {
    d <- subset(d, d[, ncol(d)] > 0)
  }

  # If aggregated data are supplied, but all possible combinations are desired:
  if (raw == FALSE & all == TRUE) {
    d <- d[rep(1:nrow(d[,-ncol(d)]), d[,ncol(d)]), -ncol(d)]
    d <- as.data.frame(table(d))
  }

  # Define variable that is otherwise undefined globally
  OBS <- NULL

  # Determine number of observations
  n <- sum(d[, ncol(d)])

  # Determine number of configurations (-1 b/c frequency in last column)
  configs <- ncol(d) - 1

  ### EXPECTED VALUES
  ## snipped taken from scfa() function in package {cfa}

  # Create vector of frequency counts
  counts <- d[, ncol(d)]

  # Creates an empty matrix with n number of rows and as many cols as there are configs
  sums <- matrix(data = NA, nrow = nrow(d), ncol = configs)

  # Creates an empty vector with as many elements as there are configs
  n.levels <- vector(mode = "numeric", length = configs)

  # Loop to populate sums with the counts of each config element
  for (i in 1:configs) {
    if (!is.factor(d[, i]))
      d[, i] <- as.factor(d[, i])
    tbl <- aggregate(counts, d[i], sum)
    n.levels[i] <- nrow(tbl)
    dimnames(tbl)[1] <- tbl[1]
    if (dimnames(tbl)[1][[1]][1] == "")
      dimnames(tbl)[1][[1]][1] <- "NA"
    cidx <- as.character(d[, i])
    cidx[cidx == ""] <- "NA"
    sums[, i] <- tbl[cidx, 2]
  }

  # A vector with expected values
  exp <- apply(sums, 1, prod)/n^(configs - 1)
  d$exp <- exp

  if (am == "locmi") {
    d$am <- apply(d[, (ncol(d)-1):ncol(d)], 1, local.mi)
  }
  if (am == "poiss") {
    d$am <- apply(d[, (ncol(d)-1):ncol(d)], 1, poiss.stirl)
  }
  if (am == "pmi") {
    d$am <- apply(d[, (ncol(d)-1):ncol(d)], 1, pmi)
  }
  if (am == "tmi") {
    d$am <- apply(d[, (ncol(d)-1):ncol(d)], 1, tmi, n = n)
  }
  if (am == "t") {
    d$am <- apply(d[, (ncol(d)-1):ncol(d)], 1, t.val)
  }
  if (am == "z") {
    d$am <- apply(d[, (ncol(d)-1):ncol(d)], 1, z.val)
  }

  # Rename
  names(d)[(ncol(d)-2):ncol(d)] <- c("OBS", "EXP", toupper(am))
  # Sort by am & reset rownames
  if (reverse == FALSE) {
    d <- d[with(d, order(d[, ncol(d)], decreasing = TRUE)), ]
  } else {
    d <- d[with(d, order(d[, ncol(d)])), ]
  }

  rownames(d) <- NULL

  # Round
  d$EXP <- round(d$EXP, decimals)
  d[,ncol(d)] <- round(d[,ncol(d)], decimals)

  # Now subset by threshold
  if (threshold > 1) {
    d <- droplevels(subset(d, OBS >= threshold))
  }

  d

}

