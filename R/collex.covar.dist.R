## in goes data frame with...
# if raw = T: 1 obs  per line w/ 3 cols for 3 conditions
# if raw = F: 1 type per line w/ 3 cols for 3 conditions + freqs in col 4
# all = F removes configurations not present
collex.covar.mult <- function(x, am = "poiss", raw = TRUE, all = FALSE, configs = 3) {

  # Check if data frame and convert to one
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }

  # If raw obs are inserted, sum the counts
  if (raw == TRUE) {
    d <- as.data.frame(table(x))
  }

  # If only attested combinations are to be calculated, reduce
  if (all == FALSE) {
    d <- subset(d, d[, ncol(d)] > 0)
  }

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
  if (am == "poiss") {
    d$am <- apply(d[, (ncol(d)-1):ncol(d)], 1, poiss.stirl)
  }
  if (am == "pmi") {
    d$am <- apply(d[, (ncol(d)-1):ncol(d)], 1, pmi)
  }
  if (am == "tmi") {
    d$am <- apply(d[, (ncol(d)-1):ncol(d)], 1, tmi, n = n)
  }

  names(d)[(ncol(d)-2):ncol(d)] <- c("OBS", "EXP", toupper(am))
  d <- d[with(d, order(d[ncol(d)], decreasing = T)), ]
  rownames(d) <- NULL
  d

}

logw0 <- function(x, base = exp(1)) {
  ifelse (x == 0, 0, log(x, base))
}

poiss.stirl <- function(x) {
  x[1] * ( logw0(x[1]) - logw0(x[2]) - 1)
}

pmi <- function(x) {
  log2(x[1]/x[2])
}

# http://dspace.uib.no/bitstream/handle/1956/11033/lyse-andersen-mwe-final.pdf?sequence=1&isAllowed=y
# True Mutual Information
tmi <- function(x, n) {
  x[1] <- ifelse(x[1] == 0, 0.001, x[1])
  (x[1]/n) * (log2(x[1]/x[2]))
}

