#' @export
#' @importFrom stats aggregate
#' @importFrom stats dhyper
#' @importFrom stats runif
collex.covar <- function(x, am = "logl", raw = TRUE, all = FALSE,
                         reverse = FALSE, decimals = 5,
                         str.dir = FALSE, p.fye = FALSE, delta.p = FALSE) {

  if (class(x) == "list") {
    z <- as.data.frame(x)
    names(z) <- names(x[[1]])
    rownames(z) <- NULL
    x <- z
  }
  # Define variables that are otherwise undefined globally
  OBS=EXP=ASSOC=COLL.STR=P.FYE=SIGNIF=NULL

  # Checks if input df is the result of some tidyverse something-something
  if (any(class(x) %in% c("tbl_df", "tbl", "grouped_df"))) {
    x <- as.data.frame(x)
  }

  if (ncol(x) >= 3 & raw == TRUE) {
    stop("\nYour 'x' contains too many columns for raw observations ('raw = TRUE').\nIf you have an aggregated list, set 'raw = FALSE'.")
  }

  # Checks if input has three colums, raw = F, and only col 3 is numeric
  if (ncol(x) == 3 & raw == FALSE) {
    if((is.numeric(x[,2])) & (is.numeric(x[,3]))) {
      stop("\nThe input does not contain two character/factor columns (cols 1 and 2) and one numeric vector (col 3); collex.covar() probably not the right function for the data type.")
    }
  }

  # checks if a two-column in put is correctly set to raw = FALSE
  if (ncol(x) == 2 & raw == FALSE) {
    if (!("numeric") %in% sapply(x, class)) {
      stop("\nThe input contains only two factor/character columns, suggesting raw observations. Try setting 'raw = TRUE' (the default).")
    } else if (any(sapply(x, class) %in% c("numeric", "integer"))) {
    stop("\nThe input contains only two columns, one of which is numeric, suggesting collex.covar() is the wrong choice for this data type.")
    }
  }

  # Define variable that is otherwise undefined globally
  STR.DIR <- NULL
  fCOVAR <- NULL

  if (is.data.frame(x) == F)
    stop("invalid 'df' argument -- must be dataframe")

  # checks if sensible argument for association measure has been entered
  na.am <- match(am, c("logl", "chisq", "cramersV", "dice", "fye", "fye.ln", "gmean",
                       "jaccard", "liddell", "mi", "ms", "mi3", "odds",
                       "pois", "t", "z", "z.cor", "random"))
  if (is.na(na.am))
    stop("Invalid 'am' argument. See ?collex.covar() for available association measures.")

  if (raw == FALSE) {
    if (is.integer(x[ ,3]) == T)
    x[ ,3] <- as.numeric(x[ ,3])
  }

  # all options
  tS1 <- unique(x[, 1]) # types in slot 1
  tS2 <- unique(x[, 2]) # types in slot 2
  attested.types <- x[!duplicated(x), c(1,2)] # types of slot1~slot2 combinations
  names(attested.types) <- c("SLOT1", "SLOT2")

  # Input df can either contain list of raw slot1~slot2 combinations (default)
  # or a frequency list of slot1~slot2 combs (three columns (raw = FALSE))
  if (raw == TRUE) {

    orig.slot.names <- names(x)
    names(x) <- c("SLOT1", "SLOT2")
    N <- length(x[, 1]) # corpus size

    # frequency lists of items in slots
    fS1 <- as.data.frame(table(x[, 1]))
    names(fS1) <- c("SLOT1", "Freq")
    fS2 <- as.data.frame(table(x[, 2]))
    names(fS2) <- c("SLOT2", "Freq")

    # determine frequency list of covar.cxn
    x1 <- x
    x1$FREQ <- 1
    covar.types <- aggregate(FREQ ~ SLOT1 + SLOT2, data = x1, FUN = length)

  } else if (raw == FALSE) {

    orig.slot.names <- names(x)
    names(x) <- c("SLOT1", "SLOT2", "FREQ")
    N <- sum(x[, 3]) # corpus size
    y1 <- x[, c(1,3)]
    y2 <- x[, c(2,3)]
    fS1 <- aggregate(FREQ ~ SLOT1, data = y1, FUN = sum)
    fS2 <- aggregate(FREQ ~ SLOT2, data = y2, FUN = sum)
    covar.types <- x

  } else {

    stop("Invalid specification for 'raw'. Must be TRUE or FALSE.")

  }

  # for negative evidence calculation, determine unattested type combinations
  if (all == TRUE) {
    poss.comb <- expand.grid(tS1, tS2)
    names(poss.comb) <- c("SLOT1", "SLOT2")
    unattested <- setdiff.df(poss.comb, attested.types)
    unattested$FREQ <- 0
    covar.types <- rbind(covar.types, unattested)
    }

  ## obs matrix
  # ads frequency in Slot1
  covar.types <- merge(covar.types, fS1, by.x = "SLOT1", by.y = "SLOT1")
  # adds frequency in Slot2
  covar.types <- merge(covar.types, fS2, by.x = "SLOT2", by.y = "SLOT2")
  names(covar.types) <- c("SLOT2", "SLOT1", "fCOVAR", "fS1", "fS2")
  covar.types <- subset(covar.types, select = c(SLOT1, SLOT2, fCOVAR, fS1, fS2))
  SLOT1 <- covar.types[,1]
  SLOT2 <- covar.types[,2]

  obs <- as.matrix(covar.types[, -c(1:2)])

  obs.val.covar <- function(x) {
    obs <- matrix(c(x[1], x[2] - x[1], x[3] - x[1], N - sum(x[1], x[2] - x[1], x[3] - x[1])), ncol = 4, byrow = F)
  }
  obs <- t(apply(obs, 1, obs.val.covar))
  exp <- t(apply(obs, 1, exp.val))

  res <- cbind(obs, exp)

  if (am == "logl") { res <- cbind(res, apply(res, 1, logl.am)) }
  if (am == "chisq") { res <- cbind(res, apply(res, 1, chi.am)) }
  if (am == "fye") { res <- cbind(res, apply(res, 1, fye.am)) }
  if (am == "fye.ln") { res <- cbind(res, apply(res, 1, fye.ln.am)) }
  if (am == "mi") { res <- cbind(res, apply(res, 1, mi.am)) }
  if (am == "pois") { res <- cbind(res, apply(res, 1, pois.am)) }
  if (am == "z") { res <- cbind(res, apply(res, 1, z.score.am)) }
  if (am == "z.cor") { res <- cbind(res, apply(res, 1, z.score.cor.am)) }
  if (am == "t") { res <- cbind(res, apply(res, 1, t.score.am)) }
  if (am == "odds") { res <- cbind(res, apply(res, 1, odds.am)) }
  if (am == "liddell") { res <- cbind(res, apply(res, 1, liddell.am)) }
  if (am == "ms") { res <- cbind(res, apply(res, 1, ms.am)) }
  if (am == "gmean") { res <- cbind(res, apply(res, 1, gmean.am)) }
  if (am == "dice") { res <- cbind(res, apply(res, 1, dice.am)) }
  if (am == "jaccard") { res <- cbind(res, apply(res, 1, jaccard.am)) }
  if (am == "mi3") { res <- cbind(res, apply(res, 1, mi3.am)) }
  if (am == "cramersV") { res <- cbind(res, apply(res, 1, cramer.am)) }
  if (am == "random") { res <- cbind(res, runif(nrow(res), 0, 1)) }

  res <- cbind(res[,1], covar.types[,4:5], res[, c(5,9)])
  names(res) <- c("OBS", "fS1", "fS2", "EXP", "COLL.STR")


  res$ASSOC <- apply(res[,1:5], 1, association.covar)
  res$STR.DIR <- apply(res[,1:5], 1, direction.covar, am)
  res$SIGNIF <- apply(res[,1:5], 1, sig.level.covar, am)
  res$SLOT1 <- SLOT1
  res$SLOT2 <- SLOT2

  # reorder (for some reason this needs)
  res <- subset(res, select = c(9, 10, 2:3, 1, 4, 6, 5, 7, 8))

  if (am %in% c("fye") & p.fye == TRUE) {
    res$P.FYE <- signif(10^-res$COLL.STR, decimals)
    res <- subset(res, select = c(SLOT1, SLOT2, fS1, fS2, OBS, EXP, ASSOC, COLL.STR, P.FYE, STR.DIR, SIGNIF))
  }

  if (delta.p == TRUE) {
    res$DP1 <- apply(obs, 1, dp1)
    res$DP2 <- apply(obs, 1, dp2)
  }

  # where association measure range from 0 to +, but changes at 0
  if (am %in% c("mi", "pois", "t", "mi3", "random",
                "chisq", "cramersV", "dice", "fye", "fye.ln", "gmean", "jaccard", "logl", "ms")) {
    if (reverse == FALSE) {
      res <- res[with(res, order(-STR.DIR)), ]
    } else if (reverse == TRUE) {
      res <- res[with(res, order(STR.DIR)), ]
    }
  }

  # where am ranges from - to + so it can be ordered
  if (am %in% c("z", "z.cor", "odds", "liddell")) {
    if (reverse == FALSE) {
      res <- res[with(res, order(COLL.STR, decreasing = TRUE)), ]
    } else {
      res <- res[with(res, order(COLL.STR)), ]
    }
  }

  # display direction strength
  if (str.dir == FALSE) {
      res <- subset(res, select = -STR.DIR)
    } else {
      res$STR.DIR <- round(res$STR.DIR, decimals)
  }

  res$EXP <- round(res$EXP, 1)
  res$COLL.STR <- round(res$COLL.STR, decimals)
  rownames(res) <- NULL
  res$ASSOC <- as.factor(res$ASSOC)
  res$SIGNIF <- as.factor(res$SIGNIF)

  colnames(res)[colnames(res) == "COLL.STR"] <- paste("COLL.STR", toupper(am), sep=".")

  if (am == "fye") {
    message("Note to users of 'fye': packages versions up to 0.0.10 used the
            negative natural logarithm for p-value transformation. Versions >= 0.1.0 use
            the negative decadic logarithm. If you want to continue using the natural logarithm
            transformation, use 'fye.ln' as an association measure and repeat procedure.
            NB I: If you see this message, you are using a version >= 0.1.0.
            NB II: I'm a message, *not* an error. I will disappear in CRAN, >= 1.0.")
  }

  res

}




