#' @export
collex.dist <- function(x, am = "logl", raw = FALSE,
                        reverse = FALSE, decimals = 5,
                        threshold = 1, cxn.freqs = NULL, str.dir = FALSE, p.fye = FALSE, delta.p = FALSE) {

  if (class(x) == "list") {
    x <- as.data.frame(x)
    names(x) <- names(x[[1]])
    rownames(x) <- NULL
  }

  # Checks if input df is the result of some tidyverse something-something
  if (any(class(x) %in% c("tbl_df", "tbl", "grouped_df"))) {
    x <- as.data.frame(x)
  }

  # checks if the correct raw argument has been set for raw input
  if (raw == FALSE & ncol(x) == 2 & !(is.numeric(x[,2]))) {
    stop("It appears your input is raw (one obs per line), not aggregated. Try setting raw = TRUE.")
  }

  # checks number of levels in first column to check if there are more than 2 levels in raw file:
  if (raw == TRUE & length(levels(as.factor(x[,1]))) > 2) {
    #print(levels(as.factor(x[,1])))
    stop("Your cxn column has more than two types, probably because of a case-sensitivity or missing values problem (you might want to run input.check(x)). If you want to run Multiple Distinctive Collexeme Analysis, use collex.covar().")
  }

  # turns raw input into aggregate:
  if (raw == TRUE) {
    x <- as.data.frame.matrix(table(x[,2], x[,1]))
    x$COLLEX <- rownames(x)
    x <- subset(x, select = c(3, 1, 2))
  }

  # checks if there are incomplete cases:
  if (raw == FALSE & any(is.na(x))) {
      print(x[!(complete.cases(x)),])
      stop("Your input contains the above lines with incomplete cases.")
    }
  # checks if input contains items where freq in both cxns is null.
  if (raw == FALSE & any(x[, 2] + x[, 3] == 0) == TRUE) {
      print(x[x[, 2] + x[, 3] == 0,])
      stop("Your input contains the above item(s), where both cxn and corpus frequencies are null.")
    }



  # Define variables that are otherwise undefined globally
  COLLEX=O.CXN1=E.CXN1=O.CXN2=E.CXN2=ASSOC=COLL.STR=P.FYE=SIGNIF=NULL
  STR.DIR <- NULL

  # checks if sensible argument for association measure has been entered
  na.am <- match(am, c("logl", "chisq", "cramersV", "dice", "fye", "fye.ln", "gmean",
                       "jaccard", "liddell", "mi", "ms", "mi3", "odds",
                       "pois", "t", "z", "z.cor", "random"))
  if (is.na(na.am))
    stop("Invalid 'am' argument. See ?collex.dist() for available association measures.")

  if (is.integer(x[ ,2])) {
    x[ ,2] <- as.numeric(x[,2])
  }
  if (is.integer(x[ ,3])) {
    x[ ,3] <- as.numeric(x[ ,3])
  }

  logw0 <- function(x, base = exp(1)) {
    ifelse (x == 0, 0, log(x, base))
  }

  if (!is.null(cxn.freqs)) {
    if (!is.vector(cxn.freqs) || length(cxn.freqs) != 2 || !is.numeric(unlist(cxn.freqs))) {
      stop("Invalid 'cxn.freqs' argument. Must be numeric vector or list of length two.")
    } else {
      cxn1.freq <- unlist(cxn.freqs[1])
      cxn2.freq <- unlist(cxn.freqs[2])
    }
  } else {
    cxn1.freq <- sum(x[, 2])
    cxn2.freq <- sum(x[, 3])
  }

  corp.freq <- cxn1.freq + cxn2.freq
  cxn1 <- colnames(x)[2]
  cxn2 <- colnames(x)[3]

  if(sum(x[, 2]) > cxn1.freq || sum(x[, 3]) > cxn2.freq) {
    stop("One or both construction frequencies you provided are too small for the data in 'x'")
  }

  # checks validity of threshold argument and reduces set if given
  if (threshold <= 0 || !is.numeric(threshold)) {
    stop("Invalid 'threshold' argument. Must be a positive number.")
  } else if (threshold != 1) {
    x <- x[x[, 2] + x[, 3] >= threshold, ]
  }

  # matrix with observed frequencies
  types <- x[,1]
  x <- as.matrix(x[,2:3])
    obs.val <- function(x) {
    obs <- matrix(c(x[1], x[2], cxn1.freq - x[1], cxn2.freq - x[2]), ncol = 4, byrow = F)
  }
  obs <- t(apply(x, 1, obs.val))
  row.names(obs) <- types

  # matrix with expected frequencies
  exp <- t(apply(obs, 1, exp.val))

  # combine
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

  res <- as.data.frame(subset.matrix(res, select = c(1:2, 5:6, 9)))
  names(res) <- c("O.CXN1", "O.CXN2", "E.CXN1", "E.CXN2", "COLL.STR")

  association.dist <- function(x) {
    if (x[1] > x[3]) { cxn1 }
    else if (x[2] > x[4]) { cxn2 }
    else { "none" }
  }

  res$ASSOC <- apply(res[,1:5], 1, association.dist)
  res$STR.DIR <- apply(res[,1:5], 1, direction.dist, am)
  res$SIGNIF <- apply(res[1:5], 1, sig.level.dist, am)
  res$COLLEX <- types

  res <- subset(res, select = c(9, 1, 3, 2, 4, 6, 5, 7, 8))
  res$E.CXN1 <- round(res$E.CXN1, 1)
  res$E.CXN2 <- round(res$E.CXN2, 1)

  # Undos the log transformation of FYE p-values
  if (am %in% c("fye") & p.fye == TRUE) {
    res$P.FYE <- signif(10^-res$COLL.STR, decimals)
    res <- subset(res, select = c(COLLEX, O.CXN1, E.CXN1, O.CXN2, E.CXN2, ASSOC, COLL.STR, P.FYE, STR.DIR, SIGNIF))
  }


  if (delta.p == TRUE) {
    res$DP1 <- apply(obs, 1, dp1)
    res$DP2 <- apply(obs, 1, dp2)
  }

  # where association measure range from 0 to +, but changes at 0
  if (am %in% c("chisq", "cramersV", "dice", "fye", "fye.ln", "gmean", "jaccard", "logl", "random", "mi3", "ms")) {
    if (reverse == FALSE) {
      res <- res[with(res, order(-STR.DIR)), ]
    } else if (reverse == TRUE) {
      res <- res[with(res, order(STR.DIR)), ]
    }
  } else if (am == "pois") {
    res <- res[with(res, order(STR.DIR)), ]
  }

  # where am ranges from - to + so it can be ordered
  # or where changing +/- does not make sense
  # note that in collex(), "jaccard" is in this section, but here it makes more sense
  # to include it above to sort by association
  if (am %in% c("liddell", "odds", "pois", "t", "z", "z.cor", "mi")) {
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


  res$COLL.STR <- round(res$COLL.STR, decimals)
  rownames(res) <- NULL

  SHARED <- ifelse(res[, 2] == 0 | res[ ,4] == 0, "N", "Y")
  res$SHARED <- as.factor(SHARED)
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

