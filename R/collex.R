#' @export
#' @importFrom stats complete.cases
collex <- function(x, corpsize = 1e+08L, am = "logl", reverse = FALSE,
                   decimals = 5, threshold = 1, cxn.freq = NULL, str.dir = FALSE, p.fye = FALSE,
                   delta.p = FALSE, p.adj = "none") {

  # if looping over list of dataframes, input data is of class list
  # which needs to be converted to class dataframe
  if (class(x) == "list") {
    x <- as.data.frame(x)
    names(x) <- names(x[[1]])
    rownames(x) <- NULL
  }

  # Checks if input df is the result of some tidyverse something-something
  if (any(class(x) %in% c("tbl_df", "tbl", "grouped_df"))) {
    x <- as.data.frame(x)
  }

  # sort all entries in descending order of frequency to avoid
  # weird sorting for non-rankable items
  x <- x[order(x[,2], decreasing = T),]

  # Define variables that are otherwise undefined globally
  COLLEX=CORP.FREQ=OBS=EXP=ASSOC=COLL.STR=STR.DIR=SIGNIF=P.FYE=NULL

  # checks if input is of class dataframe
  if (is.data.frame(x) == FALSE) {
    stop("\nInvalid 'x' argument. Must be a dataframe.")
  }

  # checks if input contains incomplete cases
  if (any(is.na(x))) {
    print(x[!(complete.cases(x)),])
    stop("\nYour input contains the above lines with incomplete cases.")
  }

  # checks if input contains items where cxn freq is larger than its corp freq
  if (any(x[, 2] > x[, 3]) == TRUE) {
    print(subset(x, x[, 2] > x[, 3]))
    stop("\nFor the above item(s), the frequency in the construction is larger than the frequency in the corpus.")
  }

  # checks if input contains items where cxn freq and corpus frequency is null.
  if (any(x[, 2] + x[, 3] == 0) == TRUE) {
    print(subset(x, x[, 2] + x[, 3] == 0))
    stop("\nFor the above item(s), both their cxn frequency and their corpus frequency is null.")
  }

  # checks if sensible argument for association measure has been entered
  na.am <- match(am, c("logl", "chisq", "cramersV", "dice", "fye", "fye.ln", "gmean",
                       "jaccard", "liddell", "mi", "ms", "mi3", "odds",
                       "pois", "t", "z", "z.cor", "random"))
  if (is.na(na.am)) {
    stop("\nInvalid 'am' argument. See ?collex() for available association measures.")
  }

  # checks if sensible argument for p-value adjustment has been entered
  na.padj <- match(p.adj, c("holm", "hochberg", "hommel", "bonferroni",
                            "BH", "BY", "fdr", "none"))
  if(is.na(na.padj)) {
    stop("\nInvalid 'p.adj' argument. See ?collex() or ?p.adjust() for available methods for p-value adjustments.")
  }

  # input dataframe must contain data of type 'numeric'
  # converts integer vectors to numeric vectors for computation
  if (is.integer(x[ ,2]) == TRUE) {
    x[ ,2] <- as.numeric(x[ ,2])
  }
  if (is.integer(x[ ,3]) == TRUE) {
    x[ ,3] <- as.numeric(x[ ,3])
  }

  # catch construction frequency before potentially reducing set by threshold
  # checks first if the argument is empty, i.e., if cxn.freq has been provided
  if (is.null(cxn.freq)) {
    cxn.freq <- sum(x[, 2])
  }

  # checks validity of threshold argument and reduces set if given
  if (threshold < 0 || !is.numeric(threshold)) {
    stop("\nInvalid 'threshold' argument. Must be a positive number.")
  } else if (threshold != 1) {
    x <- x[x[, 2] >= threshold, ]
  }

  types <- x[, 1]
  x <- as.matrix(x[, 2:3])
  obs.val.collex <- function(x) {
    obs <- matrix(c(x[1], x[2] - x[1], cxn.freq - x[1], corpsize - (sum(c(x[1], x[2]-x[1], cxn.freq - x[1])))), ncol = 4, byrow = F)
  }
  obs <- t(apply(x, 1, obs.val.collex))
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

  res <- as.data.frame(cbind(x[,2], res[, c(1,5,9)]))

  names(res) <- c("CORP.FREQ", "OBS", "EXP", "COLL.STR")
  res$ASSOC <- apply(res, 1, association.collex)
  res$STR.DIR <- apply(res[,1:4], 1, direction.collex, am)
  res$SIGNIF <- apply(res[,1:4], 1, sig.level.collex, am)
  res$COLLEX <- types

  res <- subset(res, select = c(COLLEX, CORP.FREQ, OBS, EXP, ASSOC, COLL.STR, STR.DIR, SIGNIF))

  if (am %in% c("fye") & p.fye == TRUE) {
    res$P.FYE <- signif(10^-res$COLL.STR, decimals)
    res <- subset(res, select = c(COLLEX, CORP.FREQ, OBS, EXP, ASSOC, COLL.STR, STR.DIR, P.FYE, SIGNIF))
  }

  if(p.adj != "none") {
    res$P.ADJ <- p.adjust.ca(res$COLL.STR, method = p.adj, am = am)
    y <- res[ncol(res)]
    res$SIG.ADJ <- as.factor(apply(y, 1, p.adjust.siglev))
  } else {
    res
  }

  res$EXP <- round(res$EXP, 1)

  if (delta.p == TRUE) {
    res$DP1 <- round(apply(obs, 1, dp1), decimals)
    res$DP2 <- round(apply(obs, 1, dp2), decimals)
  }

  # where association measure range from 0 to +, but changes at 0
  if (am %in% c("chisq", "cramersV", "dice", "fye", "fye.ln", "gmean", "logl", "mi", "mi3", "random")) {
    if (reverse == FALSE) {
      res <- res[with(res, order(-STR.DIR)), ]
    } else if (reverse == TRUE) {
      res <- res[with(res, order(STR.DIR)), ]
    }
  }

  # where am ranges from - to + so it can be ordered
  # or where changing +/- does not make sense
  if (am %in% c("jaccard", "liddell", "ms", "odds", "pois", "t", "z", "z.cor")) {
    if (reverse == FALSE) {
      res <- res[with(res, order(COLL.STR, decreasing = TRUE)), ]
    } else {
      res <- res[with(res, order(COLL.STR)), ]
    }
  }

  # Should str.dir be included in the output?
  if (str.dir == FALSE) {
        res <- subset(res, select = -STR.DIR)
      } else {
        res$STR.DIR <- round(res$STR.DIR, decimals)
      }

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




