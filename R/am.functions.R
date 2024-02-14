# Functions for calculations of association measures
logl.am <- function(x) {
  logl <- sum(logw0(x[1]/x[5]) * x[1],
              logw0(x[2]/x[6]) * x[2],
              logw0(x[3]/x[7]) * x[3],
              logw0(x[4]/x[8]) * x[4]) * 2
}
chi.am <- function(x) {
  chi <- sum(((x[1] - x[5])^2 / x[5]),
             ((x[2] - x[6])^2 / x[6]),
             ((x[3] - x[7])^2 / x[7]),
             ((x[4] - x[8])^2 / x[8]))
}

cramer.am <- function(x) {
  chi2 <- sum(((x[1] - x[5])^2 / x[5]),
             ((x[2] - x[6])^2 / x[6]),
             ((x[3] - x[7])^2 / x[7]),
             ((x[4] - x[8])^2 / x[8]))
  V = sqrt(chi2 / (sum(x[1:4])))
}

fye.am <- function(x) {
  m.obs <- matrix(c(x[1], x[2], x[3], x[4]), nrow = 2)
  sr <- rowSums(m.obs)
  sc <- colSums(m.obs)
  n <- sum(m.obs)
  m.exp <- outer(sr, sc, "*")/n

  if (x[1] < m.exp[1]) {
    res <- -log(sum(dhyper(0:x[1], sc[1], sc[2], sr[1])), 10)
  } else if (x[1] > m.exp[1]) {
    res <- -log(sum(dhyper(x[1]:sc[1], sc[1], sc[2], sr[1])), 10)
  } else {
    res <- 0
  }
}
fye.ln.am <- function(x) {
  m.obs <- matrix(c(x[1], x[2], x[3], x[4]), nrow = 2)
  sr <- rowSums(m.obs)
  sc <- colSums(m.obs)
  n <- sum(m.obs)
  m.exp <- outer(sr, sc, "*")/n

  if (x[1] < m.exp[1]) {
    res <- -log(sum(dhyper(0:x[1], sc[1], sc[2], sr[1])))
  } else if (x[1] > m.exp[1]) {
    res <- -log(sum(dhyper(x[1]:sc[1], sc[1], sc[2], sr[1])))
  } else {
    res <- 0
  }
}

mi.am <- function(x) {
  res <- logw0(x[1] / x[5], 10)
}
pois.am <- function(x) {
  res <- (x[1] * (logw0(x[1]) - logw0(x[5]) - 1))
}
z.score.am <- function(x) {
  res = (x[1] - x[5]) / sqrt(x[5])
}
z.score.cor.am <- function(x) {
  if (x[1] >= x[5]) {
    x[1] <- x[1] - 0.5
  } else {
    x[1] <- x[1] + 0.5
  }
  res = (x[1] - x[5]) / sqrt(x[5])
}
t.score.am <- function(x) {
  res = (x[1] - x[5]) / sqrt(x[1])
}
odds.am <- function(x) {
  res <- log( ((x[1] + 0.5) * (x[4] + 0.5)) / ((x[2] + 0.5) * (x[3] + 0.5)) , 10)
}
liddell.am <- function(x) {
  res = ((x[1] * x[4]) - (x[2] * x[3])) / ((x[1] + x[3]) * (x[2] + x[4]))
}
ms.am <- function(x) {
  res = min( x[1] / (x[1] + x[2]), x[1] / (x[1] + x[3]) )
}
gmean.am <- function(x) {
  res = (x[1] / sqrt((x[1] + x[2]) * (x[1] + x[3])))
}
dice.am <- function(x) {
  res <-  (2 * x[1]) / ((x[1] + x[2]) + (x[1] + x[3]))
}
jaccard.am <- function(x) {
  res <- x[1] / (x[1] + x[2] + x[3])
}
mi3.am <- function(x) {
  res <- log(x[1]^3 / x[5], 10)
}

## Functions for collex.covar.mult
logw0 <- function(x, base = 10) {
  ifelse (x == 0, 0, log(x, base))
}

# Adds 0.5 to the observed value of unattested combinations to avoid NaN
poiss.stirl <- function(x) {
  #x[1] * ( logw0(x[1]) - logw0(x[2]) - 1)
  #x[1] * ( logw0(x[1], 10) - logw0(x[2], 10) - 1)
  if(x[1] > 0) {
    x[1] * ( log10(x[1]) - log10(x[2]) - 1)
  } else {
    x[1] * ( log10(x[1] + 0.5) - log10(x[2]) - 1)
  }
}

# local-MI
# O * log(O/E)
# Adds 0.5 to the observed value of unattested combinations to avoid NaN
local.mi <- function(x) {
  if(x[1] > 0) {
    x[1] * (log10(x[1]/x[2]))
  } else {
    (x[1] + 0.5) * (log10((x[1]+0.5)/x[2]))

  }
}

pmi <- function(x) {
  log2(x[1]/x[2])
}

# Adds 0.5 to the observed value of unattested combinations to avoid -Inf
t.val <- function(x) {
  if (x[1] > 0) {
    (x[1]-x[2])/sqrt(x[1])
  } else {
    (x[1]-x[2])/sqrt(x[1] + 0.5)
  }
}
z.val <- function(x) {
  (x[1]-x[2])/sqrt(x[2])
}

# http://dspace.uib.no/bitstream/handle/1956/11033/lyse-andersen-mwe-final.pdf?sequence=1&isAllowed=y
# True Mutual Information
tmi <- function(x, n) {
  x[1] <- ifelse(x[1] == 0, 0.001, x[1])
  (x[1]/n) * (log2(x[1]/x[2]))
}

p.adjust.ca <- function(x, method = p.adj, am = am) {
  if(am %in% c("logl", "chisq")) {
    # convert the AM to p-values
    pvals <- stats::pchisq(x, df=1, lower.tail=FALSE)
    pvals <- stats::p.adjust(pvals, method = method)
    pvals
  } else if (am == "fye") {
    # convert the AM to p-values
    pvals <- 10^-x
    pvals <- stats::p.adjust(pvals, method = method)
    pvals
  } else if (am == "fye.ln") {
    pvals <- exp(1)^-x
    pvals <- stats::p.adjust(pvals, method = method)
    pvals
  } else if (am %in% c("z", "z.cor")) {
    pvals <- 2*stats::pnorm(-abs(x))
    pvals <- stats::p.adjust(pvals, method = method)
    pvals
  }
}


# Function to determine significance levels for adjusted p-values
p.adjust.siglev <- function(x) {
  if (x < 0.00001) { "*****" }
  else if (x < 0.0001) { "****" }
  else if (x < 0.001) { "***" }
  else if (x < 0.01) { "**" }
  else if (x < 0.05) { "*" }
  else { "ns" }
}

logw0 <- function(x, base = exp(1)) {
  ifelse (x == 0, 0, log(x, base))
}
exp.val <-  function(x) {
  n <- sum(x)
  x <- matrix(c(x[1], x[2], x[3], x[4]), nrow = 2, byrow = F)
  nr <- as.integer(nrow(x))
  nc <- as.integer(ncol(x))
  if (is.na(nr) || is.na(nc) || is.na(nr * nc))
    stop("invalid nrow(x) or ncol(x)", domain = NA)
  sr <- rowSums(x)
  sc <- colSums(x)
  E <- outer(sr, sc, "*")/n
  v <- function(r, c, n) c * r * (n - r) * (n - c)/n^3
  V <- outer(sr, sc, v, n)
  dimnames(E) <- dimnames(x)
  E
}

### Delta Ps
# Function to determine delta p, y given x (ΔP(y|x))
dp1 <- function(x) {
  (x[1] / (x[1] + x[2]) ) - (x[3] / (x[3] + x[4]))
}

# Function to determine delta p, x given y (ΔP(x|y))
dp2 <- function(x) {
  (x[1] / (x[1] + x[3]) ) - (x[2] / (x[2] + x[4]))
}

# Define global variables
am <- NULL
p.adj <- NULL

# Functions specific to collex.covar()
association.covar <- function(x) {
  if (x[1] > x[4]) { "attr" }
  else if (x[4] > x[1]) { "rep" }
  else { "none" }
}
direction.covar <- function(x, am) {
  if (am %in% c("mi", "mi3", "liddell", "pois", "t")) {
    x[5]
  } else if (am %in% c("na")) {
    -x[5]
  } else {
    if (x[1] >= x[4]) {
      x[5]
    } else {
      -x[5]
    }
  }
}
sig.level.covar <- function(x, am) {
  if (am %in% c("logl", "chisq")) {
    if (x[5] < 3.841) { "ns" }
    else if (x[5] < 6.635) { "*" }
    else if (x[5] < 10.828) { "**" }
    else if (x[5] < 15.137) { "***" }
    else if (x[5] < 19.511) { "****" }
    else { "*****" }
  } else if (am == "fye") {
    if (x[5] < 1.30103) { "ns" }
    else if (x[5] < 2) { "*" }
    else if (x[5] < 3) { "**" }
    else if (x[5] < 4) { "***" }
    else if (x[5] < 5) { "****" }
    else { "*****" }
  } else if (am == "fye.ln") {
    if (x[5] < 2.995732) { "ns" }
    else if (x[5] < 4.605170) { "*" }
    else if (x[5] < 6.907755) { "**" }
    else if (x[5] < 9.210340) { "***" }
    else if (x[5] < 11.512925) { "****" }
    else { "*****" }
  } else if (am %in% c("z", "z.cor")) {
    if (abs(as.numeric(x[5])) < 1.959964 ) { "ns" }
    else if (abs(as.numeric(x[5])) < 2.575829) { "*" }
    else if (abs(as.numeric(x[5])) < 3.290527) { "**" }
    else if (abs(as.numeric(x[5])) < 3.890592) { "***" }
    else if (abs(as.numeric(x[5])) < 4.417173) { "****" }
    else { "*****" }
  # } else if (am %in% c("t")) {
  #   if (abs(as.numeric(x[5])) < 2.021075 ) { "ns" }
  #   else if (abs(as.numeric(x[5])) < 2.704459) { "*" }
  #   else if (abs(as.numeric(x[5])) < 3.550966) { "**" }
  #   else if (abs(as.numeric(x[5])) < 4.3207) { "***" }
  #   else if (abs(as.numeric(x[5])) < 5.052953) { "****" }
  #   else { "*****" }
  } else { "na"  }
}

setdiff.df <- function(x, y) {
  colnames(x) <- colnames(y)
  z <- x[!duplicated(rbind(y, x))[-seq_len(nrow(y))], ]
  if (nrow(z) == 0) {numeric(0)} else {z}
}

# Functions specific to collex()
association.collex <- function(x) {
  if (x[2] > x[3]) { "attr" }
  else if (x[3] > x[2]) { "rep" }
  else { "none" }
}
direction.collex <- function(x, am) {
  if (am %in% c("mi", "pois", "z", "z.cor", "t", "liddell", "odds", "ms", "mi3")) { x[4] }
  else {
    if (x[2] >= x[3]) { x[4] }
    else { -x[4] }
  }
}
sig.level.collex <- function(x, am) {
  if (am %in% c("logl", "chisq")) {
    if (x[4] < 3.841) { "ns" }
    else if (x[4] < 6.635) { "*" }
    else if (x[4] < 10.828) { "**" }
    else if (x[4] < 15.137) { "***" }
    else if (x[4] < 19.511) { "****" }
    else { "*****" }
  } else if (am == "fye") {
    if (x[4] < 1.30103) { "ns" }
    else if (x[4] < 2) { "*" }
    else if (x[4] < 3) { "**" }
    else if (x[4] < 4) { "***" }
    else if (x[4] < 5) { "****" }
    else { "*****" }
  } else if (am == "fye.ln") {
    if (x[4] < 2.995732) { "ns" }
    else if (x[4] < 4.605170) { "*" }
    else if (x[4] < 6.907755) { "**" }
    else if (x[4] < 9.210340) { "***" }
    else if (x[4] < 11.512925) { "****" }
    else { "*****" }
  } else if (am %in% c("z", "z.cor")) {
    if (abs(as.numeric(x[4])) < 1.959964 ) { "ns" }
    else if (abs(as.numeric(x[4])) < 2.575829) { "*" }
    else if (abs(as.numeric(x[4])) < 3.290527) { "**" }
    else if (abs(as.numeric(x[4])) < 3.890592) { "***" }
    else if (abs(as.numeric(x[4])) < 4.417173) { "****" }
    else { "*****" }
  # } else if (am %in% c("t")) {
  #   if (abs(as.numeric(x[4])) < 2.021075 ) { "ns" }
  #   else if (abs(as.numeric(x[4])) < 2.704459) { "*" }
  #   else if (abs(as.numeric(x[4])) < 3.550966) { "**" }
  #   else if (abs(as.numeric(x[4])) < 4.3207) { "***" }
  #   else if (abs(as.numeric(x[4])) < 5.052953) { "****" }
  #   else { "*****" }
  } else { "na"  }
}

# Functions specific to collex.dist()
# because the coll.str/str.dir is in a different column
# and because the O vs. E relationship can be reversed
direction.dist <- function(x, am) {
  if (am %in% c("pois", "mi3", "chisq", "cramersV", "dice", "fye", "fye.ln", "gmean", "jaccard", "logl", "ms")) {
    if (x[1] >= x[3]) {
      x[5]
    } else {
      -x[5]
    }
  } else {
    if (x[2] >= x[3]) {
      -x[5]
    } else {
      x[5]
    }
  }
}
sig.level.dist <- function(x, am) {
  if (am %in% c("logl", "chisq")) {
    if (x[5] < 3.841) { "ns" }
    else if (x[5] < 6.635) { "*" }
    else if (x[5] < 10.828) { "**" }
    else if (x[5] < 15.137) { "***" }
    else if (x[5] < 19.511) { "****" }
    else { "*****" }
  } else if (am == "fye") {
    if (x[5] < 1.30103) { "ns" }
    else if (x[5] < 2) { "*" }
    else if (x[5] < 3) { "**" }
    else if (x[5] < 4) { "***" }
    else if (x[5] < 5) { "****" }
    else { "*****" }
  } else if (am == "fye.ln") {
    if (x[5] < 2.995732) { "ns" }
    else if (x[5] < 4.605170) { "*" }
    else if (x[5] < 6.907755) { "**" }
    else if (x[5] < 9.210340) { "***" }
    else if (x[5] < 11.512925) { "****" }
    else { "*****" }
  } else if (am %in% c("z", "z.cor")) {
    if (abs(as.numeric(x[5])) < 1.959964 ) { "ns" }
    else if (abs(as.numeric(x[5])) < 2.575829) { "*" }
    else if (abs(as.numeric(x[5])) < 3.290527) { "**" }
    else if (abs(as.numeric(x[5])) < 3.890592) { "***" }
    else if (abs(as.numeric(x[5])) < 4.417173) { "****" }
    else { "*****" }
  # } else if (am %in% c("t")) {
  #   if (abs(as.numeric(x[5])) < 2.021075 ) { "ns" }
  #   else if (abs(as.numeric(x[5])) < 2.704459) { "*" }
  #   else if (abs(as.numeric(x[5])) < 3.550966) { "**" }
  #   else if (abs(as.numeric(x[5])) < 4.3207) { "***" }
  #   else if (abs(as.numeric(x[5])) < 5.052953) { "****" }
  #   else { "*****" }
  } else { "na"  }
}
