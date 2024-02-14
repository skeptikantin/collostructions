\name{collex.dist}
\alias{collex.dist}
\title{Function for Distinctive Collexeme Analysis}
\description{Implementation of Distinctive Collexeme Analysis (Gries & Stefanowitsch 2004) to be run over a data frame with frequencies of items in two alternating constructions OR, more generally, for keyword analysis over a data frame with frequency lists of two corpora OR, most generally, over a data frame with comparison of frequencies of two conditions. Note: if you want to perform Multiple Distinctive Collexeme Analysis, use \code{collex.covar()}, see \bold{Notes} below.}
\usage{
collex.dist(x, am = "logl", raw = FALSE, reverse = FALSE, decimals = 5,
            threshold = 1, cxn.freqs = NULL, str.dir = FALSE, p.fye = FALSE,
            delta.p = FALSE)}
\arguments{
  \item{x}{Input, a data frame. Two options: EITHER as aggregated frequencies with the types in column A (\code{WORD}), and the frequencies of \code{WORD} in the first construction in column 2 and in the frequencies of \code{WORD} in the second construction in column 3, OR as raw data, i.e., one observation per line, where column 1 must contain the construction types and column 2 must contain the collexeme (see \code{head(future)} or \code{?future} for an example of this format). The names of the columns is up to you, but you must set \code{raw = TRUE} if the input data frame has one observation per line, otherwise the function will abort with an error message to that effect (see below).}
  \item{am}{Association measure to be calculated. Currently available, tested, and conventional in the collostruction literature:\cr
  \code{"logl"} (log likelihood, the default)\cr
  \code{"fye"} (negative decadic log transformed p-value of FYE, the original), can be used in conjunction witht the argument \code{p.fye} to calculate the original Fisher-Yates p-value.\cr
  \code{"fye.ln"} (negative natural log transformed p-value of FYE, variant of \code{fye}).\cr\cr
  Experimental, so use with caution (implementing Evert 2004). They are primarily collocation measures, and most make not much sense in the context of collostructions, but may be useful for some purposes:\cr\code{"chisq", "dice", "gmean", "jaccard", "liddell", "mi", "mi3", "ms", "odds", "pois", "t", "z", "z.cor", "random"}. Based on a chi-square test, there is also \code{"cramersV"}, which is an effect size for chi-squared.}
  \item{raw}{Does input data frame contain a raw list of occurrences? Leave default \code{raw = FALSE}, if your input contains aggregated frequencies with collexeme types in column 1 and the frequencies in the alternating constructions in columns 2 and 3, respectively. Set this argument to \code{raw = TRUE}, if you have one observation per line, where the type of construction is in column 1 and the collexeme in column 2 (see \code{?future} data set for an example). The function will then produce the combined frequency list before performing the distinctive collexeme analysis.}
  \item{reverse}{\code{FALSE} (default) or \code{TRUE}. If \code{FALSE}, output will be ordered in descending order of collostruction strength relative to condition A (i.e. attraction to condition A). Use \code{reverse = TRUE} for ordering in descending order of collostruction strength relative to condition B (i.e. repulsion to condition A).}
  \item{decimals}{Number of decimals in the output. Default is 5 decimal places (except for \code{EXP}, which is always 1).}
  \item{threshold}{Frequency threshold for items you want to calculate association measures for. By default, this is 1 (= calculated for all items). Any number higher than 1 will exclude items that have a combined frequency lower than the threshold. For instance, \code{threshold = 2} will exclude hapaxes, \code{threshold = 3} will exclude items occurring only once or twice, etc. Note: the threshold refers to overall frequency, i.e., if an item occurs four times in one condition, but never in the other, this item is excluded at \code{threshold = 5}, whereas items that occur twice in condition A and three times in condition B (total of 5) is included, so make sure this is what you want.}
  \item{cxn.freqs}{A numeric vector or list of length 2. This option lets you enter construction frequencies 'manually' if \code{x} is a reduced dataset (e.g., from \code{join.freqs()} with a threshold higher than 1). For full data in \code{x}, the function will extract the required construction frequencies for calulations of association directly from \code{x}, but if \code{x} contains a reduced data set, you *must* provide corpus frequencies or else your results will be wrong (although they may look reasonable). Note that \code{cxn.freqs} is independent of the setting for \code{threshold} (which refers to the output).}
  \item{str.dir}{Do you want a "directed" association measure in the output? For measures that are positive for attracted and repelled items, \code{str.dir = TRUE} returns negative values for repelled items and allows differentiation of attraction based on this measure. Default is \code{FALSE}.}
  \item{p.fye}{If \code{am = "fye"}, should traditional pFYE-value be calculated? This feature is experimental, and will not have an effect unless \code{am = "fye"}.}
  \item{delta.p}{Should delta\emph{P} be calculated? If yes, both types of delta\emph{P} will be calculated (see below).}
}
\details{
\bold{FYE}: Note to users of \code{am = "fye"}: packages versions up to 0.0.10 used the negative natural logarithm for p-value transformation. Versions >= 0.1.0 use the negative decadic logarithm. If you want to continue using the natural logarithm transformation, use \code{am = "fye.ln"} as an association measure and repeat procedure. (If you see this message, you are using a version >= 0.1.0. It will disappear in versions >= 1.0 once on CRAN).\cr\cr
\bold{Association measures}: The default \code{"logl"} as an association measure is due to the fact that for larger datasets from larger corpora the original \code{"fye"} easily returns \code{Inf} for the most strongly associated and/or dissociated items, which are then non-rankable (they are ranked by frequency of occurrence if this happens).\cr\cr
\bold{Thresholds}: The \code{threshold} argument default of 1 does not remove non-occurring items (although the logic of this argument implies as much). This is a "bug" that I decided to keep for historical reasons (and to avoid problems with \code{collex()}. Use primarily to leave out low-frequency items in the \emph{output}.
}
\value{
  OUTPUT ordered in descending order of association to cxn/condition A (unless \code{reverse = TRUE}):
  \item{COLLEX}{type/item (e.g. string, lemma...)}
  \item{O.CXN1}{Observed frequency in cxn/condition A}
  \item{E.CXN1}{Expected frequency in cxn/condition A}
  \item{O.CXN2}{Observed frequency in cxn/condition B}
  \item{E.CXN2}{Expected frequency in cxn/condition B}
  \item{ASSOC}{Name of cxn/condition with which the COLLEX type is associated.}
  \item{COLL.STR./AM/}{Association measure used. \code{/AM/}, i.e., type of association score should always be reported along with values.}
  \item{STR.DIR}{"Directed" collostruction strength. \code{STR.DIR} is positive if item is attracted to condition A and negative if attracted to condition B (i.e. 'negatively attracted' to condition A). Only displayed if \code{str.dir = TRUE}.}
    \item{DP1}{Association (word-to-cxn), i.e., delta\emph{P}(w|cxn), see Gries & Ellis (2015: 240)}
    \item{DP2}{Association (cxn-to-word), i.e., delta\emph{P}(cxn|w), see Gries & Ellis (2015: 240)}
  \item{SIGNIF}{Significance level.  \code{*****} = significant at \emph{p} < .00001, \code{****} = significant at \emph{p} < .0001, \code{***} at \emph{p} < .001, \code{**} at \emph{p} < .01, \code{*} at \emph{p} < .05, and \code{ns} is not significant (but tendency is given by the difference between \code{OBS} and \code{EXP}, i.e. as returned in \code{ASSOC}). Association measures without significance tests have \code{SIGNIF == na}.}
  \item{SHARED}{Is \code{COLLEX} a shared type between cxn/condition A and B? Returns \code{Y} (yes) or \code{N} (no).}
}
\references{
Evert, Stefan. 2004. \emph{The statistics of word cooccurrences: Word pairs and collocations}. U Stuttgart Dissertation. {\href{http://www.collocations.de/AM/}{http://www.collocations.de/AM/}}\cr\cr
Gries, Stefan Th. & Nick C. Ellis. 2015. Statistical measures for usage-based linguistics. \emph{Language Learning} 65(S1). 228–255. doi:10.1111/lang.12119.\cr\cr
Gries, Stefan Th. & Anatol Stefanowitsch. 2004. Extending collostructional analysis: A corpus-based perspective on "alternations." \emph{International Journal of Corpus Linguistics} 9(1). 97-129.\cr\cr
Hilpert, Martin. 2011. Diachronic collostructional analysis: How to use it and how to deal with confounding factors. In Kathryn Allan & Justyna A. Robinson (eds.), \emph{Current methods in historical semantics}, 133–160. Berlin & Boston: De Gruyter.\cr

Johannsen, Berit & Susanne Flach. 2015. Systematicity beyond obligatoriness in the history of the English progressive. Paper presented at ICAME 36, 27–31 May 2015, Universität Trier.
}
\note{
\bold{Multiple Distinctive Collexeme Analysis}: If you want to perform a Multiple Distinctive Collexeme Analysis (MDCA) for more than two levels of your construction, you cannot use \code{collex.dist()}, as it can only handle two levels of cxn. Instead, use \code{collex.covar()} for Co-Varying Collexeme Analysis (CCA), which can handle more than two levels of the first condition. The main difference between MDCA and CCA is conceptual (and arguably historical), but they are mathematically the same thing; the association scores of CCA and approximation-based MDCA correlate highly (e.g., for the \code{modadv} data: Pearson \emph{r} = .9987841; Spearman's rho = .9999993).
}

\author{
Susanne Flach, susanne.flach@es.uzh.ch

Thanks to Anatol Stefanowitsch, Berit Johannsen, Kirsten Middeke and Volodymyr Dekalo for discussion, suggestions, debugging, and constructive complaining.
}
\seealso{
See \code{freq.list()} for an easy function to create frequency lists from vectors of characters in preparation for Distinctive Collexeme Analysis.For use with incomplete data, see example in \code{\link{ditrdat_pub}}.
}
\examples{\dontrun{
##### Calculate Distinctive Collexeme Analysis
## This is a little lengthy, because there are multiple ways to provide
## input to DCA (Case 1, Case 2, and Case 3). There are also use cases
## to run multiple DCA over a list of files (see below).

### Case 1: An aggregated frequency list
## The easiest use case: Words in col1, and their frequencies
## in cxns A and B in col2 and col3, respectively:
# load data
data(beginStart)
# perform DCA (no settings necessary, more as required):
beginStart.dca1 <- collex.dist(beginStart)
beginStart.dca2 <- collex.dist(beginStart, am = "fye")
beginStart.dca3 <- collex.dist(beginStart, am = "fye", str.dir = TRUE)

# inspect:
head(beginStart.dca1, 15)  # 15 most strongly attracted items to cxn A
tail(beginStart.dca1, 25)  # 20 most strongly attracted items to cxn B

# cleanup (remove objects from workspace):
rm(beginStart.dca1, beginStart.dca2, beginStart.dca3)

### Case 2: Two separate aggregated frequency lists
## Like Case 1, but with separate lists of cxns A and B that need to combined:
# load data
data(beginToV)
data(startToV)

# I. Merge the lists
beginStart.in <- join.freqs(beginToV, startToV)
head(beginStart.in, 12)

# II. Calculate association
beginStart.out <- collex.dist(beginStart.in)

# III. Inspect
head(beginStart.out, 15)   # 15 most strongly attracted items to cxn A
tail(beginStart.out, 20)   # 20 most strongly attracted items to cxn B

# cleanup (remove objects from workspace):
rm(beginToV, startToV, beginStart.in, beginStart.out)

### Case 3: A list with one observation per line (i.e. raw = TRUE)
# where the cxns are in col1 and the collexemes are in col2:
# load & inspect the will/going-to-V alternation:
data(future)
head(future, 12)

# Calculate:
future.out <- collex.dist(future, raw = TRUE)
head(future.out, 6)
tail(future.out, 6)

# cleanup (remove objects from workspace):
rm(future, future.out)

##### IF YOU HAVE INCOMPLETE DATA SETS
## Illustrate the application of the cxn.freq argument if you do not have all
## types; this is *not* a sample of a larger data set, but rather as if lines
## from an aggregate frequency lists were unavailable. To illustrate, we'll
## recreate the dative alternation from Gries & Stefanowitsch (2004).
data(ditrdat_pub)

# The data is from Gries & Stefanowitsch (2004: 106), ie. the top collexemes for
# the ditransitive vs. the to-dative. That is, the low-frequent items are not
# in their results list in the table.
# So the following would lead to (linguistically) wrong results, because
# collex.dist() calculates the cxn frequencies from an incomplete data set,
# when in fact they are much higher:

collex.dist(ditrdat_pub, am = "fye")

# However, you can recreate the results table by specifying the cxn frequencies
# provided in the publication (as cxn.freq), as a vector, where the first
# element contains the total of the first cxn and the second contains the total
# of the second cxn. You can also get the traditional Fisher-Yates p-value as
# in the original publication:

ditrdat.dca <- collex.dist(ditrdat_pub, am = "fye", p.fye = TRUE,
                           cxn.freqs = c(1035, 1919), decimals = 3)
# Inspect:
head(ditrdat.dca, 20)  # left side of Table 2 (Gries & Stefanowitsch 2004: 106)
tail(ditrdat.dca, 19)  # right side of Table 2 (Gries & Stefanowitsch 2004: 106)
# the right side of Table 2 is "upside down", because collex.dist() orders by
# collostructional continuum. Run the above again with reverse = T.


# NB: If you have a raw input file, make sure you pass a vector with the correct
# frequencies for whatever R recognizes as the first element (usually
# alphabetically if column 1 is a factor); this behavior has to be checked
# carefully in R4.x, as R now standardly reads in characters as characters.

##### IN USE OVER LISTS, e.g.,
## We performed several Distictive Collexeme Analyses for (present) progressive
## vs. (simple) present over 10 25-yr periods in CLMET (Johannsen & Flach 2015).
## Note that although using historical data, this is quite different to
## Diachronic Distinctive Collexeme Analysis (e.g., Hilpert 2011), where periods
## are conditions in *one* DCA and thus mathematically *not* independent of
## each other. The sample data below runs one DCA *per period*, so the DCAs are
## mathematically independent of each other. The conditions are still
## two alternating constructions as in 'ordinary' DCA.

## So this means 'multiple' DCAs in the sense of 'several' DCAs, not in the
## sense of 'Multiple Distinctive Collexeme Analysis' (MDCA).

## load data
data(CLMETprog.qc)
data(CLMETsimple.qc)

### I. Prepare

# split data by time period and drop redudant period column
prog <- split(CLMETprog.qc[, c(1,3)], CLMETprog.qc$QUARTCENT)
simple <- split(CLMETsimple.qc[, c(1,3)], CLMETsimple.qc$QUARTCENT)

# combine frequencies for progressive & simple as input to collex.dist()
dist.in <- join.lists(prog, simple)
dist.in <- lapply(dist.in, droplevels)

### II. Perform several DCA (returns a list of DCA output for list)
# if only defaults of arguments are used, use lapply() like so:
dist.out <- lapply(dist.in, collex.dist)

# if you want to override default arguments, use Map() like so:
dist.out <- Map(collex.dist, dist.in, am = "fye")
dist.out <- Map(collex.dist, dist.in, decimals = 7)

### III. Inspect output
str(dist.out)     # structure of list
str(dist.out[1])  # structure of first item in list

## VI. Export (works if you have installed package 'openxlsx')
# Will write each DCA in a separate Excel worksheet
openxlsx::write.xlsx(dist.out, "ProgSimpleDistColl_CLMET.xlsx")
}
}
