\name{collex}
\alias{collex}
\title{
Function for Simple Collexeme Analysis
}
\description{
Implementation of Simple Collexeme Analysis (Stefanowitsch & Gries 2003) over a data frame with frequencies of verbs in a construction and their total frequencies in a corpus.
}
\usage{
collex(x, corpsize = 1e+08L, am = "logl", reverse = FALSE, decimals = 5,
       threshold = 1, cxn.freq = NULL, str.dir = FALSE, p.fye = FALSE,
       delta.p = FALSE, p.adj = "none")
}
\arguments{
  \item{x}{
A data frame with types in a construction in column 1 (\code{WORD}), construction frequencies in column 2 (\code{CXN.FREQ}) and corpus frequencies in column 3 (\code{CORP.FREQ}). The name of the columns in your input file is up to you.}
  \item{corpsize}{The size of the corpus in number of tokens (e.g., verb frequencies for all verb constructions, or total corpus size etc.). If not given, default is 100 million words, roughly the size of the BNC, but you should always provide the appropriate number. (Note: \code{corpsize} is different to the argument \code{cxn.freq}, which refers to the total number of tokens of the construction, see below.)}
  \item{am}{Association measure to be calculated. Currently available, tested, and conventional in the collostruction literature:\cr
  \code{"logl"} (log likelihood, the default)\cr
  \code{"fye"} (negative decadic log transformed p-value of FYE, the original), can be used in conjunction witht the argument \code{p.fye} to calculate the original Fisher-Yates p-value.\cr
  \code{"fye.ln"} (negative natural log transformed p-value of FYE, variant of \code{fye}).\cr\cr
  Experimental, so use with caution (implementing Evert 2004). They are association measures for the investigation of collocations, and most do not make much sense in the context of collostructions, but they may be useful for some purposes:\cr\code{"chisq", "dice", "gmean", "jaccard", "liddell", "mi", "mi3", "ms", "odds", "pois", "t", "z", "z.cor", "random"}. Based on a chi-square test, there is also \code{"cramersV"}, which is an effect size for chi-squared.}
  \item{reverse}{If \code{FALSE} (default), output will be ordered in descending order of collostruction strength (i.e. descending order of attraction). Use \code{reverse = TRUE} for ordering in ascending order of collostruction strength (i.e. descending order of repulsion).}
  \item{decimals}{Number of decimals in the output. Default is 5 decimal places (except for \code{EXP}, which is always 1).}
  \item{threshold}{Frequency threshold of items for which collostruction strength is to be calculated, i.e., if you want to exclude hapaxes from the \emph{output} (they are not excluded from the calculation).}
  \item{cxn.freq}{Frequency of construction. Use \emph{only} if \code{x} does not contain all instances of the construction (e.g., if hapaxes have been removed beforehand or if the corpus you queried does not give you access to full frequency information on all tokens/types). The default case is that \code{collex()} automatically calculates \code{cxn.freq} from the input in \code{x}. Note that \code{cxn.freq} is independent of the setting for \code{threshold} (which determines how many of your \emph{input} types you want in the \emph{output}).}
  \item{str.dir}{Do you want a "directed" association measure in the output? For measures that are positive for attracted and repelled items, \code{str.dir = TRUE} returns negative values for repelled items and allows differentiation of attraction based on this measure. Default is \code{FALSE}.}
  \item{p.fye}{Should the traditional Fisher-Yates p-value be calculated? This will not have any effect unless \code{am = "fye"}.}
  \item{delta.p}{Should delta\emph{P} be calculated? If yes, both types of delta\emph{P} will be calculated (see below).}
  \item{p.adj}{If an association measure is chosen that provides significance levels (\code{"logl", "chisq", "fye", "fye.ln"}), the significance levels in the output are adjusted for multiple hypothesis testing. The default is the Bonferroni correction. You can use any of the adjustment methods (\code{"holm", "hochberg", "hommel", "BH", "BY", "fdr"}) as used in \code{\link[stats]{p.adjust}} (follow link for further details). If you don't want a correction, use \code{p.adjust = "none"}.}
}
\details{
\bold{Corpus size}: It's highly recommended to specify \code{corpsize} to a conceptually sensible value. If you do not, the default may give you anything from an error for data from much larger corpora than the BNC to highly inaccuarate results for small corpora or infrequent phenomena. For phenomena from the BNC, this will probably not distort results too much (cf. Gries in multiple discussions of the "Fourth Cell(TM)").\cr\cr
\bold{FYE}: Note to users of \code{am = "fye"}: packages versions up to 0.0.10 used the negative natural logarithm for p-value transformation. Versions >= 0.1.0 use the negative decadic logarithm. If you want to continue using the natural logarithm transformation, use \code{am = "fye.ln"} as an association measure and repeat procedure. (If you see this message, you are using a version >= 0.1.0. It will disappear in versions >= 1.0 once on CRAN).\cr\cr
\bold{Association measures}: The default \code{"logl"} as an association measure is due to the fact that for larger datasets from larger corpora the original \code{"fye"} easily returns \code{Inf} for the most strongly associated and/or dissociated items, which are then non-rankable (they are ranked by frequency of occurrence if this happens).\cr\cr
\bold{Thresholds}: The \code{threshold} argument default of 1 does not remove non-occurring items (although the logic of this argument implies as much). This is a "bug" that I decided to keep for historical reasons. If you do not want to calculate the repulsion for non-occurring items, you need to enter a frequency list that contains only the occurring items.
}

\value{
The output of \code{collex()} is sorted by collostructional strength, with most attracted to least attracted; for ties, ordering is by frequency (further columns depent on argument settings):
  \item{COLLEX}{The collexemes.}
  \item{CORP.FREQ}{Frequency of collexeme in corpus.}
  \item{OBS}{Observed frequency of collexeme in construction.}
  \item{EXP}{Expected frequency in construction}
  \item{ASSOC}{Association of collexeme: \code{attr} (attracted) or \code{rep} (repelled), based on the difference between observed and expected.}
  \item{COLL.STR./AM/}{Association measure used. \code{/AM/}, i.e., type of association score should always be reported along with values.}
  \item{STR.DIR}{Same as collostruction strength (\code{COLL.STR}), but "directed", i.e. positive for attraction and negative for repulsion. Only displayed if \code{str.dir = TRUE}.}
    \item{DP1}{Association (word-to-cxn), i.e., delta\emph{P}(w|cxn), see Gries & Ellis (2015: 240)}
    \item{DP2}{Association (cxn-to-word), i.e., delta\emph{P}(cxn|w), see Gries & Ellis (2015: 240)}
  \item{SIGNIF}{Significance level.  \code{*****} = significant at \emph{p} < .00001, \code{****} = significant at \emph{p} < .0001, \code{***} at \emph{p} < .001, \code{**} at \emph{p} < .01, \code{*} at \emph{p} < .05, and \code{ns} is not significant (but tendency is given by the difference between \code{OBS} and \code{EXP}, i.e. as returned in \code{ASSOC}). Association measures without significance tests have \code{SIGNIF == na}.}
}
\note{
The function will abort if your input data frame has items with 'non-sensical' data, that is, if a collexeme has a higher frequency in the \emph{construction} than it has in the \emph{corpus} (which is logically impossible of course). This is a surprisingly common problem with untidy corpus frequency lists derived from messy annotation, especially when the collexemes have been manually cleaned from a rather inclusive query, but the corpus frequencies have different/erroneous part-of-speech tagging (cf. Flach 2015), where a syntactically quirky constructions in \emph{Let's \bold{go party}} was "hand-cleaned", but \emph{party} did not have any frequency as a verb, because it was always tagged as a noun. As of package version 0.2.0, the function aborts with a descriptive error message, and prints a list of the items with non-sensical frequencies. For further input checks, see \code{input.check()}.
}
\references{
Evert, Stefan. 2004. \emph{The statistics of word cooccurrences: Word pairs and collocations}. U Stuttgart Dissertation. {\href{http://www.collocations.de/AM/}{http://www.collocations.de/AM/}}\cr\cr
Flach, Susanne. 2015. Let's go look at usage: A constructional approach to formal constraints on go-VERB. In Thomas Herbst & Peter Uhrig (eds.), \emph{Yearbook of the German Cognitive Linguistics Association} (Volume 3), 231-252. Berlin: De Gruyter Mouton. doi:10.1515/gcla-2015-0013.\cr\cr
Gries, Stefan Th. & Nick C. Ellis. 2015. Statistical measures for usage-based linguistics. \emph{Language Learning} 65(S1). 228–255. doi:10.1111/lang.12119.\cr\cr
Stefanowitsch, Anatol & Stefan Th. Gries. 2003. Collostructions: Investigating the interaction of words and constructions. \emph{International Journal of Corpus Linguistics} 8(2). 209-243.
}

\author{
Susanne Flach, susanne.flach@es.uzh.ch

Thanks to Anatol Stefanowitsch, Berit Johannsen, Kirsten Middeke, Volodymyr Dekalo and Robert Daugs for suggestions, debugging, and constructive complaining, and to Stefan Hartmann, who doesn't know what a complaint is, but who provided invaluable feedback when asked how the package could be improved.
}

\examples{\dontrun{

#### Calculate Simple Collexeme Analysis
## Example 1: goVerb (cf. Flach 2015)

# load data
data(goVerb)

# inspect data (optional)
head(goVerb, 15) # displays first 15 lines of data frame

# perform collex
goV.out <- collex(goVerb, 616336708)  # total words in corpus (excl. punct)
goV.out <- collex(goVerb, 93993713)   # total verbs in corpus

# inspect output
head(goV.out, 15)  # first 15 items (strongly attracted)
tail(goV.out, 15)  # last 15 items (strongly repelled)

# clear workspace (remove objects)
rm(goVerb, goV.out)

## Example 2: beginToV (also see help file for ?beginToV)
data(beginToV)     # load data for begin-to-V
data(BNCverbL)     # load a frequency list for verb string frequencies

# merge frequency lists (see ?join.freqs):
beginToV.in <- join.freqs(beginToV, BNCverbL, all = FALSE)

# perform collex
beginToV.out <- collex(beginToV.in, sum(BNCverbL$CORP.FREQ)) # using logL

# inspect output
head(beginToV.out, 30)   # first 30 most strongly associated types
tail(beginToV.out, 20)   # last 20 items least strongly associated types

# clear workspace (remove objects)
rm(beginToV, BNCverbL, beginToV.in, beginToV.out)

##### SPECIAL: IN USE OVER LISTS
## collex() can be used to perform several Simple Collexeme Analyses
## in one function call. See ?collex.dist for an example of multiple
## analyses across time periods (simple vs. progressive). The procedure with
## collex() is almost identical, except that you should use Map(...),
## because you have to provide corpus frequencies for each period
## (i.e., for each iteration of collex()):

# 1. Create a numeric vector of corpus frequencies:
corpfreqs <- c(corpFreqPeriod1, corpFreqPeriod2, ...)

# 2. Pass 'corpfreqs' vector as an argument to Map() like so:
myList.out <- Map(collex, myList.in, corpsize = corpfreqs, ...)
}
}

