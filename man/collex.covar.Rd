\name{collex.covar}
\alias{collex.covar}
\title{Function for Covarying Collexeme Analysis}
\description{Implementation of Covarying Collexeme Analysis (Gries & Stefanowitsch 2004, Stefanowitsch & Gries 2005) to investigate the collostructional interaction between two slots of a construction.
}
\usage{
collex.covar(x, am = "logl", raw = TRUE, all = FALSE, reverse = FALSE,
            decimals = 5, str.dir = FALSE, p.fye = FALSE, delta.p = FALSE)
}
\arguments{
  \item{x}{Input, a data frame. Two options: EITHER as raw, with one observation per line, and with collexeme 1 in column 1 and collexeme 2 in column 2 (in which case \code{raw = TRUE}, the default) OR as an aggregated frequency list, which must contain a third column with the frequency of the combinations in columns 1 and 2 (in which case you must set \code{raw = FALSE}).}
  \item{am}{Association measure to be calculated. Currently available, tested, and conventional in the collostruction literature:\cr
  \code{"logl"} (log likelihood, the default)\cr
  \code{"fye"} (negative decadic log transformed p-value of FYE, the original), can be used in conjunction witht the argument \code{p.fye} to calculate the original Fisher-Yates p-value.\cr
  \code{"fye.ln"} (negative natural log transformed p-value of FYE, variant of \code{fye}).\cr\cr
  Experimental, so use with caution (implementing Evert 2004). They are primarily collocation measures, and most make not much sense in the context of collostructions, but may be useful for some purposes:\cr\code{"chisq", "dice", "gmean", "jaccard", "liddell", "mi", "mi3", "ms", "odds", "pois", "t", "z", "z.cor", "random"}. }
  \item{raw}{\code{TRUE} (default) or \code{FALSE}. If you have one observation per line (i.e. a 'raw' list), use default. If you enter an aggregated frequency list, set \code{raw = FALSE}.}
  \item{all}{\code{FALSE} (default) or \code{TRUE}. The default calculates attested combinations. If all possible combinations of slots 1 and 2 are to be included, set \code{all = TRUE}, but this can take a long time for constructions with high type frequencies.
}
  \item{reverse}{\code{FALSE} (default) or \code{TRUE}. If \code{FALSE}, output will be ordered in descending order of attraction. Use \code{TRUE} for ordering in descending order of repulsion.}
  \item{decimals}{Number of decimals in the output. Default is 5 decimal places (except for \code{EXP}, which is always 1).}
  \item{str.dir}{Do you want a "directed" collostruction strength in the output? For measures that are positive for both attracted and repelled items, \code{str.dir = TRUE} returns negative values for repelled items. Default is \code{FALSE}.}
  \item{p.fye}{If \code{am = "fye"}, should traditional pFYE-value be calculated? This feature is experimental, and will not have an effect unless \code{am = "fye"}.}
  \item{delta.p}{Should delta\emph{P} be calculated? If yes, both types of delta\emph{P} will be calculated (see below).}
}
\value{
  Output is ordered in descending order of association (unless \code{reverse = TRUE}):
  \item{SLOT1}{type/item in slot 1 (e.g. string or lemma...)}
  \item{SLOT2}{type/item in slot 2 (e.g. string or lemma...)}
  \item{fS1}{Total frequency of item 1 in slot 1 of cxn}
  \item{fS2}{Total frequency of item 2 in slot 2 of cxn}
  \item{OBS}{Observed frequency of combination}
  \item{EXP}{Expected frequency of combination}
  \item{ASSOC}{Association of combination (\code{attr} or \code{rep}).}
  \item{COLL.STR./AM/}{Value of association measure used.}
  \item{STR.DIR}{"Directed" collostruction strength. \code{STR.DIR} is positive if item is attracted and negative if repelled. (Only displayed if \code{str.dir = TRUE}.)}
    \item{DP1}{Association (slot1-to-slot2), i.e., delta\emph{P}(s2|s1), how predictive is slot 1 of slot 2? See Gries & Ellis (2015: 240)} in the constructional context.
    \item{DP2}{Association (slot2-to-slot1), i.e., delta\emph{P}(s1|s2), how predictive is slot 2 of slot 1? See Gries & Ellis (2015: 240)} in the constructional context.
  \item{SIGNIF}{Significance level.  \code{*****} = significant at \emph{p} < .00001, \code{****} = significant at \emph{p} < .0001, \code{***} at \emph{p} < .001, \code{**} at \emph{p} < .01, \code{*} at \emph{p} < .05, and \code{ns} is not significant (but tendency is given by the difference between \code{OBS} and \code{EXP}, i.e. as returned in \code{ASSOC}). Association measures without significance tests have \code{SIGNIF == na}.}
}
\note{
If you use the function on constructions with a high type frequency, be patient when setting \code{all = TRUE}. The function needs to perform \code{Types.In.A x Types.In.B} number of tests. Even on a fairly powerful computer it can take about half an hour to perform ~500,000+ tests.\cr\cr
For Multiple Distinctive Collexeme Analysis (MDCA), where you have more than two conditions/constructions, you can use \code{collex.covar()}. The association scores of CCA (\code{collex.covar()}) and approximation-based MDCA correlate highly (e.g., for the \code{\link{modadv}} data: Pearson \emph{r} = .9987841; Spearman's rho = .9999993), suggesting \code{collex.covar()} is a workable alternative to approximation.
\cr\cr
Note to users of \code{am = "fye"}: packages versions up to 0.0.10 used the negative natural logarithm for p-value transformation. Versions >= 0.1.0 use the negative decadic logarithm. If you want to continue using the natural logarithm transformation, use \code{am = "fye.ln"} as an association measure and repeat procedure. (If you see this message, you are using a version >= 0.1.0. It will disappear in versions >= 1.0 once on CRAN).
}
\author{
Susanne Flach, susanne.flach@es.uzh.ch

Thanks to Anatol Stefanowitsch, Berit Johannsen, Kirsten Middeke and Volodymyr Dekalo for suggestions, debugging, and constructive complaining.
}
\seealso{
Use \code{\link{reshape.cca}} to transform the output of \code{\link{collex.covar}} from the 'long' format to a 'wide' format (i.e., cross-tabulate association scores.)
}
\references{
Evert, Stefan. 2004. The statistics of word cooccurrences. Word pairs and collocations. Stuttgart: Universität Stuttgart Doctoral Dissertation. http://www.stefan-evert.de/PUB/Evert2004phd.pdf.\cr\cr
Gries, Stefan Th. & Nick C. Ellis. 2015. Statistical measures for usage-based linguistics. \emph{Language Learning} 65(S1). 228–255. doi:10.1111/lang.12119.\cr\cr
Gries, Stefan Th. & Anatol Stefanowitsch. 2004. Covarying collexemes in the into-causative. In Michel Archard & Suzanne Kemmer (eds.), \emph{Language, culture, and mind}, 225–236. Stanford, CA: CSLI.\cr\cr
Stefanowitsch, Anatol & Stefan Th. Gries. 2005. Covarying collexemes. \emph{Corpus Linguistics and Linguistic Theory} 1(1). 1–43. \href{http://dx.doi.org/doi:10.1515/cllt.2005.1.1.1}{doi:10.1515/cllt.2005.1.1.1}.
}
\examples{\dontrun{

### Example I: Attested combinations (only)
data(caus.into)

# inspect
head(caus.into)

# subset, because caus.into contains too many variables
into.vrbs <- caus.into[, c(2,3)]    # CCA between V1 and V2
into.voice <- caus.into[, c(1,2)]   # CCA between VOICE and V1

# perform Co-Varying Collexeme Analysis
into.vrbs.cca <- collex.covar(into.vrbs)
into.voice.cca <- collex.covar(into.voice)

# clear workspace (remove objects)
rm(into.voice, into.vrbs.cca, into.voice.cca)

### Example 2: If you want to test all possible combinations
# Depending on your machine, this may take a while, because it needs
# to perform 199*426 = 84,774 tests for the 199 slot 1 types and the
# 426 slot 2 types (rather than the 1,076 tests for attested combinations).
into.vrbs.cca.all <- collex.covar(into.vrbs, all = TRUE)

### Example 3: An aggregated list
# set raw = FALSE, otherwise the function will abort
# (output wouldn't make any sense):
data(modadv)
head(modadv, 12)
modadv.cca <- collex.covar(modadv, raw = FALSE)
}}
