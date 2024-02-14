\name{ditrdat_pub}
\alias{ditrdat_pub}
\docType{data}
\title{Data set ditransitive-dative (subset)}
\description{
A subset of the ditransitive-dative alternation in ICE-GB, taken from Gries & Stefanowitsch (2004: 106).}
\usage{data("ditrdat_pub")}
\format{
  A data frame with 39 observations on the following 3 variables.
  \describe{
    \item{\code{VERB}}{a factor variable with the top attracted collexemes to the ditransitive and to-dative, respectively}
    \item{\code{DITR}}{a numeric variable, containing the frequency of \code{VERB} in the ditransitive}
    \item{\code{DAT}}{a numeric variable, containing the frequency of \code{VERB} in the to-dative}
  }
}
\details{
Data to illustrate the use of \code{collex.dist()} with incomplete data sets. For full datasets, \code{collex.dist()} determines the construction totals from the files directly, but this is unavailable if you have incomplete data sets (i.e., with types and their frequencies missing. For example, in their publication, Gries & Stefanowitsch list only data that amounts to 957 (ditransitive) and 813 (to-dative) data points, while the construction totals are 1,035 and 1,919, respectively. In cases of incomplete data, but where construction totals are known, the construction totals need to be passed to the function, see below for an example.
}
\source{
Recreated from Table 2 in Gries & Stefanowitsch (2004: 106).
}
\references{
Gries, Stefan Th. & Anatol Stefanowitsch. 2004. Extending collostructional analysis: A corpus-based perspective on “alternations.” \emph{International Journal of Corpus Linguistics} 9(1). 97–129.
}
\examples{\dontrun{

## 1 Inspect the data: a data frame with 3 columns and the aggregated
## frequencies of verbs in the ditransitive (DITR) and to-dative (DAT).
head(ditrdat_pub)

## 2 Recreate the results in Gries & Stefanowitsch (2004: 106), with
#    the construction frequencies given as the cxn.freqs argument:

# with standards in collex.dist(), i.e., log-likelihood (G2):
collex.dist(ditrdat_pub, cxn.freqs = c(1035, 1919), decimals = 3)

# with p.fye = TRUE, to recreate the traditional p-value (Fisher-Yates)
# (note that p.fye = TRUE will only make sense if am = "fye"):
collex.dist(ditrdat_pub, am = "fye", p.fye = TRUE, cxn.freqs = c(1035, 1919))

## See ?collex.dist() for further examples.
}
}
\keyword{datasets}
