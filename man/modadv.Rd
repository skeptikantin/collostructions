\name{modadv}
\alias{modadv}
\docType{data}
\title{
Data set (modal-adverb combinations)
}
\description{
Data set of 792 modal-adverb pairs in the BNC-BABY, such as \emph{would possibly}, \emph{may well} or \emph{can hardly}.\cr
}
\usage{data("modadv")}
\format{
  A data frame with 792 observations on the following 3 variables.
  \describe{
    \item{\code{MODAL}}{A factor with 11 levels of the core modals and contractions, i.e., \code{\'d}, \code{\'ll}, \code{can}, \code{could}, \code{may}, \code{might}, \code{must}, \code{shall}, \code{should} \code{will}, \code{would}}
    \item{\code{ADVERB}}{A factor with 280 levels, for each adverb types following modal verbs, e.g., \code{certainly}, \code{essentially}, \code{even}, \code{lawfully}, \code{publicly}, \code{quickly}, and \code{well}}
    \item{\code{FREQ}}{The frequency of the combination.}
  }
}
\source{
BNC-BABY; [hw="will|would|can|could|may|might|must|shall|should" & pos="VM0"] [pos="AV0"];
cf. Flach (2020) with COCA data.
}
\references{
Flach, Susanne. 2020. Beyond modal idioms and modal harmony: A corpus-based analysis of gradient idiomaticity in MOD+ADV collocations. \emph{English Language and Linguistics}. aop.
}
\examples{\dontrun{
data(modadv)

## Inspect:
# This is an aggregated frequency list:
head(modadv, 12)

### Perform co-varying collexeme analysis
## ?collex.covar()
# since it's aggregated, you must set raw = FALSE, or it will make no sense.
cca.att <- collex.covar(modadv, am="fye", raw = FALSE, all = FALSE) # only attested combinations
cca.all <- collex.covar(modadv, am="fye", raw = FALSE, all = TRUE) # all combinations


## Reshape the cca output by association measure:
# ?reshape.cca
cca.wide.att <- reshape.cca(cca.att)
cca.wide.all <- reshape.cca(cca.all)

}}

