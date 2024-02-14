\name{causInto}
\alias{causInto}
\docType{data}
\title{
Data set (\emph{into}-causative)}
\description{
A dataset of the \emph{into}-causative (e.g., \emph{they forced us into thinking...}) from the BNC for the illustration of functions. Contains one observation/token per line.}
\usage{data("causInto")}
\format{
  A data frame with 1,426 observations on the following 3 variables.
  \describe{
    \item{\code{VOICE}}{a factor with annotation for voice (levels \code{active}, \code{passive}, and \code{reflexive})}
    \item{\code{V1}}{a factor with levels for the matrix verb, lemmatised (e.g., \emph{they \bold{force}d us into thinking...})}
    \item{\code{V2}}{a factor with levels for the content verb, lemmatised (e.g., \emph{they forced us into \bold{think}ing...})}
  }
}
\references{
Flach, Susanne. 2018. “What’s that passive doing in my argument structure construction?” A note on constructional interaction, verb classes and related issues. Talk at the Workshop on constructions, verb classes, and frame semantics, IDS Mannheim.
}
\examples{\dontrun{

## E.g., in Co-Varying Collexeme Analysis:

# load data
data(causInto)

# inspect, contains more variables than needed for collex.covar()
head(causInto)

# subset
into.cca <- subset(causInto, select = c(V1, V2))

# perform CCA
into.cca <- collex.covar(into.cca)
}}
