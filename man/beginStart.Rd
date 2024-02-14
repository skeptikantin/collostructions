\name{beginStart}
\alias{beginStart}
\docType{data}
\title{
Data set (\emph{begin/start-to}-VERB)
}
\description{
Data set of the \emph{begin/start-to-VERB} pattern in the British National Corpus (BNC), with the frequencies of the verbs in the open slot. Aggregate of \code{beginToV} and \code{startToV}, to illustrate the easiest use of \code{collex.dist()} with aggregated frequency lists, e.g., imported from outside R.
}
\usage{data("beginStart")}
\format{
  A data frame with 2163 observations on the following 3 variables.
  \describe{
    \item{\code{WORD}}{a factor with 2,163 verb types either in \emph{begin to} V and/or \emph{start to} V.}
    \item{\code{beginToV}}{numeric vector, frequencies of verbs with \emph{begin to} V.}
    \item{\code{startToV}}{numeric vector, frequencies of verbs with \emph{start to} V.}
  }
}

\examples{\dontrun{
## Distinctive Collexeme Analysis

# load data
data("beginStart")

# perform Distinctive Collexeme Analysis (with defaults)
# see ?collex.dist() for more use cases:
x <- collex.dist(beginStart)

}
}
