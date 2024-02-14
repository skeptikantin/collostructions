\name{future}
\alias{future}
\docType{data}
\title{
Data set (\emph{will/going-to}-VERB)
}
\description{
Data set (sample) of 10,000 random tokens of the \emph{will}- vs. \emph{going to}-VERB alternation in the spoken section of the BNC2014, for illustration in \code{collex.dist()} with raw input.
}
\usage{data("future")}
\format{
  A data frame with 10000 observations on the following 2 variables.
  \describe{
    \item{\code{CXN}}{a factor with levels \code{going.to} and \code{will}}
    \item{\code{COLLEXEME}}{a factor with 599 verb types in either the \emph{going to V} or \emph{will} cxn.}
  }
}

\examples{\dontrun{
## Distinctive Collexeme Analysis

# load data
data("future")

# perform Distinctive Co-Varying Collexeme Analysis (with defaults)
# see ?collex.dist() for more use cases:
x <- collex.dist(future, raw = TRUE)
# If you do not set raw = TRUE, function aborts:
x <- collex.dist(future)

}
}
