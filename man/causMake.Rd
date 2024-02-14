\name{causMake}
\alias{causMake}
\docType{data}
\title{
Data set (\emph{make}-causative)
}
\description{
A toy sample data set of pronominal \emph{make}-causative (\emph{this makes me feel X}) from the British National Corpus (BNC) to illustrate multiple distinctive co-varying collexeme analysis (see \code{collex.covar.mult()}).
}
\usage{data("causMake")}
\format{
  A data frame with 5000 observations of the following 3 variables.
  \describe{
    \item{\code{MAKE}}{a factor with 4 levels of different inflectional forms of \emph{make}}
    \item{\code{OBJ}}{a factor with 7 pronoun levels levels (i.e., \code{her} \code{him} \code{it} \code{me} \code{them} \code{us} \code{you})}
    \item{\code{V2}}{a factor with levels 633 levels, i.e., the verbs in the complement.}
  }
}

\examples{\dontrun{
## Multiple Distinctive Collexeme Analysis

# load data
data("causMake")

# perform Multiple Distinctive Co-Varying Collexeme Analysis (with defaults)
# see ?collex.covar.mult for more use cases:
x <- collex.covar.mult(causMake)
}
}
