\name{BNCverbL}
\alias{BNCverbL}
\docType{data}
\title{Data set (BNC verb lemma frequency list)
}
\description{A data frame with case-insensitive BNC verb lemma frequencies.
}
\usage{data("BNCverbL")}
\format{
  A data frame with 35939 observations on the following 2 variables.
  \describe{
    \item{\code{WORD}}{a factor with levels for all verb lemmas}
    \item{\code{CORP.FREQ}}{a numeric vector with verb lemma frequencies}
  }
}
\details{Lemmas starting in problematic characters have been removed (almost exclusively tokenization problems, e.g. single quotes, slashes, backslashes, hashes, asterisks or square brackets). Please make sure you have a clean file of frequencies if you want to use this data for joining frequency lists to avoid problems in collex().
}
\source{BNCxml version (CQP query: [class="VERB"])
}
\examples{\dontrun{

data(BNCverbL)
str(BNCverbL)
head(BNCverbL)

}}
\keyword{datasets}
