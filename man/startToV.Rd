\name{startToV}
\alias{startToV}
\docType{data}
\title{Data set (start-to-VERB)
}
\description{Data set of the \emph{start-to-}VERB construction in the British National Corpus (BNC), with the frequencies of the verbs in the open slot ([hw="start" & class="VERB"] [hw="to"] [pos="V.I"]).
}
\usage{data("startToV")}
\format{
  A data frame with 1168 observations on the following 2 variables.
  \describe{
    \item{\code{WORD}}{a factor with levels of types start-to-V}
    \item{\code{CXN.FREQ}}{a numeric vector of the frequencies in V2.}
  }
}
\examples{\dontrun{

data(startToV)     # load
str(startToV)      # inspect structure of object
head(startToV)     # view head of object
}
}
\keyword{datasets}
