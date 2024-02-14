\name{beginToV}
\alias{beginToV}
\docType{data}
\title{Data set (begin-to-VERB)
}
\description{Data set of the \emph{begin-to-VERB} construction in the British National Corpus (BNC), with the frequencies of the verbs in the open slot (CQP query: [hw="begin" & class="VERB"] [hw="to"] [pos="V.I"]).}
\usage{data("beginToV")}
\format{
  A data frame with 1957 observations on the following 2 variables.
  \describe{
    \item{\code{WORD}}{A factor with levels of types in open slot of \emph{begin-to-VERB}}
    \item{\code{CXN.FREQ}}{A numeric vector of the frequencies in V2.}
  }
}
\examples{\dontrun{

data(beginToV)     # load
str(beginToV)      # inspect structure of object
head(beginToV)     # view head of object

## Calculate Simple Collexeme Analysis

# load required data
data(beginToV)     # load frequency list for construction
data(BNCverbL)     # load frequency list for verb frequencies (lemmas)

# join frequency files to create input for collex()
beginToV.in <- join.freqs(beginToV, BNCverbL, all = FALSE) # only types in cxn
beginToV.in <- join.freqs(beginToV, BNCverbL) # all types, even if not in cxn

# calculate
beginToV.out <- collex(beginToV.in, sum(BNCverbL$CORP.FREQ)) # using logL
beginToV.out <- collex(beginToV.in, sum(BNCverbL$CORP.FREQ), "mi") # mi

# inspect result
head(beginToV.out, 20) # view first 20 lines of calculation
tail(beginToV.out, 20) # view last 20 lines of calculation

## Calculate Distinctive Collexeme Analysis

# load data
data(beginToV)
data(startToV)

# merge frequency lists
# the first argument to join.freqs() will be the 'default' by which output is
# sorted and Z.DIR is calculated
beginStart.in <- join.freqs(beginToV, startToV) # merge both data frames

# calculate
beginStart.out <- collex.dist(beginStart.in)

# inspect result
head(beginStart.out, 20)
head(beginStart.out, 20)
}
}
\keyword{datasets}
