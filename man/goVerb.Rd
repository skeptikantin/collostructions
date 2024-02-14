\name{goVerb}
\alias{goVerb}
\docType{data}
\title{Data set (go-VERB)}
\description{Data set from Flach (2015) containing the construction and corpus frequencies of all verbs that occur as V2 in the \emph{go}-VERB construction in the ENCOW14AX01 corpus (Schäfer & Bildhauer 2012). See Flach (2015:§4.2) for data extraction procedure.
}
\usage{data("goVerb")}
\format{
  A data frame with 752 observations on the following 3 variables.
  \describe{
    \item{\code{WORD}}{A factor with the levels for each of the 725 verbs that occurs in the construction}
    \item{\code{CXN.FREQ}}{A numeric vector, containing the observed frequency of V2 in \emph{go}-VERB}
    \item{\code{CORP.FREQ}}{A numeric vector, containing the string frequency of V2 in the corpus}
  }
}
\references{Flach, Susanne. 2015. Let’s go look at usage: A constructional approach to formal constraints on go-VERB. In Thomas Herbst & Peter Uhrig (eds.), \emph{Yearbook of the German Cognitive Linguistics Association (Volume 3)}, 231-252. Berlin: De Gruyter Mouton. \href{http://dx.doi.org/10.1515/gcla-2015-0013}{doi:10.1515/gcla-2015-0013}.

Schäfer, Roland & Felix Bildhauer. 2012. Building large corpora from the web using a new efficient tool chain. In Nicoletta Calzolari, Khalid Choukri, Thierry Declerck, Mehmet Uğur Doğan, Bente Maegaard, Joseph Mariani, Jan Odijk & Stelios Piperidis (eds.), \emph{Proceedings of the Eighth International Conference on Language Resources and Evaluation} (LREC’12), 486–493. Istanbul: ELRA. Available at \url{http://webcorpora.org}
}
\examples{
\dontrun{

data(goVerb)     # load
str(goVerb)      # inspect structure of object
head(goVerb)     # view head of object

collex(goVerb, 616113708) # used in Flach (2015), all tokens minus punctuation
collex(goVerb, 93993713)  # could/should probably be used, all verb tokens
collex(goVerb, 93993713, "chisq")  # returning chisquare statistic}
}
