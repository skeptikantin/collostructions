\name{CLMETsimple.qc}
\alias{CLMETsimple.qc}
\docType{data}
\title{Data set (Present Simple by quarter century, CLMET)
}
\description{
Data based on Johannsen & Flach (2015) with frequencies of all verbs occurring in the simple present in the \emph{Corpus of Late Modern English Texts} (CLMET-3.1). Contains only occurrences with unambiguous assignment of quarter century.
}
\usage{data("CLMETsimple.qc")}
\format{
  A data frame with 24,693 observations on the following 3 variables.
  \describe{
    \item{\code{WORD}}{a factor with levels for all lemmas occurring as types in present simple}
    \item{\code{QUARTCENT}}{a factor with levels \code{1700-1724} \code{1725-1749} ... \code{1900-1924}}
    \item{\code{CORP.FREQ}}{a numeric vector with the corpus frequencies}
  }

}
\details{
CQP version query: [pos="VB[PZ]"] -- frequencies for \emph{be} were reduced by their frequency value in \code{CLMETprog.qc} to avoid them being present in both data sets.
}
\source{
\emph{Corpus of Late Modern English Texts}, CLMET-3.1 (De Smet, Flach, Tyrkkö & Diller 2015), CQP version.
}
\references{

De Smet, Hendrik, Susanne Flach, Jukka Tyrkkö & Hans-Jürgen Diller. 2015. \emph{The Corpus of Late Modern English (CLMET), version 3.1: Improved tokenization and linguistic annotation.} KU Leuven, FU Berlin, U Tampere, RU Bochum.

Johannsen, Berit & Susanne Flach. 2015. Systematicity beyond obligatoriness in the history of the English progressive. Paper presented at ICAME 36, 27–31 May 2015, Universität Trier.
}
\examples{\dontrun{
data(CLMETsimple.qc)
str(CLMETsimple.qc)
head(CLMETsimple.qc)
}}
\keyword{datasets}
