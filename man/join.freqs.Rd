\name{join.freqs}
\alias{join.freqs}
\title{Function to merge two frequency lists}
\description{
Function to merge two data frames of frequency lists into a combined data frame of frequencies.
}
\usage{
join.freqs(x, y, all = TRUE, threshold = 1)
}
\arguments{
  \item{x}{A data frame with frequencies for condition A, with \code{WORD} in first colum and frequencies in the second. This could be (i) the frequency list of types in the construction under investigation (for \code{collex()}) or the frequency list of types in construction A when comparing two constructions (for \code{collex.dist())}.}
  \item{y}{A data frame with (i) corpus frequencies of an item in a construction (\code{collex()}) or the frequecies of types in the second construction (\code{collex.dist())}.}
  \item{all}{logical. If \code{TRUE} (the default), then all types from both \code{x} and \code{y} are merged, returning \code{0} in the condition for which a type is not attested. This is recommended especially for \code{collex.dist()}, but is useful for negative evidence calculation for items in \code{collex()}. If \code{FALSE}, the output will contain only the types in \code{x} with frequencies in both \code{x} and \code{y}, but types that only occur in \code{y} will be ignored. The latter is recommended only for joining frequency lists for use in \code{collex()}.}
  \item{threshold}{Numerical. How many times must an item occur overall to be included in the joint list? Default is 1, which means all items are included. If \code{all = FALSE}, setting \code{threshold = number} will mean that items in \code{x} are only included if they occur at least \code{number} times.}
}
\note{The behaviour of \code{join.freqs()} deviates from that of \code{merge()}: If you merge construction frequencies with a list of corpus frequencies with all = \code{FALSE}, and a construction item does not occur in the corpus frequencies list, the item will occur in the output with a corpus frequency of 0. While this will throw an error if you use such a list in \code{collex()}, it will allow you to identify faulty data rather than silently dropping it.\cr\cr
Also, you can merge a frequency list of two columns with a third frequency list, but then you will have to manually adjust the column headers afterwards.}
\details{Output suitable for \code{collex()} and \code{collex.dist()}. The column names of the output data frame (cols 2 and 3) will be identical to the names of the objects of \code{x} and \code{y}. Header for column one will be \code{WORD}.\cr\cr
The difference to \code{join.lists()} is that \code{join.freqs()} two frequency lists are joined, whereas \code{join.lists()} joins two \emph{lists} (i.e., 'list', the R object type) of frequency lists, which may contain, e.g., freuency lists from different periods.}
\examples{\dontrun{

#### Example for Simple Collexeme Analysis
# Using the verb lemma frequencies from the BNC
# begin to rain, begin to think, begin to blossom...
data(beginToV)

## I. Prepare
# merge by lemmas, only types that occur in the construction:
begin1.in <- join.freqs(beginToV, BNCverbL, all = FALSE)
# merge by lemmas, all types, including unattested for 'negative evidence'
begin2.in <- join.freqs(beginToV, BNCverbL)

## II. Perform SCA
# second argument is taken directly from the source data
begin1.out <- collex(begin1.in, sum(BNCverbL$CORP.FREQ))
begin2.out <- collex(begin2.in, sum(BNCverbL$CORP.FREQ))

### Example for Distinctive Collexeme Analysis
# Comparing begin to rain, start to go,...
data(beginToV)  # load set 1
data(startToV)  # load set 2

## I. Prepare
beginStart <- join.freqs(beginToV, startToV) # merge lists (all types)

## II Perform
beginStart.out <- collex.dist(beginStart)

}}
