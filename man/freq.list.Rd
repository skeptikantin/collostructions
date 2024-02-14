\name{freq.list}
\alias{freq.list}
\title{
Function to create frequency list from a character vector
}
\description{
Sometimes it is handy to create a frequency list if working with some annotated data rather than importing frequency lists from external programs. Many R functions can create frequency lists, this is just a handy alternative that creates lists that are directly compatible the \code{\link{join.freqs}} function. The output is a data frame and is sorted in descending order of frequency, which may be one or two steps quicker than R functions such as \code{table()}.
}
\usage{
freq.list(x, convert = TRUE, asFactor = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A factor or character vector.}
  \item{convert}{Should the elements in \code{x} be converted to lowercase? Default is \code{convert = TRUE}, meaning that case will be ignored and aggregated.}
  \item{asFactor}{Should the types be converted to \code{as.factor()}? Default is \code{asFactor = TRUE}, because R>=4x has \code{stringsAsFactors = FALSE} as a default, and several other functions in \code{collostructions} my otherwise not work properly.}
}
\value{
  \item{WORD}{The types}
  \item{FREQ}{The frequencies.}
}
\seealso{
Can be used in preparation for \code{\link{join.freqs}.}
}
\examples{\dontrun{

## From a list of raw observations:
head(caus.into, 12)
# so to get a frequency list of the V1 in the into-causative:
freq.list(caus.into$V1)

# or from the future-time expressions:
head(future, 12)
freq.list(future$COLLEXEME)

}}

