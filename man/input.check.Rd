\name{input.check}
\alias{input.check}
\title{
Function to help detect issues in input data frames
}
\description{
Input data from corpora can be messy. This function checks if your input to \code{collex()}, \code{collex.dist} and \code{collex.covar()} has one of three frequent, potential issues: (a) does it contain leading and trailing whitespace, (b) does it contain types that are essentially identical, but not case sensitive (e.g., \emph{fine, Fine}, and/or (c) does it contain missing values?\cr\cr
Beware, though, that the input types to all functions vary a lot and errors can be plentiful, so this function is (still) somewhat experimental and may show issues where they may not be any (and, sadly, vice versa). Function-specific errors are usually also reflected in error messages, but do let me know if you find bugs or have suggestions for more improved checks and error messages.
}
\usage{
input.check(x, raw = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The data frame that you would pass as \code{x} to \code{collex()}, \code{collex.dist()} and \code{collex.covar()}.}
  \item{raw}{Does \code{x} contain aggregated data, i.e., a frequency list (\code{raw = FALSE}) or one observation per line (\code{raw = TRUE})? Note that the function tries to guess: if \code{x} has exactly two colums with text in it, it will assume that it is a raw data file like you use in \code{collex.dist} and \code{collex.covar} with \code{raw = TRUE}, and will abort, unless you set \code{raw = TRUE}.}
}
\details{
DISCLAIMER: this function is quite experimental, so there is of course no guarantee that your data frame is free of technical or conceptual error if the function output does not report any issues. (In some cases, for example, case-sensitivity might even be wanted.) But function should detect the most common problems that arise from messy corpus data or data processing that I have come across in my own work, teaching, and workshops.\cr\cr
Note also that the function does not check for problems specific to a particular collex-function: for example, the check whether corpus frequencies are lower than construction frequencies for simple collexeme analysis is run during the execution of \code{collex()}. If you have suggestions of how to improve the checker function, let me know.\cr\cr
I'd recommend any errors are fixed in the source files, e.g., in Excel or csv files before import into R.}
\value{
Returns an attempt at diagnostics for leading/trailing whitespace, duplicate types (e.g., due to capitalisation inconsistency) and/or empty/missing frequencies and/or types.}
\author{
Report bugs and suggestions to:\cr
Susanne Flach, susanne.flach@es.uzh.ch}
\examples{\dontrun{
## dirty:
df <- data.frame(WORD  = c("Hello", "hello", "hello ", " hello", "Hi ", "  hi", "HI", "", "G'day"),
                  CXN = c(23, NA, 2, 1, 30, 59, NA, 10, 3),
                  CRP = c(569, 3049, 930, 394, 2930, 87, 9, 23, 40))
df
input.check(df)

## a little dirty:
df <- data.frame(WORD  = c("Hello", "Hi", "hi", "HI", "", "G'day"),
                 CXN = c(23, 12, 2, 1, 30, 59),
                 CRP = c(569, 3049, 930, 394, 2930, 28))
df
input.check(df)

## clean:
df <- data.frame(WORD = c("Hello", "Hi", "Goodmorning", "Goodafternoon"),
                 CXN = c(234, 139, 86, 74),
                 CRP = c(23402, 2892, 893, 20923))
input.check(df)

}
}
