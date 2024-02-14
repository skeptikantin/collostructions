\name{collex.covar.mult}
\alias{collex.covar.mult}
\title{Function for Multiple Covarying Collexeme Analysis}
\description{Implementation of Multiple/Distinctive Covarying Collexeme Analysis (Stefanowitsch & Flach 2020) to investigate the collostructional association between three or more slots/conditions of a construction.
}
\usage{
collex.covar.mult(x, am = "t", raw = TRUE, all = FALSE, reverse = FALSE,
                  threshold = 1, decimals = 5)
}
\arguments{
  \item{x}{A data frame with each (categorical) condition in one column. In principle, these conditions can be anything: open slots in a construction (verbs, nouns, prepositions, etc.), constructions (e.g., ditransitive/prep-dative, negation [y/n]), annotated variables, genre, time periods...). The function assumes a raw list by default, i.e., with one observation per line. If you have an aggregated list, the last column must contain the frequencies (in which case set \code{raw = FALSE}).}
  \item{am}{Association measure to be calculated. Currently available (though very experimental, see below):\cr
  \code{"t"} (\emph{t}-score, the default), \code{"pmi"} (positive mutual information), \code{"poiss"} (poisson-stirling), \code{"z"} (\emph{z}-score), \code{"locmi"} (local mutual information) and \code{"tmi"} (true mutual information). See below.}
  \item{raw}{Does the data frame contain raw observations, i.e., one observation per line? Default is \code{raw = TRUE}.}
  \item{all}{Should association be calculated for all possible combinations? The default is \code{FALSE}, which calculates association only for attested combinations. Please read the notes below for when \code{TRUE} makes sense and when it should be avoided.}
  \item{reverse}{The default sorting is in descending order by positive attraction. Set to \code{TRUE} for reverse sorting.}
  \item{threshold}{Set this to a positive integer if the \emph{output} should only contain combinations with a frequency above a certain threshold. Please read the notes below.}
  \item{decimals}{The number of decimals in the association measure.}
}
\details{
\bold{General:} This function uses a code section from the function \code{scfa} from the R package \code{cfa} (Mair & Funke 2017) for the calculation of expected values. Multiple Covarying Collexeme Analysis is conceptually essentially Configural Frequency Analysis (CFA; van Eye 1990) in a collostructional context. If you use this method, you can cite Stefanowitsch & Flach (2020).
\cr\cr
Note that this function can only provide measures based on the counts/observations in the data frame. That is, different to the other collexeme functions, there is no option to supply overall (corpus) frequencies if your data frame contains only a sample of all corpus counts. For instance, if you have removed hapax combinations from your data frame, then the frequencies of all types that occur in hapax types will not be included.
\cr\cr
\bold{All combinations:} While the calculation of association measures for all possible combinations of all conditions \code{all = TRUE} is necessary if you want to assess the relevance of the \emph{absence} of a combination (negative association value for most measures), you need to be aware of the consequences: for use cases with high type frequencies this will involve the calculation of a huge number of \emph{n}-types, which can break the R session. It is doubtful if it is linguistically relevant anyway: most high(er) frequent combinations that are linguistically interesting will have at least a few observations. \bold{Also note} that if you supply an aggregated data frame (with a frequency column, i.e., \code{raw = FALSE}), \code{all = TRUE} currently has no effect (and it's doubtful that this will be implemented). In this case, you can 'untable' your original data frame, see examples below.
\cr\cr
\bold{Threshold:} You can restrict the output such that only combinations that occur an \emph{n} number of times are shown. This might be useful if you have a large data frame.
\cr\cr
\bold{Association measures:} The implemented measures are relatively simple (as they only involve observed and expected values), but they do the trick of separating idiomatic, highly associated types from less strongly associated ones. Most measures are based on Evert (2004), with the exception of \code{tmi}. As there are, as of yet, no sufficient number of studies it is difficult to advise, but the \code{t-score} appears relatively robust (hence the default). However, since an observed value of 0 (if \code{all = TRUE}) would result in \code{-Inf}, 0.5 is added to the observed value of unattested combinations before the \code{t-value} is calculated.}

\value{
The output is sorted in decending order of attraction.
  \item{\emph{CONDITION}... }{The first \emph{n} columns are the conditions you originally entered.}
  \item{OBS}{Observed frequency of the combination}
  \item{EXP}{Expected frequency of the combination}
  \item{\emph{AM}}{The chosen association measure.}
}
\references{
Patrick Mair and Stefan Funke (2017). cfa: Configural Frequency Analysis (CFA). R package version
  0.10-0. \url{https://CRAN.R-project.org/package=cfa}
\cr\cr
Eye, A. von (1990). \emph{Introduction to configural frequency analysis. The search for types and anti-types in cross-classification}. Cambridge: CUP.
\cr\cr
Stefanowitsch, Anatol & Susanne Flach. 2020. \emph{Too big to fail} but \emph{big enough to pay for their mistakes}: A collostructional analysis of the patterns [too ADJ to V] and [ADJ enough to V]. In Gloria Corpas & Jean-Pierre Colson (eds.), \emph{Computational Phraseology}, 248–272. Amsterdam: John Benjamins.

}
\author{
Susanne Flach, susanne.flach@es.uzh.ch
}

\examples{\dontrun{
## Multiple Distinctive Collexeme Analysis

## Case 1: Raw list of observations
# load data

## Case 3: You only have an aggregated list with attested values,
## but want to measure the relevance of the absence (i.e., dissociation measures)
library(reshape)
# Untable (where df stands for your original data frame, minus the last column with the frequencies)
df_new <- untable(df[,-ncol(df)], num = df[, ncol(df)])
}
}
