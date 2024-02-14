\name{reshape.cca}
\alias{reshape.cca}
\title{
Function to transform collex.covar() long to wide format
}
\description{
The output of \code{\link{collex.covar}} is a so-called 'long format', where each row represents a slot1~slot2 pair (or, more generally, a cond1~cond2 pair). This function cross-tabulates the association measures in a 'wide format', where one condition occurs as rows and the other represents the columns and the cells contain the row~col association measure.
}
\usage{
reshape.cca(x, cond = "shorter", str.dir = TRUE, value = "COLL.STR", abs.dev = FALSE,
            max.assoc = FALSE, sorton = "abs.dev", decimals = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data frame containing the output of \code{\link{collex.covar}}.}
  \item{cond}{Which of the two conditions in \code{\link{collex.covar}} should be the columns of the reshaped format? The default \code{shorter} uses that condition which has the fewer number of types. But you can override this by setting with \code{cond = 1} for the first condition (usually slot 1) or \code{cond = 2} for the second condition (usually slot 2).}
  \item{str.dir}{Should the values in the cells indicate the direction of association? \code{TRUE} is the default and should be used if the association measure you used in \code{\link{collex.covar}} is one that contains only positive values (e.g., \code{logl}, \code{fye}, or \code{fye.ln}). For other association measures, set \code{str.dir = FALSE}, if the values range from - to +.}
  \item{value}{Which value should be cross-tabulated, i.e., put in the cells? For cond1~cond2 pairs that do not occur in \code{x}, \code{NA} are returned (see details below). The default cross-tabulates the association measure, but it can also be \code{'OBS'} for a co-occurrence matrix. Note, since for it does not make sense to calculate a directed value for \code{'OBS'}, the function will abort if you leave the default for \code{str.dir = TRUE} unchanged when using \code{'OBS'}.}
  \item{abs.dev}{The function sums all absolute association measures row-wise Should this be included in the output as an extra column? The default is \code{FALSE}. Note that if your \code{x} only contains attested values, the summed deviations ignore \code{NA} values.}
  \item{max.assoc}{Should the col-condition of maximum assocation per row-condition be included in the output? If \code{TRUE}, this column will contain the col-type with which the row-type has the highest positive or negative association (i.e., 'most deviant' relationship).}
  \item{sorton}{By default the output is sorted in descending order of \code{abs.dev}. You can also sort alphabetically (\code{sorton = "alphabetical")}.}
  \item{decimals}{Rounding of cell values. If this is set to a higher value than what was used for this argument in \code{\link{collex.covar}}, it will of course be ignored.}
}
\details{
The function makes most sense for a \code{\link{collex.covar}} that was run for all possible combinations. If association scores were only calculated for attested combinations, the output of \code{reshape.cca} contains \code{NA} in the cells of unattested combinations and it is up to the user to decide what to do with it. Both \code{abs.dev} and \code{max.assoc} can still be calculated and displayed, but they are based on observed combinations only. Since association measures for unobserved combinations can be read as 'negative evidence', the \code{abs.dev} will be and the \code{max.assoc} type may be different, depending on the strength of (potential) 'negative association'. See examples below for the case of unattested \emph{*}.
}
\value{
Returns cross-tabulated association scores or observed values.
}
\author{
Susanne Flach, susanne.flach@es.uzh.ch
}
\seealso{
\code{\link{collex.covar}}
}
\examples{\dontrun{
data(modadv)

## Inspect:
# This is an aggregated frequency list:
head(modadv, 12)

### Perform co-varying collexeme analysis
## ?collex.covar()
# since it's aggregated, you must set raw = FALSE, or it will make no sense.
cca.att <- collex.covar(modadv, am="fye", raw = FALSE, all = FALSE) # only attested combinations
cca.all <- collex.covar(modadv, am="fye", raw = FALSE, all = TRUE) # all combinations

## Reshape the cca output by association measure:
# ?reshape.cca
cca.wide.att <- reshape.cca(cca.att)
View(cca.wide.att)
cca.wide.all <- reshape.cca(cca.all)
View(cca.wide.all)

#### Co-occurrence of observations
modadv.obs <- reshape.cca(cca.att, value = "OBS", str.dir = FALSE) # you must set false in this case
# since we ran this on only the attested values, you can replace NA with 0:
modadv.obs[is.na(modadv.obs)] <- 0
View(modadv.obs)

}}
