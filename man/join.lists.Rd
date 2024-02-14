\name{join.lists}
\alias{join.lists}
\title{Function to merge two lists of frequency lists
}
\description{
Merges two lists of data frames pair-wise. Both lists need to be of equal length (i.e. identical number of data.frames), and where all data.frames have two columns (for \code{WORD} and \code{FREQ}); merging is by items in x[,1]. All pair-wise data frames must have an identical ID column name to match by (using \code{WORD} as a column name is recommended). Returns a list of the same length with data frames of length 3, named \code{WORD} (or name of ID column), \code{name.of.df1} and \code{name.of.df2}; the latter two contain the frequencies of \code{WORD}. Suitable to handle the output of \code{split()}, e.g. if a data frame was split by time period.}
\usage{
join.lists(x, y, all = TRUE, threshold = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{List 1 containing data frames of frequencies for construction A.}
  \item{y}{List 2 containing data frames of (i) corpus frequiences for \code{WORD} in construction A or (ii) frequencies of \code{WORD} in construction B.}
  \item{all}{If \code{all = TRUE} (the default), each data frame contains all joint types of both columns (esp. for \code{collex.dist()}). Set to \code{all = FALSE}, if only the types in the first condition (of the first data frame list) should be included (this is often the case for \code{collex()}).}
    \item{threshold}{Numerical. How many times must an item occur overall to be included in the joint list? Default is 1, which means all items are included (per list split, e.g., time period). If \code{all = FALSE}, setting \code{threshold = number} will mean that items in \code{x} are only included if they occur at least \code{number} times.}

}
\examples{\dontrun{

### We performed a series of Distictive Collexeme Analyses
### for (present) progressive vs. (simple) present over ten
### 25-year periods in CLMET (Johannsen & Flach 2015).
### Note that although working with historical data, this is something
### very different to Diachronic Distinctive Collexeme Analysis
### (e.g., Hilpert 2011), where periods are conditions in *one* DCA
### and thus mathematically not independent of each other.
### The sample data below runs one DCA per period, which are
### mathematically independent of each other. The conditions are still
### two alternating constructions as in 'ordinary' DCA.

### Also note that this means 'multiple' DCAs in the sense of 'several' DCAs,
### not in the sense of 'Multiple Distinctive Collexeme Analysis' (MDCA).

# Load data
data(CLMETprog.qc)
data(CLMETsimple.qc)
head(CLMETprog.qc, 10)
head(CLMETsimple.qc, 10)

### I. Prepare
## Make frequency lists by decade of class list, i.e.,
## split constructions by period,
## keep ITEM & FREQ only, droplevels
prog <- split(CLMETprog.qc[, c(1,3)], CLMETprog.qc$QUARTCENT)
prog <- lapply(prog, droplevels)
simp <- split(CLMETsimple.qc[, c(1,3)], CLMETsimple.qc$QUARTCENT)
simp <- lapply(simp, droplevels(x)

dist.in <- join.lists(prog, simp)
dist.in <- lapply(dist.in, droplevels)
# Cosmetics:
dist.in <- lapply(dist.in, setNames, c("WORD", "progressive", "simple"))

#### CALCULATE COLLEXEMES
dist.out.log <- lapply(prog.collexDist.in, function(x) collex.dist(x))
dist.out.fye <- Map(collex.dist, dist.in, am="fye")

### EXPORT
## Note: for this strategy, you need to install and load library(openxlsx)
write.xlsx(dist.out.log, "progCollexDistLL.xlsx")
write.xlsx(dist.out.log, "progCollexDistFYE.xlsx")
}
}

