#' @export
join.freqs <- function (x, y, all = TRUE, threshold = 1) {

  if (is.data.frame(x) == FALSE)
    stop("invalid 'df' argument - must be dataframe")
  if (is.data.frame(y) == FALSE)
    stop("invalid 'df' argument - must be dataframe")

  merged <- merge(x, y, by = 1, all = TRUE)
  merged[is.na(merged)] <- 0
  names(merged) <- c(names(x[1]), deparse(substitute(x)),
                     deparse(substitute(y)))
  merged <- merged[order(-merged[, 2]), ]
  rownames(merged) <- NULL

  if (all == FALSE) {
    merged <- subset(merged, merged[, 2] > 0)
    merged <- droplevels(merged)
    merged
  } else {
    merged
  }
  merged <- subset(merged, merged[, 2] + merged[, 3] >= threshold)
  rownames(merged) <- NULL
  merged
}

