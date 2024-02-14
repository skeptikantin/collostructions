#' @export
freq.list <- function(x, convert = TRUE, asFactor = TRUE) {

  if (!(class(x) %in% c("factor", "character"))) {
       stop("'x' must be charactor or factor vector.")
    }

  if (convert == TRUE) {
    x <- tolower(x)
  }

  if (is.character(x) == FALSE) {
    x <- as.character(x)
  }

  tabs <- rle(sort(x))
  counts <- data.frame(WORD = tabs$values, FREQ = tabs$lengths)
  counts <- counts[with(counts, order(FREQ, decreasing = TRUE)), ]
  rownames(counts) <- NULL

  ## convert to factor in R.4x:
  if (asFactor == TRUE) {
    counts$WORD <- as.factor(counts$WORD)
  }

  # outpus:
  counts

}
