#' @export
join.lists <- function(x, y, all = TRUE, threshold = 1) {

  if (is.data.frame(x[[1]]) == FALSE)
    stop("invalid 'df' argument -- must be dataframe")
  if (is.data.frame(y[[2]]) == FALSE)
    stop("invalid 'df' argument -- must be dataframe")

  comb.list = vector("list", length(x))

  for (i in 1:length(x)) {

    merged <- merge(x[[i]], y[[i]], by = 1, all = TRUE)
    merged[is.na(merged)] <- 0
    names(merged) <- c(names(x[[1]][1]),
                       deparse(substitute(x)),
                       deparse(substitute(y)))
    merged <- merged[order(-merged[, 2]), ]
    rownames(merged) <- NULL
    comb.list[[i]] <- merged

  }

  if (all == FALSE) {

    comb.list <- lapply(comb.list, function(x) x[x[, 2] > 0, ])

  }

  names(comb.list) <- names(x)
  comb.list <- lapply(comb.list, function(x) x[x[, 2] + x[, 3] >= threshold, ])
  comb.list <- lapply(comb.list, function(x) droplevels(x))
  comb.list

}

