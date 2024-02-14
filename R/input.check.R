#' @export
#' @importFrom utils head

input.check <- function(x, raw = FALSE) {

  # check if raw setting is correct
  if (ncol(x) == 2 & !(is.numeric(x[,2]))) {
#    stop("It appears your x is a data frame with two non-numerical variables,
#    which indicates a raw format (i.e., one observation per line. Set raw = TRUE)")
    raw = TRUE
    cs = data.frame()
  } else {
    raw = FALSE
    cs = data.frame()
  }

  # check if a three-column input that has raw = T
  if (raw == TRUE & ncol(x) > 2) {
    if (raw == TRUE & ncol(x) == 3 & is.numeric(x[,2]) & is.numeric(x[,3])) {
      print(head(x, 4))
      stop("Your data frame has one or more numeric columns, so it is probably not a raw format. Set raw = FALSE).")
      }
    }

  # some function below requires factors so convert all characters to factors:
  x[sapply(x, is.character)] <- lapply(x[sapply(x, is.character)], as.factor)

  # Determine columns to be checked:
  if ((is.factor(x[,1]) || is.character(x[,1])) & (is.factor(x[,2]) || is.character(x[,2]))) {
    cols = 2
  } else {
    cols = 1
  }

  ## Check if it contains leading whitespace
  # leading whitespace in col1?
  if (nrow(x[grepl("^ +.*", x[, 1]), ]) > 0) {
    leading1 <- as.integer(rownames(x[grepl("^ +.*", x[, 1]), ]))
  } else {
    leading1 <- NULL
  }
  # trailing whitespace in col1?
  if (nrow(x[grepl(" +$", x[, 1]), ]) > 0){
    trailing1 <- as.integer(rownames(x[grepl(" +$*", x[, 1]), ]))
  } else {
    trailing1 <- NULL
  }

  if (cols == 2) {
    if (nrow(x[grepl("^ +.*", x[, 2]), ]) > 0) {
      leading2 <- as.integer(rownames(x[grepl("^ +.*", x[, 2]), ]))
    } else {
      leading2 <- NULL
    }
    if (nrow(x[grepl(" +$", x[, 2]), ]) > 0) {
      trailing2 <- as.integer(rownames(x[grepl(" +$*", x[, 2]), ]))
    } else {
      trailing2 <- NULL
    }
  }

  ## combine all
  if (cols == 1) {
    # combine vectors of leading and trailing rows
    #if (leading1)
    rows.to.subset <- c(leading1, trailing1)
    # sort and remove doubles
  } else if (cols == 2) {
    rows.to.subset <- c(leading1, leading2, trailing1, trailing2)
  }

  # Now clean and subset the original dataframe
  rows.to.subset <- unique(sort(rows.to.subset))
  whitespace <- x[rows.to.subset,]

    if (raw == FALSE) {

    # Check if either column contains doubles b/c of case-sensitive
      y <- x # make copy
      y[, 1] <- tolower(y[,1]) # lower the types in col 1
      y[, 1] <- gsub("(^ +| +$)", "", y[, 1]) # remove leading/trailing whitespace
      y.types <- freq.list(y[,1]) # make frequency list
      y.types <- y.types[y.types$FREQ > 1,]$WORD # extract all types with freq greater 1
      rows.to.subset.cs <- as.numeric(rownames(subset(y, y[,1] %in% y.types))) # i
      cs <- x[rows.to.subset.cs,]

      if (cols == 2) {
        y <- x # make a new copy
        y[, 2] <- tolower(y[,2]) # lower the types in col 1
        y[, 2] <- gsub("(^ +| +$)", "", y[, 2]) # remove leading/trailing whitespace
        y.types <- freq.list(y[, 2]) # make frequency list
        y.types <- y.types[y.types$FREQ > 1,]$WORD # extract all types with freq greater 1
        rows.to.subset.cs2 <- as.numeric(rownames(subset(y, y[,2] %in% y.types))) # i
        rows.to.subset.cs2 <- unique(sort(c(rows.to.subset, rows.to.subset.cs2)))
        cs <- x[rows.to.subset.cs2,]
      }
    } else {

      # come up with a solution to check case-sensitivity in raw = F
      y <- x # make copy
      y1 <- freq.list(y[, 1], convert = FALSE)
      y1 <- freq.list(y1[, 1])
      y1.types <- droplevels(y1[y1$FREQ > 1, ]$WORD)
      y1.types <- droplevels(unique(subset(x, tolower(x[ ,1]) %in% y1.types)[,1]))

      y2 <- freq.list(y[, 2], convert = FALSE)
      y2 <- freq.list(y2[, 1])
      y2.types <- droplevels(y2[y2$FREQ > 1, ]$WORD)
      y2.types <- droplevels(unique(subset(x, tolower(x[ ,2]) %in% y2.types)[,2]))
  }

#    if((nrow(whitespace) == 0 & nrow(whitespace) == 0 & length(y1.types) == 0 & length(y2.types == 0)))
#    stop("No issues detected.")

  message("There is leading and/or trailing whitespace in these rows:")
  if (nrow(whitespace) == 0) {
    cat("None detected.\n")
  } else {
    print(whitespace)
  }

  if (raw == FALSE) {
    message("There is a potential problem with duplicate types (whitespace) in these rows:")
    if (nrow(whitespace) == 0) {
      duplicate_problem <- FALSE
      cat("None detected.\n")
      } else if(nrow(cs) == 0) {
        duplicate_problem <- FALSE
        cat("None detected.\n")
      } else {
        duplicate_problem <- TRUE
        print(cs)
      }
  } else {
    message("There is a potential problem with duplicate types (case-sensitivity):")
    if (length(y1.types) == 0 & length(y2.types) == 0) {
      duplicate_problem <- FALSE
      duplicate_problem2 <- FALSE
      cat("None detected.\n")
    } else {
      duplicate_problem <- TRUE
      duplicate_problem2 <- TRUE
      if (length(y1.types) > 0 & length(y2.types) == 0) {
        cat("Column 1:", paste(y1.types, collapse = ", "), "\n")
        cat("Column 2: none detected.\n")
      } else if (length(y1.types) > 0 & length(y2.types) > 0) {
        cat("Column 1:", paste(y1.types, collapse = ", "), "\n")
        cat("Column 2:", paste(y2.types, collapse = ", "), "\n")
      } else if (length(y1.types) == 0 & length(y2.types) > 0) {
        cat("Column 1: none detected.\n")
        cat("Column 2:", paste(y2.types, collapse = ", "), "\n")

      }
    }
  }

  if(raw == FALSE) {
    message("There is a potential problem with duplicate types (case-sensitivity) in these rows:")
    if(length(y.types)!=0) {
      duplicate_problem2 <- TRUE
      print(subset(x, tolower(x[,1]) %in% y.types | tolower(x[,2]) %in% y.types))
      } else {
      duplicate_problem2 <- FALSE
      cat("None detected.\n")
      }
  }

  message("There are incomplete cases (e.g. missing frequencies and/or types) in these rows:")
    if (any(is.na(x))) {
      incomplete_problem <- TRUE
      print(rbind(x[!(complete.cases(x)),], x[which(x[,c(1,2)]==""),]))
    } else if (nrow(x[x[,1]=="" | x[,2]=="",])>0) {
      print(x[x[,1]=="" | x[,2]=="",])
      incomplete_problem <- TRUE
    } else {
      cat("None detected.\n")
      incomplete_problem <- FALSE
    }

  ## If no problems are detected
  if (cols == 1 & nrow(whitespace) == 0 & nrow(cs) == 0 & !(any(is.na(x)))) {
    message("Summary: From what the function tests, things *appear* fine.")
  } else if (raw == TRUE & duplicate_problem == FALSE & duplicate_problem2 == FALSE & incomplete_problem == FALSE) {
    message("Summary: From what the function tests, things *appear* fine.")
  } else {
    message("Summary: You should fix these issues before running any collex analysis.
Also make sure that frequencies are correctly aggregated after fixing typo-/orthographic issues.")

  }

}
