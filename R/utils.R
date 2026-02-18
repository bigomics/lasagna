## ============================================================
## Internal utility functions copied from playbase for
## standalone operation. These are NOT exported.
## ============================================================

mofa.topSD <- function(xdata, ntop) {
  if (is.list(xdata)) {
    res <- lapply(xdata, function(x) {
      sdx <- matrixStats::rowSds(x, na.rm = TRUE)
      head(x[order(-sdx), , drop = FALSE], ntop)
    })
  } else if (is.matrix(xdata)) {
    if (all(grepl(":", rownames(xdata)))) {
      xdata <- mofa.split_data(xdata)
      res <- lapply(xdata, function(x) {
        sdx <- matrixStats::rowSds(x, na.rm = TRUE)
        head(x[order(-sdx), , drop = FALSE], ntop)
      })
      res <- mofa.merge_data(res)
    } else {
      sdx <- matrixStats::rowSds(xdata, na.rm = TRUE)
      res <- head(xdata[order(-sdx), , drop = FALSE], ntop)
    }
  } else {
    message("[mofa.topSD] WARNING: could not detect type")
    res <- xdata
  }
  return(res)
}

mofa.split_data <- function(X, keep.prefix = FALSE) {
  if (!all(grepl("[:]|SOURCE|SINK", rownames(X)))) {
    rownames(X) <- paste0("x:", rownames(X))
  }
  dtype <- sub(":.*", "", rownames(X))
  xx <- tapply(1:nrow(X), dtype, function(i) X[i, , drop = FALSE])
  if (!keep.prefix) {
    xx <- mofa.strip_prefix(xx)
  }
  xx
}

mofa.merge_data <- function(xx) {
  do.call(rbind, mofa.prefix(xx))
}

mofa.merge_data2 <- function(xdata, merge.rows = "prefix", merge.cols = "union") {
  n1 <- length(Reduce(intersect, lapply(xdata, rownames)))
  n2 <- length(Reduce(intersect, lapply(xdata, colnames)))
  rdim <- sapply(xdata, nrow)
  cdim <- sapply(xdata, ncol)
  if (n1 < min(rdim) && merge.rows != "prefix") {
    message("WARNING: rows do not match")
  }
  if (n2 < min(cdim) && merge.cols != "prefix") {
    message("WARNING: columns do not match")
  }
  prefix.rows <- (merge.rows == "prefix")
  prefix.cols <- (merge.cols == "prefix")
  if (prefix.cols) {
    for (i in 1:length(xdata)) {
      nn <- sub("^[A-Za-z]+:", "", colnames(xdata[[i]]))
      colnames(xdata[[i]]) <- paste0(names(xdata)[i], ":", nn)
    }
    merge.cols <- "union"
  }
  if (prefix.rows) {
    for (i in 1:length(xdata)) {
      nn <- sub("^[A-Za-z]+:", "", rownames(xdata[[i]]))
      rownames(xdata[[i]]) <- paste0(names(xdata)[i], ":", nn)
    }
    merge.rows <- "union"
  }
  if (merge.rows == "intersect") {
    allfeatures <- Reduce(intersect, lapply(xdata, rownames))
  } else {
    allfeatures <- unique(unlist(lapply(xdata, rownames)))
  }
  if (merge.cols == "intersect") {
    allsamples <- Reduce(intersect, lapply(xdata, colnames))
  } else {
    allsamples <- unique(unlist(lapply(xdata, colnames)))
  }
  D <- matrix(0, length(allfeatures), length(allsamples))
  nn <- matrix(0, length(allfeatures), length(allsamples))
  rownames(D) <- allfeatures
  colnames(D) <- allsamples
  for (i in 1:length(xdata)) {
    A <- xdata[[i]]
    ii <- match(rownames(D), rownames(A))
    jj <- match(colnames(D), colnames(A))
    A1 <- A[ii, jj]
    nn <- nn + !is.na(A1) * 1
    A1[is.na(A1)] <- 0
    D <- D + A1
  }
  D <- D / nn
  D[which(nn == 0)] <- NA
  rownames(D) <- allfeatures
  colnames(D) <- allsamples
  return(D)
}

mofa.prefix <- function(xx) {
  xx <- mofa.strip_prefix(xx)
  for (i in 1:length(xx)) {
    dt <- paste0(names(xx)[i], ":")
    if (is.null(dim(xx[[i]]))) {
      names(xx[[i]]) <- paste0(dt, names(xx[[i]]))
    } else {
      rownames(xx[[i]]) <- paste0(dt, rownames(xx[[i]]))
    }
  }
  xx
}

mofa.get_prefix <- function(x) {
  if (inherits(x, c("matrix", "data.frame")) || !is.null(dim(x))) {
    x <- rownames(x)
  }
  ifelse(grepl(":", x), sub(":.*", "", x), "")
}

mofa.strip_prefix <- function(xx) {
  if (is.character(xx)) {
    xx <- sub("^[A-Za-z0-9]+:", "", xx)
    return(xx)
  }
  if (is.matrix(xx)) {
    rownames(xx) <- sub("^[A-Za-z0-9]+:", "", rownames(xx))
    return(xx)
  }
  if (is.list(xx)) {
    for (i in 1:length(xx)) {
      dt <- paste0("^", names(xx)[i], ":")
      if (is.null(dim(xx[[i]]))) {
        names(xx[[i]]) <- sub(dt, "", names(xx[[i]]))
      } else {
        rownames(xx[[i]]) <- sub(dt, "", rownames(xx[[i]]))
      }
    }
    return(xx)
  }
  xx
}

expandPhenoMatrix <- function(M, drop.ref = TRUE, keep.numeric = FALSE, check = TRUE) {
  a1 <- tidy.dataframe(M)
  nlevel <- apply(a1, 2, function(x) length(setdiff(unique(x), NA)))
  nterms <- colSums(!is.na(a1))
  nratio <- nlevel / nterms
  if (inherits(a1, "data.frame")) {
    a1.typed <- utils::type.convert(a1, as.is = TRUE)
    y.class <- sapply(a1.typed, function(a) class(a)[1])
  } else {
    a1.typed <- utils::type.convert(a1, as.is = TRUE)
    y.class <- apply(a1.typed, 2, function(a) class(a)[1])
  }

  is.fac <- rep(FALSE, ncol(a1))
  is.int <- (y.class == "integer")
  ii <- which(is.int)
  is.fac[ii] <- apply(a1[, ii, drop = FALSE], 2, function(x) {
    nlev <- length(unique(x[!is.na(x)]))
    max(x, na.rm = TRUE) %in% c(nlev, nlev - 1)
  })
  is.fac2 <- (y.class == "integer" & nlevel <= 3 & nratio < 0.66)
  y.class[is.fac | is.fac2] <- "character"

  y.isnum <- (y.class %in% c("numeric", "integer"))
  kk <- which(y.isnum | (!y.isnum & nlevel > 1 & nratio < 0.66))
  if (length(kk) == 0) {
    kk <- which(y.isnum | (!y.isnum & nlevel > 1))
  }
  if (length(kk) == 0) {
    return(NULL)
  }
  a1 <- a1[, kk, drop = FALSE]
  a1.isnum <- y.isnum[kk]

  m1 <- list()
  for (i in 1:ncol(a1)) {
    if (a1.isnum[i]) {
      suppressWarnings(x <- as.numeric(a1[, i]))
      if (keep.numeric) {
        m0 <- matrix(x, ncol = 1)
        colnames(m0) <- "#"
      } else {
        if (drop.ref) {
          m0 <- matrix((x > stats::median(x, na.rm = TRUE)), ncol = 1)
          colnames(m0) <- "high"
        } else {
          mx <- stats::median(x, na.rm = TRUE)
          m0 <- matrix(cbind(x <= mx, x > mx), ncol = 2)
          colnames(m0) <- c("low", "high")
        }
      }
    } else if (drop.ref && nlevel[i] == 2) {
      x <- as.character(a1[, i])
      x1 <- utils::tail(sort(x), 1)
      m0 <- matrix(x == x1, ncol = 1)
      colnames(m0) <- x1
    } else {
      x <- as.character(a1[, i])
      x[is.na(x) | x == "NA" | x == " "] <- "_"
      m0 <- stats::model.matrix(~ 0 + x)
      colnames(m0) <- sub("^x", "", colnames(m0))
    }
    rownames(m0) <- rownames(a1)
    if ("_" %in% colnames(m0)) {
      m0 <- m0[, -which(colnames(m0) == "_")]
    }
    m1[[i]] <- m0
  }

  names(m1) <- colnames(a1)
  for (i in 1:length(m1)) {
    colnames(m1[[i]]) <- paste0(names(m1)[i], "=", colnames(m1[[i]]))
  }
  m1 <- do.call(cbind, m1)
  colnames(m1) <- sub("=#", "", colnames(m1))
  rownames(m1) <- rownames(M)
  return(m1)
}

tidy.dataframe <- function(Y) {
  Y <- Y[, which(colMeans(is.na(Y)) < 1), drop = FALSE]
  Y <- apply(Y, 2, function(x) sub("^NA$", NA, x))
  Y <- Y[, which(colMeans(is.na(Y)) < 1), drop = FALSE]
  Y <- apply(Y, 2, function(x) gsub("^[ ]*|[ ]*$", "", x))
  suppressWarnings(num.Y <- apply(Y, 2, function(x) as.numeric(as.character(x))))
  is.numeric <- (0.8 * colMeans(is.na(num.Y)) <= colMeans(is.na(Y)))
  nlevel <- apply(Y, 2, function(x) length(unique(x)))
  is.factor <- (!is.numeric | (is.numeric & nlevel <= 3))
  is.factor <- (is.factor | grepl("batch|replicat|type|clust|group", colnames(Y)))
  new.Y <- data.frame(Y, check.names = FALSE)
  new.Y[, which(is.numeric)] <- num.Y[, which(is.numeric), drop = FALSE]
  for (i in which(is.numeric)) new.Y[[i]] <- num.Y[, i]
  for (i in which(is.factor)) new.Y[[i]] <- factor(as.character(new.Y[, i]))
  new.Y <- data.frame(new.Y, check.names = FALSE)
  return(new.Y)
}

makeContrastsFromLabelMatrix <- function(lab.matrix) {
  if (!all(grepl("_vs_", colnames(lab.matrix)))) {
    stop("[makeContrastsFromLabelMatrix] FATAL:: all contrast names must include _vs_")
  }

  ct.names <- colnames(lab.matrix)
  main.grp <- sapply(strsplit(ct.names, split = "_vs_"), "[", 1)
  ctrl.grp <- sapply(strsplit(ct.names, split = "_vs_"), "[", 2)
  main.grp <- sub(".*:", "", main.grp)
  ctrl.grp <- sub("@.*", "", ctrl.grp)

  contr.mat <- matrix(0, nrow(lab.matrix), ncol(lab.matrix))
  rownames(contr.mat) <- rownames(lab.matrix)
  colnames(contr.mat) <- colnames(lab.matrix)
  for (i in 1:ncol(lab.matrix)) {
    lab1 <- trimws(lab.matrix[, i])
    lab1x <- setdiff(lab1, c(NA, ""))
    grps <- c(main.grp[i], ctrl.grp[i])
    if (all(lab1x %in% grps)) {
      j1 <- which(lab1 == main.grp[i])
      j0 <- which(lab1 == ctrl.grp[i])
    } else {
      j1 <- grep(paste0("^", toupper(main.grp[i])), toupper(lab1))
      j0 <- grep(paste0("^", toupper(ctrl.grp[i])), toupper(lab1))
    }
    contr.mat[j1, i] <- +1 / length(j1)
    contr.mat[j0, i] <- -1 / length(j0)
  }

  return(contr.mat)
}

reverse.AvsB <- function(comp) {
  reverse.AvsB.1 <- function(comp) {
    prefix <- postfix <- ""
    if (any(grepl("[:]", comp))) {
      after1 <- sub(".*:", "", comp)
      prefix <- sub(paste0(":", after1), "", comp, fixed = TRUE)
    }
    if (any(grepl("[@]", comp))) postfix <- sub(".*@", "", comp)
    comp0 <- gsub(".*:|@.*", "", comp)
    ab <- paste(rev(strsplit(comp0, split = "_vs_|_VS_")[[1]]), collapse = "_vs_")
    gsub("^:|@$", "", paste0(prefix, ":", ab, "@", postfix))
  }
  as.character(sapply(comp, reverse.AvsB.1))
}

uscale <- function(x, symm = FALSE) {
  uscale.func <- function(x) (x - min(x)) / (max(x) - min(x) + 1e-99)
  if (is.matrix(x)) {
    y <- apply(x, 2, uscale.func)
    if (nrow(x) == 1) {
      y <- matrix(y, nrow = 1)
      dimnames(y) <- dimnames(x)
    }
  } else {
    y <- uscale.func(x)
  }
  y[is.na(y)] <- NA
  if (symm) y <- (y - 0.5) * 2
  return(y)
}

iconv2utf8 <- function(s) {
  iconv(s, to = "UTF-8//TRANSLIT", sub = "")
}
