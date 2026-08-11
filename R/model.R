## ============================================================
## LASAGNA model creation functions
## ============================================================
#' Create LASAGNA multi-layer graph model
#' Builds a multi-partite graph from multi-omics data layers.
#' Edges weighted by correlation, optionally conditioned on phenotype.
#' @param data A list with \code{X} (named list of data matrices),
#'   \code{samples} (data frame), and optionally \code{contrasts}.
#' @param X Named list of data matrices, one per layer. Alternative to
#'   \code{data} when \code{data} is not supplied.
#' @param meta Sample phenotype data frame (or contrasts matrix, per
#'   \code{meta.type}). Alternative to \code{data} when \code{data} is
#'   not supplied.
#' @param meta.type Phenotype type: \code{"pheno"}, \code{"expanded"},
#'   or \code{"contrasts"}.
#' @param ntop Number of top-SD features per layer. Set 0 or NULL
#'   to keep all.
#' @param nc Maximum number of inter-layer connections per node.
#' @param annot Optional annotation.
#' @param mask.gmt Use GMT structure for masking.
#' @param mask.graphite Use graphite PPI structure for masking.
#' @param add.sink Add SOURCE/SINK nodes.
#' @param intra Include intra-layer edges.
#' @param fully_connect Fully connect all layers.
#' @param add.revpheno Add reversed phenotype contrasts.
#' @param condition.edges Weight edges by phenotype correlation.
#' @return A list with components:
#' \describe{
#'   \item{graph}{An igraph object with layers, edge weights, and
#'     connection types.}
#'   \item{X}{Merged feature matrix.}
#'   \item{Y}{Phenotype matrix.}
#'   \item{layers}{Character vector of layer names.}
#' }
#' @examples
#' set.seed(1)
#' gx <- matrix(rnorm(20 * 10), 20, 10,
#'   dimnames = list(paste0("g", 1:20), paste0("S", 1:10)))
#' px <- matrix(rnorm(15 * 10), 15, 10,
#'   dimnames = list(paste0("p", 1:15), paste0("S", 1:10)))
#' samples <- data.frame(group = rep(c("A", "B"), each = 5),
#'   row.names = paste0("S", 1:10))
#' data <- list(X = list(gx = gx, px = px), samples = samples)
#' model <- create_model(data, ntop = 10, nc = 5)
#' model$layers
#' @export
create_model <- function(data,
                         X = NULL,
                         meta = NULL,
                         meta.type = "pheno",
                         ntop = 1000,
                         nc = 20,
                         annot = NULL,
                         mask.gmt = TRUE,
                         mask.graphite = TRUE,
                         add.sink = FALSE,
                         intra = TRUE,
                         fully_connect = FALSE,
                         add.revpheno = TRUE,
                         condition.edges = TRUE) {
  
  if(!is.null(data)) {
    X <- data$X
    meta <- data$samples
    if (meta.type == "contrasts") {
      if (!"contrasts" %in% names(data)) {
        stop("ERROR: contrasts missing in data")
      }
      meta <- data$contrasts
    }
  }
  
  if(is.null(data) && (is.null(X) || is.null(meta)) ) {
    stop("must supply data or {X, meta}.")
  }
  
  if (meta.type %in% c("pheno","samples")) {
    Y <- expandPhenoMatrix(meta, drop.ref = FALSE)
    if (is.null(Y)) {
      stop("[create_model] no phenotype column with resolvable/varying groups found in samples/meta; ",
           "check that it has a factor or character column with at least 2 distinct values")
    }
  } else if (meta.type %in% c("expanded","traits")) {
    Y <- 1 * meta
  } else if (meta.type == "contrasts") {
    Y <- makeContrastsFromLabelMatrix(meta)
    Y <- sign(Y)
    if (any(grepl("^IA:", colnames(Y)))) {
      Y <- Y[, grep("^IA:", colnames(Y), invert = TRUE), drop = FALSE]
    }
    if (add.revpheno) {
      revY <- -Y
      colnames(revY) <- reverse.AvsB(colnames(Y))
      Y <- cbind(Y, revY)
    }
  } else {
    stop("[create_model] invalid 'meta.type': ", meta.type)
  }
  X[["PHENO"]] <- t(Y)

  ## restrict number of features (by SD) if requested
  xx <- X
  xx <- lapply(xx, as.matrix)
  if (!is.null(ntop) && ntop > 0) {
    xx <- mofa.topSD(xx, ntop)
  }

  ## merge data (handles non-matching samples)
  xx <- mofa.merge_data2(xx, merge.rows = "prefix", merge.cols = "union")

  ## diagnose sample overlap across layers: warn (don't fail) when most
  ## samples lack data in at least one layer, since this silently degrades
  ## downstream correlation/model fit
  row.layer <- sub(":.*", "", rownames(xx))
  layer.names <- unique(row.layer)
  layer.coverage <- sapply(layer.names, function(lyr) {
    rows <- which(row.layer == lyr)
    colSums(!is.na(xx[rows, , drop = FALSE])) > 0
  })
  complete.frac <- mean(rowSums(layer.coverage) == length(layer.names))
  na.frac <- mean(is.na(xx))
  if (complete.frac < 0.5) {
    warning(sprintf(
      "create_model: merged data is %.0f%% NA after combining %d layers with only partial sample overlap; only %.0f%% of samples have data in every layer; results may be unreliable",
      100 * na.frac, length(layer.names), 100 * complete.frac
    ), call. = FALSE)
  }

  kk <- intersect(colnames(xx), rownames(Y))
  xx <- xx[, kk]
  Y <- Y[kk, ]

  ## add SOURCE/SINK
  if (add.sink) xx <- rbind(xx, "SOURCE" = 1, "SINK" = 1)

  ## compute correlation matrix
  suppressWarnings(R <- stats::cor(t(xx), use = "pairwise"))

  ## SOURCE/SINK fully connected
  ii <- grep("SINK|SOURCE", rownames(R))
  if (length(ii)) {
    R[ii, ] <- 1
    R[, ii] <- 1
  }

  R0 <- R
  R[is.na(R)] <- 0.01234 ## replace NA with small constant????

  ## condition edges by phenotype correlation
  if (condition.edges) {
    message("conditioning edges...")
    rho <- stats::cor(t(xx), Y, use = "pairwise.complete.obs")
    maxrho <- apply(abs(rho), 1, max, na.rm = TRUE)
    ii <- grep("SINK|SOURCE", names(maxrho))
    if (length(ii)) maxrho[ii] <- 1
    rho.wt <- outer(maxrho, maxrho)
    R <- R * rho.wt
  }

  ## define layers
  dt <- sub(":.*", "", rownames(R))
  layers <- names(X)
  if (add.sink) layers <- c("SOURCE", layers, "SINK")

  ## mask for inter-layer connections
  if (!fully_connect) {
    layer_mask <- matrix(0, nrow(R), ncol(R))
    dimnames(layer_mask) <- dimnames(R)
    if (length(layers) >= 2) {
      for (i in seq_len(length(layers) - 1)) {
        ii <- which(dt == layers[i])
        jj <- which(dt == layers[i + 1])
        layer_mask[ii, jj] <- 1
        layer_mask[jj, ii] <- 1
      }
    }
    if (intra) {
      for (i in seq_along(layers)) {
        ii <- which(dt == layers[i])
        layer_mask[ii, ii] <- 1
      }
    }
    R <- R * layer_mask
  }

  ## reduce inter-connections to nc top most correlated edges per node
  if (!is.null(nc) && nc > 0) {
    message("reducing edges to maximum ", nc, " connections")
    xtypes <- setdiff(layers, c("PHENO", "SOURCE", "SINK"))
    reduce_mask <- matrix(1, nrow(R), ncol(R))
    if (length(xtypes) >= 2) {
      for (i in seq_len(length(xtypes) - 1)) {
        ii <- which(dt == xtypes[i])
        jj <- which(dt == xtypes[i + 1])
        R1 <- R[ii, jj, drop = FALSE]
        rii <- apply(abs(R1), 1, function(r) utils::tail(sort(r), nc)[1])
        rjj <- apply(abs(R1), 2, function(r) utils::tail(sort(r), nc)[1])
        rr <- abs(R1) >= rii | t(t(abs(R1)) >= rjj)
        reduce_mask[ii, jj] <- rr
        reduce_mask[jj, ii] <- t(rr)
      }
    }
    if (intra) {
      for (i in seq_along(xtypes)) {
        ii <- which(dt == xtypes[i])
        R1 <- R[ii, ii, drop = FALSE]
        rii <- apply(abs(R1), 1, function(r) utils::tail(sort(r), nc)[1])
        rjj <- apply(abs(R1), 2, function(r) utils::tail(sort(r), nc)[1])
        rr <- abs(R1) >= rii | t(t(abs(R1)) >= rjj)
        reduce_mask[ii, ii] <- rr
      }
    }
    R <- R * reduce_mask
  }

  ## create graph from correlation
  R0[which(R == 0)] <- 0
  R0[is.na(R0)] <- 0
  gr <- igraph::graph_from_adjacency_matrix(R0, diag = FALSE,
    weighted = TRUE, mode = "undirected")

  igraph::E(gr)$rho <- igraph::E(gr)$weight
  gr$layers <- layers

  ## add edge connection type as attribute
  igraph::V(gr)$layer <- sub(":.*", "", igraph::V(gr)$name)
  ee <- igraph::as_edgelist(gr)
  etype <- apply(ee, 2, function(e) sub(":.*", "", e))
  etype.idx <- apply(etype, 2, match, gr$layers)
  rev.etype <- etype.idx[, 2] < etype.idx[, 1]
  etype1 <- ifelse(rev.etype, etype.idx[, 2], etype.idx[, 1])
  etype2 <- ifelse(rev.etype, etype.idx[, 1], etype.idx[, 2])
  etype1 <- gr$layers[etype1]
  etype2 <- gr$layers[etype2]
  igraph::E(gr)$connection_type <- paste0(etype1, "->", etype2)
  ii <- which(etype1 == etype2)
  if (length(ii)) igraph::E(gr)$connection_type[ii] <- etype1[ii]

  return(list(graph = gr, X = xx, Y = Y, layers = layers))

}
