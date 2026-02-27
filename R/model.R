## ============================================================
## LASAGNA model creation functions
## ============================================================
#' Create LASAGNA multi-layer graph model
#' Builds a multi-partite graph from multi-omics data layers.
#' Edges weighted by correlation, optionally conditioned on phenotype.
#' @param data A list with \code{X} (named list of data matrices),
#'   \code{samples} (data frame), and optionally \code{contrasts}.
#' @param pheno Phenotype type: \code{"pheno"}, \code{"expanded"},
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
#' @export
create_model <- function(data,
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

  if (meta.type %in% c("pheno","samples")) {
    Y <- expandPhenoMatrix(data$samples, drop.ref = FALSE)
  } else if (meta.type %in% c("expanded","traits")) {
    Y <- 1 * data$samples
  } else if (meta.type == "contrasts") {
    if (!"contrasts" %in% names(data)) {
      message("ERROR: contrasts missing in data")
      return(NULL)
    }
    Y <- makeContrastsFromLabelMatrix(data$contrasts)
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
    message("[create_model] ERROR invalid meta.type type")
    return(NULL)
  }
  data$X[["PHENO"]] <- t(Y)

  ## restrict number of features (by SD) if requested
  xx <- data$X
  if (!is.null(ntop) && ntop > 0) {
    xx <- lapply(xx, function(x) head(x[order(-apply(x, 1, stats::sd)), , drop = FALSE], ntop))
    xx <- mofa.topSD(xx, ntop)
  }

  ## merge data (handles non-matching samples)
  X <- mofa.merge_data2(xx, merge.rows = "prefix", merge.cols = "union")
  kk <- intersect(colnames(X), rownames(Y))
  X <- X[, kk]
  Y <- Y[kk, ]

  ## add SOURCE/SINK
  if (add.sink) X <- rbind(X, "SOURCE" = 1, "SINK" = 1)

  ## compute correlation matrix
  suppressWarnings(R <- stats::cor(t(X), use = "pairwise"))

  ## SOURCE/SINK fully connected
  ii <- grep("SINK|SOURCE", rownames(R))
  if (length(ii)) {
    R[ii, ] <- 1
    R[, ii] <- 1
  }

  R0 <- R
  R[is.na(R)] <- 0.1234 ## replace NA with constant

  ## condition edges by phenotype correlation
  if (condition.edges) {
    message("conditioning edges...")
    rho <- stats::cor(t(X), Y, use = "pairwise.complete.obs")
    maxrho <- apply(abs(rho), 1, max, na.rm = TRUE)
    ii <- grep("SINK|SOURCE", names(maxrho))
    if (length(ii)) maxrho[ii] <- 1
    rho.wt <- outer(maxrho, maxrho)
    R <- R * rho.wt
  }

  ## define layers
  dt <- sub(":.*", "", rownames(R))
  layers <- names(data$X)
  if (add.sink) layers <- c("SOURCE", layers, "SINK")

  ## mask for inter-layer connections
  if (!fully_connect) {
    layer_mask <- matrix(0, nrow(R), ncol(R))
    dimnames(layer_mask) <- dimnames(R)
    for (i in 1:(length(layers) - 1)) {
      ii <- which(dt == layers[i])
      jj <- which(dt == layers[i + 1])
      layer_mask[ii, jj] <- 1
      layer_mask[jj, ii] <- 1
    }
    if (intra) {
      for (i in 1:length(layers)) {
        ii <- which(dt == layers[i])
        layer_mask[ii, ii] <- 1
      }
    }
    R <- R * layer_mask
  }

  ## reduce inter-connections to nc top most correlated edges per node
  if (!is.null(nc) && nc > 0) {
    message(paste("reducing edges to maximum", nc, "connections"))
    xtypes <- setdiff(layers, c("PHENO", "SOURCE", "SINK"))
    reduce_mask <- matrix(1, nrow(R), ncol(R))
    for (i in 1:(length(xtypes) - 1)) {
      ii <- which(dt == xtypes[i])
      jj <- which(dt == xtypes[i + 1])
      R1 <- R[ii, jj, drop = FALSE]
      rii <- apply(abs(R1), 1, function(r) utils::tail(sort(r), nc)[1])
      rjj <- apply(abs(R1), 2, function(r) utils::tail(sort(r), nc)[1])
      rr <- abs(R1) >= rii | t(t(abs(R1)) >= rjj)
      reduce_mask[ii, jj] <- rr
      reduce_mask[jj, ii] <- t(rr)
    }
    if (intra) {
      for (i in 1:length(xtypes)) {
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

  return(list(graph = gr, X = X, Y = Y, layers = layers))

}
