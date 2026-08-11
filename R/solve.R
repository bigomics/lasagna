## ============================================================
## LASAGNA graph solving functions
## ============================================================

#' Solve LASAGNA graph for a single phenotype
#' Given a LASAGNA model and a phenotype (column of Y), computes
#' per-node fold change and correlation, then weights edges accordingly.
#' The resulting graph is pruned to \code{max_edges} per connection type.
#' @param obj A LASAGNA model (output of \code{create_model}).
#' @param pheno Column name in \code{obj$Y} to solve for.
#' @param max_edges Maximum edges per connection type.
#' @param value.type Node value type: \code{"rho"} or \code{"fc"}.
#' @param min_rho Minimum absolute weight threshold.
#' @param prune Remove disconnected vertices.
#' @param fc.weights Weight edges by node fold change.
#' @param sp.weight Use shortest-path weighting.
#' @param graph Optional pre-existing graph to use instead of \code{obj$graph}.
#' @return An igraph object with solved weights and node values.
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
#' g <- solve(model, pheno = colnames(model$Y)[1], max_edges = 50, prune = FALSE)
#' igraph::vcount(g)
#' @export
solve <- function(obj,
                  pheno,
                  max_edges = 100,
                  value.type = c("rho","fc")[1],
                  min_rho = 0,
                  prune = TRUE,
                  fc.weights = TRUE,
                  sp.weight = FALSE,
                  graph = NULL) {

  if (!pheno %in% colnames(obj$Y)) stop("pheno not in Y")

  if (!"rho" %in% names(igraph::edge_attr(obj$graph))) {
    stop("graph edges should have rho attribute")
  }

  if (is.null(graph)) graph <- obj$graph

  ## score the nodes against the phenotype
  dat <- mask_neutral_samples(obj$X, obj$Y[, pheno])
  rho <- node_correlation(dat$X, dat$y)
  fc <- node_foldchange(dat$X, dat$y, rho, pheno)

  ## weight the edges from the node scores
  graph <- set_node_values(graph, rho, fc, value.type)
  graph <- set_edge_weights(graph, fc.weights)
  if (sp.weight) graph <- add_sp_weight(graph, obj$layers)

  ## zero the weight of every unwanted edge, then drop them in one go
  if (min_rho > 0) graph <- zero_weak_edges(graph, min_rho)
  if (max_edges > 0) graph <- zero_excess_edges(graph, max_edges)
  graph <- igraph::delete_edges(graph, which(igraph::E(graph)$weight == 0))

  ## prune vertices
  if (prune) {
    ewt <- igraph::E(graph)$weight
    graph <- igraph::subgraph_from_edges(graph, which(abs(ewt) > 0))
  }

  return(graph)

}

## A phenotype coded -1/0/1 uses the zeros for the samples that belong to
## neither side of the contrast. Mask those out of both the phenotype and its
## PHENO rows in the data, so that they weigh on neither correlation nor fold
## change.
mask_neutral_samples <- function(X, y) {

  ii <- grep("PHENO", rownames(X))
  has.min1 <- (min(X[ii, ], na.rm = TRUE) < 0)
  if (length(ii) && has.min1) {
    X[ii, ][which(X[ii, ] == 0)] <- NA
    y[y == 0] <- NA
  }

  return(list(X = X, y = y))

}

## Correlation of each node with the phenotype. Nodes that are constant or
## never observed together with the phenotype score zero rather than NA.
node_correlation <- function(X, y) {

  rho <- stats::cor(t(X), y, use = "pairwise")[, 1]
  rho[is.na(rho)] <- 0

  return(rho)

}

## Fold change of each node between the two sides of the phenotype. For the
## PHENO nodes a fold change is meaningless, so they take their correlation
## instead.
node_foldchange <- function(X, y, rho, pheno) {

  i0 <- which(y <= 0)
  i1 <- which(y > 0)
  if (length(i0) == 0 || length(i1) == 0) {
    stop("solve: phenotype '", pheno, "' has no samples on one side (degenerate/constant trait) -- cannot compute fold change")
  }
  m1 <- rowMeans(X[, i1, drop = FALSE], na.rm = TRUE)
  m0 <- rowMeans(X[, i0, drop = FALSE], na.rm = TRUE)
  fc <- m1 - m0

  ii <- grep("PHENO", names(fc))
  if (length(ii)) fc[ii] <- rho[ii]

  return(fc)

}

## Attach both node scores to the graph and pick the one that acts as the
## node value, which is what the edge weighting and the plots read.
set_node_values <- function(graph, rho, fc, value.type) {

  igraph::V(graph)$rho <- rho
  igraph::V(graph)$fc <- fc
  igraph::V(graph)$value <- if (value.type == "rho") rho else fc
  graph$value.type <- value.type

  return(graph)

}

## Weight each edge by its correlation, optionally scaled by the geometric
## mean of the node values at both ends, so that edges between two nodes that
## respond to the phenotype outweigh equally correlated but flat ones. The
## SOURCE and SINK edges are bookkeeping and always weigh 1.
set_edge_weights <- function(graph, fc.weights) {

  ww <- 1
  weight.type <- "rho"
  if (fc.weights) {
    ee <- igraph::as_edgelist(graph)
    ff <- igraph::V(graph)$value
    names(ff) <- igraph::V(graph)$name
    ww <- abs(ff[ee[, 1]] * ff[ee[, 2]])^0.5
    weight.type <- paste0(weight.type, "*vv")
  }

  ee.rho <- igraph::E(graph)$rho
  ee.rho[is.na(ee.rho)] <- 0.1234
  igraph::E(graph)$weight <- ee.rho * ww

  if (any(grepl("SINK|SOURCE", igraph::V(graph)$name))) {
    igraph::E(graph)[.to("SINK")]$weight <- 1
    igraph::E(graph)[.from("SOURCE")]$weight <- 1
  }
  graph$weight.type <- weight.type

  return(graph)

}

## Scale the edge weights by how well each edge carries a strong path from
## SOURCE to SINK, which favours edges of complete cross-layer chains over
## strong but isolated ones. Needs the helper nodes that add.sink adds.
add_sp_weight <- function(graph, layers) {

  if (!all(c("SOURCE", "SINK") %in% igraph::V(graph)$name)) {
    stop("sp.weight=TRUE requires a model built with add.sink=TRUE")
  }

  sp.wt <- sp_edge_weight(graph, layers)
  sp.wt <- (sp.wt / max(sp.wt, na.rm = TRUE))^2
  igraph::E(graph)$weight <- igraph::E(graph)$weight * sp.wt
  graph$weight.type <- paste0(graph$weight.type, "*sp")

  return(graph)

}

## Zero all but the 'max_edges' strongest edges of each connection type, so
## that no single type can crowd out the others.
zero_excess_edges <- function(graph, max_edges) {

  ewt <- igraph::E(graph)$weight
  esel <- tapply(
    seq_along(igraph::E(graph)), igraph::E(graph)$connection_type,
    function(ii) utils::head(ii[order(-abs(ewt[ii]))], max_edges)
  )
  dsel <- setdiff(seq_along(igraph::E(graph)), unlist(esel))
  igraph::E(graph)$weight[dsel] <- 0

  return(graph)

}

## Internal: compute shortest-path-based edge weights
sp_edge_weight <- function(graph, layers) {

  layers <- unique(c("SOURCE", layers, "SINK"))
  wt <- abs(igraph::E(graph)$weight)
  wt[is.na(wt)] <- 0
  ee <- igraph::as_edgelist(graph)
  v1 <- ee[, 1]
  v2 <- ee[, 2]
  l1 <- match(igraph::V(graph)[v1]$layer, layers)
  l2 <- match(igraph::V(graph)[v2]$layer, layers)
  p1 <- ifelse(l1 < l2, v1, v2)
  p2 <- ifelse(l1 < l2, v2, v1)
  wt <- wt + 1e-8

  s1 <- igraph::shortest_paths(graph,
    from = "SOURCE", to = p1, weights = 1 / wt, output = "epath"
  )

  s2 <- igraph::shortest_paths(graph,
    from = "SINK", to = p2, weights = 1 / wt, output = "epath"
  )

  sp <- mapply(c, s1$epath, s2$epath)
  sp.score <- vapply(sp, function(e) min(wt[e], na.rm = TRUE), numeric(1))

  return(sp.score)

}

#' Solve LASAGNA graph for multiple phenotypes
#' Solves the graph iteratively for multiple contrasts and returns
#' the root-mean-square (RMS) consensus graph.
#' @param obj A LASAGNA model.
#' @param min_rho Minimum absolute weight threshold for final graph.
#' @param traits Character vector of trait names (columns of Y). If
#'   NULL, all traits are used.
#' @param max_edges Maximum edges per connection type per trait.
#' @param value.type Node value type.
#' @param fc.weights Use fold-change weighting.
#' @param prune Remove disconnected vertices.
#' @param sp.weight Use shortest-path weighting.
#' @return An igraph object representing the RMS consensus graph.
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
#' g <- multisolve(model, min_rho = 0, prune = FALSE)
#' igraph::ecount(g)
#' @export
multisolve <- function(obj,
                       min_rho = 0.2,
                       traits = NULL,
                       max_edges = 100,
                       value.type = "rho",
                       fc.weights = TRUE,
                       prune = TRUE,
                       sp.weight = FALSE) {

  if (is.null(traits)) traits <- colnames(obj$Y)
  traits <- intersect(traits, colnames(obj$Y))

  M <- list()
  V <- list()
  for (ct in traits) {
    solved <- solve(obj,
      pheno = ct,
      min_rho = min_rho,
      max_edges = max_edges,
      value.type = value.type,
      fc.weights = fc.weights,
      sp.weight = sp.weight,
      prune = FALSE
    )
    adj <- igraph::as_adjacency_matrix(solved, attr = "weight")
    M[[ct]] <- adj
    V[[ct]] <- igraph::V(solved)$value
    names(V[[ct]]) <- igraph::V(solved)$name
  }

  ## RMS adjacency matrix as consensus solution
  M <- lapply(M, function(mat) as.matrix(mat^2))
  avgM <- sqrt(Reduce("+", M) / length(M))
  avgM <- avgM * (avgM > min_rho)
  rmsV <- rowMeans(do.call(cbind, V)**2)
  
  gr <- obj$graph
  igraph::V(gr)$value <- rmsV
  
  ee <- igraph::get.edges(gr, igraph::E(gr))
  igraph::E(gr)$weight <- sign(igraph::E(gr)$weight) * avgM[ee]
  del.ee <- which(abs(igraph::E(gr)$weight) < min_rho)
  gr <- igraph::delete_edges(gr, del.ee)
  
  if (prune) {
    del.vv <- which(igraph::degree(gr) == 0)
    gr <- igraph::delete_vertices(gr, del.vv)
  }

  return(gr)

}
