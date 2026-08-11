## ============================================================
## LASAGNA graph pruning functions
## ============================================================
#' Prune a solved LASAGNA graph
#' Filters vertices and edges by layer, node value, edge weight,
#' edge sign, and edge type (inter vs intra). Useful for reducing
#' graph complexity before visualization.
#' @param graph An igraph object (output of \code{solve} or \code{multisolve}).
#' @param ntop Number of top features per layer (by absolute value).
#' @param layers Character vector of layers to keep. If NULL, uses all.
#' @param normalize.edges Normalize edge weights per connection type.
#' @param min.rho Minimum absolute edge weight.
#' @param edge.sign Which edge signs to keep: \code{"both"},
#'   \code{"pos"}, \code{"neg"}, or \code{"consensus"}.
#' @param edge.type Which edge types: \code{"both"},
#'   \code{"inter"}, \code{"intra"}, or \code{"both2"}.
#' @param filter Named list for regex filtering nodes per layer.
#' @param select Character vector of node names to select.
#' @param prune Remove disconnected vertices.
#' @return A pruned igraph object.
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
#' g <- solve(model, pheno = colnames(model$Y)[1], prune = FALSE)
#' pruned <- prune_graph(g, ntop = 5, min.rho = 0)
#' igraph::vcount(pruned)
#' @export
prune_graph <- function(graph,
                        ntop = 100,
                        layers = NULL,
                        normalize.edges = FALSE,
                        min.rho = 0.3,
                        edge.sign = c("both", "pos", "neg", "consensus")[1],
                        edge.type = c("both", "inter", "intra", "both2")[1],
                        filter = NULL,
                        select = NULL,
                        prune = TRUE) {

  layers <- resolve_layers(graph, layers)
  graph <- igraph::subgraph(graph, igraph::V(graph)$layer %in% layers)

  if (!"value" %in% names(igraph::vertex_attr(graph))) {
    stop("vertex must have 'value' attribute")
  }

  if (!is.null(filter)) {
    if (!is.list(filter)) stop("filter must be a named list")
    if (is.null(names(filter))) stop("filter must be a named list")
  }

  ## reduce the vertex set
  if (!is.null(select)) graph <- keep_named_nodes(graph, select)
  if (!is.null(filter)) graph <- keep_matching_nodes(graph, filter)
  if (!is.null(ntop) && ntop > 0) graph <- keep_top_nodes(graph, ntop, layers)

  ## zero the weight of every unwanted edge, then drop them in one go
  if (normalize.edges) graph <- normalize_edge_weights(graph)
  if (min.rho > 0) graph <- zero_weak_edges(graph, min.rho)
  graph <- zero_edges_by_sign(graph, edge.sign, layers)
  graph <- zero_edges_by_type(graph, edge.type)

  return(drop_zero_edges(graph, prune))

}

## Layers the graph is to be reduced to: the caller's choice, otherwise the
## layers attached to the graph object, otherwise those found on its
## vertices. SOURCE and SINK are helper nodes and are never kept.
resolve_layers <- function(graph, layers) {

  if (is.null(layers)) layers <- graph$layers
  if (is.null(layers)) layers <- unique(igraph::V(graph)$layer)

  return(setdiff(layers, c("SOURCE", "SINK")))

}

## Keep the nodes named in 'select', matching either the full node name or
## the part before or after its layer prefix, so that a module can be
## selected by bare name as well as by qualified name.
keep_named_nodes <- function(graph, select) {

  vv <- igraph::V(graph)$name
  v1 <- (vv %in% select)
  v2 <- (sub(".*:", "", vv) %in% select)
  v3 <- (sub(":.*", "", vv) %in% select)

  return(igraph::subgraph(graph, vids = which(v1 | v2 | v3)))

}

## Keep the nodes matching the per-layer regular expressions in 'filter', a
## list of patterns named by layer. Each pattern only constrains its own
## layer; nodes of the layers not named in the list are always kept.
keep_matching_nodes <- function(graph, filter) {

  for (k in names(filter)) {
    vv <- igraph::V(graph)$name
    filt <- filter[[k]]
    vids <- igraph::V(graph)$layer != k | grepl(filt, vv, ignore.case = TRUE)
    graph <- igraph::subgraph(graph, which(vids))
  }

  return(graph)

}

## Keep the 'ntop' nodes of largest absolute value within each layer. Layers
## outside 'layers' are dropped entirely, which matters for the SOURCE/SINK
## nodes that survive an explicit layer selection.
keep_top_nodes <- function(graph, ntop, layers) {

  fc <- igraph::V(graph)$value
  ii <- tapply(
    seq_along(fc), igraph::V(graph)$layer,
    function(i) utils::head(i[order(-abs(fc[i]))], ntop)
  )
  ii <- unlist(ii[names(ii) %in% layers])

  return(igraph::subgraph(graph, igraph::V(graph)[ii]))

}

## Rescale the edge weights of each connection type onto a common scale, so
## that weights of different connection types can be compared.
normalize_edge_weights <- function(graph) {

  for (e in unique(igraph::E(graph)$connection_type)) {
    ii <- which(igraph::E(graph)$connection_type == e)
    max.wt <- max(abs(igraph::E(graph)$weight[ii]), na.rm = TRUE) + 1e-3
    igraph::E(graph)$weight[ii] <- igraph::E(graph)$weight[ii] / max.wt
  }

  return(graph)

}

## Zero the edges weaker than 'min.rho'.
zero_weak_edges <- function(graph, min.rho) {

  ii <- which(abs(igraph::E(graph)$weight) < min.rho)
  if (length(ii)) igraph::E(graph)$weight[ii] <- 0

  return(graph)

}

## Zero the edges of unwanted sign. "consensus" keeps the edges whose sign
## agrees with the expected direction of their source layer, which is
## inverted for the miRNA layers as those anticorrelate with their targets.
zero_edges_by_sign <- function(graph, edge.sign, layers) {

  ewt <- igraph::E(graph)$weight

  if (grepl("pos", edge.sign)) {
    igraph::E(graph)$weight[ewt < 0] <- 0
  } else if (grepl("neg", edge.sign)) {
    igraph::E(graph)$weight[ewt > 0] <- 0
  } else if (edge.sign == "consensus") {
    layersign <- rep(1, length(layers))
    names(layersign) <- layers
    layersign[grep("^mi|^mir", layers)] <- -1
    v1 <- igraph::as_edgelist(graph)[, 1]
    esign <- layersign[igraph::V(graph)[v1]$layer]
    igraph::E(graph)$weight <- ewt * (sign(ewt) == esign)
  }

  return(graph)

}

## Zero the edges of unwanted type, inter-layer edges being those whose
## connection type names two layers. "both2" keeps both types but drops the
## negative intra-layer edges.
zero_edges_by_type <- function(graph, edge.type) {

  ic <- grepl("->", igraph::E(graph)$connection_type)

  if (edge.type == "inter") {
    igraph::E(graph)$weight[!ic] <- 0
  } else if (edge.type == "intra") {
    igraph::E(graph)$weight[ic] <- 0
  } else if (edge.type == "both2") {
    sel <- (!ic & igraph::E(graph)$weight < 0)
    igraph::E(graph)$weight[sel] <- 0
  }

  return(graph)

}

## Drop the edges that the filters zeroed out and, with prune=TRUE, the
## vertices that are left without any edge.
drop_zero_edges <- function(graph, prune) {

  graph <- igraph::delete_edges(graph, which(igraph::E(graph)$weight == 0))
  if (prune) graph <- igraph::subgraph_from_edges(graph, igraph::E(graph))

  return(graph)

}
