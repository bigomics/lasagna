## ============================================================
## LASAGNA graph pruning functions
## ============================================================

#' Prune a solved LASAGNA graph
#'
#' Filters vertices and edges by layer, node value, edge weight,
#' edge sign, and edge type (inter vs intra). Useful for reducing
#' graph complexity before visualization.
#'
#' @param graph An igraph object (output of \code{solve}
#'   or \code{multisolve}).
#' @param ntop Number of top features per layer (by absolute
#'   value).
#' @param layers Character vector of layers to keep. If NULL, uses
#'   all.
#' @param normalize.edges Normalize edge weights per connection
#'   type.
#' @param min.rho Minimum absolute edge weight.
#' @param edge.sign Which edge signs to keep: \code{"both"},
#'   \code{"pos"}, \code{"neg"}, or \code{"consensus"}.
#' @param edge.type Which edge types: \code{"both"},
#'   \code{"inter"}, \code{"intra"}, or \code{"both2"}.
#' @param filter Named list for regex filtering nodes per layer.
#' @param select Character vector of node names to select.
#' @param prune Remove disconnected vertices.
#'
#' @return A pruned igraph object.
#' @export
prune_graph <- function(graph, ntop = 100, layers = NULL,
                        normalize.edges = FALSE, min.rho = 0.3,
                        edge.sign = c("both", "pos", "neg", "consensus")[1],
                        edge.type = c("both", "inter", "intra", "both2")[1],
                        filter = NULL, select = NULL, prune = TRUE) {
  if (is.null(layers)) layers <- graph$layers
  if (is.null(layers)) layers <- unique(igraph::V(graph)$layer)
  layers <- setdiff(layers, c("SOURCE", "SINK"))
  graph <- igraph::subgraph(graph, igraph::V(graph)$layer %in% layers)

  if (!"value" %in% names(igraph::vertex_attr(graph))) {
    stop("vertex must have 'value' attribute")
  }

  ## select nodes/modules
  if (!is.null(select)) {
    v1 <- (igraph::V(graph)$name %in% select)
    v2 <- (sub(".*:", "", igraph::V(graph)$name) %in% select)
    v3 <- (sub(":.*", "", igraph::V(graph)$name) %in% select)
    graph <- igraph::subgraph(graph, vids = which(v1 | v2 | v3))
  }

  if (!is.null(filter)) {
    if (!is.list(filter)) stop("filter must be a named list")
    if (is.null(names(filter))) stop("filter must be a named list")
    for (k in names(filter)) {
      vv <- igraph::V(graph)$name
      filt <- filter[[k]]
      vids <- igraph::V(graph)$layer != k | grepl(filt, vv, ignore.case = TRUE)
      graph <- igraph::subgraph(graph, which(vids))
    }
  }

  ## select ntop features
  fc <- igraph::V(graph)$value
  names(fc) <- igraph::V(graph)$name
  if (!is.null(ntop) && ntop > 0) {
    ii <- tapply(
      1:length(fc), igraph::V(graph)$layer,
      function(i) utils::head(i[order(-abs(fc[i]))], ntop)
    )
    ii <- unlist(ii[names(ii) %in% layers])
    fc <- fc[ii]
    graph <- igraph::subgraph(graph, igraph::V(graph)[ii])
  }

  if (normalize.edges) {
    for (e in unique(igraph::E(graph)$connection_type)) {
      ii <- which(igraph::E(graph)$connection_type == e)
      max.wt <- max(abs(igraph::E(graph)$weight[ii]), na.rm = TRUE) + 1e-3
      igraph::E(graph)$weight[ii] <- igraph::E(graph)$weight[ii] / max.wt
    }
  }

  if (min.rho > 0) {
    ii <- which(abs(igraph::E(graph)$weight) < min.rho)
    if (length(ii)) igraph::E(graph)$weight[ii] <- 0
  }

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

  ## delete intra or inter edges
  ic <- grepl("->", igraph::E(graph)$connection_type)
  if (edge.type == "inter") {
    igraph::E(graph)$weight[!ic] <- 0
  } else if (edge.type == "intra") {
    igraph::E(graph)$weight[ic] <- 0
  } else if (edge.type == "both2") {
    sel <- (!ic & igraph::E(graph)$weight < 0)
    igraph::E(graph)$weight[sel] <- 0
  }
  graph <- igraph::delete_edges(graph, which(igraph::E(graph)$weight == 0))

  if (prune) {
    graph <- igraph::subgraph_from_edges(graph, igraph::E(graph))
  }
  return(graph)
}
