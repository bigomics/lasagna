## ============================================================
## LASAGNA visualization functions
## ============================================================
#' Plot multi-partite graph using base R graphics
#' Draws a multi-layer LASAGNA graph as a parallel coordinate-style
#' layout using base R \code{igraph::plot}.
#' @param graph An igraph object (output of
#'   \code{solve}).
#' @param layers Character vector of layers to include.
#' @param min.rho Minimum absolute edge weight.
#' @param ntop Number of top features per layer.
#' @param labels Optional named character vector of node labels.
#' @param cex.label Label size multiplier.
#' @param vx.cex Vertex size multiplier.
#' @param xpos Numeric vector of x-positions for layers.
#' @param xlim X-axis limits.
#' @param justgraph Only plot graph (no labels/titles).
#' @param edge.cex Edge width multiplier.
#' @param edge.alpha Edge transparency.
#' @param xdist Distance between layers.
#' @param normalize.edges Normalize edges per connection type.
#' @param yheight Height of the y-axis layout.
#' @param edge.sign Edge sign filter.
#' @param edge.type Edge type filter.
#' @param labpos Label position (2=left, 4=right).
#' @param value.name Display name for node values.
#' @param strip.prefix Remove layer prefix from labels.
#' @param strip.prefix2 Remove prefix with regex.
#' @param prune Remove disconnected vertices.
#' @param edge.gamma Gamma exponent applied to scaled edge widths.
#' @param do.plot Logical; if \code{FALSE}, compute the layout without drawing.
#' @param color.var Vertex colouring: \code{"value"} (by sign) or
#'   \code{"color"} (use the vertex \code{color} attribute).
#' @param layout Layout type: \code{"parallel"} or \code{"hive"}.
#' @return Invisibly, a list with the pruned \code{graph} and its
#'   \code{layout} coordinate matrix. Called mainly for the plot side effect.
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
#' mp <- plot_multipartite(g, min.rho = 0, ntop = 10)
#' @export
plot_multipartite <- function(graph,
                              layers = NULL,
                              min.rho = 0.5,
                              ntop = 25,
                              labels = NULL,
                              cex.label = 1,
                              vx.cex = 1,
                              xpos = NULL,
                              xlim = NULL,
                              justgraph = FALSE,
                              edge.cex = 1,
                              edge.alpha = 0.33,
                              edge.gamma = 2,
                              edge.colors = c("blue3", "magenta4"),
                              xdist = 1,
                              normalize.edges = FALSE,
                              yheight = 0.85,
                              edge.sign = c("both", "pos", "neg", "consensus")[1],
                              edge.type = c("both", "inter", "intra", "both2")[1],
                              labpos = NULL,
                              value.name = NULL,
                              strip.prefix = FALSE,
                              strip.prefix2 = FALSE,
                              prune = FALSE,
                              do.plot = TRUE,
                              color.var = "value",
                              layout = c("parallel", "hive")[1]) {

  vattr <- igraph::vertex_attr_names(graph)
  edgeattr <- igraph::edge_attr_names(graph)
  if (!"rho" %in% edgeattr) warning("no rho in edge attributes")
  if (!"weight" %in% edgeattr) warning("no weight in edge attributes")
  if (!"value" %in% vattr) stop("no value in vertex attributes")
  if (!"layer" %in% vattr) stop("no layer in vertex attributes")

  if (!is.null(labels)) names(labels) <- igraph::V(graph)$name
  
  graph <- prune_graph(
    graph,
    ntop = ntop,
    layers = layers,
    normalize.edges = normalize.edges,
    min.rho = min.rho,
    edge.sign = edge.sign,
    edge.type = edge.type,
    filter = NULL,
    prune = prune
  )

  if (length(igraph::V(graph)) == 0 || length(igraph::E(graph)) == 0) {
    stop("plot_multipartite: no nodes/edges remain after filtering (layers=/ntop=/min.rho=) -- nothing to plot")
  }
  
  layers <- graph$layers
  layers <- setdiff(layers, c("SOURCE", "SINK"))
  fc <- igraph::V(graph)$value
  names(fc) <- igraph::V(graph)$name

  ## layout
  vlayer <- igraph::V(graph)$layer
  if (layout == "parallel") {
    layout.xy <- layout_multipartite(graph, xpos, xdist) 
    if (yheight <= 1) {
      yheight <- yheight * diff(range(layout.xy[,1]))
      layout.xy[, 2] <- yheight * layout.xy[, 2]
    }
  }

  if (layout == "hive") layout.xy <- layout_hiveplot(graph)  

  if(do.plot) {

    vv <- igraph::V(graph)$name
    ewt <- abs(igraph::E(graph)$weight)
    max.wt <- max(ewt, na.rm = TRUE) + 1e-3
    ew <- (ewt / max.wt)^edge.gamma
    
    xpos <- sort(unique(layout.xy[,1]))
    names(xpos) <- layers
    yheight <- diff(range(layout.xy[,2]))
    
    ## vertex size relative to centrality
    vx <- log(1000 * igraph::page_rank(graph, weights = ewt)$vector)
    vx <- (0.1 + abs(vx) / max(abs(vx)))^1

    ## node color by sign (or from color)
    vcol <- NULL
    if (color.var == "value") {
      vcol <- c("blue2", "red2")[1 + 1 * (fc[vv] > 0)]
    }
    if (color.var == "color") {
      vcol <- igraph::V(graph)$color
      vcol <- sub("^[A-Z]+","",vcol) ## strip away WGCNA prefix
    }
    if (is.null(vcol)) vcol <- "grey70"
    vcol <- ifelse(is.na(vcol),"grey70",vcol)
    
    ## edge color by sign
    ecol <- edge.colors[1 + 1 * (igraph::E(graph)$rho >= 0)]
    ii <- which(is.na(igraph::E(graph)$rho))
    if (length(ii)) ecol[ii] <- "grey70"
    ecol <- grDevices::adjustcolor(ecol, edge.alpha)
    
    ## curvature for intra-edges
    ecurv <- c(-0.25, 0)[1 + 1 * grepl("->", igraph::E(graph)$connection_type)]
    
    ##igraph::V(graph)$label <- ""
    if (is.null(xlim)) {
      xlim <- range(layout.xy[, 1])
      xlim <- xlim + c(-1.5, 1.5) * mean(diff(xpos))
    }
    
    graphics::plot(
      graph,
      layout = layout.xy,
      vertex.label = "",
      vertex.label.color = "black",
      vertex.label.dist = 0,
      vertex.label.degree = 0,
      vertex.label.size = 0.001,
      vertex.size = 6 * vx * vx.cex,
      vertex.color = vcol,
      edge.width = 5 * edge.cex * ew,
      edge.color = ecol,
      edge.curved = ecurv,
      xlim = xlim,
      ylim = range(layout.xy[, 2]) + c(-0.03, 0.08) * yheight,
      rescale = FALSE
    )
    
    x <- layout.xy[, 1]
    y <- layout.xy[, 2]
  
    ## titles
    for (i in seq_along(layers)) {
      grp1 <- layers[i]
      y1 <- y[which(vlayer == grp1)]
      graphics::text(xpos[i], max(y1), grp1, font = 2, cex = 1.25, pos = 3, adj = 0, offset = 1.3)
    }
    
    ## plot values
    xt <- min(xpos) + 0.5 * diff(range(xpos))
    if (is.null(labpos)) {
      labpos <- c(2, 4)[1 + 1 * (xpos > xt)]
    } else {
      labpos <- utils::head(rep(labpos, 99), length(layers))
    }
    
    if (is.null(value.name) && !is.null(graph$value.type)) {
      value.name <- graph$value.type
    }
    if (is.null(value.name)) value.name <- "value"
    for (i in seq_along(layers)) {
      tpos <- labpos[i]
      if (layers[i] == "PHENO") next
      graphics::text(xpos[i], -0.02, value.name,
        cex = 1.0 * cex.label, font = 2, pos = tpos,
        adj = 1, offset = 1)
    }
    labposx <- labpos[match(vlayer, layers)]
    graphics::text(x, y,
      cex = 0.85 * cex.label, round(fc[rownames(layout.xy)], 2),
      pos = labposx, adj = 1, offset = 0.8)
    
    ## plot labels
    if (!is.null(labels)) {
      vv <- rownames(layout.xy)
      labels <- labels[vv]
    } else if(!is.null(igraph::V(graph)$label)) {
      labels <- igraph::V(graph)$label
    } else {
      labels <- igraph::V(graph)$name
    }
    
    if (strip.prefix) labels <- mofa.strip_prefix(labels)
    if (strip.prefix2) labels <- sub("^[a-zA-Z]+:", "", labels)
    labels <- gsub("^NA \\(", "(", labels)
    graphics::text(x, y, labels, cex = cex.label, pos = labposx, adj = 1, offset = 3.2)

  } ## end of if do.plot
  
  out <- list(graph = graph, layout = layout.xy)
  invisible(out)

}


#' Plot LASAGNA graph as interactive network with visNetwork
#' @param graph An igraph object.
#' @param layers Character vector of layers to include.
#' @param ntop Number of top nodes by absolute value.
#' @param min_rho Minimum absolute edge weight.
#' @param mst Use minimum spanning tree layout.
#' @param vcex Vertex size multiplier.
#' @param ecex Edge width multiplier.
#' @param egamma Gamma exponent applied to edge widths.
#' @param color.var Vertex colouring: \code{"value"}, \code{"type"}/\code{"layer"},
#'   or \code{"color"} (use the vertex \code{color} attribute).
#' @param labcex Label size multiplier.
#' @param layout Optional layout matrix (rows named by vertex).
#' @param physics Enable physics simulation.
#' @return A visNetwork widget.
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
#' vis <- plot_visgraph(g, min_rho = 0, ntop = 20)
#' @export
plot_visgraph <- function(graph,
                          layers = NULL,
                          ntop = 100,
                          min_rho = 0.3,
                          mst = FALSE,
                          vcex = 1,
                          ecex = 1,
                          egamma = 1,
                          color.var = "value",
                          edge.colors = c("blue", "purple"),
                          labcex = 1,
                          layout = NULL,
                          physics = TRUE) {

  if (is.null(layers)) layers <- graph$layers

  sub <- igraph::subgraph(graph, igraph::V(graph)$layer %in% layers)

  if (mst) {
    ew <- 1 / (1e-8 + abs(igraph::E(sub)$weight)^2)
    sub <- igraph::mst(sub, weights = ew)
  }

  if (ntop > 0) {
    vsel <- utils::head(order(-abs(igraph::V(sub)$value)), ntop)
    sub <- igraph::subgraph(sub, vsel)
  }

  if (min_rho > 0) {
    sub <- igraph::subgraph_from_edges(sub, which(abs(igraph::E(sub)$weight) > min_rho))
  }

  vtype <- sub(":.*", "", igraph::V(sub)$name)
  ntypes <- length(unique(vtype))
  vcol <- "grey"

  if (color.var %in% c("type","layer")) {
    vcol <- grDevices::rainbow(ntypes)[as.factor(vtype)]
  }

  if (color.var == "value") {
    vv <- igraph::V(sub)$value
    vcol <- c("blue", "red")[1 + 1 * (vv > 0)]
  }

  if (color.var == "color") vcol <- igraph::V(sub)$color

  vcol <- sub("^[A-Z]+","",vcol) ## strip away WGCNA prefix
  igraph::V(sub)$color <- vcol
  igraph::V(sub)$color.border <- "black"
 
  vtype <- c("down", "up")[1 + 1 * (igraph::V(sub)$value > 0)]
  igraph::V(sub)$shape <- c("triangleDown", "triangle")[as.factor(vtype)]
  igraph::V(sub)$value <- abs(igraph::V(sub)$value)^2
  igraph::V(sub)$label.cex <- 0.8 * labcex

  ## make special nodes as large as largest size
  sel <- which(igraph::V(sub)$layer %in% c("SOURCE", "SINK", "PHENO"))
  if (length(sel)) {
    cex <- max(abs(igraph::V(sub)$value)) / max(abs(igraph::V(sub)$value[sel]))
    igraph::V(sub)$value[sel] <- igraph::V(sub)$value[sel] * cex
  }

  igraph::E(sub)$width <- 5 * ecex * abs(igraph::E(sub)$weight)^egamma
  igraph::E(sub)$color <- edge.colors[1 + 1 * (igraph::E(sub)$weight > 0)]
  
  data <- visNetwork::toVisNetworkData(sub, idToLabel=FALSE)

  vis <- visNetwork::visNetwork(
    nodes = data$nodes,
    edges = data$edges,
    height = "800px", width = "100%"
  ) %>%
    visNetwork::visNodes(
      scaling = list(min = 5 * vcex, max = 15*vcex),
      shadow = list(enable =TRUE),      
      color = list(border = "black"),
      font = list(align = "left")
    ) %>%
    visNetwork::visEdges(
      color = list(opacity = 0.2),
      scaling = list(min = 3, max = 30)
    ) %>%
    visNetwork::visInteraction(
      hover = TRUE
    ) %>%
    visNetwork::visOptions(
      highlightNearest = TRUE
    ) %>%
    visNetwork::visPhysics(
      enable = physics,
      barnesHut = list(
        springLength = 50
      )
    )

  if (is.null(layout) && !is.null(sub$layout)) layout <- sub$layout
  
  if (!is.null(layout)) {
    vv <- data$nodes$id
    M <- layout[vv,]
    M[,2] <- -M[,2]
    vis <- vis %>% visNetwork::visIgraphLayout(layout = "layout.norm", layoutMatrix = M)
  }
  
  return(vis)

}


#' Plot LASAGNA graph in 3D using plotly
#' Wrapper that creates a 3D plotly visualization from a solved
#' LASAGNA graph and precomputed 2D positions per layer.
#' @param graph An igraph object (output of \code{solve}).
#' @param layout Either a matrix/data frame with columns \code{x}, \code{y},
#'   \code{z} (one row per vertex, row names matching vertex names), or a named
#'   list of 2-column position matrices per layer.
#' @param draw_edges Logical; draw inter-layer edges.
#' @param num_edges Maximum number of edges per layer pair.
#' @param min_rho Minimum absolute weight for edges.
#' @param sign_rho Which edge signs to draw: \code{"pos"}, \code{"neg"} or
#'   \code{"both"}.
#' @param cex Point size multiplier.
#' @param cex.gamma Gamma exponent applied to point sizes.
#' @param color.by Vertex colouring: \code{"value"} or \code{"color"}.
#' @param znames Named character vector mapping layer codes to display names.
#' @return A plotly object.
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
#' ## build a simple per-layer 3D layout (random x/y, z = layer)
#' vn <- igraph::V(g)$name
#' layout <- data.frame(
#'   x = runif(length(vn)), y = runif(length(vn)),
#'   z = sub(":.*", "", vn),
#'   row.names = sub(".*:", "", vn))
#' fig <- plot_3d(g, layout = layout, num_edges = 20, min_rho = 0)
#' @export
plot_3d <- function(graph,
                    layout = NULL,
                    X = NULL,
                    draw_edges = TRUE,
                    num_edges = 40, 
                    min_rho = 0.1,
                    sign_rho = c("both","pos","neg")[1],
                    cex = 1,
                    cex.gamma = 1,
                    edge.cex = 1,
                    color.by = "value",
                    edge_colors = c("blue", "magenta"),
                    znames = NULL, ax=0) {
  edges <- NULL
  if (draw_edges) {
    edges <- data.frame(igraph::as_edgelist(graph),
      weight = igraph::E(graph)$weight)
  }

  vx <- igraph::V(graph)$name
  if(!is.null(layout) && is.data.frame(layout)) {
    message("using user layout")
  } else if(!is.null(layout) && is.character(layout)) {
    if(is.null(X)) stop("must provide X or layout")
    if(!all(vx %in% rownames(X))) stop("incomplete layout")    
    if(!layout %in% c("svd","tsne","umap")) stop("invalid layout")
    layout <- layout_multipartite_3d(graph, X, clust=layout)
  } else if(is.null(layout) && !is.null(graph$layout))  {
    message("using layout in graph object")
    layout <- graph$layout
  } else {
    if(is.null(X)) stop("must provide X or layout")
    if(!all(vx %in% rownames(X))) stop("incomplete layout")
    layout <- layout_multipartite_3d(graph, X, clust="umap")
  }

  ## remove SINK/SOURCE
  if(1) {
    vv <- intersect(c("SINK","SOURCE"), igraph::V(graph)$name)
    if(length(vv)) graph <- igraph::delete_vertices(graph, vv)
  }

  ## feature maps across datatypes
  vx <- igraph::V(graph)$name
  df <- data.frame(layout[vx,])
  dtype <- igraph::V(graph)$layer
  for (dt in unique(dtype)) {
    ii <- which(dtype == dt)
    if(length(ii)>1) {
      df[ii,1:2] <- apply(df[ii,1:2],2,uscale)
    }
  }

  if(is.null(igraph::V(graph)$value)) igraph::V(graph)$value <- 1
  vars <- igraph::V(graph)$value
  names(vars) <- igraph::V(graph)$name
  vars <- vars / max(abs(vars), na.rm = TRUE)
  df$value <- as.numeric(vars)
  df$size <- abs(as.numeric(vars))^cex.gamma
  df$color <- as.numeric(vars)
  
  if (!is.null(igraph::V(graph)$color) && color.by == "color") {
    df$color <- igraph::V(graph)$color
    df$color <- sub("^[A-Z]+","",df$color) ## remove prefix
  }

  if (!is.null(igraph::V(graph)$size)) {
    vsize <- igraph::V(graph)$size
    df$size <- (vsize / max(vsize, na.rm=TRUE))^cex.gamma
  }

  #levels <- intersect(graph$layers, names(layout))
  #df$z <- factor(df$z, levels = levels)
  df$text <- paste(rownames(df), "<br>value:", round(df$value, digits = 3))

  ## filter edges
  if (!is.null(edges) && min_rho >= 0) {
    if (sign_rho == "pos") {
      edges <- edges[ which(edges[, 3] > min_rho), ]
    } else if (sign_rho == "neg") {
      edges <- edges[ which(edges[, 3] < -min_rho), ]
    } else {
      edges <- edges[ which(abs(edges[, 3]) > min_rho), ]
    }
  }

  if (!is.null(edges) && num_edges > 0) {
    e1 <- sub(":.*","",edges[,1])
    e2 <- sub(":.*","",edges[,2])
    ee <- paste0(e1,'-',e2)
    for (this.ee in unique(ee)) {
      jj <- which(ee == this.ee)
      sel <- utils::head(jj[order(-abs(edges[jj, 3]))], num_edges)
      jj <- setdiff(jj, sel)
      edges[jj, 3] <- 0
    }
    edges <- edges[edges[, 3] != 0, ]
  }

  ## layer name mapping
  if (is.null(znames)) {
    znames <- c(
      "PHENO" = "Phenotype",
      "ph" = "Phenotype",
      "gset" = "Pathway",
      "mx" = "Metabolomics",
      "lx" = "Lipidomics",      
      "gx" = "Transcriptomics",
      "tx" = "Transcriptomics",
      "mir" = "microRNA",
      "px" = "Proteomics",
      "hx" = "Histone",
      "hptm" = "hPTM",
      "dr" = "Drug response",
      "me" = "Methylation",
      "mt" = "Mutation",
      "mut" = "Mutation",      
      "mu" = "Mutation"
    )
  }

  edge.colors <- c()
  
  plt <- plotlyLasagna(df, znames = znames, edges = edges,
    edge_colors = edge_colors,
    cex = cex, edge.cex = edge.cex, ax = ax)
  
  return(plt)
}

#' Internal plotly builder for LASAGNA 3D plot
#' Builds the actual plotly figure from a prepared data frame.
#' @param df Data frame with columns: feature, x, y, z, color, text.
#' @param znames Named character vector for layer display names.
#' @param cex Point size multiplier.
#' @param edges Optional data frame of edges.
#' @return A plotly object.
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   feature = paste0("f", 1:10),
#'   x = runif(10), y = runif(10),
#'   z = rep(c("gx", "px"), each = 5),
#'   color = rnorm(10), size = runif(10),
#'   text = paste0("f", 1:10))
#' fig <- plotlyLasagna(df)
#' @export
plotlyLasagna <- function(df,
                          znames = NULL, ax = 1,
                          cex = 1, edge.cex = 1,
                          edges = NULL,
                          edge_colors = c("blue", "magenta")
                          ) {
  zz <- sort(unique(df$z))
  min.x <- min(df$x, na.rm = TRUE)
  max.x <- max(df$x, na.rm = TRUE)
  min.y <- min(df$y, na.rm = TRUE)
  max.y <- max(df$y, na.rm = TRUE)

  edgetype1 <- edgetype2 <- NULL

  if (!is.null(edges)) {
    edgetype1 <- mofa.get_prefix(edges[, 1])
    edgetype2 <- mofa.get_prefix(edges[, 2])
  }

  if (is.null(df$text)) df$text <- rownames(df)

  fig <- plotly::plot_ly()

  for (k in seq_along(zz)) {
    z <- zz[k]
    df1 <- df[which(df$z == z), c("x", "y", "z", "color", "size", "text")]

    mx <- c(min.x, min.x, max.x, max.x)
    my <- c(min.y, max.y, min.y, max.y)
    mz <- rep(z, 4)

    fig <- fig %>%
      plotly::add_mesh(
        x = mx,
        y = my,
        z = mz,
        hoverinfo = "none",
        color = "grey",
        opacity = 0.15
      )

    ## add layer points
    fig <- fig %>%
      plotly::add_markers(
        data = df1,
        x = ~x,
        y = ~y,
        z = z,
        text = ~text,
        hoverinfo = "text",
        hovertemplate = "%{text}",
        type = "scattergl",
        marker = list(
          size = ~size * 15 * cex + 3,
          color = ~color,
          line = list(color = "#88888844", width = 0.0),
          colorscale = "Bluered",
          showscale = FALSE,
          showlegend = FALSE
        ),
        showlegend = FALSE
      )

    ## add segments
    if (k < length(zz) && !is.null(edges)) {
      sel1 <- which(edgetype1 == zz[k] & edgetype2 == zz[k + 1])
      sel2 <- which(edgetype2 == zz[k] & edgetype1 == zz[k + 1])
      df2 <- df[which(df$z == zz[k + 1]), c("x", "y", "z")]
      sel <- unique(c(sel1, sel2))

      if (length(sel)) {
        ee <- edges[sel, ]
        idx <- as.vector(t(as.matrix(ee[, 1:2])))
        dfe <- rbind(df1[, c("x", "y", "z")], df2[, c("x", "y", "z")])[idx, ]
        dfe$pair_id <- as.vector(mapply(rep, 1:nrow(ee), 2))
        cc <- edge_colors[1 + (ee[, 3] > 0)]
        dfe$col <- as.vector(mapply(rep, cc, 2))
        fig <- fig %>%
          plotly::add_trace(
            x = dfe$x,
            y = dfe$y,
            z = dfe$z,
            type = "scatter3d",
            mode = "lines",
            line = list(color = dfe$col, width = edge.cex, opacity = 0.2),
            split = dfe$pair_id,
            showlegend = FALSE,
            inherit = FALSE
          )
      }
    }

    ## add layer title
    ztext <- z
    if (!is.null(znames) && (is.integer(z) || z %in% names(znames))) {
      if (is.factor(z)) z <- as.character(z)
      ztext <- znames[z]
    }

    ax <- ax %% 4
    if(ax==0) {zx=min.x; zy=max.y}
    if(ax==1) {zx=min.x; zy=min.y}
    if(ax==2) {zx=max.x; zy=max.y}
    if(ax==3) {zx=max.x; zy=min.y}
    
    fig <- fig %>%
      plotly::add_text(
        x = zx,
        y = zy,
        z = z,
        mode = "text",
        text = ztext,
        textfont = list(size = 24),
        showlegend = FALSE,
        inherit = FALSE
      )
  }

  fig <- fig %>%
    plotly::layout(
      scene = list(
        xaxis = list(title = ""),
        yaxis = list(title = ""),
        zaxis = list(title = "", showticklabels = FALSE, tickfont = list(size = 0)),
        aspectmode = "cube",
        showlegend = FALSE
      )
    )

  fig
}

layout_multipartite <- function(graph, xpos=NULL, xdist=1) {

  layers <- graph$layers
  layers <- setdiff(layers, c("SOURCE","SINK"))
  if (is.null(xpos)) xpos <- c(0:(length(layers) - 1))
  xpos <- xpos * xdist
  xpos <- utils::head(rep(xpos, 10), length(layers))
  vlayer <- igraph::V(graph)$layer
  x <- xpos[match(vlayer, layers)]
  y <- igraph::V(graph)$value
  layout.xy <- cbind(x = x, y = y)
  rownames(layout.xy) <- igraph::V(graph)$name

  for (i in unique(layout.xy)[, 1]) {
    ii <- which(layout.xy[, 1] == i)
    layout.xy[ii, 2] <- rank(layout.xy[ii, 2], ties.method = "random") / length(ii)
    if (length(ii) == 1) layout.xy[ii, 2] <- 0.5 * layout.xy[ii, 2]
  }

  return(layout.xy)

}

layout_hiveplot <- function(graph) {

  layers <- graph$layers
  layers <- setdiff(layers, c("SOURCE","SINK"))
  nlayers <- length(layers)
  vlayer <- igraph::V(graph)$layer
  layout.xy <- matrix(NA, length(vlayer), 2)
  rownames(layout.xy) <- igraph::V(graph)$name
  fc <- igraph::V(graph)$value  
  names(fc) <- igraph::V(graph)$name
  for (i in seq_along(layers)) {
    ii <- which(vlayer == layers[i])
    vv <- igraph::V(graph)[ii]
    phi <- pi / 2 + i * 2 * pi / nlayers
    vf <- fc[vv$name]
    r <- 0.2 + rank(vf, ties.method = "random") / length(vf)
    x <- r * cos(phi)
    y <- r * sin(phi)
    layout.xy[ii, ] <- cbind(x, y)
  }
  rownames(layout.xy) <- igraph::V(graph)$name

  return(layout.xy)

}

#' @export
layout_multipartite_3d <- function(graph, X, clust=c("svd","tsne","umap")) {
  clust <- match.arg(clust)
  layers <- graph$layers
  layers <- setdiff(layers, c("SOURCE","SINK"))

  sx <- t(scale(t(X)))
  ff <- mofa.split_data(sx)[layers]
  xy <- list()
  k=names(ff)[1]
  for(k in names(ff)) {
    nn <- nrow(ff[[k]])
    if(nn > 5 && clust == 'tsne') {
      if (!requireNamespace("Rtsne", quietly = TRUE)) stop("package Rtsne required for clust='tsne'")
      px <- max(min(30, (nn-1)/3),2)
      xy[[k]] <- Rtsne::Rtsne(ff[[k]], perplexity=px)$Y
    } else if(nn > 5 && clust == 'umap') {
      if (!requireNamespace("uwot", quietly = TRUE)) stop("package uwot required for clust='umap'")
      nb <- max(min(15, nn/3),2)
      xy[[k]] <- uwot::umap(ff[[k]], n_neighbors=nb)
    } else if(nn == 1) {
      xy[[k]] <- matrix(0, nrow=1, ncol=2)
    } else {
      xy[[k]] <- svd(ff[[k]])$u[,1:2]
    }
  }
  for(i in 1:length(xy)) rownames(xy[[i]]) <- rownames(ff[[i]])
  xy <- mofa.merge_data(xy)
  colnames(xy) <- c("x","y")
  ii <- match(igraph::V(graph)$name, rownames(xy))
  xy <- xy[ii,]

  z <- factor(igraph::V(graph)$layer, levels=graph$layers)
  layout <- data.frame(xy, z)
  rownames(layout) <- igraph::V(graph)$name

  return(layout)
}
