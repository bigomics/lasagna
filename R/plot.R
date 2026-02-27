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
#' @param layout Layout type: \code{"parallel"} or \code{"hive"}.
#' @return Invisibly returns NULL. Called for side effect (plot).
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
  if (!"rho" %in% edgeattr) message("WARNING: no rho in edge attributes!")
  if (!"weight" %in% edgeattr) message("WARNING: no weight in edge attributes!")
  if (!"value" %in% vattr) stop("ERROR: no value in vertex attributes!")
  if (!"layer" %in% vattr) stop("ERROR: no layer in vertex attributes!")

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

  if (length(igraph::V(graph)) == 0) message("WARNING: graph has no nodes")
  if (length(igraph::E(graph)) == 0) message("WARNING: graph has no edges")
  
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
    ecol <- c("darkorange3", "magenta4")[1 + 1 * (igraph::E(graph)$rho >= 0)]
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
    for (i in 1:length(layers)) {
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
    for (i in 1:length(layers)) {
      tpos <- labpos[i]
      if (layers[i] == "PHENO") next
      graphics::text(xpos[i], -0.02, value.name,
        cex = 1.0 * cex.label, font = 2, pos = tpos,
        adj = 1, offset = 1)
    }
    labposx <- labpos[match(vlayer, layers)]
    graphics::text(x, y,
      cex = 0.85 * cex.label, round(fc[rownames(layout.xy)], 2),
      pos = labposx, adj = 1, offset = 1.0)
    
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
    graphics::text(x, y, labels, cex = cex.label, pos = labposx, adj = 1, offset = 2.8)
  }
  
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
#' @param physics Enable physics simulation.
#' @return A visNetwork widget.
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
  igraph::E(sub)$color <- c("orange", "purple")[1 + 1 * (igraph::E(sub)$weight > 0)]
  
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
#' @param pos Named list of 2-column position matrices per layer.
#' @param draw_edges Logical; draw inter-layer edges.
#' @param min_rho Minimum absolute weight for edges.
#' @param num_edges Maximum number of edges per layer pair.
#' @param znames Named character vector mapping layer codes to display names.
#' @return A plotly object.
#' @export
plot_3d <- function(graph,
                    layout,
                    draw_edges = TRUE,
                    num_edges = 40, 
                    min_rho = 0.1,
                    sign_rho = c("pos","both")[1],
                    cex = 1,
                    cex.gamma = 1,
                    color.by="value",
                    znames = NULL) {

  edges <- NULL
  if (draw_edges) {
    edges <- data.frame(igraph::as_edgelist(graph),
      weight = igraph::E(graph)$weight)
  }
  
  if(is.matrix(layout) || is.data.frame(layout)) {
    layout <- layout[,c("x","y","z")]
    layout <- tapply(1:nrow(layout), layout[,'z'], function(i) layout[i,,drop=FALSE])
  }

  ## feature maps across datatypes
  layout <- mofa.prefix(layout)
  df <- data.frame()
  for (i in names(layout)) {
    ipos <- layout[[i]][,c("x", "y")]
    ipos <- apply(ipos,2,uscale)
    df1 <- data.frame(feature = rownames(ipos), ipos, z = i)
    df <- rbind(df, df1)
  }
  
  vars <- igraph::V(graph)$value
  names(vars) <- igraph::V(graph)$name
  vars <- vars[df$feature]
  vars <- vars / max(abs(vars), na.rm = TRUE)
  df$value <- as.numeric(vars)
  df$size <- abs(as.numeric(vars))^cex.gamma
  df$color <- as.numeric(vars)
  
  if (!is.null(igraph::V(graph)$color) && color.by == "color") {
    df$color <- igraph::V(graph)[df$feature]$color
    df$color <- sub("^[A-Z]+","",df$color) ## remove prefix
  }

  if (!is.null(igraph::V(graph)$size)) {
    vsize <- igraph::V(graph)[df$feature]$size
    df$size <- (vsize / max(vsize, na.rm=TRUE))^cex.gamma
  }

  levels <- intersect(graph$layers, names(layout))
  df$z <- factor(df$z, levels = levels)
  df$text <- paste(df$feature, "<br>value:", round(df$value, digits = 3))

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
    for (i in 1:(length(levels) - 1)) {
      v1 <- rownames(layout[[i]])
      v2 <- rownames(layout[[i + 1]])
      jj <- which((edges[, 1] %in% v1) & (edges[, 2] %in% v2))
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
      "gx" = "Transcriptomics",
      "tx" = "Transcriptomics",
      "mir" = "micro-RNA",
      "px" = "Proteomics",
      "hx" = "Histone",
      "hptm" = "hPTM",
      "dr" = "Drug response",
      "me" = "Methylation",
      "mt" = "Mutation",
      "mu" = "Mutation"
    )
  }

  return(plotlyLasagna(df, znames = znames, edges = edges, cex = cex))

}

#' Internal plotly builder for LASAGNA 3D plot
#' Builds the actual plotly figure from a prepared data frame.
#' @param df Data frame with columns: feature, x, y, z, color, text.
#' @param znames Named character vector for layer display names.
#' @param cex Point size multiplier.
#' @param edges Optional data frame of edges.
#' @return A plotly object.
#' @export
plotlyLasagna <- function(df,
                          znames = NULL,
                          cex = 1,
                          edges = NULL) {

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

  for (k in 1:length(zz)) {
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
        dfe$col <- c("orange", "magenta")[1 + (ee[, 3] > 0)]
        
        fig <- fig %>%
          plotly::add_trace(
            x = dfe$x,
            y = dfe$y,
            z = dfe$z,
            type = "scatter3d",
            mode = "lines",
            line = list(color = dfe$col, width = 0.3, opacity = 0.2),
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

    fig <- fig %>%
      plotly::add_text(
        x = min.x,
        y = max.y,
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
  for (i in 1:length(layers)) {
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
