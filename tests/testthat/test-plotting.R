## layout_multipartite_3d ---------------------------------------------------

test_that("layout_multipartite_3d default clust argument does not crash", {
  skip_on_cran()
  out <- layout_multipartite_3d(toy_model$graph, toy_model$X)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("x", "y", "z") %in% colnames(out)))
  expect_equal(nrow(out), igraph::vcount(toy_model$graph))
})

test_that("layout_multipartite_3d runs with an explicit clust value", {
  skip_on_cran()
  out <- layout_multipartite_3d(toy_model$graph, toy_model$X, clust = "svd")
  expect_s3_class(out, "data.frame")
  expect_equal(rownames(out), igraph::V(toy_model$graph)$name)
})

test_that("layout_multipartite_3d handles a single-feature layer", {
  skip_on_cran()
  out <- layout_multipartite_3d(
    single_feature_model$graph, single_feature_model$X, clust = "svd"
  )
  expect_s3_class(out, "data.frame")
  px_row <- grep("^px:", rownames(out))
  expect_length(px_row, 1)
  expect_equal(unname(unlist(out[px_row, c("x", "y")])), c(0, 0))
})

## plot_multipartite ---------------------------------------------------------

test_that("plot_multipartite runs and returns the pruned graph and layout", {
  skip_on_cran()
  out <- plot_multipartite(toy_graph, min.rho = 0, ntop = 10, do.plot = FALSE)
  expect_true(is.list(out))
  expect_true(igraph::is_igraph(out$graph))
  expect_true(is.matrix(out$layout))
})

test_that("plot_multipartite errors clearly when the layers filter matches nothing", {
  skip_on_cran()
  expect_error(
    plot_multipartite(toy_graph, layers = "nonexistent", do.plot = TRUE),
    "no nodes/edges remain"
  )
  expect_error(
    plot_multipartite(toy_graph, layers = "nonexistent", do.plot = FALSE),
    "no nodes/edges remain"
  )
})

## plot_3d --------------------------------------------------------------------

test_that("plot_3d runs with clust = 'svd'", {
  skip_on_cran()
  fig <- plot_3d(toy_graph, layout = "svd", X = toy_model$X, num_edges = 5, min_rho = 0)
  expect_s3_class(fig, "plotly")
})

test_that("plot_3d uses a layout pre-attached to the graph", {
  skip_on_cran()
  g <- toy_graph
  g$layout <- layout_multipartite_3d(g, toy_model$X, clust = "svd")
  fig <- plot_3d(g, num_edges = 5, min_rho = 0)
  expect_s3_class(fig, "plotly")
})

test_that("plot_3d errors clearly when neither X nor layout is given", {
  skip_on_cran()
  expect_error(plot_3d(toy_graph), "must provide X or layout")
})

test_that("plot_3d errors clearly on an incomplete layout", {
  skip_on_cran()
  incomplete_X <- toy_model$X[-1, , drop = FALSE]
  expect_error(
    plot_3d(toy_graph, layout = "svd", X = incomplete_X),
    "incomplete layout"
  )
})

test_that("plot_3d errors clearly on an invalid layout string", {
  skip_on_cran()
  expect_error(
    plot_3d(toy_graph, layout = "badmethod", X = toy_model$X),
    "invalid layout"
  )
})

## plotlyLasagna ---------------------------------------------------------------

test_that("plotlyLasagna renders a manually constructed data frame", {
  skip_on_cran()
  set.seed(1)
  df <- data.frame(
    feature = paste0("f", 1:10),
    x = runif(10), y = runif(10),
    z = rep(c("gx", "px"), each = 5),
    color = rnorm(10), size = runif(10),
    text = paste0("f", 1:10)
  )
  fig <- plotlyLasagna(df)
  expect_s3_class(fig, "plotly")
})

## plot_visgraph -----------------------------------------------------------

test_that("plot_visgraph returns a visNetwork widget", {
  skip_on_cran()
  vis <- plot_visgraph(toy_graph, min_rho = 0, ntop = 20)
  expect_s3_class(vis, "visNetwork")
  expect_s3_class(vis, "htmlwidget")
})

test_that("plot_visgraph supports a non-default color.var and an explicit layout", {
  skip_on_cran()
  pos <- cbind(
    x = stats::rnorm(igraph::vcount(toy_graph)),
    y = stats::rnorm(igraph::vcount(toy_graph))
  )
  rownames(pos) <- igraph::V(toy_graph)$name
  vis <- plot_visgraph(toy_graph, min_rho = 0, ntop = 20, color.var = "layer", layout = pos)
  expect_s3_class(vis, "visNetwork")
})
