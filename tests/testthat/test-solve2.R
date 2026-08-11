test_that("multisolve runs across all traits and returns a valid igraph", {
  skip_on_cran()
  g <- multisolve(toy_model, min_rho = 0, prune = FALSE)
  expect_true(igraph::is_igraph(g))
  expect_equal(igraph::vcount(g), igraph::vcount(toy_model$graph))
  expect_true("value" %in% igraph::vertex_attr_names(g))
})
