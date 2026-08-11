test_that("prune_graph filters by layers and ntop", {
  skip_on_cran()
  pruned <- prune_graph(toy_graph, layers = "gx", ntop = 5, min.rho = 0)
  expect_true(all(igraph::V(pruned)$layer == "gx"))
  expect_lte(igraph::vcount(pruned), 5)
})

test_that("prune_graph filters edges by sign", {
  skip_on_cran()
  pruned_pos <- prune_graph(toy_graph, ntop = 10, min.rho = 0, edge.sign = "pos", prune = FALSE)
  expect_true(all(igraph::E(pruned_pos)$weight > 0))

  pruned_neg <- prune_graph(toy_graph, ntop = 10, min.rho = 0, edge.sign = "neg", prune = FALSE)
  expect_true(all(igraph::E(pruned_neg)$weight < 0))
})

test_that("prune_graph filters edges by type", {
  skip_on_cran()
  pruned_inter <- prune_graph(toy_graph, ntop = 10, min.rho = 0, edge.type = "inter", prune = FALSE)
  expect_true(all(grepl("->", igraph::E(pruned_inter)$connection_type)))

  pruned_intra <- prune_graph(toy_graph, ntop = 10, min.rho = 0, edge.type = "intra", prune = FALSE)
  expect_true(all(!grepl("->", igraph::E(pruned_intra)$connection_type)))
})

test_that("prune_graph applies a minimum rho threshold", {
  skip_on_cran()
  pruned <- prune_graph(toy_graph, ntop = 10, min.rho = 0.3, prune = FALSE)
  expect_true(all(abs(igraph::E(pruned)$weight) >= 0.3))
})
