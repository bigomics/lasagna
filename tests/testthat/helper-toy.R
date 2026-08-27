## Shared fixtures for the plotting/graph/solve tests: small toy models
## built once at suite load and reused across test_that() blocks so the
## suite stays fast.

make_toy_data <- function(n = 10, p1 = 20, p2 = 15, seed = 42) {
  set.seed(seed)
  gx <- matrix(rnorm(p1 * n), nrow = p1, ncol = n)
  rownames(gx) <- paste0("gene", 1:p1)
  colnames(gx) <- paste0("S", 1:n)

  px <- matrix(rnorm(p2 * n), nrow = p2, ncol = n)
  rownames(px) <- paste0("prot", 1:p2)
  colnames(px) <- paste0("S", 1:n)

  samples <- data.frame(
    group = factor(rep(c("A", "B"), each = n / 2)),
    row.names = paste0("S", 1:n)
  )

  list(X = list(gx = gx, px = px), samples = samples)
}

toy_model <- suppressMessages(
  create_model(make_toy_data(), meta.type = "pheno", ntop = 10, nc = 5)
)
toy_pheno <- colnames(toy_model$Y)[1]
toy_graph <- solve(toy_model, pheno = toy_pheno, max_edges = 50, prune = FALSE)

## a layer reduced to a single feature, the shape that used to crash
## layout_multipartite_3d's per-layer svd
single_feature_model <- suppressMessages(
  create_model(make_toy_data(n = 8, p1 = 20, p2 = 1, seed = 7),
    meta.type = "pheno", ntop = 10, nc = 5)
)
