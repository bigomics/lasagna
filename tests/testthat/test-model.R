test_that("create_model works with minimal data", {
  set.seed(42)
  n <- 10
  p1 <- 20
  p2 <- 15

  ## create two data layers
  gx <- matrix(rnorm(p1 * n), nrow = p1, ncol = n)
  rownames(gx) <- paste0("gene", 1:p1)
  colnames(gx) <- paste0("S", 1:n)

  px <- matrix(rnorm(p2 * n), nrow = p2, ncol = n)
  rownames(px) <- paste0("prot", 1:p2)
  colnames(px) <- paste0("S", 1:n)

  ## samples data frame
  samples <- data.frame(
    group = factor(rep(c("A", "B"), each = n / 2)),
    row.names = paste0("S", 1:n)
  )

  data <- list(
    X = list(gx = gx, px = px),
    samples = samples
  )

  model <- create_model(data, pheno = "pheno", ntop = 10, nc = 5)

  expect_true(is.list(model))
  expect_true("graph" %in% names(model))
  expect_true("X" %in% names(model))
  expect_true("Y" %in% names(model))
  expect_true("layers" %in% names(model))
  expect_true(igraph::is_igraph(model$graph))
  expect_true(all(c("gx", "px", "PHENO") %in% model$layers))
})

test_that("solve works on created model", {
  set.seed(42)
  n <- 10
  p1 <- 20
  p2 <- 15

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

  data <- list(
    X = list(gx = gx, px = px),
    samples = samples
  )

  model <- create_model(data, pheno = "pheno", ntop = 10, nc = 5)
  pheno <- colnames(model$Y)[1]

  solved <- solve(model, pheno = pheno, max_edges = 50)

  expect_true(igraph::is_igraph(solved))
  expect_true("value" %in% names(igraph::vertex_attr(solved)))
})
