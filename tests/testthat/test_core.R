test_that("rayleigh and evolvability agree for unit s", {
  set.seed(1)
  G <- matrix(c(0.5,0.2,0.2, 0.2,0.4,0.1, 0.2,0.1,0.3), 3, 3)
  s <- c(1,2,3); s <- s / sqrt(sum(s^2))
  expect_equal(evolvability(G, s), rayleigh_quotient(G, s), tolerance = 1e-10)
})

test_that("directional h2 is bounded between 0 and 1 if P = G + E (E >= 0)", {
  G <- matrix(c(0.5,0.1,0.1, 0.1,0.4,0.05, 0.1,0.05,0.3), 3, 3)
  E <- diag(c(0.2, 0.2, 0.2))
  P <- G + E
  s <- c(1,0,-1)
  h2 <- directional_heritability(G, P, s)
  expect_gt(h2, 0); expect_lt(h2, 1)
})

test_that("affine distance equals zero for identical matrices", {
  G <- matrix(c(0.5,0.2,0.1, 0.2,0.4,0.05, 0.1,0.05,0.3), 3, 3)
  expect_lt(affine_dist(G, G), 1e-10)
})