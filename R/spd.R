#' @title Ensure symmetric positive-definite (SPD) via eigenvalue floor
#' @description Symmetrises A and floors eigenvalues at \code{eps * max(1, max(eigs))}.
#' @param A Square numeric matrix.
#' @param eps Non-negative floor parameter.
#' @return SPD matrix.
#' @export
ensure_spd <- function(A, eps = 1e-8) {
  A <- symmetrize(A)
  ev <- eigen(A, symmetric = TRUE)
  lam <- ev$values
  floor <- eps * max(1, max(lam))
  lam <- pmax(lam, floor)
  ev$vectors %*% (lam * t(ev$vectors))
}

#' @title Matrix square-root for SPD matrices
#' @param A SPD matrix.
#' @return A^{1/2}.
#' @export
spd_sqrtm <- function(A) {
  A <- ensure_spd(A)
  ev <- eigen(A, symmetric = TRUE)
  ev$vectors %*% (sqrt(pmax(ev$values, 0)) * t(ev$vectors))
}

#' @title Inverse square-root for SPD matrices
#' @param A SPD matrix.
#' @return A^{-1/2}.
#' @export
spd_invsqrtm <- function(A) {
  A <- ensure_spd(A)
  ev <- eigen(A, symmetric = TRUE)
  ev$vectors %*% ((1 / sqrt(pmax(ev$values, .Machine$double.eps))) * t(ev$vectors))
}

#' @title Matrix logarithm for SPD matrices
#' @param A SPD matrix.
#' @return log(A) in the affine-invariant sense (via eigen-decomposition).
#' @export
spd_logm <- function(A) {
  A <- ensure_spd(A)
  ev <- eigen(A, symmetric = TRUE)
  ev$vectors %*% (log(pmax(ev$values, .Machine$double.eps)) * t(ev$vectors))
}

#' @title Affine-invariant geodesic distance between two G-matrices
#' @description Frobenius norm of log(G1^{-1/2} G2 G1^{-1/2}).
#' @param G1,G2 SPD matrices of equal dimension.
#' @return Non-negative scalar distance.
#' @export
affine_dist <- function(G1, G2) {
  G1 <- ensure_spd(G1); G2 <- ensure_spd(G2)
  S <- spd_invsqrtm(G1) %*% G2 %*% spd_invsqrtm(G1)
  M <- spd_logm(S)
  sqrt(sum(M * M))
}

#' @title Geodesic interpolation between SPD matrices
#' @param G1,G2 SPD matrices.
#' @param t Scalar in [0,1].
#' @return SPD matrix along the affine-invariant geodesic from G1 to G2.
#' @export
geodesic_interpolate <- function(G1, G2, t) {
  if (t < 0 || t > 1) stop("t must be in [0,1]")
  G1 <- ensure_spd(G1); G2 <- ensure_spd(G2)
  A <- spd_invsqrtm(G1) %*% G2 %*% spd_invsqrtm(G1)
  ev <- eigen(A, symmetric = TRUE)
  At <- ev$vectors %*% ((ev$values^t) * t(ev$vectors))
  spd_sqrtm(G1) %*% At %*% spd_sqrtm(G1)
}