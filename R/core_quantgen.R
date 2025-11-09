#' @title Evolvability e(s) = s^T G s
#' @param G Additive genetic (co)variance matrix (SPD).
#' @param s Direction vector (not necessarily unit length).
#' @return Scalar evolvability along s.
#' @export
evolvability <- function(G, s) {
  G <- ensure_spd(G); s <- as.numeric(s)
  if (length(s) != nrow(G)) stop("length(s) must equal nrow(G)")
  as.numeric(t(s) %*% G %*% s)
}

#' @title Conditional evolvability c(s) = 1 / (s^T G^{-1} s)
#' @description Uses a stable solve for G^{-1}. If G is not full rank, a ridge term may be required.
#' @param G SPD matrix.
#' @param s Direction vector.
#' @param ridge Optional non-negative ridge added to the diagonal for numerical stability.
#' @return Scalar conditional evolvability along s.
#' @export
conditional_evolvability <- function(G, s, ridge = 0) {
  G <- ensure_spd(G); s <- as.numeric(s)
  if (length(s) != nrow(G)) stop("length(s) must equal nrow(G)")
  if (ridge > 0) G <- G + diag(ridge, nrow(G))
  x <- as.numeric(t(s) %*% solve(G, s))
  1 / x
}

#' @title Directional heritability h^2_dir(s) = (s^T G s)/(s^T P s)
#' @param G Additive genetic covariance.
#' @param P Phenotypic covariance.
#' @param s Direction vector.
#' @return Scalar directional heritability along s.
#' @export
directional_heritability <- function(G, P, s) {
  G <- ensure_spd(G); P <- ensure_spd(P); s <- as.numeric(s)
  if (length(s) != nrow(G) || nrow(G) != nrow(P)) stop("Dimensions disagree")
  num <- as.numeric(t(s) %*% G %*% s)
  den <- as.numeric(t(s) %*% P %*% s)
  if (den <= 0) stop("Denominator non-positive; check P and s.")
  num / den
}

#' @title Lande's response to selection: Δz = G β
#' @param G Additive genetic covariance.
#' @param beta Selection gradient vector.
#' @return Numeric vector Δz.
#' @export
lande_response <- function(G, beta) {
  G <- ensure_spd(G); beta <- as.numeric(beta)
  if (length(beta) != nrow(G)) stop("length(beta) must equal nrow(G)")
  as.numeric(G %*% beta)
}

#' @title Response diagnostics given G and β
#' @description Computes Δz, response magnitude, angle(β, Δz), and optional directional h^2 along β.
#' @param G SPD matrix.
#' @param beta Selection gradient vector.
#' @param P Optional phenotypic covariance to compute directional h^2 along β.
#' @return List with elements: delta, norm_delta, angle_beta_delta (radians),
#'   beta_norm, h2_dir_beta (if P provided).
#' @export
response_stats <- function(G, beta, P = NULL) {
  delta <- lande_response(G, beta)
  out <- list(
    delta = delta,
    norm_delta = sqrt(sum(delta^2)),
    beta_norm = sqrt(sum(beta^2)),
    angle_beta_delta = angle_between(beta, delta)
  )
  if (!is.null(P)) {
    out$h2_dir_beta <- directional_heritability(G, P, beta)
  }
  out
}

#' @title Alignment of β to the leading eigenvectors of G
#' @param G SPD matrix.
#' @param beta Vector.
#' @param k Number of top eigenvectors to report.
#' @return Data frame with cosines between β and each of the first k eigenvectors.
#' @export
align_to_eigen <- function(G, beta, k = 3) {
  G <- ensure_spd(G); beta <- as.numeric(beta)
  ev <- eigen(G, symmetric = TRUE)
  k <- min(k, ncol(ev$vectors))
  cosines <- sapply(1:k, function(i) {
    vi <- ev$vectors[, i]
    c <- sum(beta * vi) / (sqrt(sum(beta^2)) * sqrt(sum(vi^2)))
    max(min(c, 1), -1)
  })
  data.frame(component = seq_len(k), cosine = cosines, eigenvalue = ev$values[seq_len(k)])
}

#' @title Plot evolvability along directions on the 2D unit circle
#' @description For 2-trait systems only. Returns a data.frame invisibly.
#' @param G 2x2 SPD matrix.
#' @param n Number of angles.
#' @param add Logical; add to existing plot?
#' @param scale Logical; scale radii to [0,1] for visualisation.
#' @export
plot_evolvability_circle2D <- function(G, n = 360, add = FALSE, scale = FALSE) {
  if (nrow(G) != 2 || ncol(G) != 2) stop("This visual is for 2x2 matrices only")
  G <- ensure_spd(G)
  theta <- seq(0, 2*pi, length.out = n+1)[-1]
  S <- cbind(cos(theta), sin(theta))
  e <- apply(S, 1, function(s) evolvability(G, s))
  r <- e
  if (scale) r <- (r - min(r)) / (max(r) - min(r) + .Machine$double.eps)
  x <- r * cos(theta); y <- r * sin(theta)
  if (!add) {
    plot(0, 0, type = "n", asp = 1, xlab = "x", ylab = "y",
         xlim = range(c(-r, r)), ylim = range(c(-r, r)))
    symbols(0, 0, circles = max(r), inches = FALSE, add = TRUE, lty = 3)
  }
  lines(x, y, lwd = 2)
  invisible(data.frame(theta = theta, evolvability = e, x = x, y = y))
}