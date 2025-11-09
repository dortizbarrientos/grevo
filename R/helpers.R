#' @title Symmetrise a square matrix
#' @description Return (A + t(A))/2, ensuring exact symmetry within numerical tolerance.
#' @param A A square numeric matrix.
#' @return A symmetric matrix.
#' @export
symmetrize <- function(A) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) stop("A must be square")
  (A + t(A)) / 2
}

#' @title Check dimensions of two matrices/vectors
#' @param A,B Matrices.
#' @param allow_vector If TRUE, allow B to be a vector conformable with ncol(A).
#' @noRd
.check_dims <- function(A, B, allow_vector = FALSE) {
  if (allow_vector && is.null(dim(B))) {
    if (length(B) != ncol(A)) stop("Length(B) must equal ncol(A)")
    return(invisible(TRUE))
  }
  if (ncol(A) != nrow(B)) stop("ncol(A) must equal nrow(B)")
  invisible(TRUE)
}

#' @title Unit-normalise a vector
#' @param x Numeric vector.
#' @return Unit-length vector (if norm > 0), otherwise x unchanged.
#' @export
unit_vector <- function(x) {
  x <- as.numeric(x)
  n <- sqrt(sum(x^2))
  if (n > 0) x / n else x
}

#' @title Angle between two vectors (radians)
#' @param u,v Numeric vectors of same length.
#' @return Angle in radians (0..pi).
#' @export
angle_between <- function(u, v) {
  u <- as.numeric(u); v <- as.numeric(v)
  if (length(u) != length(v)) stop("u and v must have same length")
  nu <- sqrt(sum(u^2)); nv <- sqrt(sum(v^2))
  if (nu == 0 || nv == 0) stop("Cannot compute angle with zero-length vector")
  c <- sum(u * v) / (nu * nv)
  c <- max(min(c, 1), -1) # clamp for numerical safety
  acos(c)
}

#' @title Rayleigh quotient x^T A x / (x^T x)
#' @param A Symmetric matrix.
#' @param x Vector.
#' @return Numeric scalar.
#' @export
rayleigh_quotient <- function(A, x) {
  A <- symmetrize(A)
  x <- as.numeric(x)
  if (length(x) != nrow(A)) stop("length(x) must equal nrow(A)")
  num <- as.numeric(t(x) %*% A %*% x)
  den <- sum(x^2)
  if (den == 0) stop("x has zero norm")
  num / den
}