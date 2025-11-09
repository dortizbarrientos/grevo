#' @title Load grevo toy dataset
#' @description Loads phenotypes and ground-truth matrices from \code{inst/extdata}.
#' @details Returns a list with elements:
#' \itemize{
#' \item \code{pheno}: data.frame with id, ecotype, trait1..trait3
#' \item \code{G}: 3x3 matrix (true additive genetic covariance)
#' \item \code{E}: 3x3 matrix (true environmental covariance)
#' \item \code{P}: 3x3 matrix (phenotypic covariance = G + E)
#' \item \code{beta}: length-3 selection gradient used in examples
#' }
#' @return List.
#' @export
load_toy <- function() {
  base <- system.file("extdata", package = "grevo")
  if (base == "") stop("extdata not found; is the package installed properly?")
  pheno <- utils::read.csv(file.path(base, "toy_pheno.csv"), stringsAsFactors = FALSE)
  G <- as.matrix(utils::read.csv(file.path(base, "toy_true_G.csv"), header = FALSE))
  E <- as.matrix(utils::read.csv(file.path(base, "toy_true_E.csv"), header = FALSE))
  P <- as.matrix(utils::read.csv(file.path(base, "toy_true_P.csv"), header = FALSE))
  beta <- as.numeric(utils::read.csv(file.path(base, "toy_beta.csv"), header = FALSE)[,1])
  list(pheno = pheno, G = G, E = E, P = P, beta = beta)
}