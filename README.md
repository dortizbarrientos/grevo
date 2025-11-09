# grevo: Genetic Response & Evolvability

Tools for response to selection and evolvability from the G-matrix.

`grevo` provides compact tools to study the response to selection (Lande's equation),
directional heritability, evolvability, and geometry on the space of genetic covariance
matrices (SPD). It includes a small toy dataset and a vignette.

## Install locally

1. Download and unzip the attached `grevo.zip`.
2. In R:

```r
install.packages("devtools")          # if needed
devtools::install("path/to/grevo")   # replace with the unzipped path
library(grevo)
```

## Quick start

```r
toy <- load_toy()
res <- response_stats(toy$G, toy$beta, toy$P)
res$angle_beta_delta    # radians
directional_heritability(toy$G, toy$P, toy$beta)  # h^2 along beta
plot_evolvability_circle2D(toy$G[1:2,1:2], scale = TRUE)
```
