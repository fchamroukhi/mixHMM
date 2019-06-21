
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

<!-- badges: start -->

<!-- badges: end -->

R code for the **clustering** and **segmentation** of time series
(including with regime changes) by mixture of gaussian Hidden Markov
Models (MixHMMs) and the EM algorithm, i.e functional data clustering
and segmentation.

## Installation

You can install the development version of mixHMM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/mixHMM")
```

To build *vignettes* for examples of usage, type the command below
instead:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/mixHMM", 
                         build_opts = c("--no-resave-data", "--no-manual"), 
                         build_vignettes = TRUE)
```

Use the following command to display vignettes:

``` r
browseVignettes("mixHMM")
```

## Usage

``` r
library(mixHMM)

data("simulatedtimeseries")

K <- 3 # Number of clusters
R <- 3 # Number of regimes (HMM states)
variance_type <- "heteroskedastic" # "heteroskedastic" or "homoskedastic" model

ordered_states <- TRUE
n_tries <- 1
max_iter <- 1000
init_kmeans <- TRUE
threshold <- 1e-6
verbose <- TRUE

mixhmm <- emMixHMM(t(simulatedtimeseries[, 2:ncol(simulatedtimeseries)]), K, R, variance_type, ordered_states, n_tries, max_iter, init_kmeans, threshold, verbose)
#> EM: Iteration : 1 || log-likelihood : -28768.651771201
#> EM: Iteration : 2 || log-likelihood : -22339.2695187269
#> EM: Iteration : 3 || log-likelihood : -21960.4137394824
#> EM: Iteration : 4 || log-likelihood : -21838.6832023488
#> EM: Iteration : 5 || log-likelihood : -21826.0254324452
#> EM: Iteration : 6 || log-likelihood : -21825.2945545122
#> EM: Iteration : 7 || log-likelihood : -21825.2614076716
#> EM: Iteration : 8 || log-likelihood : -21825.2600749497

mixhmm$summary()
#> -----------------------
#> Fitted mixHMM model
#> -----------------------
#> 
#> MixHMM model with K = 3 clusters and R = 3 regimes:
#> 
#>  log-likelihood nu       AIC       BIC
#>       -21825.26 44 -21869.26 -21911.32
#> 
#> Clustering table:
#>  1  2  3 
#> 20 15 15 
#> 
#> Mixing probabilities (cluster weights):
#>    1   2   3
#>  0.4 0.3 0.3
#> 
#> 
#> -------------------
#> Cluster 1 (K = 1):
#> 
#> Means:
#> 
#>     R = 1    R = 2    R = 3
#>  10.01848 7.002696 9.002662
#> 
#> Variances:
#> 
#>      R = 1     R = 2    R = 3
#>  0.9287081 0.9728516 1.077248
#> 
#> -------------------
#> Cluster 2 (K = 2):
#> 
#> Means:
#> 
#>    R = 1    R = 2    R = 3
#>  8.03851 11.00551 6.989432
#> 
#> Variances:
#> 
#>      R = 1     R = 2    R = 3
#>  0.9600214 0.9765353 1.034951
#> 
#> -------------------
#> Cluster 3 (K = 3):
#> 
#> Means:
#> 
#>     R = 1   R = 2    R = 3
#>  8.024457 10.9832 10.00617
#> 
#> Variances:
#> 
#>      R = 1   R = 2  R = 3
#>  0.9807032 1.01269 1.0784

mixhmm$plot()
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-2.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-3.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-4.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-5.png" style="display: block; margin: auto;" />
