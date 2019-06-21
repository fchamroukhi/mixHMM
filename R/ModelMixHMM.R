#' @export
ModelMixHMM <- setRefClass(
  "ModelMixHMM",
  fields = list(
    paramMixHMM = "ParamMixHMM",
    statMixHMM = "StatMixHMM"
  ),
  methods = list(

    plot = function() {

      # yaxislim <- c(min(paramMixHMM$fData$Y) - 2 * mean(sqrt(apply(paramMixHMM$fData$Y, 1, var))), max(paramMixHMM$fData$Y) + 2 * mean(sqrt(apply(paramMixHMM$fData$Y, 1, var))))

      matplot(t(paramMixHMM$fData$Y), type = "l", lty = "solid", col = "black", xlab = "x", ylab = "y(t)")
      title(main = "Original time series")

      colorsvec <- rainbow(paramMixHMM$K)
      matplot(t(paramMixHMM$fData$Y), type = "l", lty = "dotted", col = colorsvec[statMixHMM$klas], xlab = "x", ylab = "y(t)")
      title(main = "Clustered time series")

      for (k in 1:paramMixHMM$K) {
        matplot(t(paramMixHMM$fData$Y[statMixHMM$klas == k,]), type = "l", lty = "dotted", col = colorsvec[k], xlab = "x", ylab = "y(t)")
        title(main = sprintf("Cluster %1.1i", k))
        lines(statMixHMM$smoothed[, k], lwd = 1.5)
      }
    },

    summary = function() {
      digits = getOption("digits")

      title <- paste("Fitted mixHMM model")
      txt <- paste(rep("-", min(nchar(title) + 4, getOption("width"))), collapse = "")

      # Title
      cat(txt)
      cat("\n")
      cat(title)
      cat("\n")
      cat(txt)

      cat("\n")
      cat("\n")
      cat(paste0("MixHMM model with K = ", paramMixHMM$K,ifelse(paramMixHMM$K > 1, " clusters", " cluster"), " and R = ", paramMixHMM$R, ifelse(paramMixHMM$R > 1, " regimes", " regime"), ":"))
      cat("\n")
      cat("\n")

      tab <- data.frame("log-likelihood" = statMixHMM$loglik, "nu" = paramMixHMM$nu,
                        "AIC" = statMixHMM$AIC, "BIC" = statMixHMM$BIC,
                        row.names = "", check.names = FALSE)
      print(tab, digits = digits)

      cat("\nClustering table:")
      print(table(statMixHMM$klas))

      cat("\nMixing probabilities (cluster weights):\n")
      pro <- data.frame(t(paramMixHMM$w_k))
      colnames(pro) <- 1:paramMixHMM$K
      print(pro, digits = digits, row.names = FALSE)

      cat("\n\n")

      txt <- paste(rep("-", min(nchar(title), getOption("width"))), collapse = "")

      for (k in 1:paramMixHMM$K) {
        cat(txt)
        cat("\nCluster ", k, " (K = ", k, "):\n", sep = "")

        cat("\nMeans:\n\n")
        means <- data.frame(t(paramMixHMM$mu_kr[, k]))
        colnames(means) <- sapply(1:paramMixHMM$R, function(x) paste0("R = ", x))
        print(means, digits = digits, row.names = FALSE)

        cat(paste0(ifelse(paramMixHMM$variance_type == "homoskedastic", "\n", "\nVariances:\n\n")))
        sigma2 <- data.frame(t(paramMixHMM$sigma2_kr[, k]))
        if (paramMixHMM$variance_type == "homoskedastic") {
          colnames(sigma2) <- "Sigma2"
          print(sigma2, digits = digits, row.names = FALSE)
        } else {
          colnames(sigma2) = sapply(1:paramMixHMM$R, function(x)
            paste0("R = ", x))
          print(sigma2, digits = digits, row.names = FALSE)
        }
        cat("\n")
      }

    }
  )
)
