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
        matplot(t(paramMixHMM$fData$Y[statMixHMM$klas == k, ]), type = "l", lty = "dotted", col = colorsvec[k], xlab = "x", ylab = "y(t)")
        title(main = sprintf("Cluster %1.1i", k))
        lines(statMixHMM$smoothed[, k], lwd = 1.5)
      }
    }
  )
)
