FittedMixHMM <- setRefClass(
  "FittedMixHMM",
  fields = list(
    modelMixHMM = "ModelMixHMM",
    paramMixHMM = "ParamMixHMM",
    statMixHMM = "StatMixHMM"
  ),
  methods = list(
    plot = function() {

      # yaxislim <- c(min(modelMixHMM$Y) - 2 * mean(sqrt(apply(modelMixHMM$Y, 1, var))), max(modelMixHMM$Y) + 2 * mean(sqrt(apply(modelMixHMM$Y, 1, var))))

      matplot(t(modelMixHMM$Y), type = "l", lty = "solid", col = "black", xlab = "Time", ylab = "y(t)")
      title(main = "Original time series")


      colorsvec <- rainbow(modelMixHMM$K)
      matplot(t(modelMixHMM$Y), type = "l", lty = "dotted", col = colorsvec[statMixHMM$klas], xlab = "Time", ylab = "y(t)")
      title(main = "Clustered time series")

      for (k in 1:modelMixHMM$K) {
        matplot(t(modelMixHMM$Y[statMixHMM$klas == k, ]), type = "l", lty = "dotted", col = colorsvec[k], xlab = "Time", ylab = "y(t)")
        title(main = sprintf("Cluster %1.1i", k))
        lines(statMixHMM$smoothed[, k], lwd = 1.5)
      }
    }
  )
)

FittedMixHMM <- function(modelMixHMM, paramMixHMM, statMixHMM) {
  new("FittedMixHMM", modelMixHMM = modelMixHMM, paramMixHMM = paramMixHMM, statMixHMM = statMixHMM)
}
