source("R/FData.R")
source("R/enums.R")

ModelMixHMM <- setRefClass(
  "ModelMixHMM",
  contains = "FData",
  # Define the fields
  fields = list(
    K = "numeric", # Number of clusters
    R = "numeric", # Number of regimes (HMM states)
    variance_type = "numeric",
    nu = "numeric" # degree of freedom
  )
)

ModelMixHMM <- function(fData, K, R, variance_type) {
  if (variance_type == variance_types$homoskedastic) {
    nu <<- K * (R - 1 + R * (R - 1) + R + 1)
  }
  else{
    nu <<- K * (R - 1 + R * (R - 1) + R + R)
  }

  new(
    "ModelMixHMM",
    Y = fData$Y,
    X = fData$X,
    m = fData$m,
    n = fData$n,
    vecY = fData$vecY,
    K = K,
    R = R,
    variance_type = variance_type,
    nu = nu
  )
}
