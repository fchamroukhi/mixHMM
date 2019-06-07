#' @export
StatMixHMM <- setRefClass(
  "StatMixHMM",
  fields = list(
    tau_ik = "matrix",
    gamma_ikjr = "array",
    log_w_k_fyi = "matrix",
    exp_num_trans = "array",
    exp_num_trans_from_l = "array",
    loglik = "numeric",
    stored_loglik = "list",
    cputime = "numeric",
    klas = "matrix",
    z_ik = "matrix",
    smoothed = "matrix",
    mean_curves = "array",
    BIC = "numeric",
    AIC = "numeric",
    ICL1 = "numeric"
  ),
  methods = list(

    initialize = function(paramMixHMM = ParamMixHMM(fData = FData(numeric(1), matrix(1)), K = 2, R = 1, variance_type = 1)) {

      tau_ik <<- matrix(NA, paramMixHMM$fData$n, paramMixHMM$K)
      gamma_ikjr <<- array(NA, dim = c(paramMixHMM$fData$n * paramMixHMM$fData$m, paramMixHMM$R, paramMixHMM$K))
      log_w_k_fyi <<- matrix(NA, paramMixHMM$fData$n, paramMixHMM$K)
      exp_num_trans <<- array(NA, dim = c(paramMixHMM$R, paramMixHMM$R, paramMixHMM$fData$n, paramMixHMM$K))
      exp_num_trans_from_l <<- array(NA, dim = c(paramMixHMM$R, paramMixHMM$fData$n, paramMixHMM$K))
      loglik <<- -Inf
      stored_loglik <<- list()
      cputime <<- Inf
      klas <<- matrix(NA, paramMixHMM$fData$n, 1) # klas: [nx1 double]
      z_ik <<- matrix(NA, paramMixHMM$fData$n, paramMixHMM$K) # z_ik: [nxK]
      smoothed <<- matrix(NA, paramMixHMM$fData$m, paramMixHMM$K)
      mean_curves <<- array(NA, dim = c(paramMixHMM$fData$m, paramMixHMM$R, paramMixHMM$K))
      BIC <<- -Inf
      AIC <<- -Inf
      ICL1 <<- -Inf

    },

    MAP = function() {
      N <- nrow(tau_ik)
      K <- ncol(tau_ik)
      ikmax <- max.col(tau_ik)
      ikmax <- matrix(ikmax, ncol = 1)
      z_ik <<- ikmax %*% ones(1, K) == ones(N, 1) %*% (1:K) # partition_MAP
      klas <<- ones(N, 1)
      for (k in 1:K) {
        klas[z_ik[, k] == 1] <<- k
      }
    },

    computeStats = function(paramMixHMM, cputime_total) {
      cputime <<- mean(cputime_total)

      for (k in 1:paramMixHMM$K) {
        weighted_segments <- apply(gamma_ikjr[, , k] * (paramMixHMM$fData$vecY %*% matrix(1, 1, paramMixHMM$R)), 1, sum)

        dim(weighted_segments) <- c(paramMixHMM$fData$m, paramMixHMM$fData$n)
        weighted_clusters <- (matrix(1, paramMixHMM$fData$m, 1) %*% tau_ik[, k]) * weighted_segments
        smoothed[, k] <<- (apply(weighted_clusters, 1, sum)) / sum(tau_ik[, k])
      }

      # BIC AIC et ICL*
      BIC <<- loglik - (paramMixHMM$nu * log(paramMixHMM$fData$n) / 2)
      AIC <<- loglik - paramMixHMM$nu
      # ICL*
      # Compute the comp-log-lik
      cik_log_w_k_fyi <- (z_ik) * (log_w_k_fyi)
      comp_loglik <- sum(cik_log_w_k_fyi)
      ICL1 <<- comp_loglik - paramMixHMM$nu * log(paramMixHMM$fData$n) / 2 #n*m/2!

    },

    EStep = function(paramMixHMM) {
      exp_num_trans_ck  <- array(0, dim = c(paramMixHMM$R, paramMixHMM$R, paramMixHMM$fData$n))
      exp_num_trans_from_l_ck <- matrix(0, paramMixHMM$R, paramMixHMM$fData$n)

      w_k_fyi <- matrix(0, paramMixHMM$fData$n, paramMixHMM$K)

      for (k in 1:paramMixHMM$K) {
        # Run a hmm for each sequence
        log_fkr_yij <- matrix(0, paramMixHMM$R, paramMixHMM$fData$m)
        fkr_yij <- matrix(0, paramMixHMM$R, paramMixHMM$fData$m)

        Li <- matrix(0, paramMixHMM$fData$n, 1) # To store the loglik for each example

        mu_kr <- paramMixHMM$mu_kr[, k]
        num_log_post_prob <- matrix(0, paramMixHMM$fData$n, paramMixHMM$K)

        for (i in 1:paramMixHMM$fData$n) {
          y_i <- paramMixHMM$fData$Y[i, ]

          for (r in 1:paramMixHMM$R) {
            mukr <- mu_kr[r]

            if (paramMixHMM$variance_type == variance_types$homoskedastic) {
              sigma2_kr <- paramMixHMM$sigma2_kr[, k]
              sk <- sigma2_kr
            } else {
              sigma2_kr <- paramMixHMM$sigma2_kr[, k]
              sk <- sigma2_kr[r]
            }
            z <- ((y_i - mukr * matrix(1, 1, paramMixHMM$fData$m)) ^ 2) / sk
            log_fkr_yij[r, ] <- -0.5 * matrix(1, 1, paramMixHMM$fData$m) * (log(2 * pi) + log(sk)) - 0.5 * z# pdf cond ? c_i = g et z_i = k de yij
            fkr_yij[r, ] <- dnorm(y_i, mukr * matrix(1, 1, paramMixHMM$fData$m), sqrt(sk))

          }

          # Calcul of p(y) : forwards backwards
          fb <- forwardsBackwards(paramMixHMM$pi_k[, k], paramMixHMM$A_k[, , k], fkr_yij)

          gamma_ik <- fb$tau_tk
          xi_ik <- fb$xi_tk
          fwd_ik <- fb$alpha_tk
          backw_ik <- fb$beta_tk
          loglik_i <- fb$loglik

          Li[i] <- loglik_i # Loglik of the ith curve

          gamma_ikjr[(((i - 1) * paramMixHMM$fData$m + 1):(i * paramMixHMM$fData$m)), , k] <<- t(gamma_ik) # [n*m K G]

          exp_num_trans_ck[, , i] <- apply(xi_ik, MARGIN = c(1, 2), sum) # [K K n]
          exp_num_trans_from_l_ck[, i] <- gamma_ik[, 1] # [K x n]
        }

        exp_num_trans_from_l[, , k] <<- exp_num_trans_from_l_ck # [K n G]
        exp_num_trans[, , , k] <<- exp_num_trans_ck # [K K n G]

        # For the MAP partition:  the numerator of the cluster post probabilities
        # num_log_post_prob[,k] <- log(param$w_k[k]) + Li

        # For computing the global loglik
        w_k_fyi[, k] <- paramMixHMM$w_k[k] * exp(Li) # [nx1]

        log_w_k_fyi[, k] <<- log(paramMixHMM$w_k[k]) + Li
      }

      log_w_k_fyi <<- pmin(log_w_k_fyi, log(.Machine$double.xmax))
      log_w_k_fyi <<- pmax(log_w_k_fyi, log(.Machine$double.xmin))

      tau_ik <<- exp(log_w_k_fyi) / (apply(exp(log_w_k_fyi), 1, sum) %*% matrix(1, 1, paramMixHMM$K))

      # Log-likelihood
      loglik <<- sum(log(apply(exp(log_w_k_fyi), 1, sum)))

    }
  )
)
