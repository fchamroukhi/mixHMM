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
    # klas: [nx1 double]
    z_ik = "matrix",
    # z_ik: [nxK]
    smoothed = "matrix",
    mean_curves = "array",
    BIC = "numeric",
    AIC = "numeric",
    ICL1 = "numeric"
    # tau_tk = "matrix", # tau_tk: smoothing probs: [nxK], tau_tk(t,k) = Pr(z_i=k | y1...yn)
    # alpha_tk = "matrix", # alpha_tk: [nxK], forwards probs: Pr(y1...yt,zt=k)
    # beta_tk = "matrix", # beta_tk: [nxK], backwards probs: Pr(yt+1...yn|zt=k)
    # xi_tkl = "array", # xi_tkl: [(n-1)xKxK], joint post probs : xi_tk\elll(t,k,\ell)  = Pr(z_t=k, z_{t-1}=\ell | Y) t =2,..,n
    # f_tk = "matrix", # f_tk: [nxK] f(yt|zt=k)
    # log_f_tk = "matrix", # log_f_tk: [nxK] log(f(yt|zt=k))
    # loglik = "numeric", # loglik: log-likelihood at convergence
    # stored_loglik = "list", # stored_loglik: stored log-likelihood values during EM
    # cputime = "numeric", # cputime: for the best run
    # klas = "matrix", # klas: [nx1 double]
    # z_ik = "matrix", # z_ik: [nxK]
    # state_probs = "matrix", # state_probs: [nxK]
    # BIC = "numeric", # BIC
    # AIC = "numeric", # AIC
    # regressors = "matrix", # regressors: [nxK]
    # predict_prob = "matrix", # predict_prob: [nxK]: Pr(zt=k|y1...y_{t-1})
    # predicted = "matrix", # predicted: [nx1]
    # filter_prob = "matrix", # filter_prob: [nxK]: Pr(zt=k|y1...y_t)
    # filtered = "matrix", # filtered: [nx1]
    # smoothed_regressors = "matrix", # smoothed_regressors: [nxK]
    # smoothed = "matrix" # smoothed: [nx1]
    # #           X: [nx(p+1)] regression design matrix
    # #           nu: model complexity
    # #           parameter_vector
  ),
  methods = list(
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
    #######
    # compute the final solution stats
    #######
    computeStats = function(modelMixHMM, paramMixHMM, cputime_total) {
      cputime <<- mean(cputime_total)

      for (k in 1:modelMixHMM$K) {
        weighted_segments <- apply(gamma_ikjr[, , k] * (modelMixHMM$vecY %*% matrix(1, 1, modelMixHMM$R)), 1, sum)
        #weighted_segments <- sum(gamma_ikjr(:,:,g).*(ones(n*m,1)*param.mu_kr(:,k)'),2);

        dim(weighted_segments) <- c(modelMixHMM$m, modelMixHMM$n)
        weighted_clusters <- (matrix(1, modelMixHMM$m, 1) %*% tau_ik[, k]) * weighted_segments
        smoothed[, k] <<- (apply(weighted_clusters, 1, sum)) / sum(tau_ik[, k])
      }

      #stats.mean_gamma_ijk = mean_gamma_ijk;

      # # Segmentation of each cluster using the MAP rule
      # for k=1:K
      #     [segments_k Zjk] = MAP(stats.mean_gamma_ijk(:,:,k));#MAP segmentation of each cluster of sequences
      #     stats.segments(:,k) = segments_k;
      # end

      # BIC AIC et ICL*
      BIC <<- loglik - (modelMixHMM$nu * log(modelMixHMM$n) / 2)
      AIC <<- loglik - modelMixHMM$nu
      # ICL*
      # Compute the comp-log-lik
      cik_log_w_k_fyi <- (z_ik) * (log_w_k_fyi)
      comp_loglik <- sum(cik_log_w_k_fyi)
      ICL1 <<- comp_loglik - modelMixHMM$nu * log(modelMixHMM$n) / 2 #n*m/2!

    },
    #######
    # EStep
    #######
    EStep = function(modelMixHMM, paramMixHMM) {
      exp_num_trans_ck  <- array(0, dim = c(modelMixHMM$R, modelMixHMM$R, modelMixHMM$n))
      exp_num_trans_from_l_ck <- matrix(0, modelMixHMM$R, modelMixHMM$n)

      w_k_fyi <- matrix(0, modelMixHMM$n, modelMixHMM$K)

      for (k in 1:modelMixHMM$K) {
        # run a hmm for each sequence
        log_fkr_yij <- matrix(0, modelMixHMM$R, modelMixHMM$m)
        fkr_yij <- matrix(0, modelMixHMM$R, modelMixHMM$m)
        #
        Li <- matrix(0, modelMixHMM$n, 1)# to store the loglik for each example
        #
        mu_kr <- paramMixHMM$mu_kr[, k]
        num_log_post_prob <- matrix(0, modelMixHMM$n, modelMixHMM$K)

        for (i in 1:modelMixHMM$n) {
          y_i <- modelMixHMM$Y[i, ]

          for (r in 1:modelMixHMM$R) {
            mukr <- mu_kr[r]
            #sk <- sigma_kr(k);

            if (modelMixHMM$variance_type == variance_types$homoskedastic) {
              sigma_kr <- paramMixHMM$sigma_kr[, k]
              sk <- sigma_kr
            }
            else{
              sigma_kr <- paramMixHMM$sigma_kr[, k]
              sk <- sigma_kr[r]
            }
            z <- ((y_i - mukr * matrix(1, 1, modelMixHMM$m)) ^ 2) / sk
            log_fkr_yij[r, ] <- -0.5 * matrix(1, 1, modelMixHMM$m) * (log(2 * pi) + log(sk)) - 0.5 * z# pdf cond ? c_i = g et z_i = k de yij
            fkr_yij[r, ] <- dnorm(y_i, mukr * matrix(1, 1, modelMixHMM$m), sqrt(sk))

          }

          #log_fkr_yij  <- min(log_fkr_yij,log(realmax));
          #log_fkr_yij <- max(log_fkr_yij ,log(realmin));
          #fkr_yij <-  exp(log_fkr_yij);

          # calcul de p(y) : forwards backwards

          fb <- forwards_backwards(paramMixHMM$pi_k[, k], paramMixHMM$A_k[, , k], fkr_yij)

          gamma_ik <- fb$tau_tk
          xi_ik <- fb$xi_tk
          fwd_ik <- fb$alpha_tk
          backw_ik <- fb$beta_tk
          loglik_i <- fb$loglik

          #
          Li[i] <- loglik_i # loglik of the ith curve

          #
          gamma_ikjr[(((i - 1) * modelMixHMM$m + 1):(i * modelMixHMM$m)), , k] <<- t(gamma_ik)#[n*m K G]
          # xi_ikjrl(:,:,(i-1)*(m-1)+1:i*(m-1),g) <-  xi_ik;#[KxK n*m G]
          #

          exp_num_trans_ck[, , i] <- apply(xi_ik, MARGIN = c(1, 2), sum)# [K K n]
          exp_num_trans_from_l_ck[, i] <- gamma_ik[, 1]#[K x n]
          #
        }

        exp_num_trans_from_l[, , k] <<- exp_num_trans_from_l_ck#[K n G]
        exp_num_trans[, , , k] <<- exp_num_trans_ck#[K K n G]

        # # for the MAP partition:  the numerator of the cluster post
        # # probabilities
        # num_log_post_prob[,k] <- log(param$w_k[k]) + Li

        # for computing the global loglik
        w_k_fyi[, k] <- paramMixHMM$w_k[k] * exp(Li)#[nx1]

        log_w_k_fyi[, k] <<- log(paramMixHMM$w_k[k]) + Li
      }

      log_w_k_fyi <<- pmin(log_w_k_fyi, log(.Machine$double.xmax))
      log_w_k_fyi <<- pmax(log_w_k_fyi, log(.Machine$double.xmin))

      tau_ik <<- exp(log_w_k_fyi) / (apply(exp(log_w_k_fyi), 1, sum) %*% matrix(1, 1, modelMixHMM$K))
      # # log-likelihood
      loglik <<- sum(log(apply(exp(log_w_k_fyi), 1, sum)))

    }
  )
)


StatMixHMM <- function(modelMixHMM) {
  tau_ik <- matrix(NA, modelMixHMM$n, modelMixHMM$K)
  gamma_ikjr <- array(NA, dim = c(modelMixHMM$n * modelMixHMM$m, modelMixHMM$R, modelMixHMM$K))
  log_w_k_fyi <- matrix(NA, modelMixHMM$n, modelMixHMM$K)
  exp_num_trans <- array(NA, dim = c(modelMixHMM$R, modelMixHMM$R, modelMixHMM$n, modelMixHMM$K))
  exp_num_trans_from_l <- array(NA, dim = c(modelMixHMM$R, modelMixHMM$n, modelMixHMM$K))
  loglik <- -Inf
  stored_loglik <- list()
  cputime <- Inf
  klas <- matrix(NA, modelMixHMM$n, 1) # klas: [nx1 double]
  z_ik <- matrix(NA, modelMixHMM$n, modelMixHMM$K) # z_ik: [nxK]
  smoothed <- matrix(NA, modelMixHMM$m, modelMixHMM$K)
  mean_curves <- array(NA, dim = c(modelMixHMM$m, modelMixHMM$R, modelMixHMM$K))
  BIC <- -Inf
  AIC <- -Inf
  ICL1 <- -Inf
  # tau_tk <- matrix(NA, modelMixHMM$m, modelMixHMM$K) # tau_tk: smoothing probs: [nxK], tau_tk(t,k) = Pr(z_i=k | y1...yn)
  # alpha_tk <- matrix(NA, modelMixHMM$m, ncol = modelMixHMM$K) # alpha_tk: [nxK], forwards probs: Pr(y1...yt,zt=k)
  # beta_tk <- matrix(NA, modelMixHMM$m, modelMixHMM$K) # beta_tk: [nxK], backwards probs: Pr(yt+1...yn|zt=k)
  # xi_tkl <- array(NA, c(modelMixHMM$m - 1, modelMixHMM$K, modelMixHMM$K)) # xi_tkl: [(n-1)xKxK], joint post probs : xi_tk\elll(t,k,\ell)  = Pr(z_t=k, z_{t-1}=\ell | Y) t =2,..,n
  # f_tk <- matrix(NA, modelMixHMM$m, modelMixHMM$K) # f_tk: [nxK] f(yt|zt=k)
  # log_f_tk <- matrix(NA, modelMixHMM$m, modelMixHMM$K) # log_f_tk: [nxK] log(f(yt|zt=k))
  # loglik <- -Inf # loglik: log-likelihood at convergence
  # stored_loglik <- list() # stored_loglik: stored log-likelihood values during EM
  # cputime <- Inf # cputime: for the best run
  # klas <- matrix(NA, modelMixHMM$m, 1) # klas: [nx1 double]
  # z_ik <- matrix(NA, modelMixHMM$m, modelMixHMM$K) # z_ik: [nxK]
  # state_probs <- matrix(NA, modelMixHMM$m, modelMixHMM$K) # state_probs: [nxK]
  # BIC <- -Inf # BIC
  # AIC <- -Inf # AIC
  # regressors <- matrix(NA, modelMixHMM$m, modelMixHMM$K) # regressors: [nxK]
  # predict_prob <- matrix(NA, modelMixHMM$m, modelMixHMM$K) # predict_prob: [nxK]: Pr(zt=k|y1...y_{t-1})
  # predicted <- matrix(NA, modelMixHMM$m, 1) # predicted: [nx1]
  # filter_prob <- matrix(NA, modelMixHMM$m, modelMixHMM$K) # filter_prob: [nxK]: Pr(zt=k|y1...y_t)
  # filtered <- matrix(NA, modelMixHMM$m, 1) # filtered: [nx1]
  # smoothed_regressors <- matrix(NA, modelMixHMM$m, modelMixHMM$K) # smoothed_regressors: [nxK]
  # smoothed <- matrix(NA, modelMixHMM$m, 1) # smoothed: [nx1]

  new(
    "StatMixHMM",
    tau_ik = tau_ik,
    gamma_ikjr = gamma_ikjr,
    log_w_k_fyi = log_w_k_fyi,
    exp_num_trans = exp_num_trans,
    exp_num_trans_from_l = exp_num_trans_from_l,
    loglik = loglik,
    stored_loglik = stored_loglik,
    cputime = cputime,
    klas = klas,
    z_ik = z_ik,
    smoothed = smoothed,
    mean_curves = mean_curves,
    BIC = BIC,
    AIC = AIC,
    ICL1 = ICL1
  )
}
