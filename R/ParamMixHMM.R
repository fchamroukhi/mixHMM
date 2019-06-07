#' @export
ParamMixHMM <- setRefClass(
  "ParamMixHMM",
  fields = list(
    fData = "FData",

    K = "numeric", # Number of clusters
    R = "numeric", # Number of regimes (HMM states)
    variance_type = "numeric",
    nu = "numeric", # Degree of freedom

    w_k = "matrix", # Cluster weights
    pi_k = "matrix", # Initial distributions
    A_k = "array", # Transition matrices
    mu_kr = "matrix", # Means
    sigma2_kr = "matrix", # Variances
    mask = "matrix"
  ),
  methods = list(

    initialize = function(fData = FData(numeric(1), matrix(1)), K = 2, R = 1, variance_type = 1) {

      fData <<- fData

      K <<- K
      R <<- R
      variance_type <<- variance_type

      if (variance_type == variance_types$homoskedastic) {
        nu <<- (K - 1) + K * ((R - 1) + R * (R - 1) + R + 1)
      }
      else{
        nu <<- (K - 1) + K * ((R - 1) + R * (R - 1) + R + R)
      }

      w_k <<- matrix(NA, nrow = K)
      pi_k <<- matrix(NA, nrow = R, ncol = K)
      A_k <<- array(NA, dim = c(R, R, K))
      mu_kr <<- matrix(NA, nrow = R, ncol = K)

      if (variance_type == variance_types$homoskedastic) {
        sigma2_kr <<- matrix(NA, ncol = K)
      } else {
        sigma2_kr <<- matrix(NA, nrow = R, ncol = K)
      }
      mask <<- matrix(NA, R, R)

    },

    initMixHMM = function(order_constraint = TRUE, init_kmeans = TRUE, try_algo = 1) {
      # 1. Initialization of cluster weights
      w_k <<- 1 / K * matrix(1, K, 1)

      # Initialization of the model parameters for each cluster
      if (init_kmeans) {
        max_iter_kmeans <- 400
        n_tries_kmeans <- 20
        verbose_kmeans <- 0
        solution <- kmeans(fData$Y, K, n_tries_kmeans, max_iter_kmeans, verbose_kmeans)

        for (k in 1:K) {

          Yk <- fData$Y[solution$klas == k , ] # If kmeans
          initGaussHmm(Yk, k, R, variance_type, order_constraint, try_algo)

        }

      } else{

        ind <- sample(1:fData$n, fData$n)
        for (k in 1:K) {
          if (k < K) {
            Yk <- fData$Y[ind[((k - 1) * round(fData$n / K) + 1):(k * round(fData$n / K))], ]
          }
          else{
            Yk <- fData$Y[ind[((k - 1) * round(fData$n / K) + 1):fData$n], ]
          }

          initGaussHmm(Yk, k, R, variance_type, order_constraint, try_algo)

        }
      }
    },


    initGaussHmm = function(Y, k, R, variance_type, order_constraint = TRUE, try_algo) {
      # initGaussHmm  estime les parametres initiaux d'un hmm la loi conditionnelle des observations est une gaussienne
      #
      # Entrees :
      #
      #        Y(i,:,nsignal) = x(i) : observation a l'instant i du signal
      #        (sequence) nsignal (notez que pour la partie parametrisation des
      #        signaux les observations sont monodimentionnelles)
      #        R : nbre d'etats (classes) caches
      #
      # Sorties :
      #
      #         model : parametres initiaux du modele. structure
      #         contenant les champs: para: structrure with the fields:
      #         * le HMM initial
      #         1. initial_prob (k) = Pr(Z(1) = k) avec k=1,...,K. loi initiale de z.
      #         2. trans_mat(\ell,k) = Pr(z(i)=k | z(i-1)=\ell) : matrice des transitions
      #         *
      #         3.1. mur : moyenne de l'??tat k
      #         3.2 sigma2r(k) = variance de x(i) sachant z(i)=k; sigma2r(j) =
      #         sigma2_r.
      #         mu(:,k) = Esperance de x(i) sachant z(i) = k ;
      ################################################################################

      if (order_constraint) {
        # Initialisation en tenant compte de la contrainte:

        # Initialisation de la matrice des transitions
        maskM <- diag(R) # Mask of order 1
        for (r in 1:R - 1) {
          ind <- which(maskM[r, ] != 0)
          maskM[r, ind + 1] <- 1
        }

        # Initialisation de la loi initiale de la variable cachee
        pi_k[, k] <<- c(1, matrix(0, R - 1, 1))

        A_k[, , k] <<- normalize(maskM, 2)$M
        mask <<- maskM

      } else {

        # Initialisation de la loi initiale de la variable cachee
        pi_k[, k] <<- 1 / R * matrix(1, R, 1)
        A_k[, , k] <<- mkStochastic(matrix(runif(R), R, R))
      }

      # Initialisation des moyennes et des variances
      initGaussParamHmm(Y, k, R, variance_type, try_algo)
    },

    ###################################################################################
    initGaussParamHmm = function(Y, k, R, variance_type, try_algo) {
      # init_regression_model estime les parametres de la loi conditionnelle
      # des observations : une gaussienne d'un hmm homog??ne d'ordre 1
      #
      # Entrees :
      #
      #        Y : [nxm]
      #        nsignal (notez que pour la partie parametrisation des signaux les
      #        observations sont monodimentionnelles)
      #        R : nbre d'??tats (classes) cach??s
      # Sorties :
      #
      #
      #         para : parametres initiaux de la loi cond de chaque ??tat
      #         2. sigma2r(r) = variance de y(t) sachant z(t)=r; sigmar(j) =
      #         sigma2_r.
      #         3. mu(:,r) : E[y(t)|z(t) =r] ;
      ################################################################################

      n <- nrow(Y)
      m <- ncol(Y)

      if (variance_type == variance_types$homoskedastic) {
        s <- 0
      }

      if (try_algo == 1) {

        zi <- round(m / R) - 1
        for (r in 1:R) {
          i <- (r - 1) * zi + 1
          j <- r * zi

          Yij <- Y[, i:j]
          Yij <- matrix(t(Yij), 1, byrow = T)
          mu_kr[r, k] <<- mean(Yij)

          if (variance_type == variance_types$homoskedastic) {
            s <- s + sum((Yij - mu_kr[r, k]) ^ 2)
            sigma2_kr[, k] <<- s / (n * m)
          } else {
            m_r <- j - i + 1
            sigma2_kr[r, k] <<- sum((Yij - mu_kr[r, k]) ^ 2) / (n * m_r)
          }
        }

      } else {

        Lmin <- 2
        tr_init <- matrix(0, 1, R + 1)
        tr_init[1] <- 0
        R_1 <- R
        for (r in 2:R) {
          R_1 <- R_1 - 1
          temp <- seq(tr_init[r - 1] + Lmin, m - R_1 * Lmin)
          ind <- sample(1:length(temp), length(temp))
          tr_init[r] <- temp[ind[1]]
        }

        tr_init[R + 1] <- m
        for (r in 1:R) {
          i <- tr_init[r] + 1
          j <- tr_init[r + 1]
          Yij <- Y[, i:j]
          Yij <- matrix(t(Yij), ncol = 1, byrow = T)

          mu_kr[r, k] <<- mean(Yij)

          if (variance_type == variance_types$homoskedastic) {
            s <- s + sum((Yij - mu_kr[r, k]) ^ 2)
            sigma2_kr[, k] <<- s / (n * m)
          } else {
            m_r <- j - i + 1
            sigma2_kr[r, k] <<- sum((Yij - mu_kr[r, k]) ^ 2) / (n * m_r)
          }
        }
      }
    },

    MStep = function(statMixHMM, order_constraint = TRUE) {

      # Maximization of Q1 w.r.t w_k

      w_k <<- matrix(apply(statMixHMM$tau_ik, 2, sum)) / fData$n

      exp_num_trans_k <- array(0, dim = c(R, R, fData$n))

      for (k in 1:K) {
        if (variance_type == variance_types$homoskedastic) {
          s <- 0
        }

        weights_cluster_k <- statMixHMM$tau_ik[, k]

        # Maximization of Q2 w.r.t \pi^g
        exp_num_trans_k_from_l <- (matrix(1, R, 1) %*% t(weights_cluster_k)) * statMixHMM$exp_num_trans_from_l[, , k] # [K x n]

        pi_k[, k] <<- (1 / sum(statMixHMM$tau_ik[, k])) * apply(exp_num_trans_k_from_l, 1, sum) # Sum over i

        # Maximization of Q3 w.r.t Ag

        for (r in 1:R) {

          if (fData$n == 1) {
            exp_num_trans_k[r, , ] <- t(matrix(1, R, 1) %*% weights_cluster_k) * drop(statMixHMM$exp_num_trans[r, , , k])
          } else{
            exp_num_trans_k[r, , ] <- (matrix(1, R, 1) %*% t(weights_cluster_k)) * drop(statMixHMM$exp_num_trans[r, , , k])
          }
        }

        if (fData$n == 1) {
          temp <- exp_num_trans_k
        } else{
          temp <- apply(exp_num_trans_k, MARGIN = c(1, 2), sum) # Sum over i
        }


        A_k[, , k] <<- mkStochastic(temp)

        # If HMM with order constraints
        if (order_constraint) {
          A_k[, , k] <<- mkStochastic(mask * A_k[, , k])
        }

        # Maximisation de Q4 par rapport aux muk et sigmak
        # each sequence i (m observations) is first weighted by the cluster weights

        weights_cluster_k <- matrix(t(statMixHMM$tau_ik[, k]), nrow = fData$m, ncol = ncol(t(statMixHMM$tau_ik)), byrow = T)
        weights_cluster_k <- matrix(as.vector(weights_cluster_k), length(as.vector(weights_cluster_k)), 1)

        # secondly, the m observations of each sequance are weighted by the
        # wights of each segment k (post prob of the segments for each
        # cluster g)
        gamma_ijk <- statMixHMM$gamma_ikjr[, , k]# [n*m K]

        nm_kr <- apply(gamma_ijk, 1, sum) # Cardinal nbr of the segments k,k=1,...,K within each cluster g, at iteration q

        for (r in 1:R) {
          nmkr <- nm_kr[r] # Cardinal nbr of segment k for the cluster g

          # Maximization w.r.t muk
          weights_seg_k <- matrix(gamma_ijk[, r])

          mu_kr[r, k] <<- (1 / sum(weights_cluster_k * weights_seg_k)) %*% sum((weights_cluster_k * weights_seg_k) * fData$vecY)

          # Maximization w.r.t sigmak
          z <- sqrt(weights_cluster_k * weights_seg_k) * (fData$vecY - matrix(1, fData$n * fData$m, 1) * mu_kr[r, k])

          if (variance_type == variance_types$homoskedastic) {
            s <- s + (t(z) %*% z)
            ngm <- sum(apply((weights_cluster_k %*% matrix(1, 1, R)) * gamma_ijk, 1, sum))
            sigma2_kr[k] <<- s / ngm
          } else{
            ngmk <- sum(weights_cluster_k * weights_seg_k)
            sigma2_kr[r, k] <<- (t(z) %*% z) / (ngmk)
          }
        }

      }
    }
  )
)
