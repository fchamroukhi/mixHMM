source("R/enums.R")
source("R/utils.R")
source("R/myKmeans.R")
source("R/mk_stochastic.R")

ParamMixHMM <- setRefClass(
  "ParamMixHMM",
  fields = list(
    w_k = "matrix",
    # prior = "matrix",
    pi_k = "matrix",
    # Initial distributions
    # trans_mat = "matrix",
    A_k = "array",
    # Transition matrices
    mu_kr = "matrix",
    # Means
    # beta = "matrix",
    sigma_kr = "matrix",
    # Standard deviations
    mask = "matrix"
  ),
  methods = list(
    init_MixFHMM = function(modelMixHMM, order_constraint = TRUE, init_kmeans = TRUE, try_algo = 1) {
      # # 1. Initialization of cluster weights
      w_k <<- 1 / K * matrix(1, modelMixHMM$K, 1)

      # Initialization of the model parameters for each cluster
      if (init_kmeans) {
        max_iter_kmeans <- 400
        n_tries_kmeans <- 20
        verbose_kmeans <- 0
        solution <- myKmeans(modelMixHMM$Y, modelMixHMM$K, n_tries_kmeans, max_iter_kmeans, verbose_kmeans)

        for (k in 1:modelMixHMM$K) {
          Yk <- modelMixHMM$Y[solution$klas == k , ] #if kmeans

          init_gauss_hmm(Yk, k, modelMixHMM$R, modelMixHMM$variance_type, order_constraint, try_algo)

        }

      } else{
        ind <- sample(1:modelMixHMM$n, modelMixHMM$n)
        for (k in 1:modelMixHMM$K) {
          if (k < modelMixHMM$K) {
            Yk <- modelMixHMM$Y[ind[(k - 1) * round(modelMixHMM$n / modelMixHMM$K) + 1:k * round(modelMixHMM$n / modelMixHMM$K)], ]
          }
          else{
            Yk <- modelMixHMM$Y[ind[(k - 1) * round(modelMixHMM$n / modelMixHMM$K) + 1:modelMixHMM$n], ]
          }

          init_gauss_hmm(Yk, k, modelMixHMM$R, order_constraint, variance_type, try_algo)

        }
      }
    },


    init_gauss_hmm = function(Y, k, R, variance_type, order_constraint = TRUE, try_algo) {
      # init_gauss_hmm  estime les paramètres initiaux d'un hmm où la loi conditionnelle des observations est une gaussienne
      #
      # Entrees :
      #
      #        Y(i,:,nsignal) = x(i) : observation à l'instant i du signal
      #        (séquence) nsignal (notez que pour la partie parametrisation des
      #        signaux les observations sont monodimentionnelles)
      #        R : nbre d'états (classes) cachés
      #
      # Sorties :
      #
      #         model : parametres initiaux du modele. structure
      #         contenant les champs: para: structrure with the fields:
      #         * le HMM initial
      #         1. initial_prob (k) = Pr(Z(1) = k) avec k=1,...,K. loi initiale de z.
      #         2. trans_mat(\ell,k) = Pr(z(i)=k | z(i-1)=\ell) : matrice des transitions
      #         *
      #         3.1. mur : moyenne de l'état k
      #         3.2 sigma2r(k) = variance de x(i) sachant z(i)=k; sigma2r(j) =
      #         sigma^2_r.
      #         mu(:,k) = Esperance de x(i) sachant z(i) = k ;
      ################################################################################

      if (order_constraint) {
        # # Initialisation en tenant compte de la contrainte:

        # Initialisation de la matrice des transitions
        maskM <- diag(R)#mask d'ordre 1
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
        A_k[, , k] <<- mk_stochastic(matrix(runif(R), R, R))
      }

      #  Initialisation des moyennes et des variances
      init_gauss_param_hmm(Y, k, R, variance_type, try_algo)
    },

    ###################################################################################
    init_gauss_param_hmm = function(Y, k, R, variance_type, try_algo) {
      # init_regression_model estime les parametres de la loi conditionnelle
      # des observations : une gaussienne d'un hmm homogène d'ordre 1
      #
      # Entrees :
      #
      #        Y : [nxm]
      #        nsignal (notez que pour la partie parametrisation des signaux les
      #        observations sont monodimentionnelles)
      #        R : nbre d'états (classes) cachés
      # Sorties :
      #
      #
      #         para : parametres initiaux de la loi cond de chaque état
      #         2. sigma2r(r) = variance de y(t) sachant z(t)=r; sigmar(j) =
      #         sigma^2_r.
      #         3. mu(:,r) : E[y(t)|z(t) =r] ;
      ################################################################################

      n <- nrow(Y)
      m <- ncol(Y)

      if (variance_type == variance_types$homoskedastic) {
        s <- 0
      }

      if (try_algo == 1) {

        ##############################
        #decoupage de l'echantillon (signal) en K segments
        zi <- round(m / R) - 1
        for (r in 1:R) {
          i <- (r - 1) * zi + 1
          j <- r * zi

          Yij <- Y[, i:j]
          Yij <- matrix(t(Yij), 1, byrow = T)
          mu_kr[r, k] <<- mean(Yij)

          if (variance_type == variance_types$homoskedastic) {
            s <- s + sum((Yij - mu_kr[r, k]) ^ 2)
            sigma_kr[, k] <<- s / (n * m)
          } else {
            m_r <- j - i + 1
            sigma_kr[r, k] <<- sum((Yij - mu_kr[r, k]) ^ 2) / (n * m_r)
          }
        }

      } else {
        # initialisation aléatoire
        Lmin <- 2#round(m/(K+1));#nbr pts min dans un segments
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
            sigma_kr[, k] <<- s / (n * m)
          } else{
            m_r <- j - i + 1
            sigma_kr[r, k] <<- sum((Yij - mu_kr[r, k]) ^ 2) / (n * m_r)
          }
        }
      }
    },


    MStep = function(modelMixHMM, statMixHMM, order_constraint = TRUE) {
      # Maximization of Q1 w.r.t w_k

      w_k <<- matrix(apply(statMixHMM$tau_ik, 2, sum)) / modelMixHMM$n

      exp_num_trans_k <- array(0, dim = c(modelMixHMM$R, modelMixHMM$R, modelMixHMM$n))

      for (k in 1:modelMixHMM$K) {
        if (modelMixHMM$variance_type == variance_types$homoskedastic) {
          s <- 0
        }

        weights_cluster_k <- statMixHMM$tau_ik[, k]
        # Maximization of Q2 w.r.t \pi^g
        exp_num_trans_k_from_l <- (matrix(1, modelMixHMM$R, 1) %*% t(weights_cluster_k)) * statMixHMM$exp_num_trans_from_l[, , k]#[K x n]

        # pi_k[,k] <<- (1/sum(statMixHMM$tau_ik[,k]))*sum(exp_num_trans_k_from_l,2)# sum over i

        pi_k[, k] <<- (1 / sum(statMixHMM$tau_ik[, k])) * apply(exp_num_trans_k_from_l, 1, sum) # sum over i

        # Maximization of Q3 w.r.t A^g

        for (r in 1:modelMixHMM$R) {
          # squeeze=c()
          # for (i in 1:modelMixHMM$n){
          #   squeeze=cbind(squeeze,statMixHMM$exp_num_trans[r,,i,k])
          # }

          if (modelMixHMM$n == 1) {
            exp_num_trans_k[r, , ] <- t(matrix(1, modelMixHMM$R, 1) %*% weights_cluster_k) * drop(statMixHMM$exp_num_trans[r, , , k])
          } else{
            exp_num_trans_k[r, , ] <- (matrix(1, modelMixHMM$R, 1) %*% t(weights_cluster_k)) * drop(statMixHMM$exp_num_trans[r, , , k])
          }
        }

        if (modelMixHMM$n == 1) {
          temp <- exp_num_trans_k
        } else{
          temp <- apply(exp_num_trans_k, MARGIN = c(1, 2), sum)#sum over i
        }


        A_k[, , k] <<- mk_stochastic(temp)

        # if HMM with order constraints
        if (order_constraint) {
          A_k[, , k] <<- mk_stochastic(mask * A_k[, , k])
        }

        # Maximisation de Q4 par rapport aux muk et sigmak
        # each sequence i (m observations) is first weighted by the cluster weights

        weights_cluster_k <- matrix(t(statMixHMM$tau_ik[, k]), nrow = modelMixHMM$m, ncol = ncol(t(statMixHMM$tau_ik)), byrow = T)
        weights_cluster_k <- matrix(as.vector(weights_cluster_k), length(as.vector(weights_cluster_k)), 1)

        # secondly, the m observations of each sequance are weighted by the
        # wights of each segment k (post prob of the segments for each
        # cluster g)
        gamma_ijk <- statMixHMM$gamma_ikjr[, , k]# [n*m K]

        nm_kr <- apply(gamma_ijk, 1, sum)# cardinal nbr of the segments k,k=1,...,K within each cluster g, at iteration q

        for (r in 1:modelMixHMM$R) {
          nmkr <- nm_kr[r]#cardinal nbr of segment k for the cluster g
          # # Maximization w.r.t muk
          weights_seg_k <- matrix(gamma_ijk[, r])

          mu_kr[r, k] <<- (1 / sum(weights_cluster_k * weights_seg_k)) %*% sum((weights_cluster_k * weights_seg_k) * modelMixHMM$vecY)
          # mukr[r] = (1/sum(weights_cluster_k*weights_seg_k))%*%sum((weights_cluster_k*weights_seg_k)*modelMixHMM$vecY)

          # # Maximization w.r.t sigmak :)

          z <- sqrt(weights_cluster_k * weights_seg_k) * (modelMixHMM$vecY - matrix(1, modelMixHMM$n *
                                                                                   modelMixHMM$m, 1) * mu_kr[r, k])
          # z <- sqrt(weights_cluster_k*weights_seg_k)*(modelMixHMM$vecY-matrix(1,modelMixHMM$n*modelMixHMM$m,1)*mukr[r])

          if (modelMixHMM$variance_type == variance_types$homoskedastic) {
            s <- s + (t(z) %*% z)
            ngm <- sum(apply((weights_cluster_k %*% matrix(1, 1, modelMixHMM$R)) * gamma_ijk, 1, sum))

            # sigma_k <- s/ngm
            sigma_kr[k] <<- s / ngm
          }
          else{
            ngmk <- sum(weights_cluster_k * weights_seg_k)

            # sigmakr[r] <-  (t(z)%*%z)/(ngmk)
            sigma_kr[r, k] <<- (t(z) %*% z) / (ngmk)
          }
        }

        # mu_kr[,k] <<- mukr

        # if (modelMixHMM$variance_type==variance_types$homoskedastic){
        #   sigma_kr[k] <<- sigma_k
        # } else{
        #   sigma_kr[,k] <<- sigmakr
        # }
      }
    }
  )
)

ParamMixHMM <- function(modelMixHMM) {
  w_k <- matrix(NA, nrow = modelMixHMM$K)
  # prior <- matrix(NA, ncol = modelMixHMM$K - 1)
  pi_k <- matrix(NA, nrow = modelMixHMM$R, ncol = modelMixHMM$K)
  # trans_mat <- matrix(NA, modelMixHMM$K, modelMixHMM$K)
  A_k <- array(NA, dim = c(modelMixHMM$R, modelMixHMM$R, modelMixHMM$K))
  # beta <- matrix(NA, modelMixHMM$p + 1, modelMixHMM$K)
  mu_kr <- matrix(NA, nrow = modelMixHMM$R, ncol = modelMixHMM$K)
  if (modelMixHMM$variance_type == variance_types$homoskedastic) {
    sigma_kr <- matrix(NA, ncol = modelMixHMM$K)
  }
  else{
    sigma_kr <- matrix(NA, nrow = modelMixHMM$R, ncol = modelMixHMM$K)
  }
  mask <- matrix(NA, modelMixHMM$R, modelMixHMM$R)
  new(
    "ParamMixHMM",
    w_k = w_k,
    pi_k = pi_k,
    A_k = A_k,
    mu_kr = mu_kr,
    sigma_kr = sigma_kr,
    mask = mask
  )
}
