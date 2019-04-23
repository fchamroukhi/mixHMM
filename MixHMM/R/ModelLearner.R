source("R/ParamMixHMM.R")
source("R/StatMixHMM.R")
source("R/FittedMixHMM.R")
# source("R/enums.R")
# source("R/utils.R")
# source("R/ParamHMMR.R")
# source("R/StatHMMR.R")

EM <- function(modelMixHMM, order_constraint = TRUE, n_tries = 1, max_iter = 1000, init_kmeans = TRUE, threshold = 1e-6, verbose = TRUE) {
  #
  # The EM algorithm for parameter estimation of the mixture of Hidden Markov
  # Models for clustering and segmentation of time series with regime changes
  #
  # Inputs:
  # data: a set of n time series with m observations (dim: [n x m]
  # K: number of clusters
  # R: number of regimes (states)
  # options
  #
  #
  #
  # faicel chamroukhi (septembre 2009)
  #
  ## Please cite the following references for this code
  #
  # @InProceedings{Chamroukhi-IJCNN-2011,
  #   author = {F. Chamroukhi and A. Sam\'e  and P. Aknin and G. Govaert},
  #   title = {Model-based clustering with Hidden Markov Model regression for time series with regime changes},
  #   Booktitle = {Proceedings of the International Joint Conference on Neural Networks (IJCNN), IEEE},
  #   Pages = {2814--2821},
  #   Adress = {San Jose, California, USA},
  #   year = {2011},
  #   month = {Jul-Aug},
  #   url = {https://chamroukhi.com/papers/Chamroukhi-ijcnn-2011.pdf}
  # }
  #
  # @PhdThesis{Chamroukhi_PhD_2010,
  # author = {Chamroukhi, F.},
  # title = {Hidden process regression for curve modeling, classification and tracking},
  # school = {Universit\'e de Technologie de Compi\`egne},
  # month = {13 december},
  # year = {2010},
  # type = {Ph.D. Thesis},
  # url ={https://chamroukhi.com/papers/FChamroukhi-Thesis.pdf}
  # }

  ###############################################################################################

  try_EM = 0
  best_loglik = -Inf
  cputime_total = c()

  while (try_EM < n_tries){
    try_EM = try_EM +1
    print(paste("EM try n?",try_EM))
    start_time = Sys.time()

    ###################
    #  Initialization #
    ###################

    param = ParamMixHMM(modelMixHMM)
    param$init_MixFHMM(modelMixHMM, order_constraint, init_kmeans, try_EM)

    iter = 0
    converged = FALSE
    # loglik = 0
    prev_loglik=-Inf
    # stored_loglik=c()
    # stats=list(mask=param$mask)

    stat <- StatMixHMM(modelMixHMM)

    # main algorithm
    # # EM ####
    while ((iter <= max_iter) & !converged){

      ##########
      # E-Step #
      ##########
      stat$EStep(modelMixHMM, param)

      ##########
      # M-Step #
      ##########
      param$MStep(modelMixHMM, stat, order_constraint)

      iter = iter + 1

      if (verbose){
        print(paste("EM : Iteration :",iter,"log-likelihood : ",stat$loglik))
      }

      if (prev_loglik-stat$loglik > 1e-4) {
        paste("!!!!! EM log-lik is decreasing from", prev_loglik,"to",stat$loglik)
      }

      converged <- (abs((stat$loglik-prev_loglik)/prev_loglik)< threshold)
      if (is.na(converged)) {
        converged <- FALSE
      } # Basically for the first iteration when prev_loglik is Inf

      prev_loglik = stat$loglik
      stat$stored_loglik[iter] = stat$loglik

      }# end of EM  loop

    cputime_total <- cbind(cputime_total, Sys.time() - start_time)

    if (stat$loglik > best_loglik){
      statSolution <- stat$copy()
      paramSolution <- param$copy()

      best_loglik <- stat$loglik
    }

    if (try_EM>=1){
      print( paste('log-lik at convergence:', stat$loglik))
    }

  }

  if (try_EM>1){
    print(paste('log-lik max:', statSolution$loglik))
  }

  # Finding the curve partition by using the MAP rule
  statSolution$MAP()

  # FINISH computation of statSolution
  statSolution$computeStats(modelMixHMM, paramSolution, cputime_total)

  return(FittedMixHMM(modelMixHMM, paramSolution, statSolution))

}
