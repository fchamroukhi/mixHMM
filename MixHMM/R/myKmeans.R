myKmeans <- function(X, K, nbr_runs, nbr_iter_max, verbose){
  #   function res = myKmeans(X, K, nbr_runs, nbr_iter_max, verbose)
  #
  #   Algorithme des K-means
  #
  #
  #
  # Faicel CHAMROUKHI Septembre 2008 (mise a jour)
  #
  #
  #
  # distance euclidienne
  #
  ###########################################################################
  if (nargs()<5){verbose=0}
  if (nargs()<4){nbr_iter_max = 300}
  if (nargs()<3){nbr_runs = 20}


  n = nrow(X)
  p = ncol(X)
  # if one class
  global_mean = apply(X,1,mean)

  if (K==1){
    dmin = apply((X-matrix(c(1),n, 1)%*%global_mean)^2,1,sum)
    solution = append(solution,list(muk = global_mean))
    solution = append(solution,list(klas = matrix(c(1),n,1)))
    solution = append(solution,list(err = sum(dmin)))
    solution = append(solution,list(Zik = matrix(c(1),n,1)))
  }

  nbr_run = 0
  best_solution.err = 99999

  while (nbr_run < nbr_runs){
    nbr_run=nbr_run + 1
    if (nbr_runs>1 && verbose){
      fprintf("Kmeans run n? : #d  \n",nbr_run)
    }
    stored_err=c()
    solution=NULL
    iter = 0
    converged = 0
    previous_err = -99999
    Zik = matrix(c(0),n,K)# partition
    ## 1. Initialization of the centres
    rnd_indx = sample(1:n,n)
    centres = X[rnd_indx[1:K],]
    while (iter<nbr_iter_max && converged!=1){
      iter = iter+1
      old_centres = centres

      # The Euclidean distances
      eucld_dist = matrix(c(0),n, K)
      for (k in 1:K){
        muk = centres[k,]
        eucld_dist[,k] = apply((X-matrix(c(1),n,1)%*%muk)^2,1,sum)
      }
      ## classification step

      klas=c()
      for(i in 1:n){
        dmin = apply(eucld_dist,2,min)
        klas[i] = which.min(eucld_dist[i,])
      }


      Zik = matrix(as.numeric(c((klas%*%matrix(c(1),1,K))==(matrix(c(1),n,1)%*%(1:K)))),nrow=n,ncol=K)
      ## relocation step
      for (k in 1:K){
        ind_ck = which(klas==k)
        #if empty classes
        if (length(ind_ck)!=0){
          centres[k,]= old_centres[k,]
        }else{
          # update the centres
          centres[k,] = apply(X[ind_ck,],1,mean)
        }
      }

      # test of convergence
      current_err = sum(apply(Zik*eucld_dist,2,sum)) # the distorsion measure
      stored_err[iter] = current_err
      # print(iter)
      # print(current_err)
      if ((abs(current_err-previous_err))/previous_err <1e-6){
        converged = 1
      }
      previous_err = current_err
      if (verbose){
        fprintf("Kmeans : Iteration  #d  Objective: #6f  \n", iter, current_err)
      }

    }# one run
    ##
    solution = append(solution,list(stored_err=stored_err))
    solution = append(solution,list(muk = centres))
    solution = append(solution,list(Zik = Zik))
    solution = append(solution,list(klas = klas))
    solution = append(solution,list(err = current_err))
    #
    if (current_err < best_solution.err){
      best_solution = solution
    }
    solution = best_solution

  } #en of the Kmeans runs
  return(solution)
}



