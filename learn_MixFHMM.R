#function mixFHMM =  learn_MixFHMM(data, K, R, ...

learn_MixFHMM<- function(data, K, R,variance_type,order_constraint, total_EM_tries, max_iter_EM, init_kmeans, threshold, verbose){
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
  #faicel chamroukhi (septembre 2009)
  #
  ##Please cite the following references for this code
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
  source("init_MixFHMM.R")
  source("forwards_backwards.R")
  source("normalize.R")
  source("mk_stochastic.R")
  source("MAP.R")
  source("myKmeans.R")
  
  options(warn=-1)
  n=nrow(data)
  m=ncol(data)

  #n  nbre de signaux (individus); m: nbre de points pour chaque signal
  #
  Y=matrix(t(data),ncol=1,byrow=T)
  
  #
  try_EM = 0
  best_loglik = -Inf
  cputime_total = c()
  
  while (try_EM < total_EM_tries){
    try_EM = try_EM +1
    print(paste('EM try n°',try_EM))
    start_time = Sys.time()
    
    ###################
    #  Initialization #
    ###################
    
    param = init_MixFHMM(data, K, R,variance_type, order_constraint, init_kmeans, try_EM)
    
    #     Psi = zeros(nu,1);# vecteur parametre
    iter = 0
    converge = 0
    loglik = 0
    prev_loglik=-Inf
    stored_loglik=c()
    stats=list(mask=param$mask)
    # main algorithm
    # # EM ####
    while ((iter <= max_iter_EM) & (converge!=1)){
      
      #
      exp_num_trans_ck  = array(c(0),dim=c(R,R,n))
      exp_num_trans_from_l_ck = matrix(c(0),R,n)
      #
      exp_num_trans = array(c(0),dim=c(R,R,n,K))
      exp_num_trans_from_l = array(c(0),dim=c(R,n,K))
      #
      w_k_fyi = matrix(c(0),n,K)
      log_w_k_fyi = matrix(c(0),n,K)
      fkr_yij=matrix(c(0),K,m)
      
      ##########
      # E-Step #
      ##########
      gamma_ikjr = array(c(0),dim=c(n*m,R,K))
      for (k in 1:K){
        
        # run a hmm for each sequence
        log_fkr_yij =matrix(c(0),R,m)
        #
        Li = matrix(c(0),n,1)# to store the loglik for each example
        #
        mu_kr = param$mu_kr[,k]
        num_log_post_prob = matrix(c(0),n,K)
        
        for (i in 1:n){
          y_i = data[i,]
          
          for (r in 1:R){
            mukr = mu_kr[r]
            #sk = sigma_kr(k);  
            
            if (variance_type=='common'){
              sigma_kr = param$sigma_kr[k]
              sk = sigma_kr
            }
            else{
              sigma_kr = param$sigma_kr[,k]
              sk = sigma_kr[r]
            }
            z=((y_i-mukr*matrix(c(1),1,m))^2)/sk
            log_fkr_yij[r,] = -0.5*matrix(c(1),1,m)*(log(2*pi)+log(sk)) - 0.5*z# pdf cond à c_i = g et z_i = k de yij
            fkr_yij[r,] = dnorm(y_i, mukr*matrix(c(1),1,m), sqrt(sk))
            
          }
          
          #log_fkr_yij  = min(log_fkr_yij,log(realmax));
          #log_fkr_yij = max(log_fkr_yij ,log(realmin));
          #fkr_yij =  exp(log_fkr_yij);
          
          # calcul de p(y) : forwards backwards
          
          fb = forwards_backwards(param$pi_k[,k], param$A_k[,,k], fkr_yij)
          
          gamma_ik= fb$tau_tk
          xi_ik= fb$xi_tk
          fwd_ik= fb$alpha_tk
          backw_ik= fb$beta_tk
          loglik_i= fb$loglik
          
          #
          Li[i] = loglik_i # loglik of the ith curve
          
          #
          gamma_ikjr[(((i-1)*m+1):(i*m)),,k] = t(gamma_ik)#[n*m K G]
          # xi_ikjrl(:,:,(i-1)*(m-1)+1:i*(m-1),g) =  xi_ik;#[KxK n*m G]
          #
          
          exp_num_trans_ck[,,i] = apply(xi_ik, MARGIN=c(1, 2), sum)# [K K n]
          exp_num_trans_from_l_ck[,i] = gamma_ik[,1]#[K x n]
          #
        }
        
        exp_num_trans_from_l[,,k] = exp_num_trans_from_l_ck#[K n G]
        exp_num_trans[,,,k] = exp_num_trans_ck#[K K n G]
        
        # for the MAP partition:  the numerator of the cluster post
        # probabilities
        num_log_post_prob[,k] = log(param$w_k[k]) + Li
        
        # for computing the global loglik
        w_k_fyi[,k] = param$w_k[k]*exp(Li)#[nx1]
        
        log_w_k_fyi[,k] = log(param$w_k[k]) + Li
      }
      for (k in 1:K){
        for (i in 1:nrow(log_w_k_fyi)){
          log_w_k_fyi[i,k]  = min(log_w_k_fyi[i,k],log(.Machine$double.xmax))
          log_w_k_fyi[i,k] = max(log_w_k_fyi[i,k] ,log(.Machine$double.xmin))
        }
      }
      tau_ik = round(exp(log_w_k_fyi)/(apply(exp(log_w_k_fyi),1,sum)%*%matrix(c(1),1,K)),5)
      # # log-likelihood
      loglik = sum(log(apply(exp(log_w_k_fyi),1,sum)))

      ##########
      # M-Step #
      ##########
      
      # Maximization of Q1 w.r.t w_k
      w_k = t(apply(tau_ik,2,sum))/n
      pi_k=param$pi_k
      A_k=param$A_k
      mukr=matrix(c(0),R,K)
      sigmakr=matrix(c(0),R,K)
      exp_num_trans_k=array(c(0),dim=c(R,K,n))
                     
      for (k in 1:K){
        
        if (variance_type=='common'){
          s=0
        }
        
        weights_cluster_k = tau_ik[,k]
        # Maximization of Q2 w.r.t \pi^g
        exp_num_trans_k_from_l =   (matrix(c(1),R,1)%*%t(weights_cluster_k))*exp_num_trans_from_l[,,k]#[K x n]
        pi_k[,k] = (1/sum(tau_ik[,k]))*sum(exp_num_trans_k_from_l,2)# sum over i
        # Maximization of Q3 w.r.t A^g
        
        for (r in 1:R){
          squeeze=c()
          for (i in 1:n){
            squeeze=cbind(squeeze,exp_num_trans[r,,i,k])
          }  
          
          if (n==1){
            exp_num_trans_k[r,,] = t(matrix(c(1),R,1)%*%weights_cluster_k)*(squeeze)#### ACHTUNG peut etre erreur
          }
          else{
            exp_num_trans_k[r,,] = (matrix(c(1),R,1)%*%t(weights_cluster_k))*(squeeze)
          }
        }
        
        if (n==1){
          temp = exp_num_trans_k
        }else{
          temp = apply(exp_num_trans_k,MARGIN=c(1,2),sum)#sum over i
        }
        
        A_k[,,k] = mk_stochastic(temp)
        # if HMM with order constraints
        if (order_constraint==1){
          A_k[,,k] = mk_stochastic(stats$mask %*% A_k[,,k])
        }
        
        # Maximisation de Q4 par rapport aux muk et sigmak
        # each sequence i (m observations) is first weighted by the cluster weights
        weights_cluster_k =  matrix(t(tau_ik[,k]),nrow=m,ncol=ncol(t(tau_ik)),byrow = T)
        weights_cluster_k = matrix(c(as.vector(weights_cluster_k)),length(as.vector(weights_cluster_k)),1)
        
        # secondly, the m observations of each sequance are weighted by the
        # wights of each segment k (post prob of the segments for each
        # cluster g)
        gamma_ijk = gamma_ikjr[,,k]# [n*m K]
        
        nm_kr=apply(gamma_ijk,1,sum)# cardinal nbr of the segments k,k=1,...,K within each cluster g, at iteration q
        
        sigma_kr = matrix(c(0),R,1)
        sigmak=c()
        
        for (r in 1:R){
          nmkr = nm_kr[r]#cardinal nbr of segment k for the cluster g
          # # Maximization w.r.t muk
          weights_seg_k = gamma_ijk[,r]
          mu_kr[r] = (1/sum(weights_cluster_k*weights_seg_k))%*%sum((weights_cluster_k*weights_seg_k)*Y)
          # # Maximization w.r.t sigmak :)
          z = sqrt(weights_cluster_k*weights_seg_k)*(Y-matrix(c(1),n*m,1)*mu_kr[r])

          if (variance_type=='common'){
            
            s = s + (t(z)%*%z)
            ngm = sum(apply((weights_cluster_k%*%matrix(c(1),1,R))*gamma_ijk,1,sum))
            sigma_k = s/ngm
          }
          else{
            ngmk = sum(weights_cluster_k*weights_seg_k)
            sigma_kr[r]=  (t(z)%*%z)/(ngmk)
          }
          
          mukr[,k] = mu_kr
          
          if (variance_type=='common'){
            sigmak[k] = sigma_k
          }
          
          else{
            sigmakr[,k] = sigma_kr
          }
          
        }
        
        iter=iter+1
        
        if (round(prev_loglik-loglik,4) > 1e-3){
          paste('!!!!! EM log-lik is decreasing from', prev_loglik,'to',loglik)
        } 
        if (verbose==1){
          print(paste('EM : Iteration :',iter,'log-likelihood : ',loglik))
        }
        
        converge =  (abs((loglik-prev_loglik)/prev_loglik)< threshold)
        prev_loglik = loglik
        stored_loglik[iter] = loglik
        end_time=Sys.time()
        cputime_total[iter]=c(end_time-start_time)
      }# end of EM  loop
      
      param = NULL
      param=list(w_k=w_k,pi_k=pi_k,A_k=A_k,mu_kr=mukr,sigma_k=sigmak,sigma_kr=sigmakr)
      
      if (variance_type=='common'){
        stats=append(stats,list(Psi = rbind(as.vector(w_k),as.vector(A_k),as.vector(pi_k),as.vector(mukr),as.vector(sigmak))))
      }else{
        stats=append(stats,list(Psi = rbind(as.vector(w_k),as.vector(A_k),as.vector(pi_k),as.vector(mukr),as.vector(sigmakr))))
      }
      
      stats=append(stats,list(tau_ik = tau_ik))
      stats=append(stats,list(gamma_ikjr = gamma_ikjr))
      stats=append(stats,list(loglik = loglik))
      stats=append(stats,list(stored_loglik = stored_loglik))
      stats=append(stats,list(log_w_k_fyi = log_w_k_fyi))
      
      if (stats$loglik > best_loglik){
        best_loglik = stats$loglik
        best_stats = stats
      }
      if (try_EM>=1){
       print( paste('log-lik at convergence:', stats$loglik))
      }
    }
    stats = append(stats,list(loglik = best_loglik))
    #
    if (try_EM>1){
      print(paste('log-lik max:', stats$loglik))
    }
    
    stats = best_stats
    # Finding the curve partition by using the MAP rule
    stats=append(stats,MAP(stats$tau_ik))# MAP partition of the n sequences
    
    # cas ou on prend la moyenne des gamma ijkr
    smoothed = matrix(c(0),m,K)
    mean_curves = array(c(0),dim=c(m,R,K))
    #mean_gamma_ijk = zeros(m,R,K)
    
    for (k in 1:K){
      weighted_segments = apply(gamma_ikjr[,,k]*(Y%*%matrix(c(1),1,R)),1,sum)
      #weighted_segments = sum(gamma_ikjr(:,:,g).*(ones(n*m,1)*param.mu_kr(:,k)'),2);
      
      #
      weighted_segments = matrix(t(weighted_segments),nrow=m,ncol=n,byrow=T)
      weighted_clusters = (matrix(c(1),m,1)%*%stats$tau_ik[,k])* weighted_segments
      smoothed[,k] = (apply(weighted_clusters,1,sum))/sum(stats$tau_ik[,k])
      print(stats$tau_ik[,k])
    }
    
    stats = append(stats, list(smoothed = smoothed))
    stats = append(stats, list(mean_curves = mean_curves))
    #stats.mean_gamma_ijk = mean_gamma_ijk;
    
    stats = append(stats, list(cputime = mean(cputime_total)))
    
    # # Segmentation of each cluster using the MAP rule
    # for k=1:K
    #     [segments_k Zjk] = MAP(stats.mean_gamma_ijk(:,:,k));#MAP segmentation of each cluster of sequences
    #     stats.segments(:,k) = segments_k;
    # end
    
    nu = length(stats$Psi)
    # BIC AIC et ICL*
    stats = append(stats, list(BIC = stats$loglik - (nu*log(n)/2)))#n*m/2!
    stats = append(stats, list(AIC = stats$loglik - nu))
    # ICL*
    # Compute the comp-log-lik
    cik_log_w_k_fyi = (stats$Zik)*(stats$log_w_k_fyi)
    comp_loglik = sum(sum(cik_log_w_k_fyi,2))
    stats = append(stats, list(ICL1 = comp_loglik - nu*log(n)/2)) #n*m/2!
    mixFHMM=list(param=param,stats=stats)
  }
  return(mixFHMM)
}


