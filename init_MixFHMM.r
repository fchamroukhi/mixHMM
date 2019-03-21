init_MixFHMM<-function(Y, K, R, variance_type, order_constraint, init_kmeans, try_algo){
  
  #
  #
  #
  #
  #
  #
  ###################### FC ##############
  source("myKmeans.R")
  n=nrow(Y)
  m=ncol(Y)
  
  # # 1. Initialization of cluster weights
  w_k=1/K*matrix(c(1),K,1)
  # Initialization of the model parameters for each cluster
  if (init_kmeans==1){
    max_iter_kmeans = 400
    n_tries_kmeans = 20
    verbose_kmeans = 0
    solution = myKmeans(Y, K, n_tries_kmeans, max_iter_kmeans,verbose_kmeans)
    pi_k=matrix(c(0),nrow=R,ncol=K)
    A_k=array(c(0),dim=c(R,K,K))
    mu_kr=matrix(c(0),nrow=R,ncol=K)
    sigma_k=c()
    sigma_kr=matrix(c(0),nrow=K,ncol=K)
    for (k in 1:K){
      Yk = Y[solution$klas==k ,] #if kmeans
      
      param =  init_gauss_hmm(Yk, R, order_constraint, variance_type, try_algo)
      
      # 2. Initialisation de \pi_k
      pi_k[,k] = t(param$initial_prob)#[1;zeros(R-1,1)];
      # 3. Initialisation de la matrice des transitions
      A_k[,,k]  =  param$trans_mat$M
      
      if (order_constraint==1){
        mask = param$mask
      }
      # 4. Initialisation des moyennes
      mu_kr[,k] = param$mur
      if (variance_type=='common'){
        sigma_k[k] = param$sigma
      }
      
      else{
        sigma_kr[,k] = param$sigma2r
      }
    }
  }
  else{
    ind = sample(1:n,n)
    pi_k=NULL
    A_k=NULL
    for (k in 1:K){
      if (k<K){
        Yk = Y[ind[(k-1)*round(n/K) +1 : k*round(n/K)],]
      }
      else{
        Yk = Y[ind[(k-1)*round(n/K) +1 : n],]
      }
      
      param =  init_gauss_hmm(Yk, R, order_constraint, variance_type, try_algo)
      # 2. Initialisation de \pi_k
      pi_k[,k] = param$initial_prob#[1;zeros(R-1,1)];
      
      # 3. Initialisation de la matrice des transitions
      A_k[,,k]  =  param$trans_mat 
      
      if (order_constraint==1){
        mask = param$mask
      }
      
      # 4. Initialisation des moyennes
      mu_kr[,k] = param$mur
      if (variance_type=='common'){
        sigma_k[k] = param$sigma
      }
      else{
        sigma_kr[,k] = param$sigma2r
      }
    }
  }
  param=list(w_k=w_k,pi_k=pi_k,A_k=A_k,mu_kr=mu_kr,sigma_k=sigma_k,sigma_kr=sigma_kr,mask=mask)
  
}

##################################################################################################################################################
init_gauss_hmm<-function(Y, R, order_constraint, variance_type, try_EM){
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
  
  param=list()
  if (order_constraint==1){
    # # Tnitialisation en tenant compte de la contrainte:
    
    # Initialisation de la matrice des transitions
    mask = diag(R)#mask d'ordre 1
    for (r in 1:R-1){
      ind = which(mask[r,] != 0)
      mask[r,ind+1] = 1
    }
    # Initialisation de la loi initiale de la variable cachee
    param=append(param,list(initial_prob = c(1,matrix(c(0),R-1,1))))
    param=append(param,list(trans_mat = normalize(mask,2)))
    param=append(param,list(mask = mask))
  }
  else{
    # Initialisation de la loi initiale de la variable cachee
    param=append(param,list(initial_prob = 1/R*matrix(c(1),R,1)))
    param=append(param,list(trans_mat = mk_stochastic(matrix(c(runif(R)),R,R))))   ####### PAS SUR
  }
  
  #  Initialisation des moyennes et des variances.
  param_gauss = init_gauss_param_hmm(Y, R, variance_type, try_EM)
  
  param = append(param,list(mur = param_gauss$mur))
  if (variance_type=='common'){
    param = append(param,list(sigma = param_gauss$sigma))
  }
  else{
    param = append(param,list(sigma2r = param_gauss$sigma2r))
  }
  return(param)
}

###################################################################################
init_gauss_param_hmm<-function(Y, R, variance_type, try_EM){
  
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
  
  n=nrow(Y)
  m=ncol(Y)
  mur=NULL
  sigma2=NULL
  sigma2r=NULL
  
  if (variance_type=='common'){
    s=0
  }
  
  if (try_EM==1){
    ##############################
    #decoupage de l'echantillon (signal) en K segments
    zi = round(m/R)-1
    for (r in 1:R){
      i = (r-1)*zi+1
      j = r*zi
      
      Yij = Y [,i:j]
      Yij = matrix(t(Yij),1,byrow=T)
      mur[r] = mean(Yij)
      
      if (variance_type=='common'){
        s=s+ sum((Yij-mur[r])^2)
        sigma = s/(n*m)
      }
      
      else{
        m_r = j-i+1 
        sigma2r[r] = sum((Yij-mur[r])^2)/(n*m_r)
      }
    }
    
  }
  else {
    # initialisation aléatoire
    Lmin= 2#round(m/(K+1));#nbr pts min dans un segments
    tr_init = matrix(c(0),1,R+1)
    tr_init[1] = 0
    R_1=R
    for (r in 2:R){
      R_1 = R_1-1
      temp = seq(tr_init[r-1]+Lmin,m-R_1*Lmin)
      ind = sample(1:length(temp),length(temp))
      tr_init[r]= temp[ind[1]]
    }
    
    tr_init[R+1] = m
    for (r in 1:R) {
      i = tr_init[r]+1
      j = tr_init[r+1]
      Yij = Y[,i:j]
      Yij = matrix(t(Yij),ncol=1,byrow=T) 
      
      mur(r) = mean(Yij)
      
      if(variance_type=='common'){
        s=s+ sum((Yij-mur[r])^2)
        sigma = s/(n*m)
      }
      
      else{
        m_r = j-i+1 
        sigma2r[r] = sum((Yij-mur[r])^2)/(n*m_r)
      }
      
    }
    
  }
  param_gauss=list(mur=mur,sigma2=sigma2,sigma2r=sigma2r) 
  return(param_gauss)
}             