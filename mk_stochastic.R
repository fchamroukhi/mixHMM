mk_stochastic <- function(T){
  # MK_STOCHASTIC Ensure the argument is a stochastic matrix, i.e., the sum over the last dimension is 1.
  # T = mk_stochastic(T)
  #
  # If T is a vector, it will sum to 1.
  # If T is a matrix, each row will sum to 1.
  # If T is a 3D array, then sum_k T(i,j,k) = 1 for all i,j.
  
  # Set zeros to 1 before dividing
  # This is valid since S(j) = 0 iff T(i,j) = 0 for all j
#################################################################################################################  
  if ((length(dim(T))==2) & (nrow(T)==1 || ncol(T)==1)){ # isvector
    T = normalize(T)
  }else if (is.matrix(T)==TRUE){ # matrix
    S = apply(T,1,sum)
    S = S + (S==0)
    norm=matrix(rep(S,ncol(T)),nrow=length(S))
    T = T / norm
  }else{  # multi-dimensional array
    ns = dim(T)
    T = matrix(T, prod(ns[1]:ns[length(ns)-1]), ns[length(ns)])
    S = apply(T,1,sum)
    S = S + (S==0)
    norm = matrix(rep(S,ncol(T)),nrow=length(S),ncol=ns[length(ns)])
    T = T / norm
    T =array(T, ns)}
  return(T)
}
