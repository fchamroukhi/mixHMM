normalize <- function(A, dim){
# NORMALISE Make the entries of a (multidimensional) array sum to 1
# [M, c] = normalise(A)
# c is the normalizing constant
#
# [M, c] = normalise(A, dim)
# If dim is specified, we normalise the specified dimension only,
# otherwise we normalise the whole array.

if (nargs() < 2){
  z = sum(A)
  # Set any zeros to one before dividing
  # This is valid, since c=0 => all i. A(i)=0 => the answer should be 0/1=0
  s = z + (z==0)
  M = A / s
}else if(dim==1){#normalize each column
  z = apply(A,1,sum)
  s = z + (z==0)
  # M = A ./ (d'*ones(1,size(A,1)))';
  M = A %/% matrix(rep(s,L),nrow=length(s))
}else{
  # Keith Battocchi - v. slow because of repmat
  z=apply(A,1,sum)
  s = z + (z==0)
  L=ncol(A)  
  d=length(dim(A))
  v=matrix(c(1),d,1)
  v[dim]=L
  #c=repmat(s,v)
  c=matrix(rep(s,L),nrow=length(s),ncol=L)
  M=A/c
}

return(list(M=M,z=z))
}


