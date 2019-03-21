show_mixFHMM_results <-function(Y,mixFHMM)
#
#
#
###############################################

n=nrow(Y)
m=ncol(Y)

t=(0:(m-1)) 
x11()  
matplot(t,Y[1,],main="original time series",xlab ="t",ylab="y(t)",type="l",col=1)
for (k in 2:n){
  matlines(t,Y[k,],type="l",col=k)
}

K = length(mixFHMM$param$w_k)
#print(Y[,mixFHMM$stats$klas==1])
x11()
matplot(t,t(Y[mixFHMM$stats$klas==1,]),main="Clustered time series",xlab ="t",ylab="y(t)",col=1,type="l")
for (k in 2:K){
 matlines(t,t(Y[mixFHMM$stats$klas==k,]),col=k,type="l")
}

K = ncol(mixFHMM$stats$tau_ik)
#print(mixFHMM$stats$smoothed)
for (k in 1:K){
  x11()
  matplot(t,t(Y[mixFHMM$stats$klas==k,]),main=paste("Cluster ",k),xlab ="t",ylab="y(t)",type="l",col=k)
  matlines(t,mixFHMM$stats$smoothed[,k],type="l",lwd=2.5,col=(k+1))
}


