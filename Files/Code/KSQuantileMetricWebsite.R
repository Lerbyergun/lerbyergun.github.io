# This script contains the KS-quantile metric calculations and outputs the optimal k.
# The method is discussion in Danielsson, Ergun and De vries "Tail Index Estimation: Quantile-Driven Threshold Selection"

HillInv=function(Xs,TN){
  Xs=sort(Xs, decreasing = T)
  TN=min(TN,(sum(Xs>0)-1))
  LogXs=log(Xs[1:(TN+1)])
  Alpha = vector()
  Alpha[1:TN]=1/(cumsum(LogXs[1:TN])/(seq(1:TN))-LogXs[2:(TN+1)])
  return(Alpha)
}

#KS-Distance Metric
#Data is vector of data
#TN is the number of observations you consider a possible search area/ You could use the number of positive observations
KHillKSDistance = function(Data, TN){
  
  #order the data in desending order
  Data=Data[order(Data, decreasing = T)]
  
  vTN=1:TN
  
  Alpha=HillInv(Data,TN)
  
  # create distance matrix the rows show the distance for each quantile
  # where the columns give the distances for a differently chosen threshold
  a=matrix(NA,nrow=TN,ncol=TN)
  for (j in 1:TN){
    a[1:TN,j]=abs(((j/vTN)^(1/Alpha[j])*Data[j+1])-Data[1:TN])
  }
  
  Max_a=matrix(NA,nrow=dim(a)[2],ncol=1)
  
  # maximum distance per column/k
  for (j in 1:dim(a)[2])  Max_a[j,1]=max(a[,j])
  #The k with the smallest maximum distance.
  K_KS=which.min(Max_a)
  
  return(K_KS)
}
