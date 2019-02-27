#=============================================================================================
# R code for computing functional diversity in the paper 
# Chiu, C.-H. and Chao, A. (2014). Distance-based functional diversity measures 
# and their decomposition: a framework based on Hill numbers.PLoS ONE 
# 9(7):e100014. 
# This code includes three parts:
# (1) Empirical functional diversity 
# (2) Empirical functional dissimilarity
# (3) Example
#
# Note: The package "Rcpp" must be installed and loaded before running the scripts
#=============================================================================================
library(Rcpp)

#===============================================
#
# (1). Empirical functional diversity 
#
#===============================================
#' FD_MLE(q, data, Dij) computes FD of order q based on abundance data.
#' @param q a numeric or a vector of diversity order. The suggested range for q is [0, 3].
#' @param data a vector of species sample frequencies.
#' @param Dij a matrix of distance matrix.
#' @return a numerical vector of FD of orderq .
FD_MLE = function(q, data ,Dij){
  Xi <- data[data!=0]
  distance <- Dij[data!=0, data!=0]
  a <- Xi/sum(Xi)
  
  Q = sum(distance*(a %*% t(a)))
  
  Emp <- function(q){
    if(q==1){
      Empirical = exp(-sum(empirical_q1(a, as.matrix(distance), Q)))
    }else{
      Empirical = sum(empirical(a, as.matrix(distance), q, Q))^(1/(1-q))
    }
    Empirical
  }
  sapply(q, Emp)
}


#===============================================
# (2). Empirical functional dissimilarity
#===============================================
#' FD_Beta(q, data, dij, CU, method) is a function which computes functional dissimilarity measure of order q.
#' @param q a numeric or a vector of diversity order; The suggested range for q is [0, 3].
#' @param data a S*N matrix of species sample frequencies with S species (rows), N communities (columns).
#' @param dij a matrix of distance matrix.
#' @param CU a character to choose method, "C" for 1-CqN ; "U" for 1-UqN.
#' @param method a character to choose method, "relative" or "absolute".
#' @return a numerical vector of functional dissimilarity of order q.
FD_Beta = function(q, data, dij, CU, method){
  if(method == "relative") data <- sapply(1:ncol(data), FUN = function(i) data[,i]/sum(data[,i]))
  N = ncol(data)
  S = nrow(dij)
  m = dij
  
  M=matrix(rep(t(m), N), ncol = ncol(m), byrow = TRUE)
  M=matrix(rep(M, N), nrow = nrow(M), byrow = F)
  
  alpha = (1/(N)^2)*FD_MLE(q, unlist(data), M)
  gamma = FD_MLE(q, rowSums(data), dij)
  beta = gamma/alpha
  beta[beta<1] = 1
  
  MAX.dis = matrix(1, N*S, N*S)
  for(i in 1:N){
    MAX.dis[(1:S)+((i-1)*S), (1:S)+((i-1)*S)] = dij
  }
  
  MAX.gamma =  FD_MLE(q, unlist(data), MAX.dis)
  MAX.beta = MAX.gamma/alpha
  #MAX.alpha =  FD_MLE(q, unlist(data), MAX.dis)
  #MAX.beta = gamma/MAX.alpha
  if(CU=="C"){
    out = (1-beta^(1-q))/(1-MAX.beta^(1-q))
  }else{
    out = (1-(beta)^(q-1))/(1-MAX.beta^(q-1))
  }
  out[q==1] = log(beta[q==1])/log(MAX.beta[q==1])
  out
}


#=========================================================================
#The following two cpp functions are used in the main function FD_MLE.
#=========================================================================
cppFunction(
  "NumericMatrix empirical(NumericVector ai,NumericMatrix dij,float q,float Q){
  const int S = ai.size();
  NumericMatrix temp(S,S);
  for(int i=0;i<S;i++){
  for(int j = 0;j<S;j++){
  temp(i,j) = dij(i,j)*pow((ai[i]*ai[j]/Q),q);
  }
  }
  return(temp);
  }")
cppFunction(
  "NumericMatrix empirical_q1(NumericVector ai,NumericMatrix dij,float Q){
  const int S = ai.size();
  NumericMatrix temp(S,S);
  for(int i=0;i<S;i++){
  for(int j = 0;j<S;j++){
  temp(i,j) = dij(i,j)*(ai[i]*ai[j]/Q)*log(ai[i]*ai[j]/Q);
  }
  }
  return(temp);
  }")
#===============================================
#
# (3). Example 
#
#===============================================
data <- read.table("Dunes_data.txt")
dij <- read.table("Dunes_distance.txt")
dij <- as.matrix(dij)
FD_MLE(q = c(0,1,2), data = data[,1], Dij = dij)
FD_MLE(q = c(0,1,2), data = data[,2], Dij = dij)
FD_MLE(q = c(0,1,2), data = data[,3], Dij = dij)
FD_Beta(q = c(0,1,2), data = data, dij = dij, CU = "U", method = "relative")
