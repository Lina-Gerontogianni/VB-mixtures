#####----------------- Discriminative measures by Lin Lin (2016)-----------------#####
#### Beta Mixture (independent across dimensions)
## Beta function for each component
# input: 
  ## x: the data matrix 
  ## m: index of the specific component
  ## u: DxM matrix of the first variational parameter of the Beta component
  ## v: DxM matrix of the second variational parameter of the Beta component
# output: 
  ## Nx1 matrix of the Beta function values for the particular component
f <- function(x, m, u, v)
{
  f_d <- matrix(NA, nrow=dim(x)[1], ncol=dim(x)[2])
  fprod <- matrix(NA, nrow=dim(x)[1], ncol=1)
  for(n in 1:dim(x)[1])
  {
    for(d in 1:dim(x)[2])
    {
      f_d[n,d] <- dbeta(x[n,d], u[d,m], v[d,m])
    }
    fprod[n, ] <- prod(f_d[n, ])
  }
return(fprod)
}
##############----------------------------------------------------------------------------------------------
## Beta function for all components
# input: 
  ## x: the data matrix 
  ## M: the number of components
  ## u: DxM matrix of the first variational parameter of the Beta component
  ## v: DxM matrix of the second variational parameter of the Beta component
# output: 
  ## NxM matrix of the Beta function values 
f_all <- function(x, M, u, v)
{
  Func <- matrix(NA, nrow=dim(x)[1], ncol=M)
  for (j in 1:M)
  {
    Func[ ,j] <- f(x, j, u, v)
  }
  return(Func)
}
##############----------------------------------------------------------------------------------------------
# function of Ac(h) (here c in Lin's paper is m)
  # input: 
    ## x: data matrix
    ## m: the index of the specific component for which A will be computed
    ## h: the vector of dimensions (variables) to be tested
    ## M: number of components
    ## pi_M: variational vector of the mixing coefficients
  # output:
    ## Am(h) value (A value for the mth component on the h variable vector )
A <- function(x, m, h, M, pi_M, u, v)
{
  h<-sort(h)
  x <- as.matrix(x[ ,h])
  # mixture function of the form: p1*f1 + p2*f2 + ...
  if(length(pi_M)==1){
    g<-as.matrix(rowSums(f_all(x, M, u, v) %*% pi_M))
  } else{
  g <- as.matrix(rowSums(f_all(x, M, u, v) %*% diag(pi_M)))}
  # numerator of d1 (Lin's thesis)
  delta <- sum(f(x, m, u, v)*(g - pi_M[m]*f(x,m,u,v)))
  # denominator of d1 (Lin's thesis)
  Delta <- sum(f(x,m,u,v)*g)
  # d1 measure
  d <- delta/Delta
  # t+ measure (Lin's paper)
  tauplus <- 1 - d
  # integral contained in the calculation of t-
  intg2 <- sum(g^2)
  # t- measure (Lin's paper)
  tauminus <- pi_M[m]*delta / (intg2 - pi_M[m]*Delta)
  # Am(h)
  return (pi_M[m]*tauplus + (1 - pi_M[m])*(1-tauminus) )
}
##############----------------------------------------------------------------------------------------------

#### examples for the Am(h) calculation
#### posterior means of u and v in Beta densities (note that u and v follow gamma distributions)
# model is the VB BetaMix algorithm (either AVB/VB or Finite/DP)
E.z <- model_dpbeta$r
clusters <- max.col(E.z, ties.method="first")

# the variational weights to be passed into the Am(h) function
#pi_M <- model_dpbeta$pi_M[model_dpbeta$pi_M > 0.005]
pi_M<-sort(table(clusters)/dim(data)[1])

uni<-as.numeric(names(pi_M))

# posterior means of u and v to be passed into the f function included into the Am(h)
u <- model_dpbeta$alpha[ ,uni] / model_dpbeta$beta[ ,uni]
v <- model_dpbeta$mu[ ,uni] / model_dpbeta$eta[ ,uni]
