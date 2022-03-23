##############----------------------------------------------------------------------------------------------
#### initialization of the DP Gaussian hyperparameters (independent Gaussians across dimensions)
# inputs: 
  ## X: NxD data matrix
  ## M: number of initial components
# output: a list of:
  ## phi_0: an 1xM vector for the stick-breaking point Beta parameter 
  ## alpha_0/beta_0: DxM Gamma hyperparameter matrices for the variance of the Gaussians
  ## m_0/s2_0: DxM Normal hyperparameter matrices for the mean of the Gaussians

init.hyperparameters <- function(X, M) 
{
  
  X <- as.matrix(X)
  N <- dim(X)[1]
  D <- dim(X)[2]
  
  phi_0 <- rep(0.01,M)
  alpha_0 <- matrix(data=0.1, nrow=D, ncol=M)
  beta_0 <- matrix(data=0.1, nrow=D, ncol=M)
  m_0 <- matrix(data=0, nrow=D, ncol=M)
  s2_0 <- matrix(data=1, nrow=D, ncol=M)
  
  retList <- list("phi_0"=phi_0, "alpha_0"=alpha_0, "beta_0"=beta_0, "m_0"=m_0, "s2_0"=s2_0)
  
  return(retList)
}
##############----------------------------------------------------------------------------------------------

#### save the initial hypeparameters for the specific dataset (to be passed into the main VB function later)

# fix number of initial components
M <- 20

# run the hyperparameter function
init.h <- init.hyperparameters(data, M)

# the saved items
phi_0 <- init.h$phi_0
alpha_0 <- init.h$alpha_0
beta_0 <- init.h$beta_0
m_0 <- init.h$m_0
s2_0 <- init.h$s2_0


##############----------------------------------------------------------------------------------------------
#### initialization of the DP Beta variational parameters
# inputs: 
  ## X: NxD data matrix 
  ## M: number of initial components
  ## phi_0: an 1xM vector for the stick-breaking point Beta parameter 
  ## alpha_0/beta_0: DxM Gamma hyperparameter matrices for the variance of the Gaussians
  ## m_0/s2_0: DxM Normal hyperparameter matrices for the mean of the Gaussians
# output: a list of:
  ## alpha/beta: DxM Gamma initial variational matrices for the variance of the Gaussians
  ## m/s2: DxM Gamma initial variational matrices for the mean of the Gaussians
  ## p/q: 1xM initial variational vectors for the stick-breaking point Beta parameter 
  ## r: NxD initial variational matrix for the latent allocation z

init.parameters <- function(X, M, phi_0, alpha_0, beta_0, m_0, s2_0)
{
  
  X <- as.matrix(X)
  N <- dim(X)[1]
  D <- dim(X)[2]
  
  # the variational Gamma parameters of the Beta parameters take the hyperparameter values
  alpha <- alpha_0
  beta <- beta_0
  s2 <- s2_0
  
  ###### set seed for reproduction of same Kmeans
  set.seed(1234567)
  
  # Kmeans algorithm 
  kmeans.out <- kmeans(X, M, nstart=20, iter.max="20")
  kmeans.clusters <- kmeans.out$cluster
  
  # the initial mean vectors by Kmeans
  m <- kmeans.out$centers
  m <- t(m)
  
  # the variational parameter of the latent allocation z
  r <- matrix(data=0, nrow=N, ncol=M)
  
  for (n in 1:N) 
  {
    # Kmeans results are used to compute the initial r matrix
    r[n, kmeans.clusters[n]] <- 1
  }
  
  Ns <- colSums(r, na.rm=TRUE)
  
  # for the stick-breaking point variational parameters their closed form equations are used 
  # along with the Kmeans results (through Ns) to calculate their initial values
  p <- 1 + Ns
  q <- phi_0 + rev(cumsum(rev(Ns))) - Ns
  
  retList <- list("alpha"=alpha, "beta"=beta, "m"=m, "s2"=s2, "p"=p, "q"=q, "r"=r)
  
  return(retList)
}
##############----------------------------------------------------------------------------------------------

#### save the initial variational parameters for the specific dataset (to be passed into the main VB function later)

init.p <- init.parameters(data, M, phi_0, alpha_0, beta_0, m_0, s2_0)

## the saved items
alpha <- init.p$alpha
beta <- init.p$beta
s2 <- init.p$s2 

# Kmeans results
p <- init.p$p; q <- init.p$q
m <- init.p$m
