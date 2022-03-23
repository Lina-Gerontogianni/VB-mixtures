##############----------------------------------------------------------------------------------------------
#### initialization of the DP GaussianMix (multivariate) hyperparameters
# inputs: 
  ## X: NxD data matrix
  ## M: number of initial components
# output: a list of:
  ## phi_0: an 1xM vector for the stick-breaking point Beta parameter 
  ## theta_0: DxM Gaussian hyperparameter matrix for the mean of the Gaussians 
  ## beta_0: 1XM Gaussian hyperparameter vector related to the covariance matrices of the Gaussians (see Bishop 2006)
  ## W_0: DxDxM Wishart hyperparameter array for the covariance matrices of the Gaussians
  ## nu_0: 1xM Wishart hyperparameter vector for the covariance matrices of the Gaussians

init.hyperparameters <- function(X, M) 
{
  
  X <- as.matrix(X)
  N <- dim(X)[1]
  D <- dim(X)[2]
  
  phi_0 <- rep(0.01,M)
  theta_0 <- matrix(colMeans(X), nrow=D, ncol=M, byrow = TRUE)
  beta_0 <- rep(1, M)
  W_0 <- array(diag(1.2, D, D), c(D, D, M))
  nu_0 <- rep((D + 20), M)
  
  retList <- list("phi_0"=phi_0, "theta_0"=theta_0, "beta_0"=beta_0, "W_0"=W_0, "nu_0"=nu_0)
  
  return(retList)
}
##############----------------------------------------------------------------------------------------------

#### save the initial hypeparameters for the specific dataset (to be passed into the main VB function later)

# fix number of initial components
M <- 30

# run the hyperparamter function
init.h <- init.hyperparameters(data, M)

# the saved items
phi_0 <- init.h$phi_0
theta_0 <- init.h$theta_0
beta_0 <- init.h$beta_0
W_0 <- init.h$W_0
nu_0 <- init.h$nu_0


##############----------------------------------------------------------------------------------------------
#### initialization of the DP GaussianMix (multivariate) variational parameters
# inputs: 
  ## X: NxD data matrix 
  ## M: number of initial components
  ## phi_0: an 1xM vector for the stick-breaking point Beta parameter 
  ## theta_0: DxM Gaussian hyperparameter matrix for the mean of the Gaussians 
  ## beta_0: 1XM Gaussian hyperparameter vector related to the covariance matrices of the Gaussians (see Bishop 2006)
  ## W_0: DxDxM Wishart hyperparameter array for the covariance matrices of the Gaussians
  ## nu_0: 1xM Wishart hyperparameter vector for the covariance matrices of the Gaussians
# output: a list of:
  ## p/q: 1xM initial variational vectors for the stick-breaking point Beta parameter 
  ## theta: DxM Gaussian initial variational matrix for the mean of the Gaussians 
  ## beta: 1XM Gaussian initial variational vector related to the covariance matrices of the Gaussians (see Bishop 2006)
  ## W: DxDxM Wishart initial variational array for the covariance matrices of the Gaussians
  ## nu: 1xM Wishart initial variational vector for the covariance matrices of the Gaussians 

init.parameters <- function(X, M, phi_0, theta_0, beta_0, W_0, nu_0)
{
  
  X <- as.matrix(X)
  N <- dim(X)[1]
  D <- dim(X)[2]
  
  # the variational Dirichlet and Wishart parameters take the hyperparameter values
  
  nu <- nu_0 +runif(M,0,1)/100
  W <- W_0
  
  # the variational Gaussian-Wishart parameters (kmeans and hyperparameter respectively)
  
  # Kmeans algorithm 
  ###### set seed for reproduction of same Kmeans
  #set.seed(1234)
  kmeans.out <- kmeans(X, M, nstart = 25)
  kmeans.centers <- kmeans.out$centers
  kmeans.clusters <- kmeans.out$cluster
  
  theta <- t(kmeans.centers)
  beta <- beta_0
  
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
  
  retList <- list("p"=p, "q"=q, "theta"=theta, "beta"=beta, "W"=W, "nu"=nu,"Ns"=Ns)
  
  return(retList)
}
##############----------------------------------------------------------------------------------------------

#### save the initial variational parameters for the specific dataset (to be passed into the main VB function later)

init.p <- init.parameters(data, M, phi_0, theta_0, beta_0, W_0, nu_0)

## the saved items
p <- init.p$p
q <- init.p$q
theta <- init.p$theta
beta <- init.p$beta
W <- init.p$W
nu <- init.p$nu
