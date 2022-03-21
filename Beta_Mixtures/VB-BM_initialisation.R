##############----------------------------------------------------------------------------------------------
#### initialization of the Beta hyperparameters
# inputs: 
  ## X: NxD data matrix
  ## M: number of initial components
# output: a list of:
  ## a_0: 1xM Dirichlet hyperparameter matrices for the weights 
  ## alpha_0/beta_0: DxM Gamma hyperparameter matrices for the first Beta parameter
  ## mu_0/eta_0: DxM Gamma hyperparameter matrices for the second Beta parameter
init.hyperparameters <- function(X, M) 
{
  
  X <- as.matrix(X)
  N <- dim(X)[1]
  D <- dim(X)[2]
  
  a_0 <- rep(1/M,M)
  alpha_0 <- matrix(data=1, nrow=D, ncol=M)
  beta_0 <- matrix(data=0.01, nrow=D, ncol=M)
  mu_0 <- alpha_0
  eta_0 <- beta_0
  
  retList <- list("a_0"=a_0, "alpha_0"=alpha_0, "beta_0"=beta_0, "mu_0"=mu_0, "eta_0"=eta_0)
  
  return(retList)
}
##############----------------------------------------------------------------------------------------------
#### save the initial hypeparameters for the specific dataset (to be passed into the main VB function later)
# fix number of initial components
M <- 10
# run the hyperparamter function
init.h <- init.hyperparameters(data, M)
# the saved items
a_0 <- init.h$a_0
alpha_0 <- init.h$alpha_0
beta_0 <- init.h$beta_0
mu_0 <- init.h$mu_0
eta_0 <- init.h$eta_0
##############----------------------------------------------------------------------------------------------
#### initialization of the Beta variational parameters
# inputs: 
  ## X: NxD data matrix 
  ## M: number of initial components
  ## a_0: 1xM Dirichlet hyperparameter matrices for the weights 
  ## alpha_0/beta_0: DxM Gamma hyperparameter matrices for the first Beta parameter (from the function above)
  ## mu_0/eta_0: DxM Gamma hyperparameter matrices for the second Beta parameter (from the function above)
# output: a list of:
  ## alpha/beta: DxM Gamma initial variational matrices for the first Beta parameter
  ## mu/eta: DxM Gamma initial variational matrices for the second Beta parameter
  ## a: 1xM Dirichlet variational matrix for the weights   
  ## r: NxD initial variational matrix for the latent allocation z
init.parameters <- function(X, M, a_0, alpha_0, beta_0, mu_0, eta_0)
{
  X <- as.matrix(X)
  N <- dim(X)[1]
  D <- dim(X)[2]
  # the variational Gamma parameters of the Beta parameters take the hyperparameter values
  alpha <- alpha_0
  beta <- beta_0
  mu <- mu_0
  eta <- eta_0
  ###### set seed for reproduction of same Kmeans
  set.seed(1234567)
  # Kmeans algorithm 
  kmeans.out <- kmeans(X, M, nstart=50, iter.max="20")
  kmeans.clusters <- kmeans.out$cluster
  # the variational parameter of the latent allocation z
  r <- matrix(data=0, nrow=N, ncol=M)
  for (n in 1:N) 
  {
    # Kmeans results are used to compute the initial r matrix
    r[n, kmeans.clusters[n]] <- 1
  }
  
  Ns <- colSums(r, na.rm=TRUE)
  a <- a_0 + Ns
  retList <- list("alpha"=alpha, "beta"=beta, "mu"=mu, "eta"=eta, "a"=a, "r"=r)
  return(retList)
}
##############----------------------------------------------------------------------------------------------
#### save the initial variational parameters for the specific dataset (to be passed into the main VB function later)
init.p <- init.parameters(data, M, a_0, alpha_0, beta_0, mu_0, eta_0)
## the saved items
alpha <- init.p$alpha
beta <- init.p$beta
mu <- init.p$mu
eta <- init.p$eta
# Kmeans results
a <- init.p$a