##############----------------------------------------------------------------------------------------------
#### main AVB function for the DP Gaussian Mixture Model (multivariate) + extra helpful functions
##############----------------------------------------------------------------------------------------------

##############----------------------------------------------------------------------------------------------
#### safe computation of logsumexp (included in the main function)
# inputs: 
  ## y: scalar
# output: 
  ## lse: safe log sum exp value
library(matrixcalc)

logsumexp <- function(y) 
  {
  
  # Computes log(sum(exp(x))
  a <- max(y)
  lse <- log(sum(exp(y-a))) + a
  j <- which(!is.finite(lse))
  if (length(j) > 0) {lse[j] <- a}
  
  return(lse)
  }
##############----------------------------------------------------------------------------------------------

##############----------------------------------------------------------------------------------------------
#### main AVB function for DP GaussianMix (multivariate)
# inputs: 
  ## X: NxD data matrix
  ## M: number of initial components
  ## a: 1xM Dirichlet initial variational matrix for the weights 
  ## theta: DxM Gaussian initial variational matrix for the mean of the Gaussians 
  ## beta: 1XM Gaussian initial variational vector related to the covariance matrices of the Gaussians (see Bishop 2006)
  ## W: DxDxM Wishart initial variational array for the covariance matrices of the Gaussians
  ## nu: 1xM Wishart initial variational vector for the covariance matrices of the Gaussians 
  ## phi_0: an 1xM vector for the stick-breaking point Beta parameter 
  ## theta_0: DxM Gaussian hyperparameter matrix for the mean of the Gaussians 
  ## beta_0: 1XM Gaussian hyperparameter vector related to the covariance matrices of the Gaussians (see Bishop 2006)
  ## W_0: DxDxM Wishart hyperparameter array for the covariance matrices of the Gaussians
  ## nu_0: 1xM Wishart hyperparameter vector for the covariance matrices of the Gaussians
  ## T: the temperature vector for the annealing part (pre-define)
  ## max_iterations: maximum number of VB iterations
  ## epsilon: threshold to achieve convergence
# output: 
  ## p/q: 1xM variational vectors for the stick-breaking point Beta parameter
  ## theta: DxM Gaussian variational matrix for the mean of the Gaussians 
  ## beta: 1XM Gaussian variational vector related to the covariance matrices of the Gaussians (see Bishop 2006)
  ## W: DxDxM Wishart variational array for the covariance matrices of the Gaussians
  ## nu: 1xM Wishart variational vector for the covariance matrices of the Gaussians 
  ## r: NxD variational matrix for the latent allocation z
  ## L: ELBO values
  ## printL: print the ELBO values and difference to the previous one

avb_dpgm_mult <- function(X, M, p, q, theta, beta, W, nu, phi_0, theta_0, 
                          beta_0, W_0, nu_0, Temp, max_iterations=iter, epsilon=1e-4, printL=FALSE)
{
  # define the dimensions according to the dataset
  X <- as.matrix(X)
  D <- ncol(X)
  N <- nrow(X)
  Tinv<-1/Temp
  # ELBO storage variable
  L <- rep(-Inf, max_iterations)
  
  # weight storage variable
  pi <- rep(NA, M)
  
  # create the expectation found in r
  E.log_detLambda <- rep(NA, M)
  
  # create another expectation of r
  E.log_x_mu_Lambda <- matrix(NA, nrow=N, ncol=M)
  
  # create variables for terms included into the variationbal parameters 
  xbar <- matrix(NA, nrow=D, ncol=M)
  S <- array(NA, c(D, D, M))
  log_rho_lambda <- matrix(NA, nrow = N, ncol = M)
  
  ### VB scheme
  for (i in 2:max_iterations) 
  {
    
    # calculation of expectations contained into the variational parameter r of the allocation z
    E.loglambda <- digamma(p) - digamma(p + q)
    E.log1_lambda <-  digamma(q) - digamma(p + q)
    
    E.log1_lambdaj <- c(0, cumsum(E.log1_lambda)[1:(M-1)])
    N.E.log1_lambdaj <- matrix(E.log1_lambdaj, nrow=N, ncol=M, byrow=TRUE)
    N.E.loglambda <- matrix(E.loglambda, nrow=N, ncol=M, byrow=TRUE)
    
    for(m in 1:M)
    {
      E.log_detLambda[m] <- sum(digamma((nu[m] + 1-c(1:D))/2)) + D*log(2) + log(det(W[ , ,m]))

      x_mu <- sweep(X, MARGIN = 2, STATS = theta[, m], FUN = "-")
      E.log_x_mu_Lambda[ ,m] <- D*(beta[m])^(-1) + nu[m] * diag(x_mu %*% W[ , ,m] %*% t(x_mu))
      
      # the variational equation of log-rho
      log_rho_lambda[ ,m] <- - (D/2)*log(2*355/113) + (1/2)*E.log_detLambda[m] - (1/2)*E.log_x_mu_Lambda[ ,m]
    }
    
    log_rho <- log_rho_lambda + N.E.loglambda + N.E.log1_lambdaj

    # the usefulness of logsumexp function
    Q <- apply(log_rho, 1, logsumexp)
    log_r <- log_rho - Q
    
    # the variational parameter r of z
    r <- apply(log_r, 2, exp)
    # trick to avoid zero values
    r <- (r + 10^-9)^Tinv[i]

    # useful expectations wrt r (to be used for the computation of the rest of the variational parameters)
    Ns <- colSums(r, na.rm=TRUE)
    
    for (m in 1:M) 
    {
      xbar[, m] <- (1/ Ns[m]) * (r[ ,m] %*% X)    
      x_xbar <- sweep(X, MARGIN = 2, STATS = xbar[, m], FUN = "-")
      S[, , m]   <- (1/ Ns[m]) * t(x_xbar) %*% (x_xbar * r[, m])  
    }
    
    # Normal-Wishart variational parameters for the mean vectors and covariance matrices
    nu <- nu_0 + Tinv[i]*Ns +1
    beta <- beta_0 + Tinv[i]*Ns
    
    for(m in 1:M)
    {
      theta[ ,m] <- (1/beta[m]) * (beta_0[m] *theta_0[ ,m] + Tinv[i]*Ns[m] * xbar[ ,m])
      W[ , ,m] <- solve(W_0[ , , m]) + Tinv[i]*(Ns[m]*S[ , ,m] + ((beta_0[m]*Ns[m])/(beta_0[m] + Ns[m])) * tcrossprod((xbar[, m] - theta_0[ ,m])))
      W[ , ,m] <- solve(W[ , ,m])
    }
    
    # the stick-breaking point variational parameters
    p <- 1 + Tinv[i] * Ns
    q <- phi_0 + Tinv[i] * (rev(cumsum(rev(Ns))) - Ns)
    
    # the stick-breaking point parameter
    lambda <- head(p, M-1) / (head(p, M-1) + head(q, M-1))
    lambda <- c(lambda, 1)
    
    # the variational weights after the stick-breaking point computation
    for (k in 1:M)
    {
      pi[k] <- lambda[k] * prod(head(1-lambda, k-1))
    }
    
    # update the expectations contained into the ELBO
    E.loglambda <- digamma(p) - digamma(p + q)
    E.log1_lambda <-  digamma(q) - digamma(p + q)
    E.log1_lambdaj <- c(0, cumsum(E.log1_lambda)[1:(M-1)])
    N.E.log1_lambdaj <- matrix(E.log1_lambdaj, nrow=N, ncol=M, byrow=TRUE)
    N.E.loglambda <- matrix(E.loglambda, nrow=N, ncol=M, byrow=TRUE)
    
    for(m in 1:M)
    {
      E.log_detLambda[m] <- sum(digamma((nu[m] + 1-c(1:D))/2)) + D*log(2) + log(det(W[ , ,m]))
    }
    
    # functions included into the ELBO terms
    # Wishart constant
    logBeta <- function(B,v)
    {
      return(-(v/2) * log(det(B)) - (v*D/2) * log(2) - (D*(D-1)/4)*log(355/113) - sum(lgamma((v +1 -(1:D))/2)) )
    }
    
    # Entropy
    H <- function(B,v,R)
    {
      return(-logBeta(B,v) - ((v-D-1)/2) * R + (v*D/2))
    }
    
    ## ELBO 
    # each function into the ELBO has been calculated individually (l1,l2,l3,l4,l5,l6,and l7) and finally all are summed up
    l5_m<-l6_m<-l7_m<-rep(NA, M)
    
    l1 <- sum(Tinv[i]*r * (N.E.loglambda + E.log1_lambdaj ))
    l2 <- sum(log( phi_0) + (phi_0 - 1) * E.log1_lambda )
    l3 <- sum(r * log_r)
      
    log_g.p.q <- lgamma( p + q)
    log_p <- lgamma(p)
    log_q <- lgamma(q)
    
    l4_1 <-  - sum(log_g.p.q) + sum(log_p) + sum(log_q)
    l4_2 <-  - sum((p -1) * E.loglambda ) - sum((q -1) * E.log1_lambda)
    l4 <- l4_1 + l4_2
    
    for( m in 1:M)
    {
      l5_m[m] <- (1/2) * Tinv[i]*Ns[m] * ( E.log_detLambda[m] - D/beta[m] - nu[m] * matrix.trace(S[ , ,m] %*% W[ , ,m]) 
                                           - nu[m] * t(xbar[ ,m] - theta[ ,m]) %*% W[ , ,m] %*% (xbar[ ,m] - theta[ ,m]) - D *log(2*355/113) )
      
      l6_m[m] <- (1/2) * (D *log(beta_0[m]/(2*355/113)) + E.log_detLambda[m] 
                          - (D * beta_0[m])/beta[m] - beta_0[m] * nu[m] * t(theta[ ,m] - theta_0[ ,m]) %*% W[ , ,m] %*% (theta[ ,m] - theta_0[ ,m]))
      + logBeta(W_0[ , ,m], nu_0[m]) + ((nu_0[m]-D-1)/2) * E.log_detLambda[m] -(1/2) * nu[m] * matrix.trace(solve(W_0[ , ,m]) %*% W[ , ,m])
      
      l7_m[m] <-  (1/2) * E.log_detLambda[m] + (D/2) * log(beta[m]/(2*355/113)) - (D/2) - H(W[ , ,m], nu[m], E.log_detLambda[m])
    }
    
    l5 <- sum(l5_m)
    l6 <- sum(l6_m)
    l7 <- sum(l7_m)
    
    #print(c(l1,l2,-l3,-l4,l5,l6,-l7))
    
    # Total ELBO calculation
    L[i] <- l1 + l2 - l3 - l4 + l5 + l6 - l7
    #print(L[i])
    # print ELBO value and difference with the previous one
    if (printL) { cat("Iter:\t", i, "\tELBO:\t", L[i], "\tELBO_diff:\t", L[i] - L[i-1], "\n")}
    
    # test if ELBO decreases
    if (L[i] < L[i - 1]) { message("Warning: ELBO decreases\n"); }
    
    # test convergence with epsilon threshold
    if (abs(L[i] - L[i - 1]) < epsilon) { break }
    
    # test VB needs more iteration to converge
    if (i == max_iterations) {warning("VB did not converge\n")}
  }
  
  object <- structure(list(p=p, q=q, pi=pi, theta=theta, beta=beta, nu=nu, W=W, r=r, L=L[2:i]))
  
  return(object)
} 
##############----------------------------------------------------------------------------------------------

### example model for data from gmm_mult_generator
# pre-define max iterations for AVB
iter <- 10000

# pre-define T: a vector with K elements and its inverse value Tinv
length <- 100
Temp <- c(rev(seq(1, 10, length.out = length )),rep(1,iter-length ))

model <- avb_dpgm_mult(data, M, p, q, theta, beta, W, nu, phi_0, theta_0, beta_0, W_0, nu_0, Temp, max_iterations=iter)

## ELBO plot

# truncate the initial elbo iterations 
o<-10
Elbo <- data.frame(elbo=model$L[o:length(model$L)], iter=o:length(model$L))

pel <- ggplot(data=Elbo,aes(x=iter,y=elbo)) +xlab("iterations")+ylab("ELBO") + labs(title="ELBO of AVB-multivariate GaussianMix") 
pel <- pel + geom_line(size=1.1,color="royalblue2")  

pel <- pel + theme(
  panel.background = element_rect(colour="gray91", fill="gray91"), axis.text.x=element_text(size=12, face="bold"),
  axis.text.y=element_text(size=12, face="bold"),
  axis.title.x=element_text(size=14,face="bold", colour="slategrey"),
  axis.title.y=element_text(size=16,face="bold", colour = "royalblue"),
  plot.title = element_text(size=14,face="bold") 
)
print(pel)


# final component weights (discard those with weight < 0.001)
var_weights <- model$pi[model$pi>0.007]
cat("Variational weights", var_weights, "\n")
round(var_weights,3) 

