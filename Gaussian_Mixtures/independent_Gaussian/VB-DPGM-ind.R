##############----------------------------------------------------------------------------------------------
#### main AVB function for the DP Gaussian Mixture Model (independent across dimensions) + extra helpful functions
##############----------------------------------------------------------------------------------------------

##############----------------------------------------------------------------------------------------------
#### safe computation of logsumexp (included in the main function)
# inputs: 
  ## y: scalar
# output: 
  ## lse: safe log sum exp value

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
#### main AVB function for DP GaussianMix (independent across dimensions)
# inputs: 
  ## X: NxD data matrix
  ## M: number of initial components
  ## alpha/beta: DxM Gamma initial variational matrices for the variance of the Gaussians
  ## m/s2: DxM Gamma initial variational matrices for the mean of the Gaussians
  ## p/q: 1xM initial variational vectors for the stick-breaking point Beta parameter 
  ## phi_0: an 1xM vector for the stick-breaking point Beta parameter 
  ## alpha_0/beta_0: DxM Gamma hyperparameter matrices for the variance of the Gaussians
  ## m_0/s2_0: DxM Normal hyperparameter matrices for the mean of the Gaussians
  ## T: the temperature vector for the annealing part (pre-define)
  ## max_iterations: maximum number of VB iterations
  ## epsilon: threshold to achieve convergence
# output: 
  ## alpha/beta: DxM Gamma variational matrices for the variance of the Gaussians
  ## m/s2: DxM Gamma variational matrices for the mean of the Gaussians
  ## p/q: 1xM variational vectors for the stick-breaking point Beta parameter
  ## r: NxD variational matrix for the latent allocation z
  ## L: ELBO values
  ## w: weight values in each iteration (for evolution purposes)
  ## printL: print the ELBO values and difference to the previous one

avb_dpgm <- function(X, M, alpha, beta, m, s2, p, q, alpha_0, beta_0, 
                    m_0, s2_0, phi_0, T, max_iterations=iter, epsilon=1e-4, printL=FALSE)
{
  # define the dimensions according to the dataset
  X <- as.matrix(X)
  D <- ncol(X)
  N <- nrow(X)
  
  # initial objects to receive the variational results
  # initial ELBO
  L <- rep(-Inf, max_iterations)

  # initial weights in each iteration
  w <- matrix(1/M, ncol= M, nrow = max_iterations)
  
  # initial weights in final iteration
  pi <- rep(0, M)
  
    ### AVB scheme
    for (i in 2:max_iterations) 
      {
      
      # calculation of expectations contained into the variational parameters
      E.mu <- m
      V.mu <- s2 
      E.logsigma2 <- log(beta) - digamma(alpha)
      E.inv_sigma2 <- alpha/beta 
      E.loglambda <- digamma(p) - digamma(p + q)
      E.log1_lambda <-  digamma(q) - digamma(p + q)

      # calculation of terms found in the variational equations
      N.E.logsigma2 <- matrix(colSums(E.logsigma2, na.rm=TRUE),  nrow=N, ncol=M, byrow=TRUE)
      N.E.inv_sigma2 <- matrix(colSums(E.inv_sigma2, na.rm=TRUE),  nrow=N, ncol=M, byrow=TRUE)
      E.log1_lambdaj <- c(0, cumsum(E.log1_lambda)[1:(M-1)])
      N.E.log1_lambdaj <- matrix(E.log1_lambdaj, nrow=N, ncol=M, byrow=TRUE)
      N.E.loglambda <- matrix(E.loglambda, nrow=N, ncol=M, byrow=TRUE)
      log.2pi <- matrix(log(2*355/113), nrow=N, ncol=M)
      
      inv_s2.x2 <- (X^(2)) %*% E.inv_sigma2
      inv_s2.x.m <- X %*% (E.inv_sigma2 * E.mu)
      inv_s2.s2 <- matrix(colSums(E.inv_sigma2 * V.mu), ncol=M, nrow=N, byrow = TRUE)
      inv_s2.m2 <- matrix( colSums(E.inv_sigma2 * (E.mu)^(2)), ncol=M, nrow=N, byrow = TRUE)
      
      D.X_m <-  inv_s2.x2 - 2 * inv_s2.x.m + inv_s2.s2 + inv_s2.m2

      log_rho <- N.E.loglambda + N.E.log1_lambdaj - (D/2) * log.2pi - (1/2) * N.E.logsigma2 - (1/2) * D.X_m
      
      # the usefulness of logsumexp function
      S <- apply(log_rho, 1, logsumexp)
      log_r <- log_rho - S
      
      # the variational parameter of z
      r <- apply(log_r, 2, exp)
      # trick to avoid zero values
      r <- (r + 10^-9)^ Tinv[i] 
      
      # term into p,q variational parameters
      Ns <- colSums(r, na.rm=TRUE)
      
      # the stick-breaking point variational parameters
      p <- 1 + Tinv[i] * Ns
      q <- phi_0 + Tinv[i] * (rev(cumsum(rev(Ns))) - Ns)
      
      # term into alpha,mu variational parameters
      r.colSums <- matrix(Ns, nrow=D, ncol=M, byrow=TRUE)
      
      # first Gamma variational parameter for the variances of the Gaussian mixture
      alpha <- alpha_0 + (1/2) * r.colSums *Tinv[i]
      
      # calculation of terms conatined into the second Gamma variational parameter
      x2.r <- t(X^(2)) %*% r
      x.r.m <- (t(X) %*% r) * E.mu
      Ns.m2 <- r.colSums * (E.mu)^2
      Ns.s2 <- r.colSums * V.mu
      N.x_m <- x2.r - 2 * x.r.m + Ns.m2 + Ns.s2
      # second Gamma variational parameter for the variances
      beta <- beta_0 + (1/2) * N.x_m *Tinv[i]
      
      # variational Gaussian variance for the means of the Gaussian mixture
      s2 <- ((1/s2_0) + Tinv[i]* E.inv_sigma2 * r.colSums)^(-1)
      
      # variational Gaussian mean for the means of the Gaussian mixture
      m <- s2 * ((m_0/s2_0) +  Tinv[i]* E.inv_sigma2 * (t(X) %*% r))

      # the stick-breaking point parameter
      lambda <- head(p, M-1) / (head(p, M-1) + head(q, M-1))
      lambda <- c(lambda, 1)
      
      # the variational weights after the stick-breaking point computation
      for (k in 1:M)
      {
        pi[k] <- lambda[k] * prod(head(1-lambda, k-1))
      }
      
      # the variational weights in each iteration
      w[i, ] <- round(pi,3)
      
      # update the expectations contained into the ELBO
      E.mu <- m
      V.mu <- s2 
      E.logsigma2 <- log(beta) - digamma(alpha)
      E.inv_sigma2 <- alpha/beta 
      E.loglambda <- digamma(p) - digamma(p + q)
      E.log1_lambda <-  digamma(q) - digamma(p + q)
      
      ## ELBO 
      # each term into the ELBO has been calculated individually 
      l1 <- sum (r * Tinv[i] *log_rho)
      l2 <- sum(log( phi_0) + (phi_0 - 1) * E.log1_lambda )
      l3 <- - (D * M /2) * log(2 * 355/113) - (1/2) * (sum( log(s2_0)) + sum(s2_0^(-1) * ((E.mu - m_0)^(2) + V.mu)))
      l4 <- sum(alpha_0 * log(beta_0)) - sum(lgamma(alpha_0)) - sum((alpha_0 + 1) * E.logsigma2) - sum(beta_0 * E.inv_sigma2)
      l5 <- - sum(r * log(r))
      
      # useful names for terms in l6
      log_g.p.q <- lgamma( p + q)
      log_p <- lgamma(p)
      log_q <- lgamma(q)
      
      l6_1 <-  - sum(log_g.p.q) + sum(log_p) + sum(log_q)
      l6_2 <-  - sum((p -1) * E.loglambda ) - sum((q -1) * E.log1_lambda)
      l6 <- l6_1 + l6_2
      
      l7 <- (D * M /2) * log(2 * 355/113) + (1/2) * sum(log(V.mu)) + M*D/2
      l8 <- - sum(alpha * log(beta)) + sum(lgamma(alpha)) + sum((alpha + 1) * E.logsigma2) + sum(beta * E.inv_sigma2)
      
      # Total ELBO calculation
      L[i] <- l1 + l2 + l3 + l4 + l5 + l6 +l7 +l8

      # print ELBO value and difference with the previous one
      if (printL) { cat("Iter:\t", i, "\tELBO:\t", L[i], "\tELBO_diff:\t", L[i] - L[i-1], "\n")}
      
      # test if ELBO decreases
      if (L[i] < L[i - 1]) { message("Warning: ELBO decreases\n"); }
      
      # test convergence with epsilon threshold
      if (abs(L[i] - L[i - 1]) < epsilon) { break }
      
      # test VB needs more iteration to converge
      if (i == max_iterations) {warning("VB did not converge\n")}
    }
  
  object <- structure(list(alpha=alpha, beta=beta, m=m, s2=s2, p=p, q=q, r=r, L=L[2:i], w=w))
  
  return(object)
} 
##############----------------------------------------------------------------------------------------------

### example model for data from gmm_generator
# pre-define max iterations for AVB
iter <- 10000
# pre-define T: a vector with K elements and its inverse value Tinv
length <- 100
#set.seed(123)
T <- c(rev(seq(1, 10, length.out = length)),rep(1,iter-length))
Tinv <- 1/T

model <- avb_dpgm(data, M, alpha, beta, m, s2, p, q, alpha_0, beta_0, m_0, s2_0, phi_0, T=T, max_iterations = iter)

## ELBO plot

# truncate the initial elbo iterations 
o<-5
Elbo <- data.frame(elbo=model$L[o:length(model$L)], iter=o:length(model$L))

library(ggplot2)

pel <- ggplot(data=Elbo,aes(x=iter,y=elbo)) +xlab("iterations")+ylab("ELBO") + labs(title="ELBO of AVB-DP-GaussianMix") 
pel <- pel + geom_line(size=1.1,color="royalblue2")  

pel <- pel + theme(
  panel.background = element_rect(colour="gray91", fill="gray91"), axis.text.x=element_text(size=12, face="bold"),
  axis.text.y=element_text(size=12, face="bold"),
  axis.title.x=element_text(size=14,face="bold", colour="slategrey"),
  axis.title.y=element_text(size=16,face="bold", colour = "royalblue"),
  plot.title = element_text(size=14,face="bold") 
)
print(pel)

# calculate the final variational weights
lambda <- head(model$p, M-1) / (head(model$p, M-1) + head(model$q, M-1))
lambda <- c(lambda, 1)
pi_M_f <- rep(0, M)

for (i in 1:M)
{
  pi_M_f[i] <- lambda[i] * prod(head(1-lambda, i-1))
} 
pi_M<-round(pi_M_f[pi_M_f>1e-2], 3)

# the final variational weights
print(pi_M)
