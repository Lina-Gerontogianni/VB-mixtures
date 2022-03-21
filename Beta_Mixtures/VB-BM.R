##############----------------------------------------------------------------------------------------------
#### main AVB function for the Beta Mixture Model + extra helpful functions
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
#### main AVB function for BetaMix
# inputs: 
  ## X: NxD data matrix
  ## M: number of initial components
  ## alpha/beta: DxM Gamma initial variational matrices for the first Beta parameter
  ## mu/eta: DxM Gamma initial variational matrices for the second Beta parameter
  ## p/q: 1xM initial variational vectors for the stick-breaking point Beta parameter 
  ## a0: 1xM Dirichlet hyperparameter matrices for the weights 
  ## alpha_0/beta_0: DxM Gamma hyperparameter matrices for the first Beta parameter
  ## mu_0/eta_0: DxM Gamma hyperparameter matrices for the second Beta parameter
  ## Temp: the temperature vector for the annealing part (pre-define)
  ## max_iterations: maximum number of VB iterations
  ## epsilon: threshold to achieve convergence
# output: 
  ## alpha/beta: DxM Gamma variational matrices for the first Beta parameter
  ## mu/eta: DxM Gamma variational matrices for the second Beta parameter 
  ## a: 1xM Dirichlet variational matrix for the weights 
  ## r: NxD variational matrix for the latent allocation z
  ## L: ELBO values
  ## w: weight values in each iteration (for evolution purposes)
  ## printL: print the ELBO values and difference to the previous one
avb_bm <- function(X, M, alpha, beta, mu, eta, a, alpha_0, beta_0, 
                    mu_0, eta_0, a_0, Temp, max_iterations=iter, epsilon=1e-4, printL=FALSE)
{
  # define the dimensions according to the dataset
  X <- as.matrix(X)
  D <- ncol(X)
  N <- nrow(X)
  # initial objects to receive the variational results
  # initial ELBO
  L <- rep(-Inf, max_iterations)
  # Temperature constant
  Tinv <- 1/Temp
    ### AVB scheme
    for (i in 2:max_iterations) 
      {
      # calculation of expectations contained into the variational parameters
      ubar <- alpha / beta
      vbar <- mu / eta
      E.logu <- digamma(alpha) - log(beta)
      E.logv <- digamma(mu) - log(eta) 
      E.log_pi <- digamma(a) - digamma(sum(a))
      # helpful names for terms in the variational equations
      log_g.u.v <- lgamma(ubar + vbar)
      log_g.u <- lgamma(ubar)
      log_g.v <- lgamma(vbar)
      dig.u.v <- digamma(ubar + vbar)
      dig.u <- digamma(ubar)
      dig.v <- digamma(vbar)
      log.u <- log(ubar)
      log.v <- log(vbar)
    
      R <- log_g.u.v - log_g.u - log_g.v + (dig.u.v - dig.u) * ubar * (E.logu -  log.u) + (dig.u.v - dig.v) * vbar * (E.logv -  log.v)
      # calculation of expectations contained into the variational equation of r for the allocation z
      N.E.log_pi <- matrix(E.log_pi, nrow=N, ncol=M, byrow=TRUE)
      R.colsums <- matrix(colSums(R, na.rm=TRUE), nrow=N, ncol=M, byrow=TRUE )
      log_rho <- N.E.log_pi + R.colsums + log(X) %*% (ubar - 1)  +  log(1 - X) %*% (vbar - 1)
      # the usefulness of logsumexp function
      S <- apply(log_rho, 1, logsumexp)
      log_r <- log_rho - S
      # the variational parameter of z
      r <- apply(log_r, 2, exp)
      # trick to avoid zero values
      r <- (r + 10^-9)^ Tinv[i] 
      # term into p,q variational parameters
      Ns <- colSums(r, na.rm=TRUE)
      # the Dirichlet variational parameter
      a <- a_0 + Tinv[i]*Ns
      # term into alpha,mu variational parameters
      r.colSums <- matrix(Ns, nrow=D, ncol=M, byrow=TRUE)
      # Gamma variational parameters for the Beta mixture parameters
      alpha <- alpha_0 + Tinv[i] * r.colSums * ubar * (dig.u.v - dig.u) 
      beta <- beta_0 - Tinv[i] * t(t(r) %*% log(X))
      mu <- mu_0 + Tinv[i] * r.colSums * vbar * (dig.u.v - dig.v) 
      eta <- eta_0 - Tinv[i] * t(t(r) %*% log(1-X))
      # the variational weights
      pi_M <- (1/(M*a_0 + N)) *a
      # update the expectations contained into the ELBO
      ubar <- alpha / beta
      vbar <- mu / eta
      E.logu <- digamma(alpha) - log(beta)
      E.logv <- digamma(mu) - log(eta) 
      E.log_pi <- digamma(a) - digamma(sum(a))
      # Dirichlet constant
      logC <- function(u)
      {
        return(lgamma(sum(u)) - sum(lgamma(u)))
      }
      ## ELBO 
      # each function into the ELBO has been calculated individually (l1,l2,l3,l4,l5 and l6) and finally all are summed up
      l1 <- sum(r * log_rho)
      l2 <- sum((a_0 - 1) * E.log_pi) + logC(a_0)
      # helpful names for terms in l3
      log_g.alpha_0 <- lgamma(alpha_0)
      log_g.mu_0 <- lgamma(mu_0)
      l3_1 <- sum(alpha_0 * log(beta_0)) - sum(log_g.alpha_0) + sum((alpha_0 - 1) * E.logu)  - sum(beta_0 * ubar)
      l3_2 <- sum (mu_0 * log(eta_0)) - sum(log_g.mu_0) + sum((mu_0 - 1) * E.logv)  - sum(eta_0 * vbar)
      l3 <- l3_1 + l3_2
      l4 <- - sum(r * log_r)
      # helpful names for terms in l5
      log_g.alpha <- lgamma(alpha)
      log_g.mu <- lgamma(mu)
      l5_1 <- - sum(alpha * log(beta)) + sum( log_g.alpha) - sum((alpha - 1) * E.logu )  + sum(beta * ubar)
      l5_2 <- - sum(mu * log(eta)) + sum(log_g.mu) - sum((mu - 1) * E.logv)  + sum(eta * vbar)
      l5 <- l5_1 + l5_2
      # helpful names for terms in l6
      l6 <- sum((a - 1) * E.log_pi) + logC(a)
      
      # Total ELBO calculation
      L[i] <- l1 + l2 + l3 + l4 + l5 + l6

      # print ELBO value and difference with the previous one
      if (printL) { cat("Iter:\t", i, "\tELBO:\t", L[i], "\tELBO_diff:\t", L[i] - L[i-1], "\n")}
      # test if ELBO decreases
      if (L[i] < L[i - 1]) { message("Warning: ELBO decreases\n"); }
      # test convergence with epsilon threshold
      if (abs(L[i] - L[i - 1]) < epsilon) { break }
      # test VB needs more iteration to converge
      if (i == max_iterations) {warning("VB did not converge\n")}
    }
  
  object <- structure(list(alpha=alpha, beta=beta, mu=mu, eta=eta, a=a, r=r, pi_M=pi_M,ubar=ubar,vbar=vbar,L=L[2:i]))
  return(object)
} 
##############----------------------------------------------------------------------------------------------
### example model for data from bmm_generator
# pre-define max iterations for AVB
iter <- 10000
# pre-define T: a vector with K elements and its inverse value Tinv
length <-100
Temp <- c(rev(seq(1, 10, length.out = length)),rep(1,iter-length))
# without annealing replace Temp with 1
Temp <- rep(1,iter)
# the model
model_beta <- avb_bm(data, M, alpha, beta, mu, eta, a, alpha_0, beta_0, mu_0, eta_0, a_0, Temp=Temp, max_iterations = iter)

## ELBO plot
# truncate the initial elbo iterations 
o<-5
Elbo <- data.frame(elbo=model_beta$L[o:length(model_beta$L)], iter=o:length(model_beta$L))

library(ggplot2)

pel <- ggplot(data=Elbo,aes(x=iter,y=elbo)) +xlab("iterations")+ylab("ELBO") + labs(title="ELBO of AVB-BetaMix") 
pel <- pel + geom_line(size=1.1,color="royalblue2")  
pel <- pel + theme(
  panel.background = element_rect(colour="gray91", fill="gray91"), axis.text.x=element_text(size=12, face="bold"),
  axis.text.y=element_text(size=12, face="bold"),
  axis.title.x=element_text(size=14,face="bold", colour="slategrey"),
  axis.title.y=element_text(size=16,face="bold", colour = "royalblue"),
  plot.title = element_text(size=14,face="bold") 
)
print(pel)

# final component weights (discard those with weight < 0.003)
var_weights <- model_beta$pi_M[model_beta$pi_M>0.003]
cat("Variational weights:", var_weights, "\n")

var_par_u<- as.matrix(model_beta$ubar[ ,model_beta$pi_M>0.003])
cat("Variational parameters:\n");print(var_par_u)

var_par_v<- as.matrix(model_beta$vbar[ ,model_beta$pi_M>0.003])
cat("Variational parameters:\n");print(var_par_v)