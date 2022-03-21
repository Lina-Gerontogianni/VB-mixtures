##############----------------------------------------------------------------------------------------------
#### main AVB function for the DP Beta Mixture Model + extra helpful functions
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
#### main AVB function for DP BetaMix
# inputs: 
  ## X: NxD data matrix
  ## M: number of initial components
  ## alpha/beta: DxM Gamma initial variational matrices for the first Beta parameter
  ## mu/eta: DxM Gamma initial variational matrices for the second Beta parameter
  ## p/q: 1xM initial variational vectors for the stick-breaking point Beta parameter 
  ## phi_0: an 1xM vector for the stick-breaking point Beta parameter 
  ## alpha_0/beta_0: DxM Gamma hyperparameter matrices for the first Beta parameter
  ## mu_0/eta_0: DxM Gamma hyperparameter matrices for the second Beta parameter
  ## Temp: the temperature vector for the annealing part (pre-define)
  ## max_iterations: maximum number of VB iterations
  ## epsilon: threshold to achieve convergence
# output: 
  ## alpha/beta: DxM Gamma variational matrices for the first Beta parameter
  ## mu/eta: DxM Gamma variational matrices for the second Beta parameter 
  ## p/q: 1xM variational vectors for the stick-breaking point Beta parameter 
  ## r: NxD variational matrix for the latent allocation z
  ## L: ELBO values
  ## w: weight values in each iteration (for evolution purposes)
  ## printL: print the ELBO values and difference to the previous one

avb_dpbm <- function(X, M, alpha, beta, mu, eta, p, q, alpha_0, beta_0, 
                    mu_0, eta_0, phi_0, Temp, max_iterations=iter, epsilon=1e-4, printL=FALSE)
{
  # define the dimensions according to the dataset
  X <- as.matrix(X)
  D <- ncol(X)
  N <- nrow(X)
  # Temperature constant
  Tinv <- 1/Temp
  # initial ELBO
  L <- rep(-Inf, max_iterations)
  # initial weights in final iteration
  pi_M <- rep(0, M)
    ### AVB scheme
    for (i in 2:max_iterations) 
      {
      # calculation of expectations contained into the variational parameters
      ubar <- alpha / beta
      vbar <- mu / eta
      E.logu <- digamma(alpha) - log(beta)
      E.logv <- digamma(mu) - log(eta) 
      E.loglambda <- digamma(p) - digamma(p + q)
      E.log1_lambda <-  digamma(q) - digamma(p + q)
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
      E.log1_lambdaj <- c(0,cumsum(E.log1_lambda)[1:(M-1)])
      N.E.log1_lambdaj <- matrix(E.log1_lambdaj, nrow=N, ncol=M, byrow=TRUE)
      N.E.loglambda <- matrix(E.loglambda, nrow=N, ncol=M, byrow=TRUE)
      R.colsums <- matrix(colSums(R, na.rm=TRUE), nrow=N, ncol=M, byrow=TRUE )

      log_rho <- N.E.loglambda + N.E.log1_lambdaj +  R.colsums + log(X) %*% (ubar - 1)  +  log(1 - X) %*% (vbar - 1)
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
      # Gamma variational parameters for the Beta mixture parameters
      alpha <- alpha_0 + Tinv[i] * r.colSums * ubar * (dig.u.v - dig.u) 
      beta <- beta_0 - Tinv[i] * t(t(r) %*% log(X))
      mu <- mu_0 + Tinv[i] * r.colSums * vbar * (dig.u.v - dig.v) 
      eta <- eta_0 - Tinv[i] * t(t(r) %*% log(1-X))
      # the stick-breaking point parameter
      lambda <- head(p, M-1) / (head(p, M-1) + head(q, M-1))
      lambda <- c(lambda, 1)
      # the variational weights after the stick-breaking point computation
      for (k in 1:M)
      {
        pi_M[k] <- lambda[k] * prod(head(1-lambda, k-1))
      }
      # update the expectations contained into the ELBO
      ubar <- alpha / beta
      vbar <- mu / eta
      E.logu <- digamma(alpha) - log(beta)
      E.logv <- digamma(mu) - log(eta) 
      E.loglambda <- digamma(p) - digamma(p + q)
      E.log1_lambda <-  digamma(q) - digamma(p + q)
      
      ## ELBO 
      # each function into the ELBO has been calculated individually (l1,l2,l3,l4,l5 and l6) and finally all are summed up
      l1 <- sum(r * log_rho)
      l2 <- sum(log( phi_0 ) + (phi_0 - 1) * E.log1_lambda)
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
      log_g.p.q <- lgamma(p + q)
      log_p <- lgamma(p)
      log_q <- lgamma(q)
      l6_1 <-  - sum(log_g.p.q) + sum(log_p) + sum(log_q)
      l6_2 <-  - sum((p -1) * E.loglambda) - sum((q -1) * E.log1_lambda)
      l6 <- l6_1 + l6_2
      
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
  object <- structure(list(alpha=alpha, beta=beta, mu=mu, eta=eta, p=p, q=q, r=r, pi_M=pi_M, ubar=ubar, vbar=vbar, L=L[2:i]))
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
model_dpbeta <- avb_dpbm(data, M, alpha, beta, mu, eta, p, q, alpha_0, beta_0, mu_0, eta_0, phi_0, Temp=Temp, max_iterations = iter)

## ELBO plot
# truncate the initial elbo iterations 
o<-5
Elbo <- data.frame(elbo=model_dpbeta$L[o:length(model_dpbeta$L)], iter=o:length(model_dpbeta$L))

library(ggplot2)

pel <- ggplot(data=Elbo,aes(x=iter,y=elbo)) +xlab("iterations")+ylab("ELBO") + labs(title="ELBO of AVB-DP-BetaMix") 
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
var_weights <- model_dpbeta$pi_M[model_dpbeta$pi_M>0.003]
cat("Variational weights:", var_weights, "\n")

var_par_u<- as.matrix(model_dpbeta$ubar[ ,model_dpbeta$pi_M>0.003])
cat("Variational parameters:\n");print(var_par_u)

var_par_v<- as.matrix(model_dpbeta$vbar[ ,model_dpbeta$pi_M>0.003])
cat("Variational parameters:\n");print(var_par_v)