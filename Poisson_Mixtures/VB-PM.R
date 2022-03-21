##############----------------------------------------------------------------------------------------------
#### main AVB function for the Poisson Mixture Model + extra helpful functions
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
#### main AVB function for PoisMix
# inputs: 
## Y: NxD data matrix
## M: number of initial components
## phi_0: 1xM Dirichlet hyperparameter matrix for the weights
## a_0/b_0: DxM Gamma hyperparameter matrices for the Poisson parameter
## phi: 1xM Dirichlet initial variational matrix for the weights 
## a/b: DxM Gamma initial variational matrices for the Poisson parameter
## Temp: the temperature vector for the annealing part (pre-defined)
## max_iterations: maximum number of VB iterations
## epsilon: threshold to achieve convergence
# output: 
## phi: 1xM Dirichlet variational matrix for the weights 
## a/b: DxM Gamma variational matrices for the Poisson parameter
## r: NxD variational matrix for the latent allocation z
## pi_M: variational weights
## E.lambda: variational Poisson parameters (mean lambda values for each cluster)
## L: ELBO values
## printL: print the ELBO values and difference to the previous one

vb_poismix <- function(Y, M, phi_0, a_0, b_0, phi, a, b, Temp, max_iterations, epsilon=1e-4, printL=FALSE)
{
  # define the dimensions according to the dataset
  Y <- as.matrix(Y)
  D <- ncol(Y)
  N <- nrow(Y)
  # Temperature constant
  Tinv <- 1/Temp
  # initial ELBO
  L <- rep(-Inf, max_iterations)
  
  ### VB scheme
  for (i in 2:max_iterations) 
  {
    # calculation of expectations contained into the variational parameters
    E.lambda <- a/b
    E.loglambda <- digamma(a) - log(b) 
    E.log_pi <- digamma(phi) - digamma(sum(phi))
    N.E.log_pi <- matrix(E.log_pi, nrow=N, ncol=M, byrow=TRUE)
    N.E.lambda <- matrix(colSums(E.lambda, na.rm=TRUE),nrow=N,ncol=M,byrow = TRUE)
    N.Y <- matrix(rowSums(log(factorial(Y)),na.rm=TRUE), nrow=N, ncol=M, byrow = FALSE)
  
    log_rho <- N.E.log_pi + Y %*% E.loglambda - N.E.lambda  
    # the usefulness of logsumexp function
    S <- apply(log_rho, 1, logsumexp)
    log_r <- log_rho - S
    # the variational parameter of z
    r <- apply(log_r, 2, exp)
    # trick to avoid zero values
    r <- (r + 10^-9)^ Tinv[i]
    # terms into the variational parameters
    Ns <- colSums(r, na.rm=TRUE)
    D.Ns <- matrix(Ns,nc=M,nr=D,byrow=TRUE)
    # the Gamma variational parameters
    a <- a_0 + Tinv[i]*t(Y) %*% r
    b <- b_0 + Tinv[i]*D.Ns 
    # the Dirichlet variational parameter
    phi <- phi_0 + Tinv[i]*Ns
    # the variational weights
    pi_M <- (phi) / (M * phi_0 + N)
    # update expectations
    E.lambda <- a/b
    E.loglambda <- digamma(a) - log(b) 
    E.log_pi <- digamma(phi) - digamma(sum(phi))
    N.E.log_pi <- matrix(E.log_pi, nrow=N, ncol=M, byrow=TRUE)
    N.E.lambda <- matrix(colSums(E.lambda, na.rm=TRUE),nrow=N,ncol=M,byrow = TRUE)
    # Dirichlet constant
    logC <- function(u)
    {
      return(lgamma(sum(u)) - sum(lgamma(u)))
    }
    ## ELBO 
    l1 <- sum(a_0*log(b_0) - a*log(b) -lgamma(a_0) + lgamma(a) + (a_0-a)*E.loglambda +(b - b_0)*E.lambda )
    l2 <- sum(r * (log_rho - log_r - N.Y))
    l3 <- sum((phi_0 - phi) %*% (E.log_pi) + logC(phi_0) - logC(phi))
    
    # Total ELBO calculation
    L[i] <-  l1 + l2 + l3
    
    # print ELBO value and difference with the previous one
    if (printL) { cat("Iter:\t", i, "\tELBO:\t", L[i], "\tELBO_diff:\t", L[i] - L[i-1], "\n")}
    # test if ELBO decreases
    if (L[i] < L[i - 1]) { message("Warning: ELBO decreases\n"); }
    # test convergence with epsilon threshold
    if (abs(L[i] - L[i - 1]) < epsilon) { break }
    # test VB needs more iteration to converge
    if (i == max_iterations) {warning("VB did not converge\n")}
  }
  object <- structure(list(phi=phi, a=a, b=b, r=r, pi_M=pi_M, E.lambda=E.lambda, L=L[2:i]))
  return(object)
} 
##############----------------------------------------------------------------------------------------------

### example model 
# pre-define max iterations for AVB
iter <- 5000
# pre-define Temp: a vector with K elements and its inverse value Tinv
length <-100
Temp <- c(rev(seq(1, 10, length.out = length)),rep(1,iter-length))
# without annealing replace Temp with 1
Temp <- rep(1,iter)
# the model
model_pois <- vb_poismix(Y, M, phi_0, a_0, b_0, phi, a, b, Temp, max_iterations=iter, epsilon=1e-10, printL=T)

## ELBO plot
# truncate the initial elbo iterations 
o<-1
Elbo <- data.frame(elbo=model_pois$L[o:length(model_pois$L)], iter=o:length(model_pois$L))
library(ggplot2)

pel <- ggplot(data=Elbo,aes(x=iter,y=elbo)) +xlab("iterations")+ylab("ELBO") + labs(title="ELBO of AVB-PoisMix") 
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
var_weights <- model_pois$pi_M[model_pois$pi_M>0.003]
cat("Variational weights:", var_weights, "\n")

var_par<- as.matrix(model_pois$E.lambda[ ,model_pois$pi_M>0.003])
cat("Variational parameters:\n");print(var_par)