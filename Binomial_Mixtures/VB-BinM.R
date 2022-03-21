##############----------------------------------------------------------------------------------------------
#### main AVB function for the Binomial Mixture Model + extra helpful functions
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
#### main AVB function for BinMix
# inputs: 
## Y: NxD data matrix
## S: NxD matrix of the Bernoulli trials
## M: number of initial components
## phi_0: 1xM Dirichlet hyperparameter matrix for the weights
## a_0/b_0: DxM Beta hyperparameter matrices for the Binomial probability parameter
## phi: 1xM Dirichlet initial variational matrix for the weights 
## a/b: DxM Beta initial variational matrices for the Binomial probability parameter
## Temp: the temperature vector for the annealing part (pre-defined)
## max_iterations: maximum number of VB iterations
## epsilon: threshold to achieve convergence
# output: 
## phi: 1xM Dirichlet variational matrix for the weights 
## a/b: DxM Beta variational matrices for the Binomial probability parameter
## r: NxD variational matrix for the latent allocation z
## pi_M: variational weights
## par: variational Binomial parameters (mean binomial probability values for each cluster)
## L: ELBO values
## printL: print the ELBO values and difference to the previous one

vb_binmix <- function(Y, S, M, phi_0, a_0, b_0, phi, a, b, Temp, max_iterations, epsilon=1e-4, printL=FALSE)
{
  # define the dimensions according to the dataset
  Y <- as.matrix(Y)
  S<-as.matrix(S)
  D <- ncol(Y)
  N <- nrow(Y)
  S_Y <- S - Y
  # Temperature constant
  Tinv <- 1/Temp
  # initial ELBO to receive the variational results
  L <- rep(-Inf, max_iterations)
  
  ### VB scheme
  for (i in 2:max_iterations) 
  {
    # calculation of expectations involved into the variational parameters
    E.logp <- digamma(a) - digamma(a+b)
    E.log1_p <- digamma(b) - digamma(a+b) 
    E.log_pi <- digamma(phi) - digamma(sum(phi))
    # calculation of expectations contained into the variational equation of r for the allocation z
    N.E.log_pi <- matrix(E.log_pi, nrow=N, ncol=M, byrow=TRUE)
    
    log_rho <- N.E.log_pi + Y %*% E.logp + S_Y %*% E.log1_p 
    # the usefulness of logsumexp function
    Sr <- apply(log_rho, 1, logsumexp)
    log_r <- log_rho - Sr
    # the variational parameter of z
    r <- apply(log_r, 2, exp)
    r <- (r + 10^-9)^ Tinv[i] 
    # term into a,b variational parameters
    Ns <- colSums(r, na.rm=TRUE)+ 1e-10
    # the Beta variational parameters
    a <- a_0 + Tinv[i]*t(Y) %*% r
    b <- b_0 + Tinv[i]*t(S) %*% r - Tinv[i]*t(Y) %*% r 
    # the Dirichlet variational parameter
    phi <- phi_0 + Tinv[i]*Ns
    # the variational weights
    pi_M <- (phi) / (M * phi_0 + N)
    # the variation mean values of the Binomial parameter
    par <- a / (a + b)
    # update of expectations included in the variational parameters
    E.log_pi <- digamma(phi) - digamma(sum(phi))
    N.E.log_pi <- matrix(E.log_pi, nrow=N, ncol=M, byrow=TRUE)
    # term inside the ELBO
    log.ch = lchoose(S, Y)
    lch <- rowSums(log.ch,na.rm=TRUE)
    N.lch <- matrix(lch, nrow=N, ncol=M, byrow=FALSE)
    # Dirichlet constant
    logC <- function(u)
    {
      return(lgamma(sum(u)) - sum(lgamma(u)))
    }
    ## ELBO 
    l1 <- sum( lbeta(a,b) - lbeta(a_0,b_0) )
    l2 <- sum((a_0 - a) * digamma(a) + (b_0 - b) * digamma(b) + (a - a_0 + b - b_0) * digamma(a+b))
    l3 <- sum(r * (Y %*% (digamma(a)-digamma(b)) + S %*% (digamma(b)-digamma(a+b)) + N.E.log_pi + N.lch -log_r))
    l4 <- sum((phi_0 - phi) %*% (E.log_pi) + logC(phi_0) - logC(phi))
    
    # Total ELBO calculation
    L[i] <- l1 + l2 + l3 + l4
    # print ELBO value and difference with the previous one
    if (printL) { cat("Iter:\t", i, "\tELBO:\t", L[i], "\tELBO_diff:\t", L[i] - L[i-1], "\n")}
    # test if ELBO decreases
    if (L[i] < L[i - 1]) { message("Warning: ELBO decreases\n"); }
    # test convergence with epsilon threshold
    if (abs(L[i] - L[i - 1]) < epsilon) { break }
    # test VB needs more iteration to converge
    if (i == max_iterations) {warning("VB did not converge\n")}
  }
  object <- structure(list(phi=phi, a=a, b=b, r=r, pi_M=pi_M, par=par, L=L[2:i]))
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
model_bin <- vb_binmix(Y, S, M, phi_0, a_0, b_0, phi, a, b, Temp, max_iterations=iter, epsilon=1e-10, printL=FALSE)

## ELBO plot
# truncate the initial elbo iterations (for graphical reasons only)
o<-2
Elbo <- data.frame(elbo=model_bin$L[o:length(model_bin$L)], iter=o:length(model_bin$L))

library(ggplot2)
pel <- ggplot(data=Elbo,aes(x=iter,y=elbo)) +xlab("iterations")+ylab("ELBO") + labs(title="ELBO of VB-BinMix") 
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
var_weights <- model_bin$pi_M[model_bin$pi_M>0.003]
cat("Variational weights:", var_weights, "\n")

var_par<- as.matrix(model_bin$par[ ,model_bin$pi_M>0.003])
cat("Variational parameters:\n");print(var_par)