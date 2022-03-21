##############----------------------------------------------------------------------------------------------
#### main AVB function for the DP Binomial Mixture Model + extra helpful functions
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
##############----------------------------------------------------------------------------------------------
#### main AVB function for DP BinMix
# inputs: 
## Y: NxD data matrix
## S: NxD matrix of the Bernoulli trials
## M: number of initial components
## phi_0: 1xM Beta hyperparameter matrix for the stick-breaking point process
## a_0/b_0: DxM Beta hyperparameter matrices for the Binomial probability parameter
## phi/delta: 1xM initial Beta variational matrices for the stick-breaking point process 
## a/b: DxM Beta initial variational matrices for the Binomial probability parameter
## Temp: the temperature vector for the annealing part (pre-defined)
## max_iterations: maximum number of VB iterations
## epsilon: threshold to achieve convergence
# output: 
## phi/delta: 1xM Beta variational matrices for stick-breaking point process 
## a/b: DxM Beta variational matrices for the Binomial probability parameter
## r: NxD variational matrix for the latent allocation z
## pi_M: variational weights
## par: variational Binomial parameters (mean binomial probability values for each cluster)
## L: ELBO values
## printL: print the ELBO values and difference to the previous one

vb_dpbinmix <- function(Y, S, M, phi_0, a_0, b_0, delta, phi, a, b, max_iterations, epsilon=1e-4, printL=FALSE)
{
  # define the dimensions according to the dataset
  Y <- as.matrix(Y)
  S <-as.matrix(S)
  D <- ncol(Y)
  N <- nrow(Y)
  S_Y <- S - Y
  # Temperarure constant
  Tinv <- 1/Temp
  # initial ELBO to receive the variational results
  L <- rep(-Inf, max_iterations)
  
  ### VB scheme
  for (i in 2:max_iterations) 
  {
    # calculation of expectations involved into the variational parameters
    E.logp <- digamma(a) - digamma(a+b)
    E.log1_p <- digamma(b) - digamma(a+b) 
    E.logw <- digamma(delta) - digamma(delta + phi)
    E.log1_w <-  digamma(phi) - digamma(delta + phi)
    # calculation of expectations contained into the variational equation of r for the allocation z
    N.E.logw <- matrix(E.logw, nrow=N, ncol=M, byrow=TRUE)
    E.log1_wj <- c(0,cumsum(E.log1_w)[1:(M-1)])
    N.E.log1_wj <- matrix(E.log1_wj, nrow=N, ncol=M, byrow=TRUE)
    
    log_rho <- N.E.logw + N.E.log1_wj + Y %*% E.logp + S_Y %*% E.log1_p 
    # the usefulness of logsumexp function
    Sr <- apply(log_rho, 1, logsumexp)
    log_r <- log_rho - Sr
    # the variational parameter of z
    r <- apply(log_r, 2, exp)
    r <- (r + 10^-9)^ Tinv[i] 
    # term into a,b variational parameters
    Ns <- colSums(r, na.rm=TRUE)+ 1e-10
    # the Beta variational parametres
    a <- a_0 + Tinv[i]*t(Y) %*% r
    b <- b_0 + Tinv[i]*t(S) %*% r - Tinv[i]*t(Y) %*% r 
    # the Beta variational parameters for the DP
    delta <- 1 + Tinv[i]*Ns
    phi <- phi_0 + Tinv[i]*(rev(cumsum(rev(Ns))) - Ns)
    # the variational weights according to DP
    w <- head(delta, M-1) / (head(delta, M-1) + head(phi, M-1))
    w <- c(w, 1)
    pi_M <- rep(0, M)
    
    for (j in 1:M)
    {
      pi_M[j] <- w[j] * prod(head(1-w, j-1))
    }
    # the variation mean values of the Binomial parameter
    par <- a/(a+b)
    # update of expectations included in the variational parameters
    E.log1_w <-  digamma(phi) - digamma(delta + phi)
    E.logw <- digamma(delta) - digamma(delta + phi)
    N.E.logw <- matrix(E.logw, nrow=N, ncol=M, byrow=TRUE)
    E.log1_wj <- c(0,cumsum(E.log1_w)[1:(M-1)])
    N.E.log1_wj <- matrix(E.log1_wj, nrow=N, ncol=M, byrow=TRUE)
    # term inside the ELBO
    log.ch = lchoose(S, Y)
    lch <- rowSums(log.ch,na.rm=TRUE)
    N.lch <- matrix(lch, nrow=N, ncol=M, byrow=FALSE)
    
    ## ELBO 
    l1 <- sum( lbeta(a,b) - lbeta(a_0,b_0) )
    l2 <- sum((a_0 - a) * digamma(a) + (b_0 - b) * digamma(b) + (a - a_0 + b - b_0) * digamma(a+b))
    l3 <- sum(r * (Y %*% (digamma(a)-digamma(b)) + S %*% (digamma(b)-digamma(a+b)) + N.E.logw + N.E.log1_wj - log_r + N.lch ))
    l4 <- sum(log( phi_0) + (phi_0 ) * E.log1_w + lbeta(delta,phi) - (delta -1) * E.logw - (phi)*E.log1_w)
    
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
  
  object <- structure(list(delta=delta, phi=phi,  a=a, b=b, r=r, pi_M=pi_M, par=par, L=L[2:i]))
  return(object)
} 
##############----------------------------------------------------------------------------------------------

### example model 
# pre-define max iterations for VB
iter <- 7000
# pre-define T: a vector with K elements and its inverse value Tinv
length <-100
Temp <- c(rev(seq(1, 10, length.out = length)),rep(1,iter-length))
# without annealing replace Temp with 1
Temp <- rep(1,iter)
# the model
model_dpbin <- vb_dpbinmix(Y, S, M, phi_0, a_0, b_0, delta, phi, a, b, max_iterations=iter, epsilon=1e-10, printL=FALSE)

## ELBO plot
# truncate the initial elbo iterations 
o<-2
Elbo <- data.frame(elbo=model_dpbin$L[o:length(model_dpbin$L)], iter=o:length(model_dpbin$L))

library(ggplot2)
pel <- ggplot(data=Elbo,aes(x=iter,y=elbo)) +xlab("iterations")+ylab("ELBO") + labs(title="ELBO of VB-DP-BinMix") 
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
var_weights <- model_dpbin$pi_M[model_dpbin$pi_M>0.003]
cat("Variational weights:", var_weights, "\n")

var_par<- as.matrix(model_dpbin$par[ ,model_dpbin$pi_M>0.003])
cat("Variational parameters:\n");print(var_par)