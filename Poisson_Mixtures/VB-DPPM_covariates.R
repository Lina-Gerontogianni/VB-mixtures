##############----------------------------------------------------------------------------------------------
#### main VB function for the DP Poisson Mixture with covariates model + extra helpful functions
##############----------------------------------------------------------------------------------------------
library(psych)
library(nlme)
library(mvtnorm)
library(base)
############## useful finction for determinant calculation safely
# Compute the log-determinant of a matrix
ldet <- function(l) 
{
  if(!is.matrix(l)) return(log(l))
  determinant(l,logarithm = TRUE)$modulus
}
logsumexp <- function(y) 
{
  # Computes log(sum(exp(x))
  l <- max(y)
  lse <- log(sum(exp(y-l))) + l
  j <- which(!is.finite(lse))
  if (length(j) > 0) {lse[j] <- l}
  
  return(lse)
}
##############----------------------------------------------------------------------------------------------
#### main VB function for DP Poisson Mix with covariates
# inputs: 
## Y: Nx1 data vector 
## X: NxD design matrix
## mu_0/S_0: Dx1 and DxD hyperparameters of the Gaussian upon beta coefficients (mean vector and covariance matrix respectively)
## phi_0: 1xM Beta hyperparameter matrix for the stick-breaking point process
## phi/delta: 1xM initial Beta variational matrices for the stick-breaking point process 
## mu/S: Dx1 and DxD initial variational parameters of the Gaussian upon beta coefficients 
## max_iterations: maximum number of VB iterations
## epsilon: threshold to achieve convergence
# output: 
## mu/S: Dx1 and DxD variational parameters of the Gaussian upon beta coefficients (mean vector and covariance matrix respectively)
## phi/delta: 1xM Beta variational matrices for the stick-breaking point process 
## pi_M: variational weights
## r: NxD variational matrix for the latent allocation z
## L: ELBO values
## printL: print the ELBO values and difference with the previous one

vb_dppoismixreg <- function(Y, X, M, mu, S, delta, phi, mu_0, S_0, phi_0,
                         max_iterations=iter, epsilon=1e-4, printL=FALSE)
{
  # define the dimensions according to the dataset
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  D <- ncol(Y)
  t <- ncol(X)
  N <- nrow(X)
  # initial ELBO and empty useful terms for the inference 
  L <- rep(-Inf, max_iterations)
  l2_m <- matrix(NA, nrow=M, ncol=D)
  NM.logy <- matrix(rowSums(log(Y+10^-10)), nrow=N, ncol=M, byrow=FALSE)
  NM.ylogy <- matrix(rowSums(Y*(log(Y+10^-10))^2), nrow=N, ncol=M, byrow=FALSE)
  E.xbeta <- logy.E.xbeta <- t.XS <- E.Xbeta2 <- array(NA, c(N,D,M))
  D.logy.E.xbeta <- y.Exbeta2 <- matrix(NA, nc=M, nr=N)
  XS <- array(NA,c(N,t,D,M))
  S_inv <- array(NA, c(t,t,D,M))
  pi_M <- rep(0, M)
  
  ### VB scheme
  for (l in 2:max_iterations) 
  {
    # terms inside the variational parameters
    E.logw <- digamma(delta) - digamma(delta + phi)
    E.log1_w <-  digamma(phi) - digamma(delta + phi)
    E.log1_wj <- c(0,cumsum(E.log1_w)[1:(M-1)])
    N.E.log1_wj <- matrix(E.log1_wj, nrow=N, ncol=M, byrow=TRUE)
    N.E.logw <- matrix(E.logw, nrow=N, ncol=M, byrow=TRUE)
    logar.2 <- matrix(log(2),nr=N, nc=M)
    
    for (i in 1:M)
    {
      E.xbeta[,,i] <- X %*% mu[,,i]
      logy.E.xbeta[,,i] <- E.xbeta[,,i] * log(Y+10^-10) *Y
      D.logy.E.xbeta[,i] <- rowSums(logy.E.xbeta[,,i])
      for (d in 1:D)
      {
        XS[,,d,i] <- X %*% S[,,d,i] * X
        t.XS[,d,i] <- rowSums(XS[,,d,i],na.rm=TRUE)
      }  
      E.Xbeta2[ ,,i] <- t.XS[,,i] + E.xbeta[ ,,i]^(2)
      y.Exbeta2[,i] <- rowSums( Y * E.Xbeta2[,,i])
    }
    log_R <- N.E.logw + N.E.log1_wj + D.logy.E.xbeta - (1/2) *y.Exbeta2
    
    # the usefulness of logsumexp function
    K <- apply(log_R, 1, logsumexp)
    log_r <- log_R - K
    # the variational parameter of z
    r <- apply(log_r, 2, exp)
    # trick to avoid zero values
    r <- r + 10^-9 
    # variational parameters of beta 
    for (j in 1:M)
    {
      for (f in 1:D)
      {
        S_inv[,,f,j] <- solve(S_0[,,f,j]) + crossprod(X * r[,j] * Y[,f] ,X)
        S[,,f,j] <- solve(S_inv[,,f,j])
        mu[,f,j] <- S[,,f,j] %*% (crossprod(X,r[,j]*log(Y[,f]+10^-10)*Y[,f]) + solve(S_0[,,f,j]) %*% mu_0[,f,j])
      }
    }
    # term into delta,phi variational parameters
    Ns <- colSums(r, na.rm=TRUE)
    # the stick-breaking point variational parameters
    delta <- 1 + Ns
    phi <- phi_0 + rev(cumsum(rev(Ns))) - Ns
    # calculate the final variational weights
    w <- head(delta, M-1) / (head(delta, M-1) + head(phi, M-1))
    w <- c(w, 1)
    for (i in 1:M)
    {
      pi_M[i] <- w[i] * prod(head(1-w, i-1))
    }
    # update the variational equations
    for (g in 1:M)
    {
      E.xbeta[,,g] <- X %*% mu[,,g]
      for (d in 1:D)
      {
        XS[,,d,g] <- X %*% S[,,d,g] * X
        t.XS[,d,g] <- rowSums(XS[,,d,g],na.rm=TRUE)
      }  
      E.Xbeta2[ ,,g] <- t.XS[,,g] + E.xbeta[ ,,g]^(2)
    }
    E.logw <- digamma(delta) - digamma(delta + phi)
    E.log1_w <-  digamma(phi) - digamma(delta + phi)
    E.log1_wj <- c(0,cumsum(E.log1_w)[1:(M-1)])
    N.E.log1_wj <- matrix(E.log1_wj, nrow=N, ncol=M, byrow=TRUE)
    N.E.logw <- matrix(E.logw, nrow=N, ncol=M, byrow=TRUE)
    
    ## ELBO 
    for (k in 1:M)
    {
    for (c in 1:D)
      {
      l2_m[k,c] <- -(1/2) * ldet(S_0[,,c,k]) -(1/2) * t(mu[,c,k] - mu_0[,c,k]) %*% solve(S_0[,,c,k]) %*% (mu[,c,k] - mu_0[,c,k])
      -(1/2) * tr(solve(S_0[,,c,k]) %*% S[,,c,k]) + (1/2)*ldet(S[,,c,k])
      }
    }
    l1 <- sum(r * (-(1/2)*NM.logy -(1/2)*NM.ylogy + log_R )) -(D/2)*log(2*pi)
    l2 <- sum(l2_m)
    l3 <- sum(log(phi_0) + (phi_0 -1) * E.log1_w)
    l4 <- -sum(r * log_r) + D*M*t/2
    l5_1 <- sum( lbeta(delta, phi) )
    l5_2 <-  - sum((delta -1) * E.logw - (phi -1) * E.log1_w)
    l5 <- l5_1 + l5_2

    # Total ELBO calculation
    L[l] <- l1 + l2  + l3 + l4 + l5 
    # print ELBO value and difference with the previous one
    if (printL) { cat("Iter:\t", l, "\tELBO:\t", L[l], "\tELBO_diff:\t", L[l] - L[l-1], "\n")}
    # test if ELBO decreases
    if (L[l] < L[l - 1]) { message("Warning: ELBO decreases\n"); }
    # test convergence with epsilon threshold
    if (abs(L[l] - L[l - 1]) < epsilon) { break }
    # test VB needs more iterations to converge
    if (l == max_iterations) {warning("VB did not converge\n")}
  }
  object <- structure(list(mu=mu, S=S, delta=delta, phi=phi, r=r, pi_M=pi_M, L=L[2:l]))
  return(object)
} 
# pre-define max iterations for VB
iter <- 6000
model_dppoisreg <- vb_dppoismixreg(Y, X, M, mu, S, delta, phi, mu_0, S_0, phi_0,
                      max_iterations=iter, epsilon=1e-4, printL=FALSE)

## ELBO plot
# truncate the initial elbo iterations 
o<-2
Elbo <- data.frame(elbo=model_dppoisreg$L[o:length(model_dppoisreg$L)], iter=o:length(model_dppoisreg$L))

library(ggplot2)

pel <- ggplot(data=Elbo,aes(x=iter,y=elbo)) +xlab("iterations")+ylab("ELBO") + labs(title="ELBO of VB-DP-PoissonMix with Covariates") 
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
var_weights <- model_dppoisreg$pi_M[model_dppoisreg$pi_M>0.003]
cat("Variational weights:", var_weights, "\n")

var_par<- model_dppoisreg$mu[ ,,model_dppoisreg$pi_M>0.003]
cat("Variational parameters:\n");print(var_par)