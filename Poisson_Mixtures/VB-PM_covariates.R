##############----------------------------------------------------------------------------------------------
#### main VB function for the Poisson Mixture with covariates model + extra helpful functions
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
#### main VB function for DP Poisson Mix with covariates
# inputs: 
## Y: Nx1 data vector 
## X: NxD design matrix
## mu_0/S_0: Dx1 and DxD hyperparameters of the Gaussian upon beta coefficients (mean vector and covariance matrix respectively)
## phi_0: 1xM Dirichlet hyperparameter matrices
## phi: 1xM initial Dirichlet variational matrices
## mu/S: Dx1 and DxD initial variational parameters of the Gaussian upon beta coefficients 
## max_iterations: maximum number of VB iterations
## epsilon: threshold to achieve convergence
# output: 
## mu/S: Dx1 and DxD variational parameters of the Gaussian upon beta coefficients (mean vector and covariance matrix respectively)
## phi: 1xM Dirichlet variational matrices
## pi_M: variational weights
## r: NxD variational matrix for the latent allocation z
## L: ELBO values
## printL: print the ELBO values and difference with the previous one

vb_poismixreg <- function(Y, X, M, mu, S, phi, mu_0, S_0, phi_0,
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
    E.log_pi <- digamma(phi) - digamma(sum(phi))
    N.E.log_pi <- matrix(E.log_pi, nrow=N, ncol=M, byrow=TRUE)
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
    log_R <- N.E.log_pi + D.logy.E.xbeta - (1/2) *y.Exbeta2
    
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
    # term into phi variational parameter
    Ns <- colSums(r, na.rm=TRUE)
    # the Dirichlet variational parameter
    phi <- phi_0 + Ns
    # the variational weights
    pi_M <- (phi) / (M * phi_0 + N)
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
    E.log_pi <- digamma(phi) - digamma(sum(phi))
    N.E.log_pi <- matrix(E.log_pi, nrow=N, ncol=M, byrow=TRUE)
    # Dirichlet constant function
    logC <- function(u)
    {
      return(lgamma(sum(u)) - sum(lgamma(u)))
    }
    ## ELBO 
    for (k in 1:M)
    {
    for (c in 1:D)
      {
      l2_m[k,c] <- -(1/2) * ldet(S_0[,,c,k]) -(1/2) * t(mu[,c,k] - mu_0[,c,k]) %*% solve(S_0[,,c,k]) %*% (mu[,c,k] - mu_0[,c,k])
      -(1/2) * tr(solve(S_0[,,c,k]) %*% S[,,c,k]) + (1/2)*ldet(S[,,c,k])
      }
    }
    l1 <- sum( r * (-(1/2)*NM.logy -(1/2)*NM.ylogy + log_R )) -(D/2)*log(2*pi)
    l2 <- sum(l2_m)
    l3 <- -sum(r * log_r) + D*M*t/2
    l4 <- sum((phi_0 - phi) %*% (E.log_pi) + logC(phi_0) - logC(phi))

    # Total ELBO calculation
    L[l] <- l1 + l2  + l3 + l4 
    # print ELBO value and difference with the previous one
    if (printL) { cat("Iter:\t", l, "\tELBO:\t", L[l], "\tELBO_diff:\t", L[l] - L[l-1], "\n")}
    # test if ELBO decreases
    if (L[l] < L[l - 1]) { message("Warning: ELBO decreases\n"); }
    # test convergence with epsilon threshold
    if (abs(L[l] - L[l - 1]) < epsilon) { break }
    # test VB needs more iterations to converge
    if (l == max_iterations) {warning("VB did not converge\n")}
  }
  object <- structure(list(mu=mu, S=S, phi=phi, r=r, pi_M=pi_M, L=L[2:l]))
  return(object)
} 
# pre-define max iterations for VB
iter <- 6000
model_poisreg <- vb_poismixreg(Y, X, M, mu, S, phi, mu_0, S_0, phi_0,
                      max_iterations=iter, epsilon=1e-4, printL=TRUE)

## ELBO plot
# truncate the initial elbo iterations 
o<-2
Elbo <- data.frame(elbo=model_poisreg$L[o:length(model_poisreg$L)], iter=o:length(model_poisreg$L))

library(ggplot2)

pel <- ggplot(data=Elbo,aes(x=iter,y=elbo)) +xlab("iterations")+ylab("ELBO") + labs(title="ELBO of VB-PoissonMix with Covariates") 
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
var_weights <- model_poisreg$pi_M[model_poisreg$pi_M>0.003]
cat("Variational weights:", var_weights, "\n")

var_par<- model_poisreg$mu[ ,,model_poisreg$pi_M>0.003]
cat("Variational parameters:\n");print(var_par)