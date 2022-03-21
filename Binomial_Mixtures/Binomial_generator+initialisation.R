############################################### binomial simulation #######################################################
############################################################################################################################
##### general function to simulate binomial data
binm_generator <- function(N, prob, D, trial)
{
  M <- length(prob)
  data <- matrix(0, nrow=N, ncol=D)
  prob <- c(0,prob)
  i <- prob[1]
  
  for(m in 1:(M))
  {
    # index to assign the magnitude of observations within a group according to the weight values
    i <- i+prob[m+1]*N
    
    for(d in 1:D)
    {
      # Generate a random center for this cluster in this dimension
      u <- runif(1, min=0.001, max=0.99)
      data[(i-prob[m+1]*N +1):i, d] <- rbinom(prob[m+1]*N, trial, u)
    }
  }
  # 0 values return error in the inferential scheme due to log calculations
  # little trick to avoid
  #if (sum(data==0)!=0) {data[data==0]<-0.001}
  return(data)
}
# examples using this function: cases of different weight vectors
case2 <- c(0.125,0.1,0.275,0.175,0.15,0.05,0.125)
case1 <- c(0.125,0.250,0.375,0.15,0.1)
case3 <- c(0.3,0.3,0.3,0.1)
library(hitandrun)
case4 <-  simplex.sample(50, 1)$samples
# set seed to reproduce the dataset
set.seed(1234)
Y <- binm_generator(1000, sort(case1), D=5,60)
S <- matrix(60,nrow=dim(Y)[1],ncol=dim(Y)[2])
##################################################################
####### first example: binomial data with M=4, N=1000, D=1 ######
## binomial mix generator
rBinomMix = function(n, a, theta, size) {
  nindex = rmultinom(1, size=n, prob=a)
  rbinom(n, size, rep(theta, nindex))
}
## setting the desired values
N=1000
numTrials = 50
trueNumComponents = 4
phi = rep(1 / trueNumComponents, trueNumComponents)
theta = seq(0.2, 0.8, length = trueNumComponents)

set.seed(123)
data= as.matrix(rBinomMix(N, phi, theta, numTrials))
## the final binomial dataset
Y=data
## the number of trials for each observation
S=as.matrix(rep(numTrials,N))
##################################################################
####### second example: binomial data with M=2, N=1000, D=2 ######
yy<-rbinom(200,15,0.5)
yyy<-rbinom(800,15,0.8)
psi <- c(yy,yyy)
y<-rbinom(200,28,0.2)
y_<-rbinom(800,28,0.4)
psi_<-c(y,y_)
X<-cbind(psi,psi_)
l<-matrix(rep(c(15,28),each=1000),ncol=2,byrow=FALSE)
## the final binomial dataset
Y=X
## the number of trials for each observation
S=l
##################################################################
####### third example: binomial data with M=3, N=1000, D=1 ######
library(gtools)
N <- 1000
K <- 3
phi <- rdirichlet(1, rep(1/K, K))
theta <- rbeta(K, 1, 1)
### different number of trials for each datapoint
Nt <- sample(10:30,size=N,replace=TRUE)
Z <- sample(1:K, N, replace=TRUE, prob=phi)
y <- rbinom(N, Nt, theta[Z])
## the final binomial dataset
Y=as.matrix(y)
## the number of trials for each observation
S=as.matrix(Nt)
############################################################################################################################

############################################### bernoulli simulation #######################################################
############################################################################################################################
####### first example: bernoulli data with M=3, N=1000, D=5 ######
N=1000
P1 <- c(0.9, 0.9, 0.9, 0.1, 0.1)
P2 <- c(0.1, 0.1, 0.9, 0.9, 0.9)
P3 <- c(0.1, 0.9, 0.1, 0.9, 0.1)

pr <- list(P1, P2, P3)
weights <- c(0.25, 0.25, 0.5)

## the final binomial dataset
Y <- t(replicate(N, {
  ## pick a random prototype
  i <- sample(1:3, size = 1, prob =  weights)
  ## sample bits from the chosen prototype
  sapply(pr[[i]], function(p) rbinom(1, 1, p))
}))
## the number of trials for each observation (here all equal to 1)
S=matrix(1,nrow=dim(Y)[1],ncol=dim(Y)[2])
##################################################################
####### second example: binomial data with M=3, N=150, D=800K ######
D=800000
N=150
set.seed(123)
P1 <- runif(n=D, min=0.01, max=.92)
P2 <- runif(n=D, min=0.01, max=.92)
P3 <- runif(n=D, min=0.01, max=.92)

pr <- list(P1, P2, P3)
weights <- c(0.25, 0.25, 0.5)

## the final binomial dataset
Y <- t(replicate(N, {
  ## pick a random prototype
  i <- sample(1:3, size = 1, prob =  weights)
  ## sample bits from the chosen prototype
  sapply(pr[[i]], function(p) rbinom(1, 1, p))
}))
## the number of trials for each observation (here all equal to 1)
S=matrix(1,nrow=dim(Y)[1],ncol=dim(Y)[2])
############################################################################################################################

############################################# initialization (same for bernoulli/binomial) #################################
############################################################################################################################
# initial number of components
M<-10
####### for DP mix:
# hyperparameters:
phi_0<-rep(0.5,M)
a_0<- rep(0.5,M*dim(Y)[2]) + runif(M*dim(Y)[2]) / 1000
b_0<- rep(0.5,M*dim(Y)[2]) + runif(M*dim(Y)[2]) / 1000
a_0<- matrix(a_0,ncol=M)
b_0<- matrix(b_0,ncol=M)
# initial variational parameters:
delta <- phi_0 +runif(M) / 1000
phi<-phi_0+runif(M) / 1000
a<-a_0
b<-b_0

###### for finite mix:
# hyperparameters:
phi_0<-rep(1e-6,M)
a_0<- rep(0.5,M*dim(Y)[2]) + runif(M*dim(Y)[2]) / 1000
b_0<- rep(0.5,M*dim(Y)[2]) + runif(M*dim(Y)[2]) / 1000
a_0<- matrix(a_0,ncol=M)
b_0<- matrix(b_0,ncol=M)
# initial variational parameters:
phi <- phi_0
a<-a_0
b<-b_0
############################################################################################################################

############################################### logistic regression simulation #############################################
############################################################################################################################
####### first example: logistic data with M=3, N=1000, D=6 ######
## set the desired values
num=1
t=2 # the number of covariates
D=6
N=1000

beta1=matrix(c(1,2,2,2,1,1,2,1,1,1,2,1),nr=t,nc=D,byrow = TRUE)
beta2=matrix(c(10,8,9,8,9,10,-10,-9,-9.5,-8.9,-9,-10),nr=t,nc=D,byrow = TRUE)
beta3=matrix(c(-2,-4,-2,-3,-1.5,-2.2,-1.5,-2.5,-2,-3,-3,-2),nr=t,nc=D)

betas=list( beta1, beta2, beta3)
weights <- c(0.3,0.2,0.5)

set.seed(12)
x <- runif(N,0, 1) # Generating the covariates
X <- cbind(1,x)
Y <- matrix(NA, nc=D,nr=N)
for (n in 1:N)
  {
  ## pick a random beta
  i <- sample(1:3, size = 1, prob =  weights)
  ## the final logistic dataset
  Y[n,] <- sapply(plogis(X[n, ]%*%betas[[i]]), function(p) rbinom(1, num, p))
}
##################################################################
####### second example: logistic data with M=2, N=1000, D=2 ######
## set the desired values
num=1 
t=2 # the number of covariates
D=2
N=1000

beta1=matrix(c(1,3,4,2),nr=t,nc=D,byrow = TRUE)
beta2=matrix(c(-2,-2,-4,-3),nr=t,nc=D,byrow = TRUE)

betas=list( beta1, beta2)
weights <- c(0.25,0.75)

set.seed(12)
x <- runif(N,0, 1) # Generating the covariates
X <- cbind(1,x)
Y <- matrix(NA, nc=D,nr=N)

for (n in 1:N)
{
  ## pick a random beta
  i <- sample(1:2, size = 1, prob =  weights)
  ## the final logistic dataset
  Y[n,] <- sapply(plogis(X[n, ]%*%betas[[i]]), function(p) rbinom(1, num, p))
}

############################################################################################################################

##################################################### initialization (for logistic) ########################################
############################################################################################################################
# initial number of components
M<-10 
# hyperparameters:
mu_0=rep(0,M*dim(X)[2]*D) + runif(M*dim(X)[2]*D) / 100
mu_0=array(mu_0,c(t,D,M))
S_0=array(diag(100,t),c(t,t,D,M))
phi_0 <- rep(0.1,M)+ runif(M) / 100
# initial variational parameters:
mu=mu_0
S=S_0
E.omega= matrix(0.25,nr=N,nc=D) +runif(N)/100
a <- rep(0.1,M)+  runif(M) / 100
b <- rep(0.1,M) +  runif(M) / 100