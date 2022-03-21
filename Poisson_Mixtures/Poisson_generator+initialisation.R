############################################### poisson simulation #######################################################
############################################################################################################################
##### general function to simulate poisson data
poismix_generator <- function(N, prob, D)
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
      u <- rgamma(1, 4,0.5)
      data[(i-prob[m+1]*N +1):i, d] <- rpois(prob[m+1]*N, u)
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
# the dataset
Y <- poismix_generator(1000, sort(case3), D=4)
##################################################################
####### first example: poisson data with M=3, N=1000, D=2 #######
N=1000
D=2
P1 <- c(10,10)
P2 <- c(18,20)
P3 <- c(4,3)

pr <- list(P1, P2, P3)
weights <- c(0.25, 0.25, 0.5)

## the final binomial dataset
Y <- t(replicate(N, {
  ## pick a random prototype
  i <- sample(1:3, size = 1, prob =  weights)
  ## sample bits from the chosen prototype
  sapply(pr[[i]], function(p) rpois(1, p))
}))
# the dataset
Y=matrix(Y,ncol=D)
##################################################################
####### second example: real data with  M=2 (0.65, 0.35), D=1 ####
library(MixtureInf)
data(earthquake)
# the dataset
Y<- earthquake
##################################################################
####### third example: poisson data with M=4, N=1000, D=1 #######
yy<-rpois(100,35)
y<-rpois(400,20)
yyy<-rpois(200,12)
g<-rpois(300,4)
x<-c(yy,y,yyy,g)
# the dataset
Y <-as.matrix(x)
##################################################################
####### fourth example: poisson data with M=2, N=1000, D=2 #######
yy<-rpois(200,5)
yyy<-rpois(800,10)
psi <- c(yy,yyy)
y<-rpois(200,7)
y_<-rpois(800,9)
psi_<-c(y,y_)
X<-cbind(psi,psi_)
# the dataset
Y=X

############################################# initialization (same for bernoulli/binomial) #################################
############################################################################################################################
# initial number of components
M<-10 
set.seed(234)
####### for DP mix:
# hyperparameters:
phi_0 = 0.1
a_0 = 1
b_0 = 1
phi_0<-rep(phi_0,M) + runif(M) / 1000
a_0<- rep(a_0,M*dim(Y)[2]) + runif(M*dim(Y)[2]) / 1000
b_0<- rep(b_0,M*dim(Y)[2]) + runif(M*dim(Y)[2]) / 1000
a_0<- matrix(a_0,ncol=M)
b_0<- matrix(b_0,ncol=M)
# initial variational parameters:
delta<-rep(1,M)+ runif(M) / 1000
phi<-rep(1,M)+ runif(M) / 1000
a <- a_0
b <-b_0
a<-rgamma(M,1,100) + a_0
b<-rgamma(M,1,100) +b_0

###### for finite mix:
# hyperparameters:
phi_0 = 0.1
a_0 = 1
b_0 = 1
phi_0<-rep(phi_0,M) + runif(M) / 1000
a_0<- rep(a_0,M*dim(Y)[2]) + runif(M*dim(Y)[2]) / 1000
b_0<- rep(b_0,M*dim(Y)[2]) + runif(M*dim(Y)[2]) / 1000
a_0<- matrix(a_0,ncol=M)
b_0<- matrix(b_0,ncol=M)
# initial variational parameters:
phi<-phi_0+ runif(M) / 1000
a <- a_0
b <-b_0
a<-rgamma(M,1,100) + a_0
b<-rgamma(M,1,100) +b_0