##############----------------------------------------------------------------------------------------------
#### generate data from a Beta mixture model
# inputs: 
  ## N: number of observations
  ## prob: vector of the weights 
  ## D: number of dimensions
# output: 
  ## an NxD matrix of Beta Mixture Model

bmm_generator <- function(N, prob, D)
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
      u <- runif(1, min=1, max=100)
      v <- runif(1, min=1, max=100)
      data[(i-prob[m+1]*N +1):i, d] <- rbeta(prob[m+1]*N, u, v)
    }
  }
  # 0 values return error in the inferential scheme due to log calculations
  # little trick to avoid
  if (sum(data==0)!=0) {data[data==0]<-0.001}
  return(data)
}
##############----------------------------------------------------------------------------------------------
### examples
# cases of different weight vectors
case2 <- c(0.125,0.1,0.275,0.175,0.15,0.05,0.125)
case1 <- c(0.125,0.250,0.375,0.15,0.1)
case3 <- c(0.3,0.3,0.3,0.1)
#library(hitandrun)
#case4 <-  simplex.sample(50, 1)$samples
# set seed to reproduce the dataset
set.seed(123)
data <- bmm_generator(1000, sort(case1), D=10)

# contamination: in case we need to check the performance of inferential algos in contaminated scenarios
#c <- 5000 # number of contaminated dimensions

#N <- 300 # number of observations
#cont_cols <- matrix(NA, nrow=N, ncol=c)

#for(i in 1:c)
#  {
# unif <- runif(N, 0.001, 1)
# cont_cols[,i] <- unif
#   }

#set.seed(123)
#data <- cbind(data,cont_cols);data <- as.matrix(data) # bind contaminated columns to original dataset