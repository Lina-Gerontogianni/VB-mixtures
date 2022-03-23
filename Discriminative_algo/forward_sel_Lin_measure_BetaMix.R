forward_selection <- function(x,m,M,pi_M,u,v,epsilon,a,hold=c()){
hnew<-1:dim(x)[2]
iterations<-rep(NA,a)

for(iter in 1:a){
  index <-c()
  A1<-c()
  for (i in hnew){
    A1<-c(A1,A(x,m,c(hold,i),M,pi_M,u,v))
    index <-c(index,i)
  }
  
  Alist=list(max_A=round(max(A1),3),features=sort(c(index[which.max(A1)],hold)))
  print(Alist)
  
  hnew <-hnew[hnew!=index[which.max(A1)]]
  hold<-c(hold,index[which.max(A1)])

  if (max(A1)>=epsilon) {break} 

}
return(Alist)
}