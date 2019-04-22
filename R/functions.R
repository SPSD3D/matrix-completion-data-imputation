hnk<- function(b,l) 
{
  
  #b=hankel
  #b=as.vector(b)
  #=50
  
  n = length(b)
  m = n-l+1;
  x <- matrix(0,m,l)
  for (i in 1:m){
  x[i,]=b[i:(i+l-1)];
  }
x
}


#Average
ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}


#Median
mmed <- function(x,n=5){runmed(x,n)} 