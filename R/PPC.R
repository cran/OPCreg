#' The perturbation principal component method can handle online data sets.
#'
#' @param data is an online data set
#' @param m is the number of principal component 
#' @param eta is the proportion of online data to total data
#'
#' @return T2,T2k,V,Vhat,lambdahat,time
#' @export
#'
#' @examples
#' library(MASS)
#' n=2000;p=20;m=9;
#' mu=t(matrix(rep(runif(p,0,1000),n),p,n))     
#' mu0=as.matrix(runif(m,0))
#' sigma0=diag(runif(m,1))
#' F=matrix(mvrnorm(n,mu0,sigma0),nrow=n)
#' A=matrix(runif(p*m,-1,1),nrow=p)
#' D=as.matrix(diag(rep(runif(p,0,1))))
#' epsilon=matrix(mvrnorm(n,rep(0,p),D),nrow=n)
#' data=mu+F%*%t(A)+epsilon
#' PPC(data=data,m=m,eta=0.8) 
PPC<-function(data,m,eta){
X<-scale(data)
S<-cov(X)   
n<-nrow(X)          
n0<-round(eta*n)  
p<-ncol(X)               
eig<-eigen(cov(X[1:n0,])) 
s=Sys.time()  
lambda<-eig$values      
V<-eig$vectors  
V1<-V[,1:m]   
T2<-matrix(rep(0,(n-n0)),nrow=1)
T2k<-matrix(rep(0,(n-n0)),nrow=1)
s=Sys.time()
for (i in (n0+1):n) {
f<-1/i 
Xcenter<-t(X[i,])                    
lambda<-(1-f)*lambda        
Q<-sqrt(f)*t(X[i,]%*%V)  
Q2<-Q*Q                      
num<-tcrossprod(Q)              
den<-matrix(lambda+Q2,p,p,byrow=T)-matrix(Q2+lambda,p,p)  
U<-num/den       
diag(U)<-1         
V<-V%*%U      
sigma2<-.colSums(V*V,p,p)
lambda<-(lambda+Q2)*sigma2
V<-V*rep.int(1/sqrt(sigma2),rep.int(p,p))
ind<-order(lambda,decreasing=T)  
lambda<-lambda[ind]
V<-V[,ind]
V2<-V[,1:m]
t<-Xcenter%*%V2
T2[,(i-n0)]<-(Xcenter%*%V2)%*%ginv(diag(lambda[1:m]))%*%t(Xcenter%*%V2)
T2k[,(i-n0)]<-((m*(i^2-1))/(i*(i-m)))*qf(0.001,m,i-m,lower.tail=FALSE)
}
e=Sys.time()
return(list(T2=T2,T2k=T2k,V=V1,Vhat=V2,lambdahat=lambda,time=(e-s)))
}

