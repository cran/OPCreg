#' The stochastic approximate component method can handle online data sets.
#'
#' @param data is a online data set
#' @param m is the number of principal component 
#' @param eta is the proportion of online data to total data
#' @param alpha is the step size
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
#' SAPC(data=data,m=m,eta=0.8,alpha=1) 
SAPC<-function(data,m,eta,alpha){
X<-scale(data)
n<-nrow(X)          
n0<-round(eta*n) 
p<-ncol(X)           
eig1<-eigen(cov(X[1:n0,])) 
lambda<-diag(eig1$values)
V<-eig1$vectors
V1<-V[,1:m]
T2<-matrix(rep(0,(n-n0)),nrow=1)
T2k<-matrix(rep(0,(n-n0)),nrow=1)
s=Sys.time()
for (i in (n0+1):n) {
Xcenter<-t(X[i,])
gamma<-i^(-alpha)
V<-V+gamma*t(Xcenter)%*%Xcenter%*%V
V<-qr.Q(qr(V))
V2<-V[,1:m]
lambda<-lambda+gamma*(t(V)%*%t(Xcenter)%*%Xcenter%*%V)
t<-Xcenter%*%V2
T2[,(i-n0)]<-(Xcenter%*%V2)%*%ginv(diag(lambda[1:m]))%*%t(Xcenter%*%V2)
T2k[,(i-n0)]<-((m*(i^2-1))/(i*(i-m)))*qf(0.001,m,i-m,lower.tail=FALSE)
}
lambda<-diag(lambda)
e=Sys.time()
return(list(T2=T2,T2k=T2k,V=V1,Vhat=V2,lambdahat=lambda,time=(e-s)))
}

