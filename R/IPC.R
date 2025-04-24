#' The incremental principal component method can handle online data sets.
#'
#' @param data is an online data set
#' @param m is the number of principal component 
#' @param eta is the proportion of online data to total data
#'
#' @return T2,T2k,V,Vhat,lambdahat,time
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
#' IPC(data=data,m=m,eta=0.8) 
#' @export
IPC<-function(data,m,eta){
X<-scale(data)
S<-cov(X)   
n<-nrow(X)          
n0<-round(eta*n)  
p<-ncol(X)               
eig1<-eigen(cov(X[1:n0,]))
lambda<-eig1$values[1:m]      
V<-eig1$vectors[,1:m] 
V1<-V[,1:m]
H<-matrix(rep(0,(m+1)*(m+1)),nrow=(m+1))
T2<-matrix(rep(0,(n-n0)),nrow=1)
T2k<-matrix(rep(0,(n-n0)),nrow=1)
s=Sys.time() 
for (i in (n0+1):n) {
Xcenter<-t(X[i,])
g<-t(V)%*%t(Xcenter)
Xhat<-t(V%*%g)
h<-t(Xcenter)-V%*%g
hmao<-norm(h,"2")
gamma<-as.numeric(t(h/hmao)%*%t(Xcenter)) 
H[1:m,]<-cbind(((i-1)/i)*diag(lambda)+(1/i)*g%*%t(g),(1/i)*gamma*g)
H[(m+1),]<-cbind((1/i)*gamma*t(g),(1/i)*gamma^2)
eig2<-eigen(H)              
lambda<-eig2$values[1:m]               
V<-(cbind(V,h/hmao)%*%eig2$vectors)[,1:m]  
V2<-V[,1:m]
t<-Xcenter%*%V2
T2[,(i-n0)]<-(Xcenter%*%V2)%*%ginv(diag(lambda[1:m]))%*%t(Xcenter%*%V2)
T2k[,(i-n0)]<-((m*(i^2-1))/(i*(i-m)))*qf(0.001,m,i-m,lower.tail=FALSE)
}
e=Sys.time()
return(list(T2=T2,T2k=T2k,V=V1,Vhat=V2,lambdahat=lambda,time=(e-s)))
}


