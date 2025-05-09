% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IPC.R
\name{IPC}
\alias{IPC}
\title{The incremental principal component method can handle online data sets.}
\usage{
IPC(data, m, eta)
}
\arguments{
\item{data}{is an online data set}

\item{m}{is the number of principal component}

\item{eta}{is the proportion of online data to total data}
}
\value{
T2,T2k,V,Vhat,lambdahat,time
}
\description{
The incremental principal component method can handle online data sets.
}
\examples{
library(MASS)
n=2000;p=20;m=9;
mu=t(matrix(rep(runif(p,0,1000),n),p,n))     
mu0=as.matrix(runif(m,0))
sigma0=diag(runif(m,1))
F=matrix(mvrnorm(n,mu0,sigma0),nrow=n)
A=matrix(runif(p*m,-1,1),nrow=p)
D=as.matrix(diag(rep(runif(p,0,1))))
epsilon=matrix(mvrnorm(n,rep(0,p),D),nrow=n)
data=mu+F\%*\%t(A)+epsilon
IPC(data=data,m=m,eta=0.8) 

}
