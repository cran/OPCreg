#' The stochastic principal  component  regression with varying learning-rate can handle online data sets.
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
#' 
#' library(MASS)
#' n=2000;p=20;m=9;
#' mu=t(matrix(rep(runif(p,0,1000),n),p,n))     
#' mu0=as.matrix(runif(m,0))
#' sigma0=diag(runif(m,1))
#'ro=as.matrix(c(runif(round(p/2),-1,-0.8),runif(p-round(p/2),0.8,1)))
#'R0=ro%*%t(ro);diag(R0)=1    
#'Sigma0=sigma0%*%t(sigma0)*R0 
#'x=mvrnorm(n,mu0,Sigma0)  
#'colnames(x)<-paste("x",1:p,sep="")
#'e=rnorm(n,0,1)
#'B=sample(1:3,(p+1),replace=T)
#'en<-matrix(rep(1,n*1),ncol=1)
#'y=cbind(en,x)%*%B+e
#'colnames(y)<-paste("y")
#'data<-data.frame(cbind(y,x))
#' spcrl(data=data,m=m,eta=0.8,alpha=0.5)
#' 
#' @importFrom stats reformulate
#' @importFrom stats sd
#' 
#' 
#' 


spcrl <- function(data, m, eta, alpha) {
  # Check if eta is within the valid range
  if (eta > 1 || eta <= 0) {
    stop("eta must be a proportion between 0 and 1.")
  }
  # Standardize predictors (excluding the response variable)
  X <- scale(data[, -1])  # Standardized predictors
  n <- nrow(X)            # Number of observations
  p <- ncol(X)            # Number of predictors
  
  # Calculate initial sample size based on eta
  n0 <- round(eta * n)
  if (n0 > n) {
    stop("n0 cannot be greater than the number of observations.")
  }
  
  # Initial covariance and eigen decomposition
  eig1 <- eigen(cov(X[1:n0, ]))
  lambda <- diag(eig1$values)  # Initial eigenvalues
  V <- eig1$vectors             # Initial eigenvectors
  
  # Record the start time
  s <- Sys.time()

# Iteratively update the principal components
for (i in (n0 + 1):n) { 
  Xcenter <- (X[i, , drop = FALSE])  # Current standardized observation
    
  # Update the principal components
  gamma <- i^(-alpha)  # Corrected learning rate
  V <- V + gamma * t(Xcenter) %*% Xcenter %*% V
  V <- qr.Q(qr(V))     # Orthogonalize
    
  # Update eigenvalues
  lambda <- lambda + gamma * (t(V) %*% t(Xcenter) %*% Xcenter %*% V)
}
  # Record the end time
  e <- Sys.time()
  
  # Post-processing
  lambda <- diag(lambda)  # Ensure diagonal form
  V2 <- V[, 1:m, drop = FALSE]  # Select first m components
  
  # Project data onto principal components
  Z <- X %*% V2
  colnames(Z) <- paste0("z", 1:m)
  data1 <- cbind(data, Z)
  
  # Fit linear model with dynamic formula
  formula <- reformulate(paste0("z", 1:m), response = "y")
  lm.sol <- lm(formula, data = data1)
  
 # Inverse transformation to obtain the regression coefficients in the original space
  A <- coef(lm.sol)
  xbar <- colMeans(data1[, 2:(p + 1)])  # Mean of the predictors
  xsd <- sqrt(diag(cov(data1[, 2:(p + 1)])))  # Standard deviation of the predictors
  b <- (A[2] * V[, 1] + A[3] * V[, 2]) / xsd  # Transform coefficients back to original space
  b0 <- A[1] - sum(xbar * b)  # Intercept term
  B <- c(b0, b)  # Combined coefficients
  
  # Calculate RMSE
  predicted <- cbind(1, as.matrix(data[, -1])) %*% B
  RMSE <- sqrt(mean((data$y - predicted)^2))
  
  # Return results
  list(
    Bhat = B,
    RMSE = RMSE,
    summary = summary(lm.sol),
    Vhat = V2,
    lambdahat = lambda[1:m],
    time = e - s,
    yhat = predicted
  )
}