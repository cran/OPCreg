#' The stochastic principal component method can handle online data sets.
#'
#' @param data A data frame containing the response variable and predictors.
#'@param eta  proportion (between 0 and 1) determining the initial sample size for PCA. 
#'@param m The number of principal components to retain. 
#' @return A list containing the following elements:
#' \item{Bhat}{The estimated regression coefficients, including the intercept.}
#' \item{RMSE}{The Root Mean Square Error of the regression model.}
#' \item{summary}{The summary of the linear regression model.}
#' \item{yhat}{The predicted values of the response variable.}
#' 
#' @examples
#' # Example data
#' library(MASS);library(car);library(carData);library(stats);library(matrixcalc)
#' set.seed(1234)
#' n <- 2000
#' p <- 10
#' mu0 <- as.matrix(runif(p, 0))
#' sigma0 <- as.matrix(runif(p, 0, 10))
#' ro <- as.matrix(c(runif(round(p / 2), -1, -0.8), runif(p - round(p / 2), 0.8, 1)))
#' R0 <- ro %*% t(ro)
#' diag(R0) <- 1
#' Sigma0 <- sigma0 %*% t(sigma0) * R0
#' x <- mvrnorm(n, mu0, Sigma0)
#' colnames(x) <- paste("x", 1:p, sep = "")
#' e <- rnorm(n, 0, 1)
#' B <- sample(1:3, (p + 1), replace = TRUE)
#' en <- matrix(rep(1, n * 1), ncol = 1)
#' y <- cbind(en, x) %*% B + e
#' colnames(y) <- paste("y")
#' data <- data.frame(cbind(y, x))
#'
#' # Call the SPCR function
#' result <- SPCR(data,  eta = 0.0035, m = 3)
#' @export
SPCR <- function(data, eta, m) {
  # Check if eta is within the valid range
  if (eta > 1 || eta <= 0) {
    stop("eta must be a proportion between 0 and 1.")
  }
  
  # Standardize predictors (excluding the response variable)
  X <- scale(data[, -1])  # x
  n <- nrow(X)  # Number of observations
  p <- ncol(X)  # Number of predictors
  
  # Calculate initial sample size based on eta
  n0 <- round(eta * n)
  if (n0 > n) {
    stop("n0 cannot be greater than the number of observations.")
  }
  
  # Initial covariance and eigen decomposition
  eig1 <- eigen(cov(X[1:n0,]))  #
  lambda <- diag(eig1$values)  #
  V <- eig1$vectors  #
  
  # Record the start time
  s <- Sys.time()
  
  # Iteratively update the principal components
  for (i in (n0 + 1):n) { 
    Xcenter <- t(X[i,])  # Centered observation
    
    # Update the principal components using the online update rule
    gamma <- 1 / i  # Learning rate
    V <- V + gamma * t(Xcenter) %*% Xcenter %*% V
    V <- qr.Q(qr(V))  # Orthogonalize the principal components
    lambda <- lambda + gamma * (t(V) %*% t(Xcenter) %*% Xcenter %*% V)
  }
  
  # Record the end time
  e <- Sys.time()
  
  # Diagonalize the eigenvalues
  lambda <- diag(lambda)
  
  # Select the first m principal components
  V2 <- V[, 1:m]
  
  # Project the standardized predictors onto the principal components
  Z <- X %*% V2  #
  
  # Create a matrix to store the principal components
  z <- matrix(rep(0, n * m), ncol = m)
  colnames(z) <- paste("z", 1:m, sep = "")
  for (i in 1:m) {
    z[, i] <- Z[, i]
  }
  
  # Combine the original data with the principal components
  data1 <- cbind(data, z)
  
  # Fit a linear regression model using the principal components as predictors
  lm.sol.1 <- lm(y ~ z1 + z2 + z3, data = data1)
  
  # Inverse transformation to obtain the regression coefficients in the original space
  A <- coef(lm.sol.1)  # Coefficients of the regression model
  xbar <- colMeans(data1[, 2:(p + 1)])  # Mean of the predictors
  xsd <- sqrt(diag(cov(data1[, 2:(p + 1)])))  # Standard deviation of the predictors
  b <- (A[2] * V[, 1] + A[3] * V[, 2]) / xsd  # Transform coefficients back to original space
  b0 <- A[1] - sum(xbar * b)  # Intercept term
  B <- c(b0, b)  # Combined coefficients
  
  # Calculate the Root Mean Squared Error (RMSE)
  en <- matrix(rep(1, n * 1), ncol = 1)  # en = (1, ..., 1)T
  RMSE <- (sum((data$y - as.matrix(cbind(en, data1[, 2:(p + 1)])) %*% B)^2) / n)^(1/2)
  
  # Calculate the computation time
  computation_time <- e - s
  
  # Return the results
  return(list(
    Bhat = B,  # Estimated regression coefficients
    RMSE = RMSE,  # Root Mean Squared Error
    summary = summary(lm.sol.1),  # Summary of the linear regression model
    Vhat = V2,  # Estimated principal components
    lambdahat = lambda,  # Estimated eigenvalues
    time = computation_time,  # Computation time
    yhat = as.matrix(cbind(en, data1[, 2:(p + 1)])) %*% B  # Predicted values
  ))
}