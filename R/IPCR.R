#' Incremental Principal Component Regression for Online Datasets
#'
#' The IPCR function implements an incremental Principal Component Regression (PCR) method designed to handle online datasets.
#' It updates the principal components recursively as new data arrives, making it suitable for real-time data processing.
#'
#' @param data A data frame where the first column is the response variable and the remaining columns are predictor variables.
#' @param eta The proportion of the initial sample size used to initialize the principal components (0 < eta < 1). Default is 0.0035.
#' @param m The number of principal components to retain. Default is 3.
#' @param alpha The significance level used for calculating critical values. Default is 0.05.
#'
#' @return A list containing the following elements:
#' \item{Bhat}{The estimated regression coefficients, including the intercept.}
#' \item{RMSE}{The Root Mean Square Error of the regression model.}
#' \item{summary}{The summary of the linear regression model.}
#' \item{yhat}{The predicted values of the response variable.}
#' 
#' @details The IPCR function performs the following steps:
#' 1. Standardizes the predictor variables.
#' 2. Initializes the principal components using the first \code{n0 = round(eta * n)} samples.
#' 3. Recursively updates the principal components as each new sample arrives.
#' 4. Fits a linear regression model using the principal component scores.
#' 5. Back-transforms the regression coefficients to the original scale.
#'
#' This method is particularly useful for datasets where new observations are continuously added, and the model needs to be updated incrementally.
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' library(MASS)
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
#' result <- IPCR(data = data, m = 3, eta = 0.0035, alpha = 0.05)
#' print(result$Bhat)
#' print(result$yhat)
#' print(result$RMSE)
#' print(result$summary)
#' }
#'
#' @seealso \code{\link{lm}}: For fitting linear models.
#' @seealso \code{\link{eigen}}: For computing eigenvalues and eigenvectors.
#'
#' @export
#' 
#' @importFrom MASS ginv
#' @importFrom stats cov qf
#' @importFrom stats lm coef
IPCR <- function(data, eta, m, alpha) {
  if (m > ncol(data) - 1) {
    stop("Number of components m cannot be greater than the number of predictors.")
  }
  
  n <- nrow(data)
  if (eta > 1 || eta <= 0) {
    stop("eta must be a proportion between 0 and 1.")
  }
  n0 <- round(eta * n)  # Calculate initial sample size based on eta
  if (n0 > n) {
    stop("n0 cannot be greater than the number of observations.")
  }
  
  X <- scale(data[, -1])  # Standardize predictors
  S <- cov(X)
  p <- ncol(X)
  
  # Initial covariance and eigen decomposition
  eig1 <- eigen(cov(X[1:n0,]))
  lambda <- eig1$values[1:m]
  V <- eig1$vectors[, 1:m]  # Initialize V
  
  H <- matrix(0, (m + 1), (m + 1))
  T2 <- matrix(0, (n - n0), 1)
  T2k <- matrix(0, (n - n0), 1)
  
  for (i in (n0 + 1):n) {
    Xcenter <- t(X[i,])  # x_i = x
    g <- t(V) %*% t(Xcenter)  # g = V_k^T X_{i+1}
    Xhat <- t(V %*% g)  # \hat{X}_{i+1}
    h <- t(Xcenter) - V %*% g  # h = X_{i+1} - \hat{X}_{i+1}
    hmao <- norm(h, "2")  # h_unit = h / ||h||_2
    gamma <- as.numeric(t(h / hmao) %*% t(Xcenter))  # \gamma = <X_{i+1}, h_unit>
    
    H[1:m, ] <- cbind(((i - 1) / i) * diag(lambda) + (1 / i) * g %*% t(g), (1 / i) * gamma * g)
    H[(m + 1), ] <- cbind((1 / i) * gamma * t(g), (1 / i) * gamma^2)
    
    eig2 <- eigen(H)  # H_k+1 = H_k + H_{i+1}
    lambda <- eig2$values[1:m]  # \lambda_{k+1}
    V <- (cbind(V, h / hmao) %*% eig2$vectors)[, 1:m]  # V_{k+1} = (V_k, h_unit) H_k+1
    t <- (t(V[, 1:m])) %*% t(t(X[i,]))
  }
  
  V2 <- V[, 1:m]
  Z <- X %*% V2
  z <- matrix(0, n * m, ncol = m)
  colnames(z) <- paste("z", 1:m, sep = "")
  for (i in 1:m) {
    z[, i] <- Z[, i]
  }
  data1 <- cbind(data, z)
  
  # Assuming 'y' is the first column of data
  lm.sol.1 <- lm(data1[, 1] ~ z1 + z2 + z3, data = data1)
  
  # Inverse transformation to original regression model
  A <- coef(lm.sol.1)
  xbar <- colMeans(data1[, 2:(p + 1)])
  xsd <- (diag(cov(data1[, 2:(p + 1)]))^(1/2))
  b <- (A[2] * V[, 1] + A[3] * V[, 2]) / xsd
  b0 <- A[1] - sum(xbar * b)
  Bhat <- c(b0, b)
  
  # Calculate RMSE
  predictions <- as.matrix(cbind(rep(1, n), data1[, 2:(p + 1)])) %*% Bhat
  RMSE <- sqrt(sum((data1[, 1] - predictions)^2) / n)
  
  return(list(Bhat = Bhat, RMSE = RMSE, summary = summary(lm.sol.1), yhat = predictions))
}