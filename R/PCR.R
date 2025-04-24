#' Principal Component Regression (PCR)
#'
#' The \code{PCR} function performs Principal Component Regression (PCR) on a given dataset.
#' It standardizes the predictor variables, determines the number of principal components to retain based on a specified threshold,
#' and fits a linear regression model using the principal component scores.
#'
#' @param data A data frame where the first column is the response variable and the remaining columns are predictor variables.
#' @param threshold The proportion of variance to retain in the principal components (default is 0.95).
#'
#' @return A list containing the following elements:
#' \item{Bhat}{The estimated regression coefficients, including the intercept.}
#' \item{RMSE}{The Root Mean Square Error of the regression model.}
#' \item{summary}{The summary of the linear regression model.}
#' \item{yhat}{The predicted values of the response variable.}
#' 
#' @details The function performs the following steps:
#' 1. Standardize the predictor variables.
#' 2. Compute the covariance matrix of the standardized predictors.
#' 3. Perform eigen decomposition on the covariance matrix to obtain principal components.
#' 4. Determine the number of principal components to retain based on the cumulative explained variance exceeding the specified threshold.
#' 5. Project the standardized predictors onto the retained principal components.
#' 6. Fit a linear regression model using the principal component scores.
#' 7. Back-transform the regression coefficients to the original scale.
#'
#' @examples
#' \dontrun{
#' # Example data
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
#' # Call the PCR function
#' result <- PCR(data, threshold = 0.9)
#'
#' # Access the estimated regression coefficients
#' print(Bhat <- result$Bhat)
#'
#' # Access the predicted values
#' print(yhat <- result$yhat)
#'
#' # Print the summary of the regression model
#' print(result$summary)
#'
#' # Print the RMSE
#' print(paste("RMSE:", result$RMSE))
#' }
#'
#' @seealso \code{\link{lm}}: For fitting linear models.
#' @seealso \code{\link{eigen}}: For computing eigenvalues and eigenvectors.
#'
#' @export
PCR <- function(data, threshold) {
  # Extract response and predictor variables
  y <- data[, 1]
  X <- data[, -1]
  
  # Standardize the independent variables
  X_scaled <- scale(X)  # Assume the first column is y
  p <- ncol(X_scaled)
  n <- nrow(X_scaled)
  S <- cov(X_scaled)
  eig <- eigen(S)
  
  # Calculate the number of components (m) based on the threshold
  sum_eig <- sum(eig$values)
  m <- 0
  for (i in 1:p) {
    if (sum(eig$values[1:i]) / sum_eig >= threshold) {
      m <- i
      break
    }
  }
  
  if (m == 0) {
    stop("No components meet the threshold criteria. Consider lowering the threshold.")
  }
  
  V <- eig$vectors[, 1:m]  # Principal components
  Z <- X_scaled %*% V      # Projected data
  
  # Fit linear model using the principal component scores
  lm.sol <- lm(y ~ Z - 1)  # Fit model without intercept in the Z space
  
  # Extract coefficients and transform back to original scale
  A <- coef(lm.sol)
  # B <- V %*% A  # Transform coefficients back to original scale
  B <- V %*% as.numeric(A)  # Ensure A is numeric
  # Calculate intercept
  xbar <- colMeans(X)
  b0 <- mean(y) - sum(xbar * B)
  Bhat <- c(b0, B)
  # Ensure Bhat is a numeric vector
Bhat <- as.numeric(Bhat)
# Ensure X is a numeric matrix
X <- as.matrix(X)
  # Calculate RMSE
  predictions <- cbind(1, X) %*% Bhat
  RMSE <- sqrt(mean((y - predictions)^2))
  
  return(list(Bhat = Bhat, RMSE = RMSE, summary = summary(lm.sol), yhat = predictions))
}