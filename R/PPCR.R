#' Perturbation-based Principal Component Regression
#'
#' This function performs Perturbation-based Principal Component Regression (PPCR) on the provided dataset.
#' It combines Principal Component Analysis (PCA) with linear regression, incorporating perturbation to enhance robustness.
#'
#' @param data A data frame containing the response variable and predictors.
#' @param eta A proportion (between 0 and 1) determining the initial sample size for PCA.
#' @param m The number of principal components to retain.
#' @param alpha Significance level (currently not used in the function).
#' @param perturbation_factor A factor controlling the magnitude of perturbation added to the principal components.
#'
#' @return A list containing the following components:
#' \item{Bhat}{Estimated regression coefficients in the original space.}
#' \item{RMSE}{Root Mean Squared Error of the regression model.}
#' \item{summary}{Summary of the linear regression model.}
#' \item{Vhat}{Estimated principal components.}
#' \item{lambdahat}{Estimated eigenvalues.}
#' \item{yhat}{Predicted values from the regression model.}
#'
#' @details The function first standardizes the predictors, then performs PCA on an initial subset of the data.
#' It iteratively updates the principal components by incorporating new observations and adding random perturbations.
#' Finally, it fits a linear regression model using the principal components as predictors and transforms the coefficients back to the original space.
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
#' # Call the PPCR function
#' result <- PPCR(data, eta = 0.0035, m = 3, alpha = 0.05, perturbation_factor = 0.1)
#'
#' # Print results
#' print(result$Bhat)  # Estimated regression coefficients
#' print(result$RMSE)  # RMSE of the model
#' print(result$summary)  # Summary of the regression model
#' }
#'
#' @seealso \code{\link{lm}}: For linear regression models.
#' @seealso \code{\link{prcomp}}: For principal component analysis.
#'
#' @keywords regression multivariate
#'
#' @export
#' 
#' @importFrom stats rnorm lm
PPCR <- function(data, eta = 0.0035, m = 3, alpha = 0.05, perturbation_factor = 0.1) {
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
  lambda <- eig1$values[1:m]  # Eigenvalues
  V <- eig1$vectors[, 1:m]  # Eigenvectors (first m components)
  
  # Initialize matrices
  H <- matrix(0, (m + 1), (m + 1))
  T2 <- matrix(0, (n - n0), 1)
  T2k <- matrix(0, (n - n0), 1)
  
  # Iteratively update the principal components with perturbation
  for (i in (n0 + 1):n) { 
    Xcenter <- t(X[i,])  # x_i = x
    g <- t(V) %*% t(Xcenter)  # g = V_k^T X_{i+1}
    Xhat <- t(V %*% g)  # \hat{X}_{i+1}
    h <- t(Xcenter) - V %*% g  # h = X_{i+1} - \hat{X}_{i+1}
    hmao <- norm(h, "2")  # h_unit = h / ||h||_2
    gamma <- as.numeric(t(h / hmao) %*% t(Xcenter))  # \gamma = <X_{i+1}, h_unit>
    
    # Introduce perturbation
    perturbation <- rnorm(length(h), mean = 0, sd = perturbation_factor * hmao)
    h <- h + perturbation  # Perturbed h
    
    H[1:m,] <- cbind(((i - 1) / i) * diag(lambda) + (1 / i) * g %*% t(g), (1 / i) * gamma * g) 
    H[(m + 1),] <- cbind((1 / i) * gamma * t(g), (1 / i) * gamma^2)
    eig2 <- eigen(H)  # H_k+1 = H_k + H_{i+1}
    lambda <- eig2$values[1:m]  # \lambda_{k+1}
    V <- (cbind(V, h / hmao) %*% eig2$vectors)[, 1:m]  # V_{k+1} = (V_k, h_unit) H_k+1
    t <- (t(V[, 1:m])) %*% t(t(X[i,]))
  }
  
  # Order the eigenvalues and corresponding eigenvectors in descending order
  ind <- order(lambda, decreasing = TRUE)
  lambda <- lambda[ind]
  V <- V[, ind]
  
  # Select the first m principal components
  V2 <- V[, 1:m]
  
  # Project the standardized predictors onto the principal components
  Z <- X %*% V2
  
  # Create a matrix to store the principal components
  z <- matrix(0, n * m, ncol = m)
  colnames(z) <- paste("z", 1:m, sep = "")
  for (i in 1:m) {
    z[, i] <- Z[, i]
  }
  
  # Combine the original data with the principal components
  data1 <- cbind(data, z)
  
  # Fit a linear regression model using the principal components as predictors
  y <- data[, 1]  # Assume the first column is y
  lm.sol.1 <- lm(y ~ ., data = data1)
  
  # Inverse transformation to obtain the regression coefficients in the original space
  A <- coef(lm.sol.1)  # Coefficients of the regression model
  xbar <- colMeans(data1[, 2:(p + 1)])  # Mean of the predictors
  xsd <- sqrt(diag(cov(data1[, 2:(p + 1)])))  # Standard deviation of the predictors
  b <- (A[2] * V[, 1] + A[3] * V[, 2]) / xsd  # Transform coefficients back to original space
  b0 <- A[1] - sum(xbar * b)  # Intercept term
  Bhat <- c(b0, b)  # Combined coefficients
  
  # Calculate the Root Mean Squared Error (RMSE)
  predictions <- as.matrix(cbind(rep(1, n), data1[, 2:(p + 1)])) %*% Bhat
  RMSE <- sqrt(sum((y - predictions)^2) / n)
  
  # Return the results
  return(list(
    Bhat = Bhat,  # Estimated regression coefficients
    RMSE = RMSE,  # Root Mean Squared Error
    summary = summary(lm.sol.1),  # Summary of the linear regression model
    Vhat = V2,  # Estimated principal components
    lambdahat = lambda,  # Estimated eigenvalues
    yhat = predictions  # Predicted values
  ))
}