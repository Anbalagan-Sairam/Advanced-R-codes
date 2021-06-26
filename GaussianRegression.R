# This project implements gaussian fit a polynomial regression and gaussian process regression with parameters
# optimised using BFGS optimization.


# Clean up the current environment
rm(list=ls())

# Make results reproducible
set.seed(12345)

library(mvtnorm)

# Define criterion to be minimised in Gaussian process regression
gp_criterion = function(p,x,y) {
  sig_sq = exp(p[1])
  rho_sq = exp(p[2])
  tau_sq = exp(p[3])
  Mu = rep(0, length(x))
  Sigma = sig_sq * exp( - rho_sq * outer(x, x, '-')^2 ) + tau_sq * diag(length(x))
  ll = dmvnorm(y, Mu, Sigma, log = TRUE)
  return(-ll)
}

# Implement the function regression_fit here
# We create a function regression fit with x_g, x, y, p and method as input parameters

regression_fit = function(x_g,x,y,p = 1,method = "BFGS"){
  #Creating the design matrix
  x_rep = matrix(rep(x, p+1), ncol = p+1, nrow = length(x))
  X = sweep(x_rep, 2, 0:p, '^')
  X_g_rep = matrix(rep(x_g, p+1), ncol = p+1, nrow = length(x_g))
  X_g = sweep(X_g_rep, 2, 0:p, '^')
  
  #Fitting the polynomial model
  mod = lsfit(X, y, intercept = FALSE)
  
  #pred_lm contains the predicted output which can be returned for plotting the fitted curve.
  pred_lm = X_g%*%mod$coefficients
  
  # The optimal values of rho_sq, sig_sq and tau_sq are obtained using optim function
  answer_BFGS = optim(c(0,0,0),gp_criterion,x = x,y =y, method = 'BFGS')
  
  rho_sq = answer_BFGS$par[1]
  sig_sq = answer_BFGS$par[2]
  tau_sq = answer_BFGS$par[3]
  C = sig_sq * exp( - rho_sq * outer(x_g, x, '-')^2 ) #= outer: outer product
  
  # Next Sigma = sig_sq * exp( - rho_sq * (x[i] - x[j])^2 ) + diag(tau_sq)
  Sigma = sig_sq * exp( - rho_sq * outer(x, x, '-')^2 ) + tau_sq * diag(length(x))
  
  # pred_gauss  contains the predicted gaussian model output which can be returned for plotting the fitted curve.
  pred_gauss = C %*% solve(Sigma, y)
  
  return(list("gaussian" = pred_gauss,"polynomial" = pred_lm))
  }

# Import the data and create the plot here

# Importing prostate dataset
prostate = read.csv('prostate.csv',header = TRUE)

# Checking first five elements of the dataset
head(prostate)

# Scaling lcavol and storing it in variable x
x = scale(prostate$lcavol,)[,1]
x

# Scaling lcavol and storing it in variable y
y = scale(prostate$lpsa,)[,1]
y

# Create a grid of new x-values for predictors
x_g = pretty(x, n = 100)
x_g

pred = regression_fit(x_g, x, y, 7, 'BFGS')

plot(x, y, xlab = "lcavol", ylab = "lpsa", main = "Regression models")
lines(x_g, pred$gaussian, col = 'blue') 
lines(x_g, pred$polynomial, col = 'red')
legend(x = "topleft",legend = c("Gaussian", "Polynomial p = 7"),bty = "n",col=c("blue","red"), lty = c(1, 1))