rm(list=ls())

library(glmnet)

data = read.csv('data.txt')
X = data[,c(2:5,7:10)]
X = scale(X)
intercept = 1
famhist = as.numeric(as.factor(data$famhist))-1
X = cbind(X,  famhist, intercept)
X = data.matrix(X)



# standardize design matrix
y = data$chd
sigmoid = function(x) return((1+exp(-x))^(-1))
prob = function(beta, x) return( sigmoid(as.vector(beta) %*% as.vector(x)) )
p = ncol(X) 
n = nrow(X) 

soft_threshold = function(rho, lambda, q){
  if (rho < -lambda/2) return( (rho + lambda/2)/q )
  else if (rho > lambda/2) return((rho - lambda/2)/q)
  else return(0)
}

threshold = 1e-5
lambda = 0

# initial beta
beta0 = rep(0,p)

# logistic regression
repeat{
  pi = rep(1,n)
  for (i in 1:n){
    pi[i] = prob(beta0, X[i,]) 
  }
  bi = pi * (1-pi)
  W = diag(bi)
  beta1 = beta0 - solve(t(X) %*% W %*% X) %*% t(X) %*% (pi-y)
  if ( sum((beta1 - beta0)^2) < threshold  ) break
  else beta0 = beta1
  print(beta0)
}

# compare result with glmnet
fit = glmnet(X[,-10], y, family="binomial", lambda=0)
coef(fit)


# logistic regression with lasso penalty using CCD
# initial beta
S = function(z,t){
  if (z > 0 & t < abs(z)) return(z-t)
  else if (z < 0 & t < abs(z)) return(z+t)
  else if (t >= abs(z)) return(0)
}

soft_threshold = function(beta, lambda) {
  ifelse(abs(beta)>lambda && beta > 0,
         beta-lambda, 
         ifelse(abs(beta) > lambda && beta < 0, 
                beta + lambda,
                0))}

beta = rep(0,ncol(X))
lambda = 0.010088
alpha = 1
epsilon = 1e-5
iter = 0
prevloglik = loglik = 1e6
result = c(iter, loglik, beta)
threshold = 0.00001
maxiter = 300

repeat{
  iter = iter + 1
  for (j in 1:ncol(X)){
    XBeta = X %*% beta
    p = 1/(1+exp(-XBeta))
    weight = p * (1-p)
    
    # to avoid divergence of coefficients
    weight = ifelse(abs(weight-0) < epsilon, epsilon, weight)
    z = XBeta + (y-p)/weight
    fitted_without_j = X[,-j] %*% beta[-j]
    
    beta[j] = S( mean(weight * X[,j] * (z - fitted_without_j)), lambda * alpha  ) / (mean(weight * X[,j]^2) + lambda*(1-alpha))
  }
  
  #calculate new log-likelihood
  loglik = 1/(2*nrow(X))*sum(weight*(z-X%*%beta)^2) + lambda*sum(abs(beta))
  
  if (abs(loglik - prevloglik) < threshold) break
  else {
    preloglik = loglik
    result = rbind(result, c(iter, loglik, beta))
  }
  if (iter > maxiter) break
}
tail(result)
sum(abs(result[nrow(result),][3:12]))

fit = cv.glmnet(X,y,family='binomial')
fit = glmnet(X,y,family='binomial',lambda=fit$lambda.min)
coef(fit)
