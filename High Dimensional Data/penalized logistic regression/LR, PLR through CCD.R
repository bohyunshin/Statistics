rm(list=ls())
#####################################################################################
############################ 1. Setup & Data Preparing ##############################
#####################################################################################
library(glmnet)
library(glue)
library(ggplot2)
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


#####################################################################################
######################### 2. Logistic Regression via CCD ############################
#####################################################################################
threshold = 1e-5
lambda = 0
# initial beta
beta0 = rep(0,p)

# logistic regression with no penalty via CCD
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

# compare result with glmnet with no penalty
fit = glmnet(X[,-10], y, family="binomial", lambda=0)
coef(fit)


#####################################################################################
#################### 3. Penalized Logistic Regression via CCD #######################
#####################################################################################
# soft-thresholding function
S = function(z,t){
  if (z > 0 & t < abs(z)) return(z-t)
  else if (z < 0 & t < abs(z)) return(z+t)
  else if (t >= abs(z)) return(0)
}

Penalized_LR = function(beta, lambda, alpha, epsilon, maxiter){
  iter = 0
  prevloglik = loglik = 1e6
  result = c(iter, loglik, beta)
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
  return(result)
}

# function to calculate l1norm
beta_l1norm = function(result){
  # first two elements are iterations and likelihood
  # last element is estimates of intercept
  # so, exclude them
  return( sum(result[nrow(result),][3:(ncol(result)-1)]) )
}

# function to extract converged coefficients
converged_coef = function(result){
  #estimates without intercept
  return(result[nrow(result),3:(ncol(result)-1)])
}

# compare with glmnet
beta = rep(0,ncol(X))
alpha = 1
lambda = 0.05
epsilon = 1e-5
threshold = 0.00001
maxiter = 300
test = converged_coef( Penalized_LR(beta, lambda, alpha, epsilon, maxiter)  )
names(test) = colnames(X)[-10]
fit = glmnet(X,y,family='binomial',lambda=lambda)
coef(fit)
test

#####################################################################################
########### 4. Store L1 norm and estimated coefficients varying lambdas #############
#####################################################################################

grid = seq(0,0.15,0.005)
alpha = 1
epsilon = 1e-5
threshold = 0.00001
maxiter = 300

result = rep(1,11)
cat(glue('Evaluated at {length(grid)} lambdas'))
for (lam in grid){
  beta = rep(0,ncol(X))
  PLR_result = Penalized_LR(beta,lam,alpha,epsilon,maxiter)
  result = rbind(result, c(lam, beta_l1norm(PLR_result), converged_coef(PLR_result)))
  cat(glue('{lam} completed'))
  cat('\n')
}
colnames(result) = c('lambda', 'l1norm', colnames(X)[-10])
result = data.frame(result[-1,])

#####################################################################################
########################## 5. Visualizing solution path #############################
#####################################################################################

ggplot(data = result, aes(x=l1norm, y=sbp, colour='black')) + geom_line() +
  geom_line(aes(x=l1norm, y=tobacco, colour = 'red')) + 
  geom_line(aes(x=l1norm, y=ldl, colour='blue')) + 
  geom_line(aes(x=l1norm, y=adiposity , colour='green')) + 
  geom_line(aes(x=l1norm, y=typea, colour='pink')) + 
  geom_line(aes(x=l1norm, y=obesity, colour='orange')) + 
  geom_line(aes(x=l1norm, y=alcohol, colour='brown')) + 
  geom_line(aes(x=l1norm, y=age, colour='purple')) + 
  geom_line(aes(x=l1norm, y=famhist, colour='gray')) + 
  labs(y='estimated coefficients', x='l1norm') + 
  ggtitle("Figure 4.13") +
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5)) + 
  scale_color_discrete(name = '',labels = colnames(X)[-10])
