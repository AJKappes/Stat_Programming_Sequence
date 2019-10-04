# solve x^2-2 = 0
# approximate the value of sqrt(2)

# newton-rhapson 
g <- function(x) {
  
  return(x^2 - 2)
  
}

gprime <- function(x) {
  
  return(2*x)
  
}

x_old <- 1
tol <- 1e-5
max_itr <- 3000
itr <- 1
error <- 1

while (error > tol & itr < max_itr) {
  
  x_new <- x_old - g(x_old)/gprime(x_old)
  error <- abs(x_new - x_old)
  itr <- itr + 1
  x_old <- x_new
  
}

for (i in 1:max_itr) {
  
  x_new <- x_old - g(x_old)/gprime(x_old)
  error <- abs(x_new - x_old)
  itr <- itr + 1
  x_old <- x_new
  
  if (error < tol) {
    
    break
    
  }
  
}

# newton rhapson for logistic regression

# generate data
n <- 300
x <- rnorm(n)
b <- c(0.25, 0.75)
X <- cbind(1, x)
invlogit <- function(X, b) {
  
  p <- exp(X %*% b)/(1 + exp(X %*% b))
  return(p)
  
}

p <- invlogit(X, b)
y <- c()
for (i in 1:n) {
  
  y[i] <- rbinom(1, 1, p[i])
  
}

# we've created the data set
# x features
# y labels
# estimate beta coefficients

tol <- 1e-4
error <- 1
itr <- 1
max_itr <- 4000
b_old <- c(0, 2)
p_mat <- matrix(0, length(p), length(p))

while(error > tol & itr < max_itr) {
  
  p_old <- invlogit(X, b_old)
  diag(p_mat) <- p_old*(1 - p_old)
  H_inv <- solve(-t(X) %*% p_mat %*% X)
  S <- t(X) %*% (y - p_old)
  b_new <- b_old  - H_inv %*% S
  error <- max(abs(b_new - b_old))
  itr <- itr + 1
  b_old <- b_new
  
}

b_new
