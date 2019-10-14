# simulated annealing

# cauchy distribution example for mle of alpha in C(alpha, beta)
# beta known, beta = 0.1

# observed data
x <- c(-4.8, -2.8, -1.35, -0.02, 0.7, 0.98, 2.92, 5.5)

loglike_cauch <- function(alpha) {
  
  y <- -sum(log(0.1^2 + (x - alpha)^2))
  return(y)
  
}

grd <- seq(-6, 6, length.out = 500)
y <- sapply(as.list(al), loglike_cauch)
plot(al, y)# showing local and global optimums

# implement simulated annealing to find the mle of alpha
set.seed(123)
r <- 0.5 # user chosen

# simulated annealing optimization

n <- 15E+3
Tt <- function(t) {
  
  y <- 1 / log(1 + t)
  return(y)
  
}

alpha_old <- 0
updates <- c()
for (k in 1:n) {

  alpha_cand <- runif(1, alpha_old - r, alpha_old + r)
  
  f_cand <- loglike_cauch(alpha_cand)
  f_old <- loglike_cauch(alpha_old)
  
  
  if (f_cand > f_old) {
    
    alpha_new <- alpha_cand
    
  } else {
    
    rho <- exp((f_cand - f_old) / Tt(k))
    b <- rbinom(1, 1, rho)
    
    if (b == 1) {
      
      alpha_new <- alpha_cand
      
    }
    
    if (b == 0) {
      
      alpha_new <- alpha_old
      
    }
    
  }
  
  updates <- c(updates, alpha_new)
  alpha_old <- alpha_new
  
}
plot(1:n, updates, type = 'l')
alpha_new

# find the optimizer of the function h(x) = (cos(50x) + sin(50x))^2
# x > 0

h <- function(x) {
  
  f <- (cos(50*x) + sin(50*x))^2
  return(f)
  
}

grd <- seq(0, 10, length.out = 100)
plot(grd, h(grd), type = 'l')

r <- 0.5
temp <- function(t) {
  
  f <- 1 / (1 + t)^2
  return(f)
}

n <- 15E+3
x_old <- 0
updates <- c()
for (k in 1:n) {
  
  x_cand <- runif(1, max(0, x_old - r), x_old + r)
  h_cand <- h(x_cand)
  h_old <- h(x_old)
  
  if (h_cand > h_old) {
    
    x_new <- x_cand
    
  } else {
    
    rho <- exp((h_cand - h_old) / temp(k))
    b <- rbinom(1, 1, rho)
    
    if (b == 1) {
      
      x_new <- x_cand
      
    }
    
    if (b == 0) {
      
      x_new <- x_old
      
    }
    
  }
  
  updates <- c(updates, x_new)
  x_old <- x_new
  
}
plot(1:n, updates, type = 'l')

# we know the global optimum of the function is 2
max_val <- 2
acc <- abs(h(x_new) - max_val)
cat('Global optimum estimator accuracy is', acc)
