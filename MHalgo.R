# MH algorithms
remove(list = objects())
nsim <- 15E+3

# simulate a gamma(3.5, 1)
# proposal is N(0,1)

# random walk and non-symmetric proposal

dgam <- function(x, alpha = 3.5, beta = 1) {
  
  y <- x^(alpha-1)*exp(-beta*x)
  
  if (x < 0){
    
    y <- 0
    
  }
  
  return(y)
  
}

proposal <- function(sym = FALSE) {
  
  if (sym == 'sym') {
    
    r <- rnorm(1)
    return(r)
    
  } else {
    
    r <- rexp(1)
    return(r)
  }
  
}

qdens <- function(y, x, sym = FALSE) {
  
  if (sym == 'sym') {
    
    # symmetric proposal cancels in test condition
    qdens <- 1
    return(qdens)
    
  } else {
    
    qdens <- exp(-x)/exp(-y)
    return(qdens)
  }
  
}

test_fun <- function(y, x, sym = FALSE) {
  
  tf <- dgam(y)/dgam(x)*qdens(y, x, sym)
  return(tf)
  
}

xn <- runif(1, 0, 3)
chain <- c()
for (i in 1:nsim) {
  
  yn <- xn + proposal(sym = 'sym')
  rho <- min(test_fun(yn, xn, sym = 'sym'), 1)
  b <- rbinom(1, 1, rho)
  chain[i] <- ifelse(b == 1, yn, xn)
  xn <- chain[i]

}
# check histogram
hist(chain)

# monte carlo approximation
mu_mh <- c()
for (i in 1:nsim) {
  
  mu_mh[i] <- mean(chain[1:i])
  
}

# convergence plot
plot(mu_mh, type = 'l')


