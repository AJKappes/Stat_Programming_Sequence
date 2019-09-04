# starting with uniform(0,1) simulation
n <- 100
x <- runif(n)

# imlementing inverse method to simulate exponential
# 1. generate U(0,1)
# 2. x = -(1/lambda)*log(U)
n <- 10^4
U <- runif(n)
E <- -log(U) # exp(1) dist
hist(U)
hist(E)

# bernoulli(p) using inverse method
# 1. generate U(0,1)
# 2. x = 0 if U < 1-p and x = 1 if U > 1-p
p <- 0.25
U <- runif(n)
x <- ifelse(U < 1 - p, 0, 1)
X <- ifelse(U < 1 - p, 0, 1)
mean(X == 1) # gives p

# transformation method
# simulate gamma(alpha, beta) alpha --> integer
n <- 10^4
alpha <- 2
beta <- 1
U <- matrix(runif(alpha*n), nrow = alpha, ncol = n)
E <- -log(U) # exp(1)
G <- beta*colSums(E)
hist(G)

# make a function to simulate gammas
gam_sim <- function(n, alpha, beta){
  if (floor(alpha) != alpha) {
    stop('alpha parameter must be an integer')
  }
  
  U <- matrix(runif(alpha*n), alpha, n)
  E <- -log(U)
  return(beta*colSums(E))
}

hist(gam_sim(1000, 1, 0.2))

# simulate beta distribution using a transformation of uniforms
beta_sim <- function(n, a, b){
  U <- matrix(runif(n*(a+b)), a+b, n)
  E <- -log(U)
  numer <- colSums(matrix(E[1:a, ], nrow=a)) # creating the ratio of n_sim cols by rows a, b
  denom <- colSums(E)
  return(numer/denom)
}

hist(beta_sim(10, 1, 1))

# accept-reject algorithm 
# simulate gamma(alpha, beta) from gamma(a, b)
alpha <- 2.5
beta <- 2

# instrumental density gamma(a, b)
a <- floor(alpha)
b <- a*beta/alpha
M <- (beta^alpha)*(b^(-a))*(((alpha-a)/(beta-b))^(alpha-a))*exp(a-alpha)

bound <- function(M, x){
  z <- dgamma(x, shape = alpha, scale = (1/beta))/(M*dgamma(x, shape = a, scale = (1/b))) # AR ratio
  return(z)
}

# implementation of AR
# 1. simulate from instrument and U(0,1)
x <- gam_sim(1, a, b)
u <- runif(1)
y <- ifelse(u <= bound(M, x), x, NaN)

AR_gam_sim <- function(n){
  AR_vec <- c()
  for (i in 1:n){
    x <- gam_sim(1, a, b)
    u <- runif(1)
    y <- ifelse(u <= bound(M, x), x, NaN)
    AR_vec[i] <- y
  }
  return(AR_vec[!is.nan(AR_vec)])
}

AR_gam_sim(10)

# for a target sample we modify with while AR_vec < n continue
AR_gam_sim <- function(n){
  AR_vec <- c()
  while (length(AR_vec) < n){
    x <- gam_sim(1, a, b)
    u <- runif(1)
    y <- ifelse(u <= bound(M, x), x, NaN)
    if (!is.nan(y) == T){
      AR_vec <- c(AR_vec, y)
    }
  }
  return(AR_vec[!is.nan(AR_vec)])
}
