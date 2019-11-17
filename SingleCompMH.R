# implementation of single component MH
remove(list = objects())

# beta-binomial example
# proposals binomial(n, 0.5) and U(0,1)
alpha <- 2
beta <- 1
n <- 5
p <- 0.5
nsim <- 2E+4

fx_y <- function(x, y) {
  
  f <- y^x*(1 - y)^(n - x)
  return(f)
  
}

fy_x <- function(y, x) {
  
  f <- y^(x + alpha - 1)*(1 - y)^(n - x + beta - 1)
  return(f)
  
}

q1 <- function(x) {
  
  q <- p^x*(1 - p)^(n - x)
  return(q)
  
}

q2 <- 1

x_chain <- c()
y_chain <- c()
x <- 1
y <- 0.5
for (i in 1:nsim) {
  
  z1 <- rbinom(1, n, p)
  z2 <- runif(1)
  
  rho1 <- min(fx_y(z1, y)/fx_y(x, y), 1)
  rho2 <- min(fy_x(z2, x)/fy_x(y, x), 1)
  
  b1 <- rbinom(1, 1, rho1)
  b2 <- rbinom(1, 1, rho2)
  
  x_chain[i] <- ifelse(b1 == 1, z1, x)
  y_chain[i] <- ifelse(b2 == 1, z2, y)
  
  x <- x_chain[i]
  y <- y_chain[i]
  
}

# hist(x_chain)
# hist(y_chain)

# change point problem
# priors gamma(a1, b1), gamma(a2, b2), discreteU()
# proposals gamma, gamma, discrete unif
# a1 = a2 = a, b1 = b2 = b

a <- 10
b <- 1
n <- 150
k <- 50
mu <- 4
lambda <- 15
x <- c(rpois(k, mu), rpois((n-k), lambda))
n <- length(x)

rdunif <- function(nsim, a, b) {
  
  f <- floor(runif(nsim, a, b + 1))
  return(f)
  
}

u <- c()
v <- c()
w <- c()
u[1] <- 1
v[1] <- 2
w[1] <- round(n/2)

for (j in 2:10000) {
  
  # first component mu update
  u[j] <- rgamma(1,
                 shape = a + sum(x[1:w[j - 1]]),
                 rate = w[j - 1] + b)
  
  # second component lambda update
  v[j] <- rgamma(1,
                 shape = a + sum(x[-(1:w[j - 1])]),
                 rate = n - w[j - 1] + b)
  
  # update for third component
  # proposal
  q <- rdunif(1,1,n)
  cr <- w[j-1]
  mu <- u[j-1]
  la <- v[j-1]
  sterm1 <- sum(x[1:q])-sum(x[1:cr]) 
  sterm2 <- sum(x[-(1:q)])-sum(x[-(1:cr)])
  sterm3 <- (q - cr)*mu + (cr - q)*la
  logp1 <- sterm1*log(mu) + sterm2*log(la) - sterm3
  rho <- min(exp(logp1), 1)
  binom <- rbinom(1, 1, rho)
  w[j] <- ifelse(binom == 1, q, w[j-1])
  
}

chain <- cbind(u,v,w)
chain_burned <- chain[-(1:5000),]

cat('Change point parameter estimes (mu, lambda, k):',
    '\n(', colMeans(chain_burned), ')')

# alternative way for change point period k
# frequintist approach
# any arbitrary starting point
sp <- round(n/2)
imu1 <- mean(x[1:sp])
imu2 <- mean(x[-(1:sp)])

t <- 10
loss <- c()
for (j in 1:n-1) {
  
  t <- j
  loss[j] <- sum((x[1:t] - imu1)^2) + sum((x[-(1:t)] - imu2)^2)

}

plot(loss, type = 'l')
which.min(loss)
