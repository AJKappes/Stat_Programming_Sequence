# EM algorith is used to approximate parameters of a normal mixture model
# (pie, mu1, mu2, sig1, sig2)

remove(list = objects())
n <- 500

# parameter target values 
p <- 0.4
mu1 <- -2
mu2 <- 2
sig1 <- 0.25
sig2 <- 0.25

# generate data
z1 <- rbinom(n, 1, p)
z2 <- 1 - z1
x1 <- rnorm(n, mu1, sig1)
x2 <- rnorm(n, mu2, sig2)
x <- z1*x1 + z2*x2

# EM algorithm
p_old <- 0.5
mu1_old <- -1
mu2_old <- 1
sig1_old <- 1
sig2_old <- 1

f <- function(x, mu, sig) {
  
  fx <- (1/(2*pi)*sig)*exp(-1/2*(x-mu)^2/sig^2)
  return(fx)
  
}

err <- 1
tol <- 1e-3
itr <- 0
max_itr <- 1000
while (err > tol) {
  
  ez1_old <- f(x, mu1_old, sig1_old)*p_old / (f(x, mu1_old, sig1_old)*p_old +
                                                f(x, mu2_old, sig2_old)*(1-p_old))
  ez2_old <- 1 - ez1_old
  
  p_new <- sum(ez1_old) / n
  
  mu1_new <- sum(ez1_old * x) / sum(ez1_old)
  mu2_new <- sum(ez2_old * x) / sum(ez2_old)
  
  sig1_new <- sum(ez1_old*(x - mu1_new)^2) / sum(ez1_old)
  sig2_new <- sum(ez2_old*(x - mu2_new)^2) / sum(ez2_old)
  
  err_vec <- abs(c(p_new - p_old, mu1_new - mu1_old, mu2_new - mu2_old,
                   sig1_new - sig1_old, sig2_new - sig2_old))
  
  err <- max(err_vec)
  cat('Iteration:', itr, 'Error:', err, '\n')
  itr <- itr + 1
  
  p_old <- p_new
  mu1_old <- mu1_new
  mu2_old <- mu2_new
  sig1_old <- sig1_new
  sig2_old <- sig2_new

}

