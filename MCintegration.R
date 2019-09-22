# lab 4

# monte carlo integration

#mc approximation of U(0,1) dist
nsim <- 500
u <- runif(nsim)
mc <- mean(u)
#gives mean of 0.5

mc <- c()
v <- c()
upper <- c()
lower <- c()
for (j in 1:nsim) {
  
  mc[j] <- mean(u[1:j])
  v[j] <- j^(-1)*var(u[1:j])
  upper[j] <- mc[j] + 1.96*sqrt(v[j])
  lower[j] <- mc[j] - 1.96*sqrt(v[j])
}
plot(mc, type = 'l')

library(ggplot2)
values <- c(mc, upper, lower)
type <- c(rep('mc', nsim), rep('upper', nsim), rep('lower', nsim))
iter <- rep(seq(1:nsim), 3)
data <- data.frame(val = values,
                   tp = type,
                   itr = iter)
ggplot(data, aes(itr, val)) +
  geom_line(aes(colour = factor(tp)))

# q1 approximate p(0<x<3), x~cauchy(0,1) using mc approximation
# result x1~n(0,1), x2~n(0,1), x1 and x2 are independent
# use ordinary monte carlo approximation sum_{j=1}^nsim 1[0<x<3]/nsim

# simulate cauchy from two standard normals and divide the standard normals for cauchy
# using box-mueller or marsaglia's method
I <- 10000
cauchy_sim <- function(I){
  
  # box-mueller n(0,1) sim
  x1_vec <- c()
  x2_vec <- c()
  
  for (i in 1:I) {
    
    u_vec <- runif(2)
    x1_vec[i] <- sqrt(-2*log(u_vec[1]))*cos(2*pi*u_vec[2])
    x2_vec[i] <- sqrt(-2*log(u_vec[1]))*sin(2*pi*u_vec[2])
                                            
  }
  cauchy_vec <- x1_vec / x2_vec
  return(cauchy_vec)
}
rv_cauchy <- cauchy_sim(I)
p_x <- length(rv_cauchy[rv_cauchy > 0 & rv_cauchy < 3]) / I

# or using an indicator function
ind <- function(x) {
  
  y <- ifelse((x > 0 & x < 3), 1, 0)
  return(y)
  
}

# showing both the length argument and indicator function provide the same results
p_x == sum(ind(rv_cauchy)) / I

# to evaluate convergence
mc <- c()
v <- c()
upper <- c()
lower <- c()
for (j in 1:nsim) {
  
  mc[j] <- mean(ind(rv_cauchy[1:j]))
  v[j] <- j^(-1)*var(ind(rv_cauchy[1:j]))
  upper[j] <- mc[j] + 1.96*sqrt(v[j])
  lower[j] <- mc[j] - 1.96*sqrt(v[j])
  
}

values <- c(mc, upper, lower)
type <- c(rep('mc', nsim), rep('upper', nsim), rep('lower', nsim))
iter <- rep(seq(1:nsim), 3)
data <- data.frame(val = values,
                   tp = type,
                   itr = iter)
ggplot(data, aes(itr, val)) +
  geom_line(aes(colour = factor(tp)))

# p = p(0<x<3) unknown but we have 10000 relizations

# 95% probability the unkown number p will lie within the CI
# make 100 different simulations, 95 of them will contain the real value p within the bounds

# q2 do the same as q1 except using importance sampling
# instrumental density is u(0,3)
# cauchy f(x) = 1/(pi*(1+x^2))
# h(x) = indcator[0<x<3]
# g(x) = 1/3
# IS estimator p approx 1/n sum_{i=1}^n h(x)f(x)/g(x) = 1/n sum_{i=1}^n 3*f(x)
# h(x_i) = 1, always between 0 and 3 from U(0,3)

cauch_dens <- function(x) {
  
  y <- 1/(pi*(1 + x^2))
  return(y)
}

nsim <- 500
x <- runif(nsim, 0, 3)
IS_mc <- mean(3*cauch_dens(x))