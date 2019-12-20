# Bayesian Addative Regression Trees #
# Chipman, George, McCulloch (1998) -> CGM98
remove(list = objects())
library(tidyverse)

#### Generate data ----------------------------------
# data generation following Hill (2010)
set.seed(77)
p <- 0.5
z <- rbinom(500, 1, p)
y <- c()
x <- c()
for (i in 1:length(z)) {
  
  if (z[i] == 1) {
    
    x[i] <- rnorm(1, 40, 10)
    y[i] <- rnorm(1, 90 + exp(0.06*x[i]), 1)
    
  } else {
    
    x[i] <- rnorm(1, 20, 10)
    y[i] <- rnorm(1, 72 + 3*sqrt(x[i]), 1)
    
  }
    
}

#set.seed(NULL)

# resale y to (-.5, .5) interval
ridx <- which(is.na(y))
y <- y[-ridx]
y <- (y - min(y))/(max(y) - min(y)) - 0.5
org_y <- y
data <- data.frame(y = y,
                   z = z[-ridx],
                   x = x[-ridx],
                   org_y = org_y)
remove(x, y, z, org_y)

#### Priors and likelihood -------------------------
# back to CGM98

# tree prior
p_split <- function(d, a = alpha, b = beta) {
  
  s <- a*(1 + d)^(-beta)
  # rv <- rbinom(1, 1, s)
  return(s)
  
}

# leaf param prior (mu)
mus_draw <- function(n_tn, s2) {
  
  mus <- rnorm(n_tn, sd = sqrt(s2))
  return(mus)
  
}

# leaf param hyper
sig2_draw <- function(v = nu, lam = lambda) {
  
  s2 <- 1/rgamma(1, v/2, scale = v*lam/2)
  return(s2)
  
}

# log-likelihood 
loglik <- function(n_vec, y_list, mu, n_tn, c = a, v = nu, lam = lambda) {
  
  if (length(which(n_vec == 0)) > 0) {
    
    j <- which(n_vec == 0)
    y_list <- y_list[-j]
    n_vec <- n_vec[-j]
    n_tn <- n_tn - 1
    
  }
  
  s_vec <- map2_dbl(y_list, n_vec,
                    function(y, n) (n - 1)*var(y))
  # 1 obs in y_list -> NA for var(y), define as o
  s_vec[is.na(s_vec)] <- 0
  
  t_vec <- map2_dbl(y_list, n_vec,
                    function(y, n) ((n*c)/(n + c))*(mean(y) - mu)^2)
  
  f <- -sum(n_vec/2)*log(pi) + n_tn*v/2*log(v*lam) + n_tn/2*log(c) -
    1/2*sum(log(n_vec + c)) + sum(log(gamma((n_vec + v)/2))) - 
    n_tn*log(gamma(v/2)) - sum((n_vec + v)/2)*sum(log(s_vec + t_vec + v*lam))
  
  return(f)
  
}

# Tree acceptance draw
rho_draw <- function(k_grow, k_prune, Lmp1, Lm, tdraw) {
  
  if (Lmp1 < 0 & Lm > 0) {
    
    l_rat <- 0
    
  } else if (Lmp1 > 0 & Lm < 0) {
    
    l_rat <- 1
    
  } else if (Lmp1 == Inf & Lm == Inf) {
    
    l_rat <- 1
    
  } else {
    
    l_rat <- Lmp1/Lm
    
  }
  
  if (tdraw == 'grow') {
    
    rho <- min(l_rat*k_prune, 1)
    
  } else if (tdraw == 'prune') {
    
    rho <- min(l_rat*1/k_grow, 1)
    
  } else {
    
    rho <- min(l_rat, 1)
    
  }
  
  bdraw <- rbinom(1, 1, rho)
  return(bdraw)
  
}

#### Tree action setup -----------------------------------

# initialize z partition nodes and split value
nodes_zl <- list()
nodes_zr <- list()
split_z <- 0.5

# initialize z_left and z_right nodes for > than conditions (go right)
nodes_zlr <- list(NA)
nodes_zrr <- list(NA)

# initialize split value vectors
splits_zl <- c(split_z)
splits_zr <- c(split_z)

# initialize leaf size constraint
min_leaf_obs <- 3

# initialize leaf parameter list
leaf_params <- list()

# initiate tree random action
tree_step <- c('grow', 'prune', 'change', 'swap')
draw_action <- function() {
  
  # kernel probabilities
  # g <- p_split(d)
  # p <- (1 - g)/3
  g_kern <- .25
  p_kern <- g_kern
  c <- .4
  s <- .1
  s <- sample(tree_step, 1, prob = c(g_kern, p_kern, c, s))
  out <- list(g_kern = g_kern, p_kern = p_kern, draw = s)
  return(out)
  
}

#### Tree action functions --------------------------------

# if action is grow
t_grow <- function() {
  
  # draw initial z node direction, < 0.5 or > 0.5
  s_direct <- sample(c('l', 'r'), 1)
  inaction <- FALSE
  
  if(length(nodes_zl) == 0 | length(nodes_zr) == 0) {
    
    # define intial partition on z 
    nodes_zl[[1]] <<- data %>%
      filter(z < split_z) %>%
      .$y
    nodes_zr[[1]] <<- data %>%
      filter(!y %in% nodes_zl[[1]]) %>%
      .$y
    
    # define initial depth
    i <- 1
    
  } else {
    
    # ifelse and if_else only handles conditions of length 1
    if (s_direct == 'l') {
      
      nodes <- nodes_zl
      nodes_r <- nodes_zlr
      splits <- splits_zl
      
    } else {
      
      nodes <- nodes_zr
      nodes_r <- nodes_zrr
      splits <- splits_zr
      
    }
    
    # new depth length to grow on
    i <- length(nodes) + 1
    
    # define uniform splitting rultes
    if (length(nodes) == 1) {
      
      splits[i] <- sample(data$x, 1)
      
    } else {
      
      splits[i] <- data$x %>%
        .[between(., min(.), splits[i - 1])] %>%
        sample(1)
      
    }
    
    split <- splits[i]
    
    # partition data at current split point for data in z <> .05
    split_y <- data %>%
      filter(x < split & y %in% nodes[[1]]) %>%
      .$y
    
    # check if minimum leaf obs is met
    if(length(split_y) >= min_leaf_obs) {
      
      # assign data to new node
      nodes[[i]] <- split_y
      
      # parition compliment data depending on node length
      if (length(nodes) <= 2) {
        
        split_comp_y <- data %>%
          filter(!y %in% nodes[[i]] &
                   y %in% nodes[[1]]) %>%
          .$y
        
      } else {
        
        split_comp_y <- data %>%
          filter(!y %in% unlist(nodes_r) &
                   !y %in% nodes[[i]] &
                   y %in% nodes[[1]]) %>%
          .$y
        
      }
      
      # assign data to compliment node
      nodes_r[[i]] <- split_comp_y
      
    } else {
      
      inaction <- TRUE
      splits <- splits[-length(splits)]
      cat('Cannot partition.',
          '\nMinimum leaf obs of', min_leaf_obs, 'not met.',
          '\n')
      
    }
    
    if (s_direct == 'l') {
      
      nodes_zl <<- nodes
      nodes_zlr <<- nodes_r
      splits_zl <<- splits
      
    } else {
      
      nodes_zr <<- nodes
      nodes_zrr <<- nodes_r
      splits_zr <<- splits
      
    }
    
  }
  
  d <- ifelse(s_direct == 'l', length(nodes_zl), length(nodes_zr))
  out <- list(left_nodes = nodes_zl, right_nodes = nodes_zr,
              left_comps = nodes_zlr, right_comps = nodes_zrr,
              left_splits = splits_zl, right_splits = splits_zr,
              depth = d, direc = s_direct, inaction = inaction)
  return(out)
  
}

# if action is prune
t_prune <- function() {
  
  s_direct <- sample(c('l', 'r'), 1)
  inaction <- FALSE
  
  if(length(nodes_zl) == 1  & s_direct == 'l') {
    
    inaction <- TRUE
    cat('Cannot prune', s_direct, 'any further.',
        '\n')
    
  } else if (length(nodes_zr) == 1  & s_direct == 'r') {
    
    inaction <- TRUE
    cat('Cannot prune', s_direct, 'any further.',
        '\n')
    
  } else {
    
    if (s_direct == 'l') {
      
      i <- length(nodes_zl)
      nodes_zl <<- nodes_zl[-i]
      nodes_zlr <<- nodes_zlr[-i]
      splits_zl <<- splits_zl[-i]
      
    } else {
      
      i <- length(nodes_zr)
      nodes_zr <<- nodes_zr[-i]
      nodes_zrr <<- nodes_zrr[-i]
      splits_zr <<- splits_zr[-i]
      
    }
    
  }
  
  d <- ifelse(s_direct == 'l', length(nodes_zl), length(nodes_zr))
  out <- list(left_nodes = nodes_zl, right_nodes = nodes_zr,
              left_comps = nodes_zlr, right_comps = nodes_zrr,
              left_splits = splits_zl, right_splits = splits_zr,
              depth = d, direc = s_direct, inaction = inaction)
  return(out)
  
}

# if action is change
t_change <- function() {
  
  s_direct <- sample(c('l', 'r'), 1)
  inaction <- FALSE
  
  if(length(nodes_zl) == 1 & s_direct == 'l') {
    
    inaction <- TRUE
    cat('Cannot change splitting rule for direction:', s_direct,
        '\nCurrent node is the initial node.',
        '\n')
    
  } else if (length(nodes_zr) == 1 & s_direct == 'r') {
    
    inaction <- TRUE
    cat('Cannot change splitting rule for direction:', s_direct,
        '\nCurrent node is the initial node.',
        '\n')
    
  } else {
    
    if (s_direct == 'l') {
      
      nodes <- nodes_zl
      nodes_r <- nodes_zlr
      splits <- splits_zl
      
    } else {
      
      nodes <- nodes_zr
      nodes_r <- nodes_zrr
      splits <- splits_zr
      
    }
    
    # random node draw conditioned on random side (z <> 0.5)
    # specify change conditions based on node length
    if (length(nodes) > 2) {
      
      i <- sample(2:length(nodes), 1)
      
    } else {
      
      i <- 2
      
    }
    
    if (i == length(nodes)) {
      
      nodes <- nodes
      nodes_r <- nodes_r
      splits <- splits
      
    } else {
      
      nodes <- nodes[-c(i + 1:length(nodes))]
      nodes_r <- nodes_r[-c(i + 1:length(nodes_r))]
      splits <- splits[-c(i + 1:length(splits))]
      
    }
    
    # store current split
    old_split <- splits[i]
    
    # check split vector length to define new split
    if (i > 2) {
      
      splits[i] <- data$x %>%
        .[between(., min(.), splits[i - 1])] %>%
        sample(1)
      
    } else {
      
      splits[i] <- sample(data$x, 1)
      
    }
    
    split <- splits[i]
    
    # partition data at current split point for data in z <> .05
    split_y <- data %>%
      filter(x < split & y %in% nodes[[1]]) %>%
      .$y
    
    # check minimum leaf obs condition
    if(length(split_y) >= min_leaf_obs) {
      
      # assign data to new random node
      nodes[[i]] <- split_y
      
      # parition compliment data depending on node length
      if (length(nodes) <= 2) {
        
        split_comp_y <- data %>%
          filter(!y %in% nodes[[i]] &
                   y %in% nodes[[1]]) %>%
          .$y
        
      } else {
        
        split_comp_y <- data %>%
          filter(!y %in% unlist(nodes_r) &
                   !y %in% nodes[[i]] &
                   y %in% nodes[[1]]) %>%
          .$y
        
      }
      
      # assign data to new random compliment node
      nodes_r[[i]] <- split_comp_y
      
    } else {
      
      # iterate through new splitting rule until meeting obs const
      itr <- 1
      max_itr <- 20
      while (length(split_y) < min_leaf_obs & itr <= max_itr) {
        
        if (i > 2) {
          
          splits[i] <- data$x %>%
            .[between(., min(.), splits[i - 1])] %>%
            sample(1)
          
        } else {
          
          splits[i] <- sample(data$x, 1)
          
        }
        
        split <- splits[i]
        
        # partition data at current split point for data in z <> .05
        split_y <- data %>%
          filter(x < split & y %in% nodes[[1]]) %>%
          .$y
        
        cat('Finding valid splitting rule, iteration:', itr,
            '\n')
        itr <- itr + 1
        
        if (itr > max_itr) {
          
          # restore old split and partitioned data
          splits[i] <- old_split
          inaction <- TRUE
          cat('Max splitting iterations reached. Original split kept.',
              '\n')
          
        }
        
      }
      
      # partition data based on splitting rule result
      split <- splits[i]
      split_y <- data %>%
        filter(x < split & y %in% nodes[[1]]) %>%
        .$y
      nodes[[i]] <- split_y
      
      # parition compliment data depending on node length
      if (length(nodes) <= 2) {
        
        split_comp_y <- data %>%
          filter(!y %in% nodes[[i]] &
                   y %in% nodes[[1]]) %>%
          .$y
        
      } else {
        
        split_comp_y <- data %>%
          filter(!y %in% unlist(nodes_r) &
                   !y %in% nodes[[i]] &
                   y %in% nodes[[1]]) %>%
          .$y
        
      }
      
    }
    
    if (s_direct == 'l') {
      
      nodes_zl <<- nodes
      nodes_zlr <<- nodes_r
      splits_zl <<- splits
      
    } else {
      
      nodes_zr <<- nodes
      nodes_zrr <<- nodes_r
      splits_zr <<- splits
      
    }
      
  }
  
  d <- ifelse(s_direct == 'l', length(nodes_zl), length(nodes_zr))
  out <- list(left_nodes = nodes_zl, right_nodes = nodes_zr,
              left_comps = nodes_zlr, right_comps = nodes_zrr,
              left_splits = splits_zl, right_splits = splits_zr,
              depth = d, direc = s_direct, inaction = inaction)
  return(out)
  
}

# if action is swap
t_swap <- function() {
  
  s_direct <- sample(c('l', 'r'), 1)
  inaction <- FALSE
  
  if (length(nodes_zl) == 1 & s_direct == 'l') {
    
    inaction <- TRUE
    cat('Cannot swap splitting rule for direction:', s_direct,
        '\nParent node is the initial node.',
        '\n')
    
  } else if (length(nodes_zr) == 1 & s_direct == 'r') {
    
    inaction <- TRUE
    cat('Cannot swap splitting rule direction:', s_direct,
        '\nParent node is the initial node.',
        '\n')
    
  } else {
    
    if (s_direct == 'l') {
      
      nodes <- nodes_zl
      nodes_r <- nodes_zlr
      splits <- splits_zl
      
    } else {
      
      nodes <- nodes_zr
      nodes_r <- nodes_zrr
      splits <- splits_zr
      
    }
    
    # conditions for random parent node draw
    if (length(nodes) > 2) {
      
      i <- sample(2:length(nodes), 1)
      cat('Sample node', i, 'in length', length(nodes),
          '\n')
      
    } else {
      
      i <- 2
      
    }
    
    if (i == length(nodes)) {
      
      inaction <- TRUE
      cat('There is no child node to swap with. Original nodes kept.',
          '\n')
      
    } else {
      
      # save child data to map to parent node (given symmetric child split rule)
      child_data <- nodes[[i + 1]]
      # subset nodes to i - 1 before parent to define child split rule on parent
      nodes <- nodes[-c(i:length(nodes))]
      nodes_r <- nodes_r[-c(i:length(nodes_r))]
      # swap splitting rule of parent i and child i + 1
      split <- splits[i + 1]
      splits <- splits[-c(i:length(splits))]
      splits[i] <- split
      
      # swap parent and child node data based on child swap rule
      nodes[[i]] <- child_data
      
      # partition new compliment data based on new parent node
      if (length(nodes) <= 2) {
        
        split_comp_y <- data %>%
          filter(!y %in% nodes[[i]] &
                   y %in% nodes[[1]]) %>%
          .$y
        
      } else {
        
        split_comp_y <- data %>%
          filter(!y %in% unlist(nodes_r) &
                   !y %in% nodes[[i]] &
                   y %in% nodes[[1]]) %>%
          .$y
        
      }
      
      # assign new compliment data based on child split
      nodes_r[[i]] <- split_comp_y
      
    }
    
    if (s_direct == 'l') {
      
      nodes_zl <<- nodes
      nodes_zlr <<- nodes_r
      splits_zl <<- splits
      
    } else {
      
      nodes_zr <<- nodes
      nodes_zrr <<- nodes_r
      splits_zr <<- splits
      
    }
    
  }
  
  d <- ifelse(s_direct == 'l', length(nodes_zl), length(nodes_zr))
  out <- list(left_nodes = nodes_zl, right_nodes = nodes_zr,
            left_comps = nodes_zlr, right_comps = nodes_zrr,
            left_splits = splits_zl, right_splits = splits_zr,
            depth = d, direc = s_direct, inaction = inaction)
  return(out)
  
}

# get tree parameters needed for MH iters
get_treedata <- function(i) {
  
  # get tree information
  all_y <- c(nodes_zl[[1]],  nodes_zr[[1]])
  y_list <- c(nodes_zl[length(nodes_zl)], nodes_zlr[-1],
              nodes_zr[length(nodes_zr)], nodes_zrr[-1])
  n_tnl <- 1 + length(nodes_zlr[-1])
  n_tnr <- 1 + length(nodes_zrr[-1])
  n_tns <- n_tnl + n_tnr
  n_vec <- sapply(y_list, length)
  l_splits <- splits_zl
  r_splits <- splits_zr
  params <- mus_draw(n_tns, sig2_mu)
  mu <- mean(params)
  
  # set summed params for resids iter
  leaf_params[i, 1:n_tns] <<- params
  J_psum <- colSums(leaf_params)
  p_sum <- J_psum[1:n_tns]
  
  # get residual information
  r_list <- map2(y_list, p_sum, function(l, p) l - p)
  
  out <- list(all_y = all_y, y_list = y_list, n_tns = n_tns,
              n_tnl = n_tnl, n_tnr = n_tnr, n_vec = n_vec,
              l_splits = l_splits, r_splits = r_splits,
              params = params, p_sum = p_sum, mu = mu,
              r_list = r_list)
  return(out)
  
}

# set residuals for tree_mp1 iteration
set_resdata <- function() {
  
  # set y data to match residual index
  yidx <- match(all_y, data$y)
  data$y[yidx] <<- unlist(r_list)
  
  # replace initial node with residual values
  nodes_zl[[1]] <<- data %>% filter(z < split_z) %>% .$y
  nodes_zr[[1]] <<- data %>% filter(z > split_z) %>% .$y
  #nodes_zl[[1]] <<- unlist(r_list[1:length(nodes_zl)])
  #nodes_zr[[1]] <<- unlist(r_list[(length(nodes_zl) + 1):length(r_list)])
  
  #replace < terminal nodes with residual values
  len_ln <- length(nodes_zl)
  if (len_ln > 1) {

    nodes_zl[[len_ln]] <<- r_list[[1]]

  }

  if (length(nodes_zr) > 1) {

    nodes_zr[[length(nodes_zr)]] <<- r_list[[len_ln + 1]]

  }
  
  # replace terminal compliment nodes with residuals
  len_lnr <- length(nodes_zlr)
  len_rnr <- length(nodes_zrr)
  if (len_lnr > 1) {

    nodes_zlr[2:len_lnr] <<- r_list[2:len_lnr]

  }

  if (len_rnr > 1) {

    nodes_zrr[2:len_rnr] <<- r_list[(len_lnr + 2):length(r_list)]

  }
  
}

#### Generate Bayesian trees --------------------------------------

# choosen parameters -> CGM98 (alpha, beta, nu, lambda, a)
# res <- resid(lm(y ~ x, data = data))
# alpha <- 0.95
# beta <- 2
# nu <- 3
# lambda <- unname(quantile(res, 0.9))
# a <- var(res)/var(data$y)
# 
# # define initial tree
# tree_m <- t_grow()
# tree_sims <- list(tree_m)
# 
# # initialize likelihood from base tree
# y_list <- get_treedata()$ys
# n_tn <- get_treedata()$n_tns
# n_vec <- get_treedata()$ns
# sig2 <- sig2_draw()
# param_draw <- mus_draw(n_tn, mean(data$y), sig2)
# mu <- mean(param_draw)
# leaf_params <- list(mu)
# loglik_m <- loglik(n_vec, y_list, mu, n_tn)
# marg_ll <- 0 
# 
# mh_n <- 5000
# for (i in 2:mh_n) {
#   
#   # draw action
#   action <- draw_action(tree_m$depth)
#   draw <- action$draw
#   
#   if (draw == 'grow') {
#     
#     tree_mp1 <- t_grow()
#     g_kern <- draw_action(tree_mp1$depth)$g_kern
#     p_kern <- draw_action(tree_mp1$depth)$p_kern
#     
#   } else if (draw == 'prune') {
#     
#     tree_mp1 <- t_prune()
#     g_kern <- draw_action(tree_mp1$depth)$g_kern
#     p_kern <- draw_action(tree_mp1$depth)$p_kern
#     
#   } else if (draw == 'change') {
#     
#     tree_mp1 <- t_change()
#     
#   } else {
#     
#     tree_mp1 <- t_swap()
#     
#   }
#   
#   y_list <- get_treedata()$ys
#   n_tn <- get_treedata()$n_tns
#   n_vec <- get_treedata()$ns
#   loglik_mp1 <- loglik(n_vec, y_list, mu, n_tn)
#   rho <- rho_draw(g_kern, p_kern, loglik_mp1, loglik_m, draw)
#   
#   if (rho == 1) {
#     
#     tree_sims[[i]] <- tree_mp1
#     diff_ll <- -1*(abs(loglik_mp1) - abs(loglik_m)) %>% round(6)
#     marg_ll <- c(marg_ll, diff_ll)
#     sig2 <- sig2_draw()
#     param_draw <- mus_draw(n_tn, mu, sig2)
#     mu <- mean(param_draw)
#     leaf_params[[i]] <- param_draw
#     loglik_m <- loglik_mp1
#     
#   } else {
#     
#     tree_sims[[i]] <- tree_m
#     leaf_params[[i]] <- param_draw
#     marg_ll <- c(marg_ll, 0)
#     nodes_zl <- tree_m$left_nodes
#     nodes_zlr <- tree_m$left_comps
#     nodes_zr <- tree_m$right_nodes
#     nodes_zrr <- tree_m$right_comps
#     splits_zl <- tree_m$left_splits
#     splits_zr <- tree_m$right_splits
#     
#   }
#   
#   tree_m <- tree_sims[[i]]
#   cat('------Bayes tree draw action:', draw,
#       '\n------Iteration', i, '/', mh_n,
#       '\n\n')
# }

#### Bayesian backfitting additive regression tree ------------------------------

# choosen parameters -> CGM98 (alpha, beta, nu, lambda, a)
alpha <- 0.95
beta <- 2

m <- 1000
k <- 2
sig2_mu <- 0.5/(k*sqrt(m)) 

nu <- 3
res <- resid(lm(y ~ x, data = data))
lambda <- unname(quantile(res, 0.9))
a <- var(res)/var(data$y)

# define initial tree
tree_m <- t_grow()

# initialize parameter matrix and rho vector
leaf_params <- matrix(0, nrow = m, ncol = m/min_leaf_obs)
rho_vec <- NA

# initialize likelihood from base tree, 
treedata <- get_treedata(1)
all_y <- treedata$all_y
y_list <- treedata$y_list
n_tn <- treedata$n_tns
n_vec <- treedata$n_vec
mu <- treedata$mu
loglik_m <- loglik(n_vec, y_list, mu, n_tn)

# get initial residuals and updated mu for tree_mp1
r_list <- treedata$r_list

# update response data to current residuals
set_resdata()

# initialize tree sim list
tree_sims <- list(tree_m)
treedata_list <- list(treedata)

# MH implementation
for (i in 2:m) {
  
  # draw action
  draw <- draw_action()$draw
  
  if (draw == 'grow') {
    
    tree_mp1 <- t_grow()
    g_kern <- draw_action()$g_kern
    p_kern <- draw_action()$p_kern
    
  } else if (draw == 'prune') {
    
    tree_mp1 <- t_prune()
    g_kern <- draw_action()$g_kern
    p_kern <- draw_action()$p_kern
    
  } else if (draw == 'change') {
    
    tree_mp1 <- t_change()
    
  } else {
    
    tree_mp1 <- t_swap()
    
  }
  
  # get new tree data
  treedata <- get_treedata(i)
  all_y <- treedata$all_y
  y_list <- treedata$y_list
  n_tn <- treedata$n_tns
  n_vec <- treedata$n_vec
  r_list <- treedata$r_list
  mu <- treedata$mu
  loglik_mp1 <- loglik(n_vec, r_list, mu, n_tn)
  
  treedata_list[[i]] <- treedata
  
  # tree_mp1 acceptance parameter 
  rho <- rho_draw(g_kern, p_kern, loglik_mp1, loglik_m, draw)
  rho_vec[i] <- rho
  
  if (rho == 1) {
    
    # T* = tree_mp1, store tree_mp1 data from above
    # compute residuals for next iter
    set_resdata()
    
    # store tree_mp1 as tree_m
    tree_sims[[i]] <- tree_mp1
    
    # store updated likelihood value as mth iter likelihood
    loglik_m <- loglik_mp1
    
    cat('T(m+1) = T*',
        '\n')
    
    } else {
      
      # T* = tree_m
      tree_sims[[i]] <- tree_m
      
      # cancel drawn parameters for unaccepted tree
      leaf_params[i, ] <- 0
      
      cat('T(m+1) = T(m)',
          '\n')
      
    }
  
  tree_m <- tree_sims[[i]]
  
  cat('------Bayes tree draw action:', draw,
      '\n------Iteration', i, '/', m,
      '\n\n')
}

#### Prediction ----------------------------

# weighted nodes expectation
vis_nodes <- data.frame(l = sapply(tree_sims, function(x) length(x$left_splits)),
                        r = sapply(tree_sims, function(x) length(x$right_splits)))

get_Enodes <- function(direc) {
  
  w_node <- sum(table(vis_nodes[[direc]])/m*sort(unique(vis_nodes[[direc]])))
  idx <- which.min(abs(w_node - unique(vis_nodes[[direc]])))
  Enode <- unique(vis_nodes[[direc]])[idx]
  
  return(Enode)
  
}

Enodes <- c(get_Enodes('l'), get_Enodes('r'))

# expected left splits
sub_lsplits_idx <- sapply(tree_sims, function(x) length(x$left_splits) == Enodes[1])
sub_lsplits <- sapply(tree_sims[c(sub_lsplits_idx)], function(x) x$left_splits)
if (Enodes[1] > 1) {
  
  # ifelse not handling return rowMeans()
  Elsplits <- rowMeans(sub_lsplits)
  
} else {
  
  Elsplits <- split_z
  
}

# expected right splits
sub_rsplits <- sapply(tree_sims, function(x) length(x$right_splits) == Enodes[2])
sub_rsplits <- sapply(tree_sims[c(sub_rsplits)], function(x) x$right_splits)
if (Enodes[2] > 1) {
  
  Ersplits <- rowMeans(sub_rsplits)
  
} else {
  
  Ersplits <- split_z
  
}

# get associated leaf parameters
get_leafparams <- function(len_sp, direc) {
  
  p_vec <- c()
  
  if (direc == 'right') {
    
    for (i in 1:len_sp) {
      
      if (len_sp == 1) {
        
        p_vec[i] <- data %>% 
          filter(z > Ersplits[1]) %>% 
          .$org_y %>% 
          mean()
          
      } else if (i == len_sp) {
        
        p_vec[i] <- data %>% 
          filter(z > Ersplits[1] &
                   x < Ersplits[i]) %>% 
          .$org_y %>% 
          mean()
        
      } else if (i == 1) {
        
        p_vec[i] <- data %>% 
          filter(z > Ersplits[1] &
                   x > Ersplits[i + 1]) %>% 
          .$org_y %>% 
          mean()
        
      } else {
        
        p_vec[i] <- data %>% 
          filter(z > Ersplits[1] &
                   x < Ersplits[i] &
                   x > Ersplits[i + 1]) %>% 
          .$org_y %>% 
          mean()
        
      }
      
    }
    
  } else {
    
    for (i in 1:len_sp) {
      
      if (len_sp == 1) {
        
        p_vec[i] <- data %>% 
          filter(z < Elsplits[1]) %>% 
          .$org_y %>% 
          mean()
        
      } else if (i == len_sp) {
        
        p_vec[i] <- data %>% 
          filter(z < Elsplits[1] &
                   x < Elsplits[i]) %>% 
          .$org_y %>% 
          mean()
        
        } else if (i == 1) {
        
        p_vec[i] <- data %>% 
          filter(z < Elsplits[1] &
                   x > Elsplits[i + 1]) %>% 
          .$org_y %>% 
          mean()
        
      } else {
        
        p_vec[i] <- data %>% 
          filter(z < Elsplits[1] &
                   x < Elsplits[i] &
                   x > Elsplits[i + 1]) %>% 
          .$org_y %>% 
          mean()
        
      }
      
    }
    
  }
  
  return(p_vec)
  
}

# retrieve parameter values
Elparams <- get_leafparams(length(Elsplits), 'left')
Erparams <- get_leafparams(length(Ersplits), 'right')

# splits may be too close together and produce no obs
if (any(is.na(Elparams))) {
  
  naidx <- which(is.na(Elparams))
  Elparams <- Elparams[-naidx]
  Elsplits <- Elsplits[-(naidx + 1)]
  
} else if (any(is.na(Erparams))) {
  
  naidx <- which(is.na(Erparams))
  Erparams <- Erparams[-naidx]
  Ersplits <- Ersplits[-(naidx + 1)]
  
}

# set values for plotting fits
data$yfit <- 0
llparams <- length(Elparams)
lrparams <- length(Erparams)
  
if (llparams == 1) {
  
  idx <- which(data$z < Elsplits[1])
  data$yfit[idx] <- Elparams[1]
  
}

if (llparams == 2) {
  
  idx1 <- which(data$z < Elsplits[1] &
                  data$x >= Elsplits[2])
  data$yfit[idx1] <- Elparams[1]
  idx2 <- which(data$z < Elsplits[1] &
                  data$x < Elsplits[2])
  data$yfit[idx2] <- Elparams[2]
  
}

for (i in 1:llparams) {
  
  if (i == 1) {
    
    idx <- which(data$z < Elsplits[1] &
                   data$x >= Elsplits[2])
    data$yfit[idx] <- Elparams[1]
    
  } else {
    
    idx <- which(data$z < Elsplits[1] &
                   data$x < Elsplits[i] & 
                   data$x >= Elsplits[i + 1])
    data$yfit[idx] <- Elparams[i]
    
  } 
  
  if (i == llparams) {
    
    idx <- which(data$z < Elsplits[1] &
                   data$x < Elsplits[i])
    data$yfit[idx] <- Elparams[i]
    
  }
  
}

if (lrparams == 1) {
  
  idx <- which(data$z > Ersplits[1])
  data$yfit[idx] <- Erparams[1]
  
}

if (lrparams == 2) {
  
  idx1 <- which(data$z > Ersplits[1] &
                  data$x >= Ersplits[2])
  data$yfit[idx1] <- Erparams[1]
  idx2 <- which(data$z > Ersplits[1] &
                  data$x < Ersplits[2])
  data$yfit[idx2] <- Erparams[2]
  
}

for (i in 1:lrparams) {
  
  if (i == 1) {
    
    idx <- which(data$z > Ersplits[1] &
                   data$x >= Ersplits[2])
    data$yfit[idx] <- Erparams[1]
    
  } else {
    
    idx <- which(data$z > Ersplits[1] &
                   data$x < Ersplits[i] & 
                   data$x >= Ersplits[i + 1])
    data$yfit[idx] <- Erparams[i]
    
  }
  
  if (i == lrparams) {
    
    idx <- which(data$z > Ersplits[1] &
                   data$x < Ersplits[i])
    data$yfit[idx] <- Erparams[i]
    
  }

}

plot(data$x, data$org_y, xlab = 'x', ylab = 'y')
points(data[data$z < 0.5, 'x'], data[data$z < 0.5, 'yfit'],
       col = 'red', pch = 19, cex = 0.75)
points(data[data$z > 0.5, 'x'], data[data$z > 0.5, 'yfit'],
       col = 'blue', pch = 19, cex = 0.75)
title(main = 'Generated Treatment Data and BART Prediction')
legend("topleft", pch = 19, col = c('blue', 'red'), box.lty = 0,
       legend = c('Yhat | Z = 1', 'Yhat | Z = 0'))







