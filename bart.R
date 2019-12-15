# Bayesian Addative Regression Trees #
# Chipman, George, McCulloch (1998) -> CGM98
remove(list = objects())
library(tidyverse)

#### Generate data ----------------------------------
# data generation following Hill (2010)
set.seed(7777)
p <- 0.5
z <- rbinom(200, 1, p)
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

set.seed(NULL)
ridx <- which(is.na(y))
data <- data.frame(y = y[-ridx],
                   z = z[-ridx],
                   x = x[-ridx])
remove(x, y, z)

#### Priors and likelihood -------------------------
# back to CGM98

# tree prior
p_split <- function(d, a = alpha, b = beta) {
  
  s <- a*(1 + d)^(-beta)
  # rv <- rbinom(1, 1, s)
  return(s)
  
}

# leaf param prior (mu)
mus_draw <- function(tn, mu, s2, c = a) {
  
  mus <- rnorm(tn, mu, sqrt(s2/c))
  return(mus)
  
}

# leaf param hyper
sig2_draw <- function(v = nu, lam = lambda) {
  
  s2 <- 1/rgamma(1, v/2, scale = v*lam/2)
  return(s2)
  
}

# log-likelihood 
loglik <- function(n_vec, y_list, mu, n_tn, c = a, v = nu, lam = lambda) {
  
  for (i in 1:length(y_list)) {
    
    j <- ifelse(length(y_list[[i]]) == 0, i, NA)
    
  }
  
  if (!is.na(j)) {
    
    y_list <- y_list[-j]
    n_vec <- n_vec[-j]
    n_tn <- n_tn - 1
    
  }
  
  s_vec <- map2_dbl(y_list, n_vec,
                    function(y, n) (n - 1)*var(y))
  t_vec <- map2_dbl(y_list, n_vec,
                    function(y, n) ((n*c)/(n + c))*(mean(y) - mu)^2)
  
  f <- -sum(n_vec/2)*log(pi) + n_tn*v/2*log(v*lam) + n_tn/2*log(c) +
    -1/2*sum(log(n_vec + c)) + sum(log(gamma((n_vec + v)/2))) - 
    n_tn*log(gamma(v/2)) - sum((n_vec + v)/2)*sum(log(s_vec + t_vec + v*lam))
  
  return(f)
  
}

# Tree acceptance draw
rho_draw <- function(k_grow, k_prune, Lmp1, Lm, tdraw) {
  
  if(tdraw == 'grow') {
    
    rho <- min(k_prune*Lmp1/Lm, 1)
    
  } else if (tdraw == 'prune') {
    
    rho <- min(Lmp1/(k_grow*Lm), 1)
    
  } else {
    
    rho <- min(Lmp1/Lm, 1)
    
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
splits_zl <- c(NA)
splits_zr <- c(NA)
# initialize leaf size constraint
min_leaf_obs <- 4
# initialize leaf parameter vectors
leaf_params <- setNames(data.frame(matrix(0, nrow = 1, ncol = 4)),
                        c('l', 'lcomp', 'r', 'rcomp'))
# initiate tree random action
tree_step <- c('grow', 'prune', 'change', 'swap')
draw_action <- function(d) {
  
  # kernel probabilities
  g <- p_split(d)
  p <- (1 - g)/3
  s <- sample(tree_step, 1, prob = c(p, rep(p, length(tree_step) - 1)))
  out <- list(g_kern = g, p_kern = p, draw = s)
  return(out)
  
}

#### Tree action functions --------------------------------

# if action is grow
t_grow <- function() {
  
  # draw initial z node direction, < 0.5 or > 0.5
  s_direct <- sample(c('l', 'r'), 1)
  
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
          filter(between(x, splits[i], splits[i - 1]) &
                   !y %in% unlist(nodes_r) &
                   y %in% nodes[[1]]) %>%
          .$y
        
      }
      
      # assign data to compliment node
      nodes_r[[i]] <- split_comp_y
      
    } else {
      
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
              depth = d, direc = s_direct)
  return(out)
  
}

# if action is prune
t_prune <- function() {
  
  s_direct <- sample(c('l', 'r'), 1)
  
  if(length(nodes_zl) == 1  & s_direct == 'l') {
    
    cat('Cannot prune', s_direct, 'any further.',
        '\n')
    
  } else if (length(nodes_zr) == 1  & s_direct == 'r') {
    
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
              depth = d, direc = s_direct)
  return(out)
  
}

# if action is change
t_change <- function() {
  
  s_direct <- sample(c('l', 'r'), 1)
  
  if(length(nodes_zl) == 1 & s_direct == 'l') {
    
    cat('Cannot change splitting rule for direction:', s_direct,
        '\nCurrent node is the initial node.',
        '\n')
    
  } else if (length(nodes_zr) == 1 & s_direct == 'r') {
    
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
          filter(between(x, splits[i], splits[i - 1]) &
                   !y %in% unlist(nodes_r) &
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
          filter(between(x, splits[i], splits[i - 1]) &
                   !y %in% unlist(nodes_r) &
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
              depth = d, direc = s_direct)
  return(out)
  
}

# if action is swap
t_swap <- function() {
  
  s_direct <- sample(c('l', 'r'), 1)
  
  if (length(nodes_zl) == 1 & s_direct == 'l') {
    
    cat('Cannot swap splitting rule for direction:', s_direct,
        '\nParent node is the initial node.',
        '\n')
    
  } else if (length(nodes_zr) == 1 & s_direct == 'r') {
    
    cat('Cannot change splitting rule direction:', s_direct,
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
          filter(between(x, splits[i], splits[i - 1]) &
                   !y %in% unlist(nodes_r) &
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
            depth = d, direc = s_direct)
  return(out)
  
}

# get tree parameters needed for MH iters
get_treedata <- function() {
  
  y_list <- c(nodes_zl, nodes_zlr[-1],
              nodes_zr, nodes_zrr[-1])
  n_termnodes <- length(y_list)
  n_vec <- sapply(y_list, length)
  
  out <- list(ys = y_list, n_tns = n_termnodes, ns = n_vec)
  return(out)
  
}

#### MH algorithm --------------------------------------

# choosen parameters -> CGM98 (alpha, beta, nu, lambda, a)
res <- resid(lm(y ~ x, data = data))
alpha <- 0.95
beta <- 2
nu <- 3
lambda <- unname(quantile(res, 0.9))
a <- var(res)/var(data$y)

# initial parameters for updating
mu <- mean(data$y)
sig2 <- resid(lm(y ~ x, data = data))^2 %>% mean 

# define initial tree
tree_m <- t_grow()
tree_sims <- list(tree_m)

# initialize likelihood from base tree
y_list <- get_treedata()$ys
n_tn <- get_treedata()$n_tns
n_vec <- get_treedata()$ns
loglik_m <- loglik(n_vec, y_list, mu, n_tn)

max_itr <- 5
for (i in 2:max_itr) {
  
  # draw action
  action <- draw_action(tree_m$depth)
  draw <- action$draw
  
  if (draw == 'grow') {
    
    tree_mp1 <- t_grow()
    g_kern <- draw_action(tree_mp1$depth)$g_kern
    p_kern <- draw_action(tree_mp1$depth)$p_kern
    
  } else if (draw == 'prune') {
    
    tree_mp1 <- t_prune()
    g_kern <- draw_action(tree_mp1$depth)$g_kern
    p_kern <- draw_action(tree_mp1$depth)$p_kern
    
  } else if (draw == 'change') {
    
    tree_mp1 <- t_change()
    
  } else {
    
    tree_mp1 <- t_swap()
    
  }
  
  y_list <- get_treedata()$ys
  n_tn <- get_treedata()$n_tns
  n_vec <- get_treedata()$ns
  sig2 <- sig2_draw()
  mu <- mean(mus_draw(n_tn, mu, sig2))
  loglik_mp1 <- loglik(n_vec, y_list, mu, n_tn)
  rho <- rho_draw(g_kern, p_kern, loglik_mp1, loglik_m, draw)
  if (rho == 1) {
    
    tree_sims[[i]] <- tree_mp1
    loglik_m <- loglik_mp1
    
  } else {
    
    tree_sims[[i]] <- tree_m
    nodes_zl <- tree_m$left_nodes
    nodes_zlr <- tree_m$left_comps
    nodes_zr <- tree_m$right_nodes
    nodes_zrr <- tree_m$right_comps
    splits_zl <- tree_m$left_splits
    splits_zr <- tree_m$right_splits
    
  }
  
  tree_m <- tree_sims[[2]]
  cat('----Bayes tree draw action:', draw,
      '\n----Iteration', i, '/', max_itr,
      '\n\n')
}





