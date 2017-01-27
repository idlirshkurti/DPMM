# -------------- NEIL ALGORITHM 8 (l=3) -------------- #

rearrange_s <- function(s){
  b <- unique(s)
  l <- length(s)
  k <- length(b)
  m <- max(s)
  
  gap <- setdiff(c(1:m),s)
  
  if(length(setdiff(c(1:m),s))>0){
  for(i in length(gap):1 ){
    s[which(s>gap[i])]<-s[which(s>gap[i])]-1
  }
  }
  return(s)
}

# Likelihood function
likelihood <- function(x1,pi1,a,b){
  return(1-pi1+pi1*dbeta(x1,a,b))}



exp_trans <- function(alpha){
  
  L_alpha = exp(-abs(alpha))
  L_alpha
}

exp_trans1 <- function(alpha){
  
  L_alpha = exp(abs(alpha))
  L_alpha
}


# Generate random data

data <- c(rbeta(10,0.1,6.1),rbeta(10,0.5,3.5))
# pi_1 = 0.4

###################################
#dpmm <- function(data){
  
  C <- length(data)
  COUNT_NUM <- C
  k <- C # number of clusters
  n.iter <- 20000
  
  s <- 1:C
  l <- 3
  tau <- 1
  mua = 1
  mub = 1
  mup = 1
  sigma2a = 1
  sigma2b = 1
  sigma2p = 0.5
  alpha_val <- rep(0, C)
  beta_val <- rep(0,C)
  for (i in 1:C){
    alpha_val[i] <- rnorm(1,mua,sigma2a)
    beta_val[i] <- rnorm(1,mub,sigma2b)
  }
  A <- matrix(0, ncol=C, nrow=n.iter)
  B <- matrix(0, ncol=C, nrow=n.iter)
  P <- NULL
  K <- NULL
  
  a <- 6
  b <-5 
  pi1 <- rnorm(1,mup,sigma2p)
  
  ptm <- proc.time()
  start.time <- Sys.time()
  #-------------- Start algorithm -----------------#
  for(t in 1:n.iter){
  #  if(t==3000){sigma2a=sigma2a*0.5}
  #  if(t==6000){sigma2a=sigma2a*0.5}
  #  if(t==10000){sigma2a=sigma2a*0.5}
    print(t)
    
    L_pi1 <- exp_trans(pi1)
    k <- length(unique(s))
    s <- rearrange_s(s)

    m_s <- table(s[-i])
    
    for (i in 1:C){
      s1 <- s
      singleton <- table(s)[s[i]] == 1
      if(singleton){
        u = runif(1)
        if(u < (k-1)/k){}else{
    #      print(" rearrange s. m=1 case")
          # rearrange s and \Phi  
          idx = which(s == k)
          tmp = s[i]
          s[idx] = tmp
          s[i] = k
          
          tmp_alpha_val = alpha_val[tmp]
          alpha_val[tmp] = alpha_val[k]
          alpha_val[k] = tmp_alpha_val
          tmp_beta_val = beta_val[tmp]
          beta_val[tmp] = beta_val[k]
          beta_val[k] = tmp_beta_val
          
          m_s <- table(s[-i])
          k <- length(m_s)
          alpha_val[c(k+2):c(k+l)] <- rnorm(l-1,mua,sigma2a)
          beta_val[c(k+2):c(k+l)] <- rnorm(l-1,mub,sigma2b)
          
          p_xi <- NULL
          for(j in 1:k){
            L_alpha_val  = exp_trans(alpha_val[j])
            L_beta_val  = exp_trans1(beta_val[j])
            p_xi[j]=prod(likelihood(data[i], L_pi1, L_alpha_val,L_beta_val))
          }
          w <- NULL
          w[1:k] <- (table(s[-i])/(C-1+tau))*p_xi[1:k]

          for(j in 1:l){
            L_alpha_val  = exp_trans(alpha_val[k+j])
            L_beta_val  = exp_trans1(beta_val[k+j])
            w[k+j] = (tau/l)/(C-1+tau)*prod(likelihood(data[i], L_pi1, L_alpha_val,L_beta_val))
          }
            population = 1:c(k+l)

          s[i] = sample(x = population ,size = 1, replace = TRUE, prob = w)
          
          k <- length(unique(s))
          m_s <- table(s)
          
          rem <- setdiff(population,s)
          if(length(rem)>0){
            alpha_val <- alpha_val[-rem]}
          if(length(rem)>0){
            beta_val <- beta_val[-rem]}
          s <- rearrange_s(s)
        }
        
      }else{ 
        # m_s(i) > 1 case
        # print(" rearrange s. m>1 case")
        
        m_s <- table(s[-i])
        k <- length(m_s)
        s1 <- s
        alpha_val[c(k+1):c(k+l)] <- rnorm(l,mua,sigma2a)
        beta_val[c(k+1):c(k+l)] <- rnorm(l,mub,sigma2b)
        p_xi <- NULL
        for(j in 1:k){
          L_alpha_val  = exp_trans(alpha_val[j])
          L_beta_val  = exp_trans1(beta_val[j])
          p_xi[j] = prod(likelihood(data[i], L_pi1, L_alpha_val, L_beta_val))
        }
        w <- NULL
        w[1:k] <- (table(s[-i])/(C-1+tau))*p_xi[1:k]
       
        temp1 <- NULL
        for(j in 1:l){
          L_alpha_val  = exp_trans(alpha_val[k+j])
          L_beta_val  = exp_trans1(beta_val[k+j])
          temp1[j] = prod(likelihood(data[i], L_pi1, L_alpha_val, L_beta_val))
        }
        
          temp2 = (tau/l)/(C-1+tau)*temp1
          w[c(k+1):c(k+l)] = temp2
          p_xi[c(k+1):c(k+l)] = temp1
          population = 1:c(k+l)

        s[i] = sample(x = population ,size = 1, replace = TRUE, prob = w)
        
        k <- length(unique(s))
        m_s <- table(s)
        
        rem <- setdiff(population,s)
        
        if(length(rem)>0){
          alpha_val <- alpha_val[-rem]}
        if(length(rem)>0){
          beta_val <- beta_val[-rem]}
        
        s <- rearrange_s(s)
      }
    }
    
    # Updates
    k <- length(unique(s))
    m_s <- table(s)
    s <- rearrange_s(s)
    
    # ----------- MCMC updates step -------------- #
    
    # Part 1 / Phi_E
    
    #   print("part 1")

    
    # Part 2 / Phi_F
    
    #   print("part 2")
    for(i in 1:k){
      Va <- sigma2a*rep(1,C)
      Vb <- sigma2b*rep(1,C)
      
      #       for(cou_ind in 1:COUNT_NUM){
      alpha_val_new = rnorm(1, alpha_val[i], sigma2a)
      beta_val_new = rnorm(1, beta_val[i], sigma2b)
      
      idx <- which(s==i)
      # new values
      L_alpha_val <- exp_trans(alpha_val_new)
      L_beta_val <- exp_trans1(beta_val_new)
      # old values 
      L_alpha_val_t <- exp_trans(alpha_val[i])
      L_beta_val_t <- exp_trans1(beta_val[i])
      
      p_xi <- NULL
      p_xi_t <- NULL
      
      for(m in 1:length(idx)){
        p_xi[m] <- sum(log(likelihood(x = data[idx[m]], L_pi1, L_alpha_val, L_beta_val))) # likelihood 
        p_xi_t[m] <- sum(log(likelihood(x = data[idx[m]], L_pi1, L_alpha_val_t, L_beta_val_t)))
      }
      
      
      p_x <- sum(p_xi)
      p_x_t <- sum(p_xi_t)
      
      fx <- dnorm(alpha_val_new, mean = 0, sd = sigma2a)*dnorm(beta_val_new, mean = 0, sd = sigma2b)
      
      fx_t <- dnorm(alpha_val[i], mean = 0, sd = sigma2a)*dnorm(beta_val[i], mean = 0, sd = sigma2b)
      
      # Metropolis-Hastings sampling
    #    print("M1")
      
      tmp <- log(fx) + p_x - log(fx_t) - p_x_t
      tmp <- exp(tmp)
      u <- runif(1)
      if(u<tmp){
        alpha_val[i]=alpha_val_new}
      if(u<tmp){
        beta_val[i]=beta_val_new}
      #       }
      
      #      if(count>=3){accept <- 1}else{
      #        sigma2a <- sigma2a*0.98
      #        if(sigma2a<0.01){sigma2a <- 0.01}
      #        if(count > 20){sigma2a <- sigma2a/0.98}
      #      }
      
      
    }

    
    # Step 3: resample \pi
    
    m_s <- table(s)
    k <- length(unique(s))
    pi1_new <- rnorm(1, pi1, sigma2p)
    
    # old values
    L_alpha_val <- exp_trans(alpha_val)
    L_beta_val <- exp_trans1(beta_val)
    L_pi1 <- exp_trans(pi1)
    
    # new values 
    L_pi1_new <- exp_trans(pi1_new)
    
    p_xi <- NULL
    p_xi_t <- NULL
    
    for(i in 1:C){
    p_xi[i] <- sum(log(likelihood(data[i], L_pi1_new, L_alpha_val[s[i]], L_beta_val[s[i]]))) # likelihood 
    p_xi_t[i] <- sum(log(likelihood(data[i], L_pi1, L_alpha_val[s[i]], L_beta_val[s[i]])))
    }
    
    p_x <- sum(p_xi)
    p_x_t <- sum(p_xi_t)
    
    fx = dnorm(pi1_new,0,sigma2p)
    fx_t = dnorm(pi1,0,sigma2p)
    
    # Metropolis-Hastings sampling
  # print("M3")
    tmp <- log(fx) + p_x - log(fx_t) - p_x_t
    tmp <- exp(tmp)
    u <- runif(1)
    if(u<tmp){
      pi1 = pi1_new
    }
    
    A[t,] <- alpha_val[s]
    B[t,] <- beta_val[s]
    P[t] <- pi1
    K[t] <- k
  #  cat(paste("No of clusters =",k))
    if (t>25000 & k == 2) break
  }  

  proc.time() - ptm
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
 
  n.burnin <- 15000
  f <- NULL
  for(j in 15000:n.iter){
    f1 <- sum(dbeta(data,exp_trans(A[j,]), exp_trans1(B[j,])))/(1+length(data))  
    f[j] <- 1 - exp_trans(P[j]) + exp_trans(P[j])*f1
  }
  
  f1<-NULL
  for(j in 15000:n.iter){
    f1[j] <- sum(dbeta(x,Alpha[j,],2))/(1+length(x))  
  }
  
  plot(f)
  
  
  # Proportion of true null hypotheses
  
  pi_bayes <- 1 - 1/(n.iter-n.burnin)*sum(exp_trans(P))
  
  # Density plots
  denplot(exp_trans(A)) 
