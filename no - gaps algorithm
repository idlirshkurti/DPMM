#-----------------------------------------------------#
#                NO GAPS ALGORITHM                    #
#-----------------------------------------------------#

# Create a function which rearranges the vector c_i 
# of parameter allocations

rearrange_s <- function(s){
  b <- unique(s); l <- length(s); k <- length(b); m <- max(s)
  gap <- setdiff(c(1:m),s)
  if(length(setdiff(c(1:m),s))>0){
    for(i in length(gap):1 ){
      s[which(s>gap[i])]<-s[which(s>gap[i])]-1
    }
  } return(s)
}


# Likelihood function
likelihood <- function(x1,pi1,a,b){
  return(1-pi1+pi1*dbeta(x1,a,b))}

# Transformation functions
exp_trans <- function(alpha){L_alpha = exp(-abs(alpha))
  return(L_alpha)
}

exp_trans1 <- function(alpha){L_alpha = exp(abs(alpha))
  return(L_alpha)
}


dpmm <- function(data, n.iter, tau, mua, mub, mup,
sigma2a, sigma2b, sigma2p, alpha_val, beta_val, s, pi1){

  C <- length(data)
  A <- matrix(0, ncol=C, nrow=n.iter)
  B <- matrix(0, ncol=C, nrow=n.iter)
  P <- NULL
  K <- NULL
  pb <- txtProgressBar(min = 0, max = n.iter, style = 3)

  #-------------- Start iterations -----------------#
  for(t in 1:n.iter){
     if(t==3000){sigma2a=sigma2a*0.75}
     if(t==6000){sigma2a=sigma2a*0.75}
     if(t==10000){sigma2a=sigma2a*0.75}
     # Do the same with sigmab and sigmap

    L_pi1 <- exp_trans(pi1); k <- length(unique(s))
    s <- rearrange_s(s); m_s <- table(s[-i])
    for (i in 1:C){
      s1 <- s
      singleton <- table(s)[s[i]] == 1
      if(singleton){
        u = runif(1)
        if(u < (k-1)/k){}else{
        
          # rearrange c and \Phi  
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
          
          m_s <- table(s[-i]); k <- length(m_s)
          
          p_xi <- NULL
          for(j in 1:k){
            L_alpha_val  = exp_trans(alpha_val[j])
            L_beta_val  = exp_trans1(beta_val[j])
            p_xi[j]=prod(likelihood(data[i], L_pi1,
            L_alpha_val,L_beta_val))
          }
          
          w <- NULL
          w[1:k] <- table(s[-i])*p_xi[1:k]
          if(k<C){
            L_alpha_val  = exp_trans(alpha_val[k+1])
            L_beta_val  = exp_trans1(beta_val[k+1])
            w[k+1] = tau/(k+1)*prod(likelihood(data[i], L_pi1,
            L_alpha_val,L_beta_val))
            population = 1:c(k+1)
          }else{population = 1:k}
          s[i] = sample(x = population ,size = 1,
          replace = TRUE, prob = w)
        }
        
        
        k <- length(unique(s))
        m_s <- table(s)
        rem <- setdiff(s1,s)
        if(length(rem)>0){
          alpha_val <- alpha_val[-rem]}
        if(length(rem)>0){
          beta_val <- beta_val[-rem]}
        s <- rearrange_s(s)
        
      }else{ 

        # m_s(i) > 1 case
        
        m_s <- table(s[-i])
      	k <- length(m_s)
        alpha_val[k+1] <- rnorm(1,mua,sigma2a)
        beta_val[k+1] <- rnorm(1,mub,sigma2b)
        p_xi <- NULL
        for(j in 1:k){
          L_alpha_val  = exp_trans(alpha_val[j])
          L_beta_val  = exp_trans1(beta_val[j])
          p_xi[j] = prod(likelihood(data[i], L_pi1,
          L_alpha_val, L_beta_val))
        }
        w <- NULL
        w[1:k] <- table(s[-i])*p_xi[1:k]
        if(k<C){
          L_alpha_val  = exp_trans(alpha_val[k+1])
          L_beta_val  = exp_trans1(beta_val[k+1])
          temp1 = prod(likelihood(data[i], L_pi1,
          L_alpha_val, L_beta_val))
          temp2 = tau/(k+1)*temp1
          w[k+1] = temp2
          p_xi[k+1] = temp1
          population = 1:c(k+1)
        }else{
          population = 1:k
        }
s[i] = sample(x = population ,size = 1, 
replace = TRUE, prob = w)
        
        k <- length(unique(s)); m_s <- table(s)
 
        rem <- setdiff(s1,s)
        if(length(rem)>0){
          alpha_val <- alpha_val[-rem]}
        if(length(rem)>0){
          beta_val <- beta_val[-rem]}
        
        if(s[i]!=k+1){alpha_val <- alpha_val[-c(k+1)]}
        if(s[i]!=k+1){beta_val <- beta_val[-c(k+1)]}
        s <- rearrange_s(s)
      }
      
      k <- length(unique(s))
      m_s <- table(s)
      rem <- setdiff(s1,s)
      if(length(rem)>0){
        alpha_val <- alpha_val[-rem]}
      if(length(rem)>0){
        beta_val <- beta_val[-rem]}
      s <- rearrange_s(s)
    }
    
    
    # ----------- MCMC updates step -------------- #
    
    # Part 1 / Phi_E
    
    if(k<C){
      for(i in c(k+1):C){
        alpha_val[i] <- rnorm(1,mua,sigma2a)
        beta_val[i] <- rnorm(1,mub,sigma2b)
      }
    }
    
    # Part 2 / Phi_F

    for(i in 1:k){
      Va <- sigma2a*rep(1,C)
      Vb <- sigma2b*rep(1,C)
      
      alpha_val_new = rnorm(1, alpha_val[i], sigma2a)
      beta_val_new = rnorm(1, beta_val[i], sigma2b)
      
      idx <- which(s==i)
      # new values
      L_alpha_val <- exp_trans(alpha_val_new)
      L_beta_val <- exp_trans1(beta_val_new)
      # old values 
      L_alpha_val_t <- exp_trans(alpha_val[i])
      L_beta_val_t <- exp_trans1(beta_val[i])
      # Note that alpha_val are actually the L_a from text
      
      p_xi <- NULL
      p_xi_t <- NULL
      for(m in 1:length(idx)){
        p_xi[m] <- sum(log(likelihood(x = data[idx[m]],
        L_pi1, L_alpha_val, L_beta_val))) # new likelihood 
        p_xi_t[m] <- sum(log(likelihood(x = data[idx[m]],
        L_pi1, L_alpha_val_t, L_beta_val_t))) # old likelihood
      }
      
      
      p_x <- sum(p_xi)
      p_x_t <- sum(p_xi_t)
      
      # Prime densities
      fx <- dnorm(alpha_val_new, mean = 0, sd = sigma2a)*dnorm(beta_val_new, mean = 0, sd = sigma2b)
      fx_t <- dnorm(alpha_val[i], mean = 0, sd = sigma2a)*dnorm(beta_val[i], mean = 0, sd = sigma2b)
      
      # Metropolis-Hastings sampling
      
      tmp <- log(fx) + p_x - log(fx_t) - p_x_t
      tmp <- exp(tmp)
      u <- runif(1)
      if(u<tmp){
        alpha_val[i]=alpha_val_new}
      if(u<tmp){
        beta_val[i]=beta_val_new} 
    }
    alpha_val <- alpha_val[-c(c(k+1):C)]
    beta_val <- beta_val[-c(c(k+1):C)]
    
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
      p_xi[i] <- sum(log(likelihood(data[i], L_pi1_new,
      L_alpha_val[s[i]], L_beta_val[s[i]]))) # new likelihood 
      p_xi_t[i] <- sum(log(likelihood(data[i], L_pi1,
      L_alpha_val[s[i]], L_beta_val[s[i]]))) # old likelihood
    }
    p_x <- sum(p_xi)
    p_x_t <- sum(p_xi_t)
    
    # Priors
    fx = dnorm(pi1_new,0,sigma2p)
    fx_t = dnorm(pi1,0,sigma2p)
    
    # Metropolis-Hastings sampling
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
    setTxtProgressBar(pb, t)
  }  
