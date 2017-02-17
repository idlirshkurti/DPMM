dpmm_neal_8 <- function(data, n.iter, tau, l, s, mua, mub, mup,
sigma2a, sigma2b, sigma2p, alpha_val, beta_val, pi1 ){
  
  # note that l is the number of temporary auxiliary parameters introduced
  
  C <- length(data)
  A <- matrix(0, ncol = C, nrow=n.iter)
  B <- matrix(0, ncol = C, nrow=n.iter)
  P <- NULL
  K <- NULL
  
  #-------------- Start iterations -----------------#
  for(t in 1:n.iter){
    if(t==3000){sigma2a=sigma2a*0.5}
    if(t==6000){sigma2a=sigma2a*0.5}
    if(t==10000){sigma2a=sigma2a*0.5}
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
            p_xi[j]=prod(likelihood(data[i], L_pi1,
            L_alpha_val,L_beta_val))
          }
          w <- NULL
          w[1:k] <- (table(s[-i])/(C-1+tau))*p_xi[1:k]

          for(j in 1:l){
            L_alpha_val  = exp_trans(alpha_val[k+j])
            L_beta_val  = exp_trans1(beta_val[k+j])
            w[k+j] = (tau/l)/(C-1+tau)*prod(likelihood(data[i], L_pi1,
            L_alpha_val,L_beta_val))
          }
            population = 1:c(k+l)

          s[i] = sample(x = population ,size = 1,
          replace = TRUE, prob = w)
          
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

        
        m_s <- table(s[-i])
        k <- length(m_s)
        s1 <- s
        alpha_val[c(k+1):c(k+l)] <- rnorm(l,mua,sigma2a)
        beta_val[c(k+1):c(k+l)] <- rnorm(l,mub,sigma2b)
        p_xi <- NULL
        for(j in 1:k){
          L_alpha_val  = exp_trans(alpha_val[j])
          L_beta_val  = exp_trans1(beta_val[j])
          p_xi[j] = prod(likelihood(data[i], L_pi1,
          L_alpha_val, L_beta_val))
        }
        w <- NULL
        w[1:k] <- (table(s[-i])/(C-1+tau))*p_xi[1:k]
       
        temp1 <- NULL
        for(j in 1:l){
          L_alpha_val  = exp_trans(alpha_val[k+j])
          L_beta_val  = exp_trans1(beta_val[k+j])
          temp1[j] = prod(likelihood(data[i], L_pi1,
          L_alpha_val, L_beta_val))
        }
        
          temp2 = (tau/l)/(C-1+tau)*temp1
          w[c(k+1):c(k+l)] = temp2
          p_xi[c(k+1):c(k+l)] = temp1
          population = 1:c(k+l)

        s[i] = sample(x = population ,size = 1, 
        replace = TRUE, prob = w)
        
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
    # PROCEED SIMILARLY TO THE PREVIOUS NO-GAPS 
    # ALGORIGHM IN ORDER TO GET THE PARAMETER
    # UPDATES USING A METROPOLIS RANDOM WALK
}
}
