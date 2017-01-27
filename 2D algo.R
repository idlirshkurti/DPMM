# -------------- NO-GAPS ALGORITHM -------------- #

smallgd<-golub[1:1000,]
classlabel<-golub.cl
# Permutation unadjusted p-values and adjusted p-values for maxT procedure
res1<-mt.maxT(smallgd,classlabel)
rawp<-res1$rawp[order(res1$index)]

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

data <- c(runif(30),rbeta(10,0.1,6.1),rbeta(10,0.5,3.5))
# data <- 1 - 0.2 + 0.2*data 
# pi_1 = 0.4


###################################
#dpmm <- function(data){
  
  C <- length(data)
  COUNT_NUM <- C
  k <- C # number of clusters
  n.iter <- 10000
  
  s <- 1:C
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
  pb <- txtProgressBar(min = 0, max = n.iter, style = 3)
  
  ptm <- proc.time()
  start.time <- Sys.time()
  #-------------- Start algorithm -----------------#
  for(t in 1:n.iter){
     if(t==3000){sigma2a=sigma2a*0.5}
     if(t==6000){sigma2a=sigma2a*0.5}
     if(t==10000){sigma2a=sigma2a*0.5}
   # print(t)
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
          
          p_xi <- NULL
          for(j in 1:k){
            L_alpha_val  = exp_trans(alpha_val[j])
            L_beta_val  = exp_trans1(beta_val[j])
            p_xi[j]=prod(likelihood(data[i], L_pi1, L_alpha_val,L_beta_val))
          }
          w <- NULL
          w[1:k] <- table(s[-i])*p_xi[1:k]
          if(k<C){
            L_alpha_val  = exp_trans(alpha_val[k+1])
            L_beta_val  = exp_trans1(beta_val[k+1])
            w[k+1] = tau/(k+1)*prod(likelihood(data[i], L_pi1, L_alpha_val,L_beta_val))
            population = 1:c(k+1)
          }else{population = 1:k}
          s[i] = sample(x = population ,size = 1, replace = TRUE, prob = w)
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
        #   print(" rearrange s. m>1 case")
        
        m_s <- table(s[-i])
        k <- length(m_s)
        s1 <- s
        alpha_val[k+1] <- rnorm(1,mua,sigma2a)
        beta_val[k+1] <- rnorm(1,mub,sigma2b)
        p_xi <- NULL
        for(j in 1:k){
          L_alpha_val  = exp_trans(alpha_val[j])
          L_beta_val  = exp_trans1(beta_val[j])
          p_xi[j] = prod(likelihood(data[i], L_pi1, L_alpha_val, L_beta_val))
        }
        w <- NULL
        w[1:k] <- table(s[-i])*p_xi[1:k]
        if(k<C){
          L_alpha_val  = exp_trans(alpha_val[k+1])
          L_beta_val  = exp_trans1(beta_val[k+1])
          temp1 = prod(likelihood(data[i], L_pi1, L_alpha_val, L_beta_val))
          temp2 = tau/(k+1)*temp1
          w[k+1] = temp2
          p_xi[k+1] = temp1
          population = 1:c(k+1)
        }else{
          population = 1:k
        }
        s[i] = sample(x = population ,size = 1, replace = TRUE, prob = w)
        
        k <- length(unique(s))
        m_s <- table(s)
        
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
    #  print("M3")
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
    #cat(paste("No of clusters =",k))
    setTxtProgressBar(pb, t)
  }  

  proc.time() - ptm
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  

library(coda);library(ggplot2);library(MCMCpack);library(mcmc);library(HDInterval)
library(gridExtra);library(mcmcplots);library(reshape2)

  
# --------- Estimators ----------- #

  n.burnin <- 15000
  M <- n.iter - n.burnin
  # f^(j)
  f <- NULL
  for(j in n.burnin:n.iter){
    f1 <- sum(dbeta(data,exp_trans(A[j,]), exp_trans1(B[j,])))/(tau+length(data))  
    f[j] <- 1 - exp_trans(P[j]) + exp_trans(P[j])*f1
  }

  # Proportion of true null hypotheses
  
  pi_bayes <- 1 - 1/(n.iter-n.burnin)*sum(exp_trans(P[-c(1:n.burnin)]))

  # F^(j)
  
  Fj <- NULL
  for(j in 1:M){
    F1 <- sum(qbeta(sig.level, exp_trans(A[c(n.burnin+j),]), exp_trans1(B[c(n.burnin+j),])))/(tau+length(data))  
    Fj[j] <- 1 - exp_trans(P[j]) + exp_trans(P[j])*F1
  }

 # pFDR_bayes

 pFDR_bayes <- 1/M*sum((exp_trans(P[-c(1:n.burnin)])*sig.level)/Fj)
 
 # pFDR frequentist 
 
 # lambda <- tunning parameter
 pi_0_freq <- (length(which(data > lambda)))/(length(data)*(1-lambda))
 F_freq <- length(which(data<sig.level))/length(data)
 pFDR_freq <- (pi_0_freq*sig.level)/F_freq

 source("https://bioconductor.org/biocLite.R")
 biocLite("qvalue")
 library(qvalue)
 
# Proportion of true null hypothesis pi_0 frequentist
qvalue(data,fdr =TRUE)$pi0

# Traceplot for pi_1
PP = data.frame(chain1,chain2) # run separately and allocated to each chain

pl <- ggplot(PP, aes(x=c(1:length(P)),y=chain2))
pl <- pl + labs(x = "Iteration", y = "pi_1") 
pl <- pl + geom_line(colour="deepskyblue4")
pl <- pl + geom_hline(yintercept = 0.4)  
pl <- pl + theme(axis.text=element_text(size=15),
                 axis.title=element_text(size=15,face="bold"))
plot(pl) 

# Traceplot for number of clusters
test_data <-
  data.frame(
    var0 = K[1:2000],
    var1 = K_gol[1:2000],
    date = c(1:2000)
  )

ggplot(test_data, aes(date)) + 
  geom_line(aes(y = var0, colour = "Neal Alg 8")) + 
  geom_line(aes(y = var1, colour = "No-gaps"))

pl1 <- ggplot(as.data.frame(K), aes(x=c(1:length(P)),y=K))
pl1 <- pl1 + labs( x = "Iteration", y = "Number of clusters") 
pl1 <- pl1 + geom_line(colour="deepskyblue4")
pl1 <- pl1 + geom_hline(yintercept = 0.4)  
pl1 <- pl1 + theme(axis.text=element_text(size=15),
                 axis.title=element_text(size=15,face="bold"))
pl1 <- pl1 + theme(plot.title = element_text(size = 25))
plot(pl1)  
  
  
# ACF plot
bacf <- acf(P_gol, lag.max = 400, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))

q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + labs(title = "No-gaps algorithm")
q

bacf1 <- acf(P, lag.max = 400 , plot = FALSE)
bacfdf1 <- with(bacf1, data.frame(lag, acf))

q1 <- ggplot(data = bacfdf1, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + labs(title = "Neal algorithm 8")
q1

grid.arrange(q,q1)


# Scatterplot of p-values
df <- data.frame(y = data, x = c(1:length(data)))  
pl2 <- ggplot(df, aes(x, y))
pl2 <- pl2 + labs( x = "Index", y = "P-value")
pl2 <- pl2 + geom_point(colour = "chartreuse4", size = 2)
pl2 <- pl2 + geom_hline(yintercept = 0.05, size = 1.5, colour = "red") 
pl2 <- pl2 + theme(axis.text=element_text(size=15),
                   axis.title=element_text(size=15,face="bold"))
plot(pl2)


# p_values histogram
his <- ggplot(df, aes(y))
his <- his + geom_histogram(bins = 10)
his <- his + labs( x = "P-values", y = "Frequency")
his <- his + theme(axis.text=element_text(size=15),
                   axis.title=element_text(size=15,face="bold"))
plot(his)

# Define a data frame PP with each column being a different chain
# HPD Interval for pi_1
hdi(exp_trans(PP),credMass = 0.90)

# Mean for each chain
apply(exp_trans(PP),2,mean)

# Visual hdi interval
caterplot(exp_trans(PP), denstrip=T, cex.labels = 0.85,
          add=TRUE, quantiles = c(0.9,0.5))

# Density plots
denplot(exp_trans(P), xlab= "Pi_1", ylab="Density", main = "Chain 1")

# Effective sample size (mean of all chains)
mean(effectiveSize(exp_trans(PP)))

# FWER
FWER <- NULL
for (i in 1:500){
  FWER[i] <- 1-(1-0.05)^i
}

# Scatterplot of FWER
Edf <- data.frame(y = FWER, x = c(1:length(FWER)))  
plE <- ggplot(Edf, aes(x, y))
plE <- plE + labs( x = "Number of hypotheses", y = "FWER")
plE <- plE + geom_line(colour = "chartreuse4", size = 3)
plE <- plE + theme(axis.text=element_text(size=15),
                   axis.title=element_text(size=15,face="bold"))
plot(plE)


# Dirichlet samples figure
alpha <- 1
draws <- 15
dimen <- 10
x <- rdirichlet(draws, rep(alpha, dimen))

dat <- data.frame(item=factor(rep(1:10,15)), 
                  draw=factor(rep(1:15,each=10)), 
                  value=as.vector(t(x)))

ggplot(dat,aes(x=item,y=value,ymin=0,ymax=value)) + 
  geom_point(colour=I("blue"))       + 
  geom_linerange(colour=I("blue"))   + 
  facet_wrap(~draw,ncol=5)           + 
  scale_y_continuous(lim=c(0,1))     +
  theme(panel.border = element_rect(fill=0, colour="black"))

