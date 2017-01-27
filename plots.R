 library(coda, ggplot2, MCMCpack, mcmc, HDInterval,
        gridExtra, mcmcplots, reshape2, stargazer)
  
# ------------- Plots ----------------- #  
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

pl1 <- ggplot(as.data.frame(K), aes(x=c(1:length(P)),y=K))
pl1 <- pl1 + labs( x = "Iteration", y = "Number of clusters") 
pl1 <- pl1 + geom_line(colour="deepskyblue4")
pl1 <- pl1 + geom_hline(yintercept = 0.4)  
pl1 <- pl1 + theme(axis.text=element_text(size=15),
                 axis.title=element_text(size=15,face="bold"))
pl1 <- pl1 + theme(plot.title = element_text(size = 25))
plot(pl1)  
  
  
# ACF plot
bacf <- acf(P1, lag.max = 100, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))

q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) + labs(title = "No-gaps algorithm")
q

bacf1 <- acf(P, lag.max = 100 , plot = FALSE)
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
pl2 <- pl2 + geom_point(colour = "chartreuse4", size = 3)
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
denplot(exp_trans(PP))

# Effective sample size (mean of all chains)
mean(effectiveSize(exp_trans(PP)))

# FWER
FWER <- NULL
for (i in 1:500){
  FWER[i] <- 1-(1-0.05)^i
}

# Scatterplot of p-values
Edf <- data.frame(y = FWER, x = c(1:length(FWER)))  
plE <- ggplot(Edf, aes(x, y))
plE <- plE + labs( x = "Number of hypotheses", y = "FWER")
plE <- plE + geom_line(colour = "chartreuse4", size = 3)
plE <- plE + theme(axis.text=element_text(size=15),
                   axis.title=element_text(size=15,face="bold"))
plot(plE)

# Dirichlet Process Plot

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
