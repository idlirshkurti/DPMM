N <- 30
N1 <- 1000
Dn = c()
CT.2 = CT.15 = CT.1 = CT.05 = CT.01 = c()

for(i in 4:N){
for(j in 1:N1){
S<-rnorm(i)
# We create two loops. One for the increase in sample size and one for the amount of samples
# needed for each sample size. Also we create a random normal sample S.
ord <- sort(S)
SD <- sd(S); m <- mean(S)

Ecdf <- ecdf(S)
F0 <- function(x) return(pnorm(x,mean=m,sd=SD))
DnStarP <- max(Ecdf(ord) - F0(ord))
DnStarM <- max(F0(ord) - c(0,Ecdf(ord)[1:(i-1)]))
DnStar <- max(DnStarP,DnStarM)
# This calculates the value of the KS statistic for the sample S.
Dn[j] <- DnStar
}

ECDF <- ecdf(Dn)
# Function cc gives the critical values of the ecdf of Dn at a specific significance level a
cc <- function(a)for(x in Dn){if(ECDF(x)==a)return(x)}

CT.2[i] <- cc(0.8)
CT.15[i] <- cc(0.85)
CT.1[i] <- cc(0.9)
CT.05[i] <- cc(0.95)
CT.01[i] <- cc(0.99)
dim(CT.2[4:30]) <- dim(CT.15[4:30]) <- dim(CT.1[4:30]) <- dim(CT.05[4:30]) <- dim(CT.01[4:30]) <- c(27,1)
CT <- cbind(CT.2[4:30],CT.15[4:30],CT.1[4:30],CT.05[4:30],CT.01[4:30])
}
CT
