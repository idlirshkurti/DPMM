# Script for pre-processing the Golub et al. (1999) ALL AML training dataset

# ---------------------------------------------------------
# Get data from Whitehead Institute website
# ---------------------------------------------------------

URL<-"http://www-genome.wi.mit.edu/mpr/publications/projects/Leukemia/data_set_ALL_AML_train.txt"	 
golub.all<-read.table(URL,sep="\t",quote="",header=TRUE,row.names=NULL,comment.char="")

#### Exploratory Min-Max plotting ####

max.var<-apply(golub.all, 2, max)
min.var<-apply(golub.all, 2, min)
plot(max.var, ylim=c(-5,5), type='l', col='red')
lines(min.var, col='blue') # Outliers, reason to normalise/Scale?

min.df<-factor(rep('min',times=length(max.var)))

max.df<-factor(rep('max',times=length(max.var)))
iterations<-seq(1:length(max.var))
max.var.df<-as.vector(max.var)
min.var.df<-as.vector(min.var)
max.var.df2<-data.frame(Function=max.df, Variable=iterations, Value=max.var.df)
min.var.df2<-data.frame(Function=min.df, Variable=iterations, Value=min.var.df)
df.it<-rbind(max.var.df2, min.var.df2)

library(ggplot2)
#### FIX THE GGPLOT!!! #####
ggplot()+geom_line(data=max.var.df2, aes(x=Variable, y=Value), colour='dodgerblue4',
                   linetype=1)+
  geom_line(data=min.var.df2,aes(x=Variable, y=Value),color='darkred') +
  theme(text = element_text(size=16))

# ---------------------------------------------------------
# Reformat and tidy up
# ---------------------------------------------------------

# Gene names and tumor class labels
golub.gnames<-cbind(dimnames(golub.all)[[1]],
                    as.character(golub.all[,1]),as.character(golub.all[,2]))
golub.cl<-c(rep(1,27),rep(2,11))

# Re-order columns
golub<-golub.all
golub<-golub[,1+2*(1:38)]
golub<-golub[,c(1:27,33:38,28:32)]
golub<-as.matrix(golub)

# Floor & ceiling
golub[golub<100]<-100
golub[golub>16000]<-16000

# Preliminary selection of genes
tmp1<-apply(golub,1,max)
tmp2<-apply(golub,1,min)
which1<-(1:7129)[(tmp1/tmp2)>5]
which2<-(1:7129)[(tmp1-tmp2)>500]
golub.sub<-intersect(which1,which2)
golub<-golub[golub.sub,]

# Log_10 transformation
golub<-log(golub,10)

# Normalization	
golub.expr<-scale(golub,TRUE,TRUE)
dimnames(golub.expr)<-list(NULL,NULL)

#export to multtest
golub<-golub.expr
golub.cl<-c(rep(0,27),rep(1,11))
golub.gnames<-golub.gnames[golub.sub,]
