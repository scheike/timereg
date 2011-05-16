data(pafah)

geno <- as.matrix(pafah[,9:18])
setup <- geno.setup(geno);

## names of haplotypes
setup$uniqueHaploNames

setup <- geno.setup(geno,haplo.baseline="C,C,R,I,A") 

out<-haplo.freqs(geno,Nit=100,geno.setup=setup,step=0.2); 
out$haplo.freq ## approximate MLE
out$score

## MLE hard to find and unstable therefore group frequencies 

# Let's group together all those haplotypes with an initial frequency
# estimate less than 0.02
guesFreq <- out$haplo.freq

sum(guesFreq[-32] < 0.02) # there are 19 of these
sum(guesFreq[-32] >= 0.02) # and 12 other haplotypes,

index.small<-(1:32)[guesFreq < 0.02] # there are 26 of these
index.small
index.large<-(1:32)[guesFreq >= 0.02][-6] # there are 5 of these
index.large

# everything compared to C,C,R,I,A, 6 paramaters

## 6 parameters, create the design matrix 
X <- matrix(0,31,6)
X[index.small,1]<-1; 
k<-0; 
for (i in index.large) 
{
k<-k+1; 
X[i,(1:length(index.large))[k]+1]<-1
}


# haplotype names going along with permutaions 
# largest frequency is baseline here 
setup$uniqueHaploNames[index.small]
setup$uniqueHaploNames[index.large]
setup$uniqueHaploNames[32]

perm.names<-c("small",setup$uniqueHaploNames[index.large])
perm.names
colnames(X)<-perm.names
rownames(X)<-setup$uniqueHaploNames[-32]
X ## design for haplotype frequencies

hapfit<-haplo.freqs(geno,Nit=100,geno.setup=setup,
		 haplo.design=X,step=0.2); 
hapfit
hapfit$score
hapfit$haplo.freq ## MLE estimates for structured haplo model 

## {{{ designs

designX<-function(x,z,h) { return(x)}

designZ<-function(x,z,h) {
h<-round(h);
vecZ<-c()
for (i in (c(1,6,8,20,21)-1))  # first component as baseline
{
vecZ<-c(vecZ,c((h[1]==i)+(h[2]==i)))
}
y<-c(vecZ)
return(y)
}
## }}}

##############################################################
## running proportional model 
##############################################################
## {{{

dummy<-rep(1,nrow(geno))
paf1<-haplo.surv(Surv(time,status)~1+prop(dummy),data=pafah,
designX,designZ,Nit=5,detail=0,start.time=0,n.sim=500,
geno.type=geno,geno.setup=setup,
haplo.freq=hapfit$haplo.freq,haplo.design=X,
step=0.1,two.stage=1,covnamesZ=colnames(X)[-1])
paf1$score
summary(paf1)

## effects of specific haplotype relative to other types
## other types is a mix of the most frequent and rare types
## somewhat strange model 

plot(paf1,xlab="Time (years)",ylab="Cumulative baseline",sim.ci=2)
par(mfrow=c(2,3))
plot(paf1,score=1,xlab="Time (years)",ylab="Score process")

##}}}

