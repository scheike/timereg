
library(mets)
set.seed(100)
### use of clayton oakes binomial additive gamma model
###########################################################
data <- simbinClaytonOakes.family.ace(10000,2,1,beta=NULL,alpha=NULL)
str(data)
head(data$cluster,n=200)
aa <- margbin <- glm(ybin~x,data=data,family=binomial())
ps <- predict(margbin,newdata=data,type="response")
margbin
head(data)
data$number <- c(1,2,3,4)
data$child <- 1*(data$number==3)
head(data)

### make ace random effects design
out <- ace.family.design(data,member="type",id="cluster")
out$pardes
head(out$des.rv)

### now specify fitting via specific pairs 
source("../R/binomial.twostage.R")
ts <- binomial.twostage(margbin,data=data,clusters=data$cluster,
          theta=c(2,1),var.link=0,step=1.0,Nit=10,detail=1,
          random.design=out$des.rv,
          theta.des=out$pardes)
summary(ts)

### first all pairs 
###cluster.index(data$cluster)
mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=20)
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
tail(pairs,n=12)
head(pairs,n=12)
## make all pairs and pair specific design and pardes 
## same as ts0 but pairs specified 

ts <- binomial.twostage(margbin,data=data,clusters=data$cluster,
               theta=c(2,1),var.link=0,step=1.0,Nit=10,detail=1,
               random.design=out$des.rv,
               theta.des=out$pardes,pairs=pairs)
ts$theta
summary(ts)

###source("../R/binomial.twostage.R"); 
### random sample of pairs 
set.seed(100)
ssid <- sort(sample(1:60000,40000))
###
### take some of all 
tsd <- binomial.twostage(aa,data=data,clusters=data$cluster,
               theta=c(2,1),var.link=0,step=1.0,
               random.design=out$des.rv,iid=1,Nit=10,
	      theta.des=out$pardes,pairs=pairs[ssid,])
summary(tsd)
tsd$score
tsd$Dscore


### same analyses but now gives only data that is used in the relevant pairs 
head(pairs[ssid,])
ids <- sort(unique(c(pairs[ssid,])))
###
pairsids <- c(pairs[ssid,])
pair.new <- matrix(fast.approx(ids,c(pairs[ssid,])),ncol=2)
head(pair.new)

## this requires that pair.new refers to id's in dataid (survival, status and so forth)
## random.design and theta.des are constructed to be the array 3 dims via individual specfication from ace.family.design
dataid <- dsort(data[ids,],"cluster")
outid <- ace.family.design(dataid,member="type",id="cluster")
outid$pardes
head(outid$des.rv)
###
tsdid <- binomial.twostage(aa,data=dataid,clusters=dataid$cluster,
               theta=c(2,1),var.link=0,step=1.0,
               random.design=outid$des.rv,Nit=10,
               theta.des=outid$pardes,pairs=pair.new)
summary(tsdid)

###
tsdid$score
tsdid$Dscore
### same as tsd 


### now direct specification of random.design and theta.design 
### rather than taking the rows of the des.rv for the relevant pairs
### can make a pair specific specification of random effects 

pair.types <-  matrix(dataid[c(t(pair.new)),"type"],byrow=T,ncol=2)
head(pair.new)
head(pair.types)

### here makes pairwise design , simpler random.design og pardes, parameters
### stil varg, varc 
### mother, child, share half rvm=c(1,1,0) rvc=c(1,0,1),
### thetadesmcf=rbind(c(0.5,0),c(0.5,0),c(0.5,0),c(0,1))
###
### father, child, share half rvf=c(1,1,0) rvc=c(1,0,1), 
### thetadescf=rbind(c(0.5,0),c(0.5,0),c(0.5,0),c(0,1))
###
### child, child,  share half rvc=c(1,1,0) rvc=c(1,0,1),
### thetadesmf=rbind(c(0.5,0),c(0.5,0),c(0.5,0),c(0,1))
###
### mother, father, share 0 rvm=c(1,0) rvf=c(0,1), 
### thetadesmf=rbind(c(1,0),c(1,0),c(0,1))

theta.des  <- array(0,c(4,2,nrow(pair.new)))
random.des <- array(0,c(2,4,nrow(pair.new)))
### random variables in each pair 
rvs <- c()
for (i in 1:nrow(pair.new))
{
	if (pair.types[i,1]=="mother" & pair.types[i,2]=="father")
	{
	theta.des[,,i] <- rbind(c(1,0),c(1,0),c(0,1),c(0,0))
       	random.des[,,i] <- rbind(c(1,0,1,0),c(0,1,1,0))
	rvs <- c(rvs,3)
	} else {
  	theta.des[,,i] <- rbind(c(0.5,0),c(0.5,0),c(0.5,0),c(0,1))
	random.des[,,i] <- rbind(c(1,1,0,1),c(1,0,1,1))
	rvs <- c(rvs,4)
	}
}
### 3 rvs here 
random.des[,,7]
theta.des[,,7]
### 4 rvs here 
random.des[,,1]
theta.des[,,1]
head(rvs)

tsdid2 <- binomial.twostage(aa,data=dataid,clusters=dataid$cluster,
               theta=c(2,1),var.link=0,step=1.0,
               random.design=random.des,
               theta.des=theta.des,pairs=pair.new,pairs.rvs=rvs)
summary(tsdid2)
tsd$theta
tsdid$theta


### simpler specification via kinship coefficient for each pair

kinship  <- c()
for (i in 1:nrow(pair.new))
{
if (pair.types[i,1]=="mother" & pair.types[i,2]=="father") pk1 <- 0 else pk1 <- 0.5
kinship <- c(kinship,pk1)
}
head(kinship,n=10)

out <- make.pairwise.design(pair.new,kinship,type="ace") 
names(out)
### 4 rvs here , here independence since shared component has variance 0 !
out$random.des[,,9]
out$theta.des[,,9]

tsdid3 <- binomial.twostage(aa,data=dataid,clusters=dataid$cluster,
               theta=c(2,1),var.link=0,step=1.0,
               random.design=out$random.design,
               theta.des=out$theta.des,pairs=pair.new,pairs.rvs=out$ant.rvs)
summary(tsdid3)

### same as above  tsdid2


bin <- binomial.twostage(margbin,data=data,clusters=data$cluster,detail=0,
 model="clayton.oakes",
 theta=-1.3,var.link=1,step=1.0)
summary(bin)

bin <- binomial.twostage(margbin,data=data,clusters=data$cluster,detail=0,
 model="clayton.oakes",
 theta=0.2,var.link=0,step=1.0)
summary(bin)



