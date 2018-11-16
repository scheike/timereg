
library(mets)
###
set.seed(1000)
data <- simClaytonOakes.family.ace(8000,2,1,0,3)
head(data)
data$number <- c(1,2,3,4)
data$child <- 1*(data$number==3)
out <- ace.family.design(data,member="type",id="cluster")
### 8 random effects with 1/4 * var.gene, and one shared environment  1 * var.env
out$pardes
head(out$des.rv)
###
aa <- phreg(Surv(time,status)~+cluster(cluster),data=data)

## {{{ additive gamma models with and without pair call 
### make ace random effects design

### simple random effects call 
ts0 <- survival.twostage(aa,data=data,clusters=data$cluster,
	detail=0,var.par=0,var.link=0,
        theta=c(2,1),
        random.design=out$des.rv,theta.des=out$pardes)
summary(ts0)

### simple random effects call 
ts1 <- survival.twostage(aa,data=data,clusters=data$cluster,
	detail=0,var.par=1,var.link=0,
        theta=c(2,1),
        random.design=out$des.rv,theta.des=out$pardes)
summary(ts1)
### parameters c(2,1)/(2+1)^2

checkderiv=1
if (checkderiv==1) {# {{{

ts0 <- twostage(aa,data=data,clusters=data$cluster,
	detail=1,numDeriv=1,Nit=0,var.par=1,
        theta=log(c(2,1)/9),var.link=1,step=1.0,
        random.design=out$des.rv,theta.des=out$pardes)
ts0$score
ts0$score1

ts0 <- twostage(aa,data=data,clusters=data$cluster,
	detail=1,numDeriv=1,Nit=0,var.par=1,
        theta=c(2,1)/9,var.link=0,step=1.0,
        random.design=out$des.rv,theta.des=out$pardes)
ts0$score
ts0$score1


ts0 <- twostage(aa,data=data,clusters=data$cluster,
	detail=1,numDeriv=1,Nit=0,var.par=0,
        theta=log(c(2,1)),var.link=1,step=1.0,
        random.design=out$des.rv,theta.des=out$pardes)
ts0$score
ts0$score1

ts0 <- twostage(aa,data=data,clusters=data$cluster,
	detail=1,numDeriv=1,Nit=0,var.par=0,
        theta=c(2,1),var.link=0,step=1.0,
        random.design=out$des.rv,theta.des=out$pardes)
ts0$score
ts0$score1

}# }}}


### now specify fitting via specific pairs 
### first construct all pairs 
mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=10)
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
tail(pairs,n=12)
## make all pairs and pair specific design and pardes 
## same as ts0 but pairs specified 
ts <- twostage(aa,data=data,clusters=data$cluster,
               theta=c(2,1)/9,var.link=0,step=1.0,
               random.design=out$des.rv,
               theta.des=out$pardes,pairs=pairs)
summary(ts)

### random sample of pairs 
ssid <- sort(sample(1:32000,20000))
###
### take some of all 
tsd <- twostage(aa,data=data,clusters=data$cluster,
               theta=c(2,1)/10,var.link=0,step=1.0,
               random.design=out$des.rv,iid=1,
	      theta.des=out$pardes,pairs=pairs[ssid,])
summary(tsd)

str(aa)
aa$id


### same analyses but now gives only data that is used in the relevant pairs 
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
tsid <- twostage(aa,data=dataid,clusters=dataid$cluster,
               theta=c(2,1)/10,var.link=0,step=1.0,
               random.design=outid$des.rv,iid=1,
               theta.des=outid$pardes,pairs=pair.new)
summary(tsdid)
coef(tsdid)
coef(tsd)
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

tsdid2 <- twostage(aa,data=dataid,clusters=dataid$cluster,
               theta=c(2,1)/10,var.link=0,step=1.0,
               random.design=random.des,
               theta.des=theta.des,pairs=pair.new,pairs.rvs=rvs)
summary(tsdid2)
tsd$theta
tsdid2$theta
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


tsdid3 <- twostage(aa,data=dataid,clusters=dataid$cluster,
               theta=c(2,1)/10,var.link=0,step=1.0,
               random.design=out$random.design,
               theta.des=out$theta.des,pairs=pair.new,pairs.rvs=out$ant.rvs)
summary(tsdid3)
coef(tsdid3)

### same as above  tsdid2

## }}} 

##### simple models, test for pairs structure ## {{{ 

library(mets)
source("../R/twostage.R")
ts0 <- twostage(aa,data=data,clusters=data$cluster,
	detail=0,numDeriv=1,Nit=10,
        theta=c(0.17),var.link=0,step=1.0)
summary(ts0)
ts0$score; ts0$score1
ts0$Dscore; ts0$hess

mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=10)
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
head(pairs,n=12)
tail(pairs,n=12)
dim(pairs)
#
cc <- cluster.index(data$cluster)
###
ts0 <- twostage(aa,data=data,clusters=data$cluster,
        detail=1,Nit=0,
        theta=ts0$theta,var.link=0,pairs=pairs)
summary(ts0)


## {{{  simple models with pair call 

library(mets)

set.seed(100)
data <- simClaytonOakes.family.ace(8000,2,1,0,3)
head(data)
data$number <- c(1,2,3,4)
data$child <- 1*(data$number==3)

### make ace random effects design
out <- ace.family.design(data,member="type",id="cluster")
out$pardes
head(out$des.rv)

### makes marginal model (same for all) 
aa <- aalen(Surv(time,status)~+1,data=data,robust=0)


mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=10)
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
head(pairs,n=12)
tail(pairs,n=12)
dim(pairs)
#

ts0 <- twostage(aa,data=data,clusters=data$cluster,
	 detail=1,Nit=10,
        theta=c(0.2),var.link=0,step=1.0)
summary(ts0)

ts0 <- twostage(aa,data=data,clusters=data$cluster,
	 detail=1,Nit=10,numDeriv=1,
        theta=c(0.2),var.link=0,step=1.0,pairs=pairs)
summary(ts0)
ts0$score
ts0$score1

ts0 <- twostage(aa,data=data,clusters=data$cluster,
	 detail=1,Nit=10,
        theta=c(0.2),var.link=0,step=1.0,model="plackett")
summary(ts0)

ts0 <- twostage(aa,data=data,clusters=data$cluster,
	 detail=1,Nit=10,
        theta=c(0.2),var.link=0,step=1.0,model="plackett",pairs=pairs)
summary(ts0)



theta.des <- model.matrix(~x1,data=data)

ts0 <- twostage(aa,data=data,clusters=data$cluster,
	 detail=1,Nit=10,theta.des=theta.des,
        theta=c(0.2),var.link=0,step=1.0)
summary(ts0)

ts0 <- twostage(aa,data=data,clusters=data$cluster,
	 detail=1,Nit=10,theta.des=theta.des,
        theta=c(0.2),var.link=0,step=1.0,pairs=pairs)
summary(ts0)

ts0 <- twostage(aa,data=data,clusters=data$cluster,
	 detail=1,Nit=10,theta.des=theta.des,
        theta=c(0.2),var.link=0,step=1.0,model="plackett")
summary(ts0)

ts0 <- twostage(aa,data=data,clusters=data$cluster,
	 detail=1,Nit=10,theta.des=theta.des,
        theta=c(0.2),var.link=0,step=1.0,model="plackett",pairs=pairs)
summary(ts0)



## }}} 


## }}} 



