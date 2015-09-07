
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
ts0 <- twostage(aa,data=data,clusters=data$cluster,
	 detail=0,
        theta=c(2,1),var.link=0,step=0.5,
        random.design=out$des.rv,theta.des=out$pardes)
summary(ts0)

### now sample some pairs for analyses 

### first all pairs 
mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=10)
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
tail(pairs,n=12)
## make all pairs and pair specific design and pardes 
## same as ts0 but pairs specified 
ts <- twostage(aa,data=data,clusters=data$cluster,
               theta=c(2,1),var.link=0,step=0.5,
               random.design=out$des.rv,
               theta.des=out$pardes,pairs=pairs)
summary(ts)


### random sample of pairs 
ssid <- sort(sample(1:48000,20000))
###
### take some of all 
tsd <- twostage(aa,data=data,clusters=data$cluster,
               theta=c(2,1),var.link=0,step=0.5,
               random.design=out$des.rv,iid=1,
	      theta.des=out$pardes,pairs=pairs[ssid,])
summary(tsd)

### same analyses but now gives only data that  is used in the relevant pairs 
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
system.time(
tsdid <- twostage(aa,data=dataid,clusters=dataid$cluster,
               theta=c(2,1),var.link=0,step=0.5,
               random.design=outid$des.rv,iid=1,
               theta.des=outid$pardes,pairs=pair.new)
)
summary(tsdid)
### same as tsd 


### now direct specification of random.design and theta.design 

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

theta.des  <- array(0,c(nrow(pair.new),4,2))
random.des <- array(0,c(nrow(pair.new),2,4))
rvs <- c()
for (i in 1:nrow(pair.new))
{
	if (pair.types[i,1]=="mother" & pair.types[i,2]=="father")
	{
	theta.des[i,,] <- rbind(c(1,0),c(1,0),c(0,1),c(0,0))
       	random.des[i,,] <- rbind(c(1,0,1,0),c(0,1,1,0))
	rvs <- c(rvs,3)
	} else {
  	theta.des[i,,] <- rbind(c(0.5,0),c(0.5,0),c(0.5,0),c(0,1))
	random.des[i,,] <- rbind(c(1,1,0,1),c(1,0,1,1))
	rvs <- c(rvs,4)
	}
}
random.des[5,,]
theta.des[5,,]
random.des[1,,]
theta.des[1,,]

tsdid2 <- twostage(aa,data=dataid,clusters=dataid$cluster,
               theta=c(2,1),var.link=0,step=0.5,
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

out <- make.pairwise.design(pair.new,kinship,type="ace") 

tsdid3 <- twostage(aa,data=dataid,clusters=dataid$cluster,
               theta=c(2,1),var.link=0,step=0.5,
               random.design=out$random.design,
               theta.des=out$theta.des,pairs=pair.new,pairs.rvs=out$ant.rvs)
summary(tsdid3)
### same as above  tsdid2



