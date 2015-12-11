library(mets)

### example 1, twin stut data 
data(twinstut)
twinstut <- subset(twinstut,zyg%in%c("mz","dz"))
twinstut$binstut <- 1*(twinstut$stutter=="yes")

b0 <- bptwin(binstut~sex,data=twinstut,id="tvparnr",zyg="zyg",DZ="dz",type="ae")
summary(b0)

out <- polygen.design(twinstut,id="tvparnr",zygname="zyg",zyg="dz",type="ae")
margbin <- glm(binstut~sex,data=twinstut,family=binomial())
bintwin <- binomial.twostage(margbin,data=twinstut,
     clusters=twinstut$tvparnr,detail=0,theta=c(0.1)/1,var.link=0,
     random.design=out$des.rv,theta.des=out$pardes)
summary(bintwin)

concordance.twin.ace(bintwin,type="ae")
summary(b0)


### example 2, simulated data  
data <- simbinClaytonOakes.twin.ace(10000,0.5,0.5,beta=NULL,alpha=NULL)
out <- polygen.design(data,id="cluster",zygname="zygosity")
out$pardes
head(out$des.rv)
tail(out$des.rv)
margbin <- glm(ybin~x,data=data,family=binomial())
p=exp(0.5)/(1+exp(0.5))
###
system.time(
bintwin <- binomial.twostage(margbin,data=data,
     clusters=data$cluster,detail=0,var.par=1,
     theta=c(0.5,0.5)/1,var.link=0,
     random.design=out$des.rv,theta.des=out$pardes)
)
summary(bintwin)

system.time(
b0 <- bptwin(ybin~x,
             data=data,
             id="cluster",zyg="zygosity",DZ="DZ",type="ace")
)
summary(b0)
summary(bintwin)

rv1 <- out$des.rv[c(1,nrow(data)-1),]
rv2 <- out$des.rv[c(2,nrow(data)),]
theta.des <- out$pardes

### concordance from additive gamma model
concordance.twostage(bintwin$theta,rep(p,2),rv1,rv2,theta.des,var.par=1)


