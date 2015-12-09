

library(mets)

 data(twinstut)
 b0 <- bptwin(stutter~sex,
              data=droplevels(subset(twinstut,zyg%in%c("mz","dz"))),
              id="tvparnr",zyg="zyg",DZ="dz",type="ae")
 summary(b0)

library(mets)
data <- simbinClaytonOakes.twin.ace(40000,2,1,beta=NULL,alpha=NULL)
out <- polygen.design(data,id="cluster",zygname="zygosity")
out$pardes
head(out$des.rv)
margbin <- glm(ybin~x,data=data,family=binomial())
###
system.time(
bintwin <- binomial.twostage(margbin,data=data,
     clusters=data$cluster,detail=0,var.par=1,
     theta=c(2,1)/9,var.link=0,
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


