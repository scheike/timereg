
data(prt)
dim(prt)
table(prt$status)
### 21000 7000 1000

library(mets)
###
set.seed(100)
prt<-simnordic(7500,cordz=3,cormz=4,cratemz=0.1,cratedz=0.1)
prt$status <-prt$cause
table(prt$status)

prt<-simnordic(75000,cordz=3.0,cormz=5.,pcensmz=0.0,pcensdz=0.0,cratemz=200.4,cratedz=100.4)
prt$status <-prt$cause
table(prt$status)
prt$cancer <- (prt$status==1)

prt<-simnordic(7500,cordz=3,cormz=4,pcensmz=0.9,pcensdz=0.9,cratemz=1.0,cratedz=1.0)
prt$status <-prt$cause
table(prt$status)


prt<-simnordic(75000,cordz=3.3,cormz=5.,pcensmz=0.0,pcensdz=0.0,cratemz=200.4,cratedz=100.4)
prt$status <-prt$cause
table(prt$status)
prt$cancer <- (prt$status==1)
###
bp3 <- bptwin(cancer~1,zyg="zyg",DZ="DZ",id="id",type="ace",data=prt)
summary(bp3)

coef(bp3)
exp(-0.1)/( exp(-0.1)+exp(-0.77)+1)
exp(-0.77)/( exp(-0.1)+exp(-0.77)+1)

gem <- c()
for (pcens in seq(0,0.95,length=5))
{
prt<-simnordic(7500,cordz=3,cormz=4,pcensmz=pcens,pcensdz=pcens,cratemz=1.0,cratedz=1.0)
prt$status <-prt$cause
tt <- table(prt$status)
if (length(tt)==2) tt <- c(0,tt)
prt$cancer <- (prt$status==1)
###
bp3 <- bptwin(cancer~1,zyg="zyg",DZ="DZ",id="id",type="ace",data=prt)
summary(bp3)
###
gem <- rbind(gem,c(pcens,tt,coef(bp3)))
}
gem


gemmzdz <- c()
cens <- gemmzdzc <- matrix(0,5,5)
gemmzdza <- matrix(0,5,5)
j <- i <- 0 
for (pcensmz in seq(0,0.95,length=5))
{
i <- i+1
j <- 0
for (pcensdz in seq(0,0.95,length=5))
{
j <- j+1
prt<-simnordic(100000,cordz=3.5,cormz=5,pcensmz=pcensmz,pcensdz=pcensdz,cratemz=1.0,cratedz=1.0)
prt$status <-prt$cause
tt <- table(prt$status)
if (length(tt)==2) tt <- c(0,tt)
prt$cancer <- (prt$status==1)
###
bp3 <- bptwin(cancer~1,zyg="zyg",DZ="DZ",id="id",type="ace",data=prt)
summary(bp3)
###
gemmzdz <- rbind(gemmzdz,c(pcensmz,pcensdz,tt,coef(bp3)))
print(c(i,j))
gemmzdza[i,j] <- coef(bp3)[2]
gemmzdzc[i,j] <- coef(bp3)[3]
cens[i,j] <- tt[1]/sum(tt)
}
}
###
gemmzdz
gemmzdza
gemmzdzc
###
h <- exp(gemmzdza)/(exp(gemmzdza)+exp(gemmzdzc)+1)/h[1,1]
c <- exp(gemmzdzc)/(exp(gemmzdza)+exp(gemmzdzc)+1)/c[1,1]
###
round(h,2)
round(c,2)
round(cens,2)


table(prt$zyg,prt$status)

