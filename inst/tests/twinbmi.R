#
# BMI data on twin pairs
#

library(mets)

data(twinbmi)

str(twinbmi)
head(twinbmi)

# restrict to data, where response is not missing
twinbmi <- twinbmi[!is.na(twinbmi$bmi),]

# install.packages("lattice")
library(lattice)
plot( histogram( ~ bmi| gender, type="density", col="red", xlab="kg/m^2", 
				main="Histogram of BMI", data=twinbmi) )

# BMI is often studied on log-scale.
# boxcox(bmi ~ age*gender, data = twinbmi)
twinbmi$logbmi <- log(twinbmi$bmi)

# Saturated model
lnbmi.sat <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="sat",control=list(method="NR"))
lnbmi.sat$estimate$opt$message
lnbmi.sat

#
lnbmi.flex <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="flex",control=list(method="NR"))
lnbmi.flex$estimate$opt$message
lnbmi.flex

compare(lnbmi.sat,lnbmi.flex)

#
lnbmi.u <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="u",control=list(method="NR"))
lnbmi.u$estimate$opt$message
lnbmi.u

compare(lnbmi.u,lnbmi.flex)

#
lnbmi.ace <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="ae",control=list(trace=1))
lnbmi.ace$estimate$opt$message
score(lnbmi.ace$estimate)
lnbmi.ace2 <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="ace",control=list(method="NR",trace=1,lambda=0.1,constrain=TRUE))
lnbmi.ace$estimate$opt$message
lnbmi.ace   

#
lnbmi.ade <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="ade",control=list(method="NR"))
lnbmi.ade$estimate$opt$message
lnbmi.ade  

AIC(lnbmi.ace,lnbmi.ade)

#
lnbmi.ae <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="ae",control=list(method="NR"))
lnbmi.ae$estimate$opt$message
lnbmi.ae  

compare(lnbmi.ace,lnbmi.ae)

#CE
lnbmi.ce <- twinlm(logbmi~age*gender, id="tvparnr", DZ="DZ", zyg="zyg",data=twinbmi,
			type="ce",control=list(method="NR"))
lnbmi.ce$estimate$opt$message
lnbmi.ce  

AIC(lnbmi.ace,lnbmi.ce)


# GOF-Table?
# mx and openmx for same data
# reshape wide -  twin-twin plot.
