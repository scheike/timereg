\name{designX}
\alias{designX}
\alias{design}
\alias{designZ}
\title{Specifies design of Cox-Aalen survival or competing risks model.}
\description{Specifies design of Cox-Aalen survival model or the competing risks regression model.

Basic intensity model is on the form
\deqn{ \lambda_{i}(t) = Y_i(t) 
( X(h)^T(t) \alpha(t) )  \exp(Z(h)^T \beta ) 
}
where \eqn{X(h)} and \eqn{Z(h)} are the design vectors
specified by the design functionsi given observed covariates 
\eqn{X} and \eqn{Z} and the unobserved haplotype \eqn{h}.
Note that Cox's regression model is when \eqn{X(h)=1}. 
 
The competing risks model is on the form :
\deqn{
P(T \leq t, cause=1 | X,Z) = P_1(t,X,Z) = 1 - \exp(- g(t,X(h),Z(h)) )
}
for the probability of dying from cause 1 in a situation with 
competing causes of death where  
\eqn{X(h)}  and \eqn{Z(h)}  
are design vectors that depends on the unobserved haplotype. 


Covariates X and Z are given:

designX(x,z,h) is an R-function that returns the X-design of the Cox-Aalen 
model, given x,z and haplotype pair \eqn{h=(h_1,h_2)}. 

designZ(x,z,h) is an R-function that returns the Z-design of the Cox-Aalen 
model, given x,z and  haplotype pair \eqn{h=(h_1,h_2)}. 

h refers to the integer code for haplotypes that can be read of in the
geno.setup() function applied to the genotype data.

To study the effect (relative risk) of observing haplotypes h* (integer code).
We want designZ(x,z,h) = \eqn{(x,I(h_1=h*)+I(h_2=h*))}
thus counting the occurrence of h* in the haplotype pair.

Let us say that the haplotype h* has code 5 in the geno.type() listing, then
the above design is given as follows : 

designX<-function(x,z,h)
\{
h<-round(h);  
return(x)
\}

designZ<-function(x,z,h)
\{
h<-round(h); 
y<-c(z,(h[1]==(5-1))+(h[2]==(5-1))); 
return(y)
\}

Note that x,z and refers to the observed risk covariates for the baseline and
the relative risk parameters, respectively. 

For the matched haplotype analyses with donor and patient the design could 
be on the form: 

designX<-function(x,z,hd,hp)
\{
hd<-round(hd); hp<-round(hp);  

return(x)
\}

designZ<-function(x,z,hd,hp)
\{
hd<-round(hd); hp<-round(hp)   

y<-c(z,(hd[1]==hp[1])*(hp[2]==hd[2])*1)

return(y)
\}

that here gives the effect of a full match in the haplotype pairs, 
with a non-parametric baseline.
}
\author{Thomas Scheike}
\keyword{survival}
