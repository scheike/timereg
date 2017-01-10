##' Fits Clayton-Oakes or bivariate Plackett (OR) models for binary data 
##' using marginals that are on logistic form. 
##' If clusters contain more than two times, the algoritm uses a compososite likelihood
##' based on all pairwise bivariate models.
##'
##' The pairwise pairwise odds ratio model provides an alternative to the alternating logistic
##' regression (ALR).
##'
##' The reported standard errors are based on a cluster corrected score equations from the 
##' pairwise likelihoods assuming that the marginals are known. This gives correct standard errors
##' in the case of the Plackett distribution (OR model for dependence), but incorrect standard
##' errors for the Clayton-Oakes types model. For the additive gamma version of the stanard errors
##' are adjusted for the uncertainty in the marginal models via an iid deomposition using the iid() function of
##' lava. For the clayton oakes model that is not speicifed via the random effects these can be 
##' fixed subsequently using the iid influence functions for the marginal model, but typically this does not
##' change much. 
##'
##' For the Clayton-Oakes version of the model, given the gamma distributed random effects it is 
##' assumed that the probabilities are indpendent, and that the marginal survival functions are on logistic form 
##' \deqn{
##' logit(P(Y=1|X)) = \alpha + x^T \beta
##' }
##' therefore conditional on the random effect the probability of the event is 
##' \deqn{
##' logit(P(Y=1|X,Z)) = exp( - Laplace^{-1}(lamtot,lamtot,P(Y=1|x)) )  
##' }
##'
##' Can also fit a structured additive gamma random effects model, such
##' the ACE, ADE model for survival data: 
##'
##' Now random.design specificies the random effects for each subject within a cluster. This is
##' a matrix of 1's and 0's with dimension n x d.  With d random effects. 
##' For a cluster with two subjects, we let the random.design rows be 
##'  \eqn{v_1} and \eqn{v_2}. 
##' Such that the random effects for subject 
##' 1 is \deqn{v_1^T (Z_1,...,Z_d)}, for d random effects. Each random effect
##' has an associated parameter \eqn{(\lambda_1,...,\lambda_d)}. By construction
##' subjects 1's random effect are Gamma distributed with 
##' mean \eqn{\lambda_j/v_1^T \lambda}
##' and variance \eqn{\lambda_j/(v_1^T \lambda)^2}. Note that the random effect 
##' \eqn{v_1^T (Z_1,...,Z_d)} has mean 1 and variance \eqn{1/(v_1^T \lambda)}.
##' It is here asssumed that  \eqn{lamtot=v_1^T \lambda} is fixed over all clusters
##' as it would be for the ACE model below.
##'
##' The DEFAULT parametrization uses the variances of the random effecs (var.par=1)
##' \deqn{
##' \theta_j  = \lambda_j/(v_1^T \lambda)^2
##' }
##'
##' For alternative parametrizations (var.par=0) one can specify how the parameters relate 
##' to \eqn{\lambda_j} with the function 
##'
##' Based on these parameters the relative contribution (the heritability, h) is 
##' equivalent to  the expected values of the random effects  \eqn{\lambda_j/v_1^T \lambda}
##'
##' Given the random effects the probabilities  are independent and on the form 
##' \deqn{
##' logit(P(Y=1|X)) = exp( - Laplace^{-1}(lamtot,lamtot,P(Y=1|x)) )  
##' }
##' with the inverse laplace of the gamma distribution with mean 1 and variance lamtot.
##'
##' The parameters \eqn{(\lambda_1,...,\lambda_d)}
##' are related to the parameters of the model
##' by a regression construction \eqn{pard} (d x k), that links the \eqn{d} 
##' \eqn{\lambda} parameters
##' with the (k) underlying \eqn{\theta} parameters 
##' \deqn{
##' \lambda = theta.des  \theta 
##' }
##' here using theta.des to specify these low-dimension association. Default is a diagonal matrix. 
##'
##' @export
##' @aliases binomial.twostage binomial.twostage.time
##' @references
##' Two-stage binomial modelling 
##' @examples
##' library("timereg")
##' data("twinstut",package="mets")
##' twinstut0 <- subset(twinstut, tvparnr<2300000)
##' twinstut <- twinstut0
##' twinstut$binstut <- (twinstut$stutter=="yes")*1
##' theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
##' margbin <- glm(binstut~factor(sex)+age,data=twinstut,family=binomial())
##' bin <- binomial.twostage(margbin,data=twinstut,var.link=1,
##' 		         clusters=twinstut$tvparnr,theta.des=theta.des,detail=0,
##' 	                 score.method="fisher.scoring")
##' summary(bin)
##' 
##' twinstut$cage <- scale(twinstut$age)
##' theta.des <- model.matrix( ~-1+factor(zyg)+cage,data=twinstut)
##' bina <- binomial.twostage(margbin,data=twinstut,var.link=1,
##' 		         clusters=twinstut$tvparnr,theta.des=theta.des)
##' summary(bina)
##' 
##' theta.des <- model.matrix( ~-1+factor(zyg)+factor(zyg)*cage,data=twinstut)
##' bina <- binomial.twostage(margbin,data=twinstut,var.link=1,
##' 		         clusters=twinstut$tvparnr,theta.des=theta.des)
##' summary(bina)
##' 
##' ## refers to zygosity of first subject in eash pair : zyg1
##' ## could also use zyg2 (since zyg2=zyg1 within twinpair's))
##' out <- easy.binomial.twostage(stutter~factor(sex)+age,data=twinstut,
##'                           response="binstut",id="tvparnr",var.link=1,
##' 	             	      theta.formula=~-1+factor(zyg1))
##' summary(out)
##' 
##' ## refers to zygosity of first subject in eash pair : zyg1
##' ## could also use zyg2 (since zyg2=zyg1 within twinpair's))
##' desfs<-function(x,num1="zyg1",num2="zyg2")
##'     c(x[num1]=="dz",x[num1]=="mz",x[num1]=="os")*1
##'
##' out3 <- easy.binomial.twostage(binstut~factor(sex)+age,
##'       data=twinstut,response="binstut",id="tvparnr",var.link=1,
##'       theta.formula=desfs,desnames=c("mz","dz","os"))
##' summary(out3)
##' 
##' ### use of clayton oakes binomial additive gamma model 
##' ###########################################################
##' \donttest{ ## Reduce Ex.Timings
##' data <- simbinClaytonOakes.family.ace(10000,2,1,beta=NULL,alpha=NULL)  
##' margbin <- glm(ybin~x,data=data,family=binomial())
##' margbin
##' 
##' head(data)
##' data$number <- c(1,2,3,4)
##' data$child <- 1*(data$number==3)
##' 
##' ### make ace random effects design
##' out <- ace.family.design(data,member="type",id="cluster")
##' out$pardes
##' head(out$des.rv)
##' 
##' bints <- binomial.twostage(margbin,data=data,
##'      clusters=data$cluster,detail=0,var.par=0,
##'      theta=c(2,1)/9,var.link=0,
##'      random.design=out$des.rv,theta.des=out$pardes)
##' summary(bints)
##' 
##' data <- simbinClaytonOakes.twin.ace(10000,2,1,beta=NULL,alpha=NULL)
##' out  <- twin.polygen.design(data,id="cluster",zygname="zygosity")
##' out$pardes
##' head(out$des.rv)
##' margbin <- glm(ybin~x,data=data,family=binomial())
##' 
##' bintwin <- binomial.twostage(margbin,data=data,
##'      clusters=data$cluster,detail=1,var.par=0,
##'      theta=c(2,1),
##'      random.design=out$des.rv,theta.des=out$pardes)
##' summary(bintwin)
##' concordance.twin.ace(bintwin)
##' }
##' 
##' @keywords binomial regression 
##' @author Thomas Scheike
##' @export
##' @param margbin Marginal binomial model 
##' @param data data frame
##' @param score.method Scoring method default is "fisher.scoring" among "fisher.scoring","nlminb","optimize","nlm"
##' @param Nit Number of iterations
##' @param detail Detail
##' @param clusters Cluster variable
##' @param silent Debug information
##' @param weights Weights for log-likelihood, can be used for each type of outcome in 2x2 tables. 
##' @param control Optimization arguments
##' @param theta Starting values for variance components
##' @param theta.des design for dependence parameters, when pairs are given this is could be a (pairs) x (numer of parameters)  x (max number random effects) matrix
##' @param var.link Link function for variance 
##' @param var.par parametrization 
##' @param var.func when alternative parametrizations are used this function can specify how the paramters are related to the \eqn{\lambda_j}'s.
##' @param iid Calculate i.i.d. decomposition when iid>=1, when iid=2 then avoids adding the uncertainty for marginal paramters for additive gamma model (default).
##' @param step Step size
##' @param notaylor Taylor expansion
##' @param model model
##' @param marginal.p vector of marginal probabilities 
##' @param beta.iid iid decomposition of marginal probability  estimates for each subject, if based on GLM model this is computed. 
##' @param Dbeta.iid derivatives of marginal model wrt marginal parameters, if based on GLM model this is computed. 
##' @param strata strata for fitting: considers only pairs where both are from same strata 
##' @param max.clust max clusters
##' @param se.clusters clusters for iid decomposition for roubst standard errors
##' @param numDeriv uses Fisher scoring aprox of second derivative if 0, otherwise numerical derivatives 
##' @param random.design random effect design for additive gamma model, when pairs are given this is a (pairs) x (2) x (max number random effects) matrix, see pairs.rvs below
##' @param pairs matrix with rows of indeces (two-columns) for the pairs considered in the pairwise composite score, useful for case-control sampling when marginal is known.
##' @param pairs.rvs for additive gamma model and random.design and theta.des are given as arrays, this specifice number of random effects for each pair. 
##' @param additive.gamma.sum this is specification of the lamtot in the models via a matrix that is multiplied onto the parameters theta (dimensions=(number random effects x number of theta parameters), when null then sums all parameters. Default is a matrix of 1's 
##' @param pair.ascertained if pairs are sampled only when there are events in the pair i.e. Y1+Y2>=1. 
##' @param case.control if data is case control data for pair call, and here 2nd column of pairs are probands (cases or controls)
##' @param twostage default twostage=1, to fit MLE use twostage=0
##' @param beta is starting value for beta for MLE version 
binomial.twostage <- function(margbin,data=sys.parent(),
     score.method="fisher.scoring",Nit=60,detail=0,clusters=NULL,silent=1,weights=NULL,
     control=list(),theta=NULL,theta.des=NULL,var.link=0,var.par=1,var.func=NULL,
     iid=1,step=1.0,notaylor=1,model="plackett",marginal.p=NULL,beta.iid=NULL,Dbeta.iid=NULL,
     strata=NULL,max.clust=NULL,se.clusters=NULL,numDeriv=0,
     random.design=NULL,pairs=NULL,pairs.rvs=NULL,additive.gamma.sum=NULL,
     pair.ascertained=0,case.control=0,
     twostage=1,beta=NULL) 
{ ## {{{
    ## {{{ seting up design and variables
    rate.sim <- 1; sym=1; 
    if (model=="clayton.oakes") dep.model <- 1 else if (model=="plackett") dep.model <- 2 else stop("Model must by either clayton.oakes or plackett \n"); 
    antpers <- NROW(data); 

### marginal prediction and binomial response, two types of calls ## {{{
    if (class(margbin)[1]=="glm") {
        ps <- predict(margbin,newdata=data,type="response")
        ### takes data to extract response and predictions, these could be different for pairs call
###     cause <- margbin$y
###     print(all.vars(margbin$formula)[1])
        cause <- data[,all.vars(margbin$formula)[1]]
	if (!is.numeric(cause)) stop(paste("response in data",margbin$formula)[1],"not numeric\n"); 
	if (is.null(beta.iid))   beta.iid <- iid(margbin,id=clusters)
	if (is.null(Dbeta.iid)) Dbeta.iid <- model.matrix(margbin$formula,data=data) * ps
	if (twostage==0)            Xbeta <- model.matrix(margbin$formula,data=data)
    }
    else if (class(margbin)[1]=="formula") {
        margbin <- glm(margbin,data=data,family=binomial())
        ps <- predict(margbin,type="response")
        cause <- margbin$y
	if (twostage==0)            Xbeta <- model.matrix(margbin$formula,data=data)
	if (is.null(Dbeta.iid)) Dbeta.iid <- model.matrix(margbin$formula,data=data) * ps
	if (is.null(beta.iid))   beta.iid <- iid(margbin,id=clusters)
    }  else if (is.null(marginal.p))
        stop("without marginal model, marginal p's must be given\n"); 

    if (!is.null(marginal.p)) {
        if (length(margbin)!=antpers) 
            stop("with marginal margbin is response \n")
        else cause <- margbin
        if (length(marginal.p)!=antpers) 
            stop("length same as data dimension  \n")
        else ps <- marginal.p
    }
    ## }}}

    notaylor <- 1
    if (is.null(weights)==TRUE) weights <- rep(1,antpers); 
    if (is.null(strata)==TRUE) strata<- rep(1,antpers); 
    if (length(strata)!=antpers) stop("Strata must have length equal to number of data points \n"); 

    out.clust <- cluster.index(clusters);  
    clusters <- out.clust$clusters
    maxclust <- out.clust$maxclust 
    antclust <- out.clust$cluster.size
    clusterindex <- out.clust$idclust
    clustsize <- out.clust$cluster.size
    call.secluster <- se.clusters

    if (is.null(se.clusters)) { se.clusters <- clusters; antiid <- nrow(clusterindex);} else  {
                                    iids <-  unique(se.clusters); 
                                    antiid <- length(iids); 
                                    if (is.numeric(se.clusters)) se.clusters <-  fast.approx(iids,se.clusters)-1
                                    else se.clusters <- as.integer(factor(se.clusters, labels = seq(antiid)))-1
                                }
    if (length(se.clusters)!=length(clusters)) stop("Length of seclusters and clusters must be same\n"); 

    if ((!is.null(max.clust))) if (max.clust< antiid) {
        coarse.clust <- TRUE
	qq <- unique(quantile(se.clusters, probs = seq(0, 1, by = 1/max.clust)))
	qqc <- cut(se.clusters, breaks = qq, include.lowest = TRUE)    
	se.clusters <- as.integer(qqc)-1
	max.clusters <- length(unique(se.clusters))
	maxclust <- max.clust    
	antiid <- max.clusters
    }                                                         

    ratesim<-rate.sim; 

   if (!is.null(random.design)) { ### different parameters for Additive random effects 
     dep.model <- 3
###     if (is.null(random.design)) random.design <- matrix(1,antpers,1); 
     dim.rv <- ncol(random.design); 
     if (is.null(theta.des)) theta.des<-diag(dim.rv);
###     ptheta <- dimpar <- ncol(theta.des); 
 
###   if (dim(theta.des)[2]!=ncol(random.design)) 
###   stop("nrow(theta.des)!= ncol(random.design),\nspecifies restrictions on paramters, if theta.des not given =diag (free)\n"); 
 } else { random.design <- matrix(0,1,1);  dim.rv <- 1; }


    if (is.null(theta.des)==TRUE) ptheta<-1; 
    if (is.null(theta.des)==TRUE) theta.des<-matrix(1,antpers,ptheta) 
###    else theta.des<-as.matrix(theta.des); 
###    ptheta<-ncol(theta.des); 
###    if (nrow(theta.des)!=antpers) stop("Theta design does not have correct dim");
  if (length(dim(theta.des))==3) ptheta<-dim(theta.des)[2] else if (length(dim(theta.des))==2) ptheta<-ncol(theta.des)
  if (nrow(theta.des)!=antpers & dep.model!=3 ) stop("Theta design does not have correct dim");

   if (length(dim(theta.des))!=3) theta.des <- as.matrix(theta.des)
###   theta.des <- as.matrix(theta.des)

    if (is.null(beta)==TRUE & twostage==0) beta <- coef(margbin) else beta <- 0
    dimbeta <- length(beta); 
    if (is.null(theta)==TRUE) {
        if (var.link==1) theta<- rep(-0.7,ptheta);  
        if (var.link==0) theta<- rep(exp(-0.7),ptheta);   
    }       
    if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta); 
    theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 

    if (maxclust==1) stop("No clusters, maxclust size=1\n"); 

  antpairs <- 1; ### to define 

  if (is.null(additive.gamma.sum)) additive.gamma.sum <- matrix(1,dim.rv,ptheta)

  if (!is.null(pairs)) { pair.structure <- 1; } else  pair.structure <- 0;  
  if (pair.structure==1 & dep.model==3) { ## {{{ 
### something with dimensions of rv.des 
### theta.des
       antpairs <- nrow(pairs); 
       if ( (length(dim(theta.des))!=3)  & (length(dim(random.design))==3) )
       {
          Ptheta.des <- array(0,c(nrow(theta.des),ncol(theta.des),antpairs))
          for (i in 1:antpairs) Ptheta.des[,,i] <- theta.des
       theta.des <- Ptheta.des
       }
       if ( (length(dim(theta.des))==3)  & (length(dim(random.design))!=3) )
       {
           rv.des <- array(0,c(2,ncol(random.design),antpairs))
           for (i in 1:antpairs) {
		   rv.des[1,,i] <- random.design[pairs[i,1],]
		   rv.des[2,,i] <- random.design[pairs[i,2],]
	   }
       random.design <- rv.des
       }
       if ( (length(dim(theta.des))!=3)  & (length(dim(random.design))!=3) )
       {
          Ptheta.des <- array(0,c(nrow(theta.des),ncol(theta.des),antpairs))
          rv.des <- array(0,c(2,ncol(random.design),antpairs))
          for (i in 1:antpairs) {
		   rv.des[1,,i] <- random.design[pairs[i,1],]
		   rv.des[2,,i] <- random.design[pairs[i,2],]
                   Ptheta.des[,,i] <- theta.des
	   }
       theta.des <- Ptheta.des
       random.design <- rv.des
       }

       if (max(pairs)>antpers) stop("Indices of pairs should refer to given data \n"); 
       if (is.null(pairs.rvs)) pairs.rvs <- rep(dim(random.design)[2],antpairs)
###       if (max(pairs.rvs)> dim(random.design)[3] | max(pairs.rvs)>ncol(theta.des[1,,])) 
###	       stop("random variables for each cluster higher than  possible, pair.rvs not consistent with random.design or theta.des\n"); 
       clusterindex <- pairs-1; 
  } ## }}} 

    if (pair.structure==1 & dep.model!=3) {
        clusterindex <- pairs-1; 
        antpairs <- nrow(pairs); 
        pairs.rvs <- 1
     }

    if (pair.structure==1) {
	    if (length(case.control)!=antpairs)         case.control <- rep(case.control[1],antpairs)
	    if (length(pair.ascertained)!=antpairs) pair.ascertained <- rep(pair.ascertained[1],antpairs)
	    if (any(case.control+pair.ascertained==2))  stop("Each pair is either case.control pair or pair.ascertained \n"); 
    }

    if (is.null(Dbeta.iid)) Dbeta.iid <- matrix(0,length(cause),1); 
    ptrunc <- rep(1,antpers); 

    ## }}}

    if (dep.model==3) model <- "clayton.oakes"

###  print(dep.model)
###  print(summary(case.control))
###  print(length(case.control))
###  print(summary(pair.ascertained))
###  print(length(pair.ascertained))
###  print(dim(random.design)); print(dim(theta.des)); 
###  print(dim(additive.gamma.sum))
###  print(head(pairs.rvs))
###  print(antpairs)

    loglike <- function(par) 
    { ## {{{
		
      if (pair.structure==0 | dep.model!=3) Xtheta <- as.matrix(theta.des) %*% matrix(c(par[seq(1,ptheta)]),nrow=ptheta,ncol=1);
      if (pair.structure==1 & dep.model==3) Xtheta <- matrix(0,antpers,1); ## not needed 
      DXtheta <- array(0,c(1,1,1));

      if (twostage==0) epar <- par[seq(1,ptheta)] else epar <- par
      if (twostage==0) { ### update, marginal.p og score for logistic model
###	     print(c(ptheta+1,ptheta+dimbeta))
	    beta <- par[seq(ptheta+1,ptheta+dimbeta)]
	    lp <- c(Xbeta %*%  beta)
###	    print(summary(ps))
            psu <- exp(lp)/(1+exp(lp)) 
###	    print(summary(psu))
            dpsu <- psu/(1+exp(lp)) 
	    ### update predictions and DbetaP
###	    print(beta); print(dim(Xbeta)); print(length(dpsu)); 
	    Dbeta.iid <- Xbeta * dpsu
###	    print(Dbeta.iid)
###	    print(apply(Dbeta.iid,2,sum))
	    ps <- psu
      } 

      if (var.link==1 & dep.model==3) epar <- c(exp(epar)) else epar <- c(epar)
      partheta <- epar

      if (var.par==1 & dep.model==3) {
       ## from variances to 
       if (is.null(var.func)) {
	    sp <- sum(epar)
	    partheta <- epar/sp^2 
         } else partheta <- epar; ## par.func(epar)
      } 

      if (pair.structure==0) 
            outl<-.Call("twostageloglikebin", ## {{{
            icause=cause,ipmargsurv=ps, 
            itheta=c(partheta),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),ithetades=theta.des,
            icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
            ivarlink=var.link,iiid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
            itrunkp=ptrunc,istrata=strata,iseclusters=se.clusters,iantiid=antiid, 
            irvdes=random.design,iags=additive.gamma.sum,ibetaiid=Dbeta.iid,pa=pair.ascertained,twostage=twostage)
            ## }}}
      else outl<-.Call("twostageloglikebinpairs", ## {{{
            icause=cause,ipmargsurv=ps, 
            itheta=c(partheta),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),ithetades=theta.des,
            icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
            ivarlink=var.link,iiid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
            itrunkp=ptrunc,istrata=strata,iseclusters=se.clusters,iantiid=antiid, 
            irvdes=random.design,
            idimthetades=dim(theta.des),idimrvdes=dim(random.design),irvs=pairs.rvs,
            iags=additive.gamma.sum,ibetaiid=Dbeta.iid,pa=pair.ascertained,twostage=twostage,
	    icasecontrol=case.control)
             ## }}} 

      if (detail==3) print(c(par,outl$loglike))

         ## variance parametrization, and inverse.link 
	 if (dep.model==3) {# {{{
	    if (var.par==1) {
		 ## from variances to and with sum for all random effects 
		 if (is.null(var.func)) {
		 if (var.link==0)  {
	###		 print(c(sp,epar))
		     mm <- matrix(-epar*2*sp,length(epar),length(epar))
		     diag(mm) <- sp^2-epar*2*sp
		 } else {
		    mm <- -c(epar) %o% c(epar)*2*sp
		    diag(mm) <- epar*sp^2-epar^2*2*sp
		 }
		    mm <- mm/sp^4
		 } else mm  <- numDeriv::hessian(var.func,par)
	      } else {
	          if (var.link==0) mm <- diag(length(epar)) else mm <- diag(c(epar))
	      }
              if (twostage==0) {  ### beta for logistic regression also part of model
		      mm0 <- diag(length(par))
		      mm0[1:nrow(mm),1:nrow(mm)] <- mm
		      mm <- mm0
	      }

	      }# }}}


	    if (dep.model==3) {# {{{
	       outl$score <-  t(mm) %*% outl$score
	       if (twostage==1) outl$DbetaDtheta <-  t(mm) %*% outl$DbetaDtheta
	       outl$Dscore <- t(mm) %*% outl$Dscore %*% mm
               if (iid==1) outl$theta.iid <- t(t(mm) %*% t(outl$theta.iid))

        ###       print(crossprod(outl$theta.iid)); print(outl$Dscore)
	###       print(c(outl$score))
	###       print(apply(outl$theta.iid,2,sum))
	    }# }}}

            attr(outl,"grad") <- attr(outl,"gradient") <-outl$score 
	    attr(outl,"hessian") <- outl$Dscore
            if (oout==0) ret <- c(-1*outl$loglike) else if (oout==1) ret <- sum(outl$score^2) else if (oout==3) ret <- outl$score else ret <- outl
            return(ret)
        } ## }}}

    if (score.method=="optimize" && ptheta!=1) {cat("optimize only works for d==1, score.mehod set to nlminb \n"); score.method <- "nlminb";}

    theta.iid <- NULL
    logl <- NULL
    p <- theta
    if (twostage==0) p <- c(p,beta); 
    theta <- p
    if (score.method=="fisher.scoring") { ## {{{
        oout <- 2;  ### output control for obj

        if (Nit>0) 
            for (i in 1:Nit)
            {
                    out <- loglike(p)
                    hess <-  -1* out$Dscore
                    if (!is.na(sum(hess))) hessi <- lava::Inverse(out$Dscore) else hessi <- hess 
                    if (detail==1) {## {{{
                        cat(paste("Fisher-Scoring ===================: it=",i,"\n")); 
                        cat("theta:");print(c(p))
                        cat("loglike:");cat(c(out$loglike),"\n"); 
                        cat("score:");cat(c(out$score),"\n"); 
                        cat("hess:\n"); cat(out$Dscore,"\n"); 
                    }## }}}
                    delta <- hessi %*% out$score *step 
	            if (Nit>0) {
                    p <- p + delta* step
                    theta <- p; 
		    }
                    if (is.nan(sum(out$score))) break; 
                    if (sum(abs(out$score))<0.00001) break; 
                    if (max(abs(theta))>20 & var.link==0) { cat("theta too large lacking convergence \n"); break; }
                }
        if (!is.nan(sum(p))) { 
            if (detail==1 && iid==1) cat("iid decomposition\n"); 
            out <- loglike(p) 
            logl <- out$loglike
            score1 <- score <- out$score
            oout <- 0; 
            hess1 <- hess <- -1*out$Dscore 
###            if (iid==1) theta.iid <- out$theta.iid
            if (detail==1 && iid==1) cat("finished iid decomposition\n"); 
        }
        if (numDeriv==1) {
            if (detail==1 ) cat("starting numDeriv for second derivative \n"); 
            oout <- 0; 
            score2 <- numDeriv::jacobian(loglike,p)
	    if (detail==1) {
		    cat("Derivative\n"); 
	            cat(c(out$score))
		    cat("\n Numerical Derivative\n"); 
	            cat(c(score2))
	            cat("\n")
	    }
	    score1 <- matrix(score2,ncol=1)
            oout <- 3
            hess <- numDeriv::jacobian(loglike,p)
            if (detail==1 ){
                  cat("finished numDeriv for second derivative \n"); 
	          cat("Numerical derivative second derivative \n"); 
	          print(hess)
	          cat("average second moment, for fisher-scoring\n"); 
	          print(out$Dscore)
	    }
        }
        if (detail==1 & Nit==0) {## {{{
            cat(paste("Fisher-Scoring ===================: final","\n")); 
            cat("theta:");print(c(p))
            cat("loglike:");cat(c(out$loglike),"\n"); 
            cat("score:");cat(c(out$score),"\n"); 
            cat("hess:\n"); cat(out$Dscore,"\n"); 
        }## }}}
        if (!is.na(sum(hess))) hessi <- lava::Inverse(hess) else hessi <- diag(nrow(hess))
        ## }}}
    } else if (score.method=="nlminb") { ## {{{ nlminb optimizer
        oout <- 0; 
        tryCatch(opt <- nlminb(theta,loglike,control=control),error=function(x) NA)
        if (detail==1) print(opt); 

        if (detail==1 && iid==1) cat("iid decomposition\n"); 
        oout <- 2
        theta <- opt$par
        out <- loglike(opt$par)
        logl <- out$loglike
        score1 <- score <- out$score
        hess1 <- hess <- - out$Dscore
        if (iid==1) theta.iid <- out$theta.iid
        if (numDeriv==1) {
            oout <- 3; 
            p <- theta
            hess <- -1 * numDeriv::jacobian(loglike,theta)
        }
        hessi <- lava::Inverse(hess); 
        ## }}}
    } else if (score.method=="optimize" && ptheta==1) { ## {{{  optimizer
        oout <- 0; 
        if (var.link==1) {mino <- -20; maxo <- 10;} else {mino <- 0.001; maxo <- 100;}
        tryCatch(opt <- optimize(loglike,c(mino,maxo)));
        if (detail==1) print(opt); 

        opt$par <- opt$minimum
        theta <- opt$par
        if (detail==1 && iid==1) cat("iid decomposition\n"); 
        oout <- 2
        out <- loglike(opt$par)
        logl <- out$loglike
        score1 <- score <- out$score
        hess1 <- hess <- - out$Dscore
        if (iid==1) theta.iid <- out$theta.iid
        if (numDeriv==1) {
            oout <- 3; 
            p <- opt$par
            hess <-  -1* numDeriv::jacobian(loglike,p)
        }
        hessi <- lava::Inverse(hess); 
        ## }}}
    } else if (score.method=="nlm") { ## {{{ nlm optimizer
        iid <- 0; oout <- 0; 
        tryCatch(opt <- nlm(loglike,theta,hessian=TRUE,print.level=detail),error=function(x) NA)
        iid <- 1; 
        hess <-  opt$hessian
        score <- opt$gradient
        if (detail==1) print(opt); 
        hessi <-  lava::Inverse(hess); 
        theta <- opt$estimate
        if (detail==1 && iid==1) cat("iid decomposition\n"); 
        oout <- 2
        out <- loglike(opt$estimate)
        logl <- out$loglike
        score1 <- out$score
        hess1 <- -1* out$Dscore
        if (iid==1) theta.iid <- out$theta.iid
        ## }}}
    }  else stop("score.methods = optimize(dim=1) nlm nlminb fisher.scoring\n"); 


    ## {{{ handling output
    iid.tot <- NULL
    var.tot <- robvar.theta <- NULL
    beta <- NULL
    if (iid>=1) {
        theta.iid <- out$theta.iid %*% hessi
        if (dep.model==3 & iid!=2  & (!is.null(beta.iid)))
	if (nrow(beta.iid)==nrow(out$theta.iid) & twostage==1) {
	     theta.beta.iid <- (beta.iid %*% t(out$DbetaDtheta) ) %*% hessi
	     theta.iid  <- theta.iid+theta.beta.iid
	     iid.tot <- cbind(theta.iid,beta.iid)
	     var.tot <- crossprod(iid.tot)
	} 
        var.theta <- robvar.theta  <- (t(theta.iid) %*% theta.iid) 
	if (is.null(var.tot)) var.tot <- var.theta
    } else var.theta <- -1* hessi 

  if (class(margbin)[1]=="glm") beta <- coef(margbin); 
  if (twostage==0) beta <- theta[seq(ptheta,ptheta+dimbeta)]

  if (iid==1) var.theta <- robvar.theta else var.theta <- -hessi
  if (!is.null(colnames(theta.des))) thetanames <- colnames(theta.des) else thetanames <- paste("dependence",1:length(theta),sep="")
### fix names !!!
###    theta <- matrix(theta,length(theta),1)
###    if (length(thetanames)==nrow(theta)) { rownames(theta) <- thetanames; rownames(var.theta) <- colnames(var.theta) <- thetanames; }

    ud <- list(theta=theta,score=score,hess=hess,hessi=hessi,var.theta=var.theta,model=model,robvar.theta=robvar.theta,
               theta.iid=theta.iid,thetanames=thetanames,
	       loglike=-logl,score1=score1,Dscore=out$Dscore,
	       margsurv=ps,iid.tot=iid.tot,var.tot=var.tot,beta=beta); 
    class(ud)<-"mets.twostage" 
    attr(ud, "binomial") <- TRUE
    attr(ud, "ptheta") <- ptheta
    attr(ud, "Formula") <- formula
    attr(ud, "Clusters") <- clusters
    attr(ud,"sym")<-sym; 
    attr(ud,"var.link")<-var.link; 
    attr(ud,"var.par")<-var.par; 
    attr(ud,"var.func")<-var.func; 
    attr(ud,"antpers")<-antpers; 
    attr(ud,"antclust")<-antclust; 
    attr(ud, "Type") <- model
    attr(ud,"DbetaDtheta")<-out$DbetaDtheta; 
    attr(ud,"ags")<- additive.gamma.sum 
    attr(ud,"twostage")<- twostage
    attr(ud,"pair.ascertained")<- pair.ascertained 
    ### to be consistent with structure for survival twostage model 
    attr(ud, "additive-gamma") <- (dep.model==3)*1
    if (dep.model==3 & pair.structure==1) attr(ud, "likepairs") <- c(out$likepairs)
    if (dep.model==3 & pair.structure==0) attr(ud, "pardes") <- theta.des
    if (dep.model==3 & pair.structure==0) attr(ud, "theta.des") <- theta.des
    if (dep.model==3 & pair.structure==1) attr(ud, "pardes") <- theta.des[,,1]
    if (dep.model==3 & pair.structure==1) attr(ud, "theta.des") <- theta.des[,,1]
    if (dep.model==3 & pair.structure==0) attr(ud, "rv1") <- random.design[1,]
    if (dep.model==3 & pair.structure==1) attr(ud, "rv1") <- random.design[1,,1]
 
    attr(ud, "response") <- "binomial"
    return(ud);
    ## }}}

} ## }}}

##' @export
p11.binomial.twostage.RV <- function(theta,rv1,rv2,p1,p2,pardes,ags=NULL,link=0,i=1,j=1) { ## {{{ 
	## computes p11 pij for additive gamma binary random effects model 
     if (is.null(ags)) ags <- matrix(1,dim(pardes))
     out <- .Call("claytonoakesbinRV",theta,i,j,p1,p2,rv1,rv2,pardes,ags,link)$like
 return(out)
} ## }}} 

##' @export
concordance.twostage<- function(theta,p,rv1,rv2,theta.des,additive.gamma.sum=NULL,link=0,var.par=0)
{# {{{

   ### takes dependence paramter from output
      ptheta <- length(theta)
###      ptheta <- attr(object,"ptheta")
###   theta <- object$theta[seq(1,ptheta)]
###   robvar.theta <- object$robvar.theta[seq(1,ptheta),seq(1,ptheta)]

   if (var.par==1) theta <- theta/sum(theta)^2

   nn <- nrow(as.matrix(rv1))
   if (is.matrix(p)==FALSE) { ll <- length(p); p <- matrix(p,ll,2); }
   if (ncol(p)!=2) p <- matrix(p,ncol=2)
   if (is.matrix(rv1)==FALSE) rv1 <- matrix(rv1,nn,length(rv1),byrow=TRUE)
   if (is.matrix(rv2)==FALSE) rv2 <- matrix(rv2,nn,length(rv1),byrow=TRUE)

   if (is.null(additive.gamma.sum)) ags <- matrix(1,ncol(rv1),ptheta) else ags <- additive.gamma.sum


   tabs <- list()
   for (i in 1:nn)
   {
        p1 <- p[i,1]
        p2 <- p[i,2]
	p11 <- p11.binomial.twostage.RV(theta,rv1[i,],rv2[i,],p1,p2,theta.des,ags=ags,link=link)
	p01 <- p[i,1]-p11
	p10 <- p[i,2]-p11
	p00 <- 1-p01-p10+p11
	tabs[[i]] <- list(pmat=matrix(c(p00,p10,p01,p11),2,2),casewise=c(p11/p1,p11/p2),marg=c(p1,p2))
   }

   return(tabs)
}	# }}}

##' @export
concordance.twin.ace<- function(object,rv1=NULL,rv2=NULL,xmarg=NULL,type="ace")
{# {{{

if (type=="ace" | type=="ade") {
if (is.null(rv1))  rv1 <- rbind(c(1,0,0,0,1),c(0,1,1,0,1))
if (is.null(rv2))  rv2 <- rbind(c(1,0,0,0,1),c(0,1,0,1,1))
}
if (type=="ae" | type=="de") {
if (is.null(rv1))  rv1 <- rbind(c(1,0,0,0),c(0,1,1,0))
if (is.null(rv2))  rv2 <- rbind(c(1,0,0,0),c(0,1,0,1))
}

var.par <- attr(object,"var.par")
var.link <- attr(object,"var.link")
theta.des <- attr(object,"theta.des")
twostage <- attr(object,"twostage")

if (twostage==0) {
	pt <- attr(object,"ptheta")
	theta <- object$theta[seq(1,pt)]
        beta <- object$theta[seq(pt+1,length(object$theta))] 
	var.theta <- object$robvar.theta[seq(1,pt),seq(1,pt)]
	var.tot <- object$var.theta
} else { 
	theta <- object$theta; 
	beta <- object$beta
	var.theta <- object$var.theta; 
	var.tot <- object$var.tot
}

if (var.link==1) theta <- exp(theta)
###print(beta)
ags <- attr(object,"ags"); 
if (is.null(xmarg)) {xmarg <- rep(0,length(beta)); xmarg[1] <- 1;} 

nn <- nrow(rv1)
if (is.matrix(rv1)==FALSE) rv1 <- matrix(rv1,nn,length(rv1),byrow=TRUE)
if (is.matrix(rv2)==FALSE) rv2 <- matrix(rv2,nn,length(rv1),byrow=TRUE)
if (is.null(ags)) ags <- matrix(1,ncol(rv1),length(theta));

   fp <- function(par) {# {{{
       pp <- par[1:length(theta)];
       beta <- par[(length(theta)+1):length(par)];
       xp  <- sum(xmarg*beta); 
       pm <- exp(xp)/(1+exp(xp)); 
       if (var.par==1) pp <- pp/sum(pp)^2
       p11 <- p11.binomial.twostage.RV(pp,rv1l,rv2l,pm,pm,theta.des,ags=ags,link=0)
       casewise <- p11/pm
       return(c(p11,casewise,pm))
   }# }}}

print(c(theta,beta)); print(var.tot)

   tabs <- list()
   for (i in 1:nn)
   {
      rv1l <- rv1[i,]
      rv2l <- rv2[i,]
      tabs[[i]] <- lava::estimate(coef=c(theta,beta),vcov=var.tot,f=function(p) fp(p))
   }

   return(tabs)
}	# }}}

##' @export
binomial.twostage.time <- function(formula,data,id,...,silent=1,fix.censweights=1,
breaks=Inf,pairsonly=TRUE,fix.marg=NULL,cens.formula,cens.model="aalen",weights="w") {
   ## {{{ 
    m <- match.call(expand.dots = FALSE)
    m <- m[match(c("","data"),names(m),nomatch = 0)]
    Terms <- terms(cens.formula,data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    M <- eval(m,envir=parent.frame())

    censtime <- model.extract(M, "response")
    status <- censtime[,2]
    time <- censtime[,1]
    timevar <- colnames(censtime)[1]
    outcome <- as.character(terms(formula)[[2]])    
    if (is.null(breaks)) breaks <-  quantile(time,c(0.25,0.5,0.75,1))

    outcome0 <- paste(outcome,"_dummy")
    res <- list()
    logor <- cif <- conc <- c()
    k <- 0
    for (tau in rev(breaks)) {
        if ((length(breaks)>1) & (silent==0)) message(tau)
        ### construct min(T_i,tau) or T_i and related censoring variable, 
        ### thus G_c(min(T_i,tau)) or G_c(T_i) as weights
        if ((fix.censweights==1 & k==0) | (fix.censweights==0)) {
        data0 <- data
        time0 <- time 
        status0 <- status 
        }
        cond0 <- (time>tau)
        if ((fix.censweights==1 & k==0) | (fix.censweights==0)) status0[cond0 & status==1] <- 3 
        if (fix.censweights==0) status0[cond0 & status==1] <- 3 
        data0[,outcome] <- data[,outcome]
        data0[cond0,outcome] <- FALSE
        if ((fix.censweights==1 & k==0) | (fix.censweights==0)) time0[cond0] <- tau
###        if (fix.censweights==0 ) time0[cond0] <- tau
        if ((fix.censweights==1 & k==0) | (fix.censweights==0)) {
		data0$S <- survival::Surv(time0,status0==1)        
	}
	if ((fix.censweights==1 & k==0) | (fix.censweights==0))
        dataw <- ipw(update(cens.formula,S~.),data=data0,cens.model=cens.model,
		     obsonly=TRUE)
        if ((fix.censweights==1)) 
		dataw[,outcome] <- (dataw[,outcome])*(dataw[,timevar]<tau)
	marg.bin <- glm(formula,data=dataw,weights=1/dataw[,weights],family="quasibinomial")
        pudz <- predict(marg.bin,newdata=dataw,type="response")
	dataw$pudz <- pudz
	datawdob <- fast.reshape(dataw,id=id)
        datawdob$minw <- pmin(datawdob$w1,datawdob$w2)
        dataw2 <- fast.reshape(datawdob)
	### removes second row of singletons 
	dataw2  <- subset(dataw2,!is.na(dataw2$minw)) 
	k <- k+1
	if (!is.null(fix.marg)) dataw2$pudz <- fix.marg[k]
        suppressWarnings( b <- binomial.twostage(dataw2[,outcome],data=dataw2,clusters=dataw2[,id],marginal.p=dataw2$pudz,weights=1/dataw2$minw,...))
        theta0 <- b$theta[1,1]
        prev <- prev0 <- exp(coef(marg.bin)[1])/(1+exp(coef(marg.bin)[1]))
	if (!is.null(fix.marg)) prev <- fix.marg[k]
        concordance <- plack.cif2(prev,prev,theta0)
	conc <- c(conc,concordance)
	cif <- c(cif,prev0)
	logor <- rbind(logor,coef(b))
###     res <- c(res,list(coef(b),concordance=concordance,cif=prev0))
    }
###    if (length(breaks)==1) return(b)
    res <- list(varname="Time",var=breaks,concordance=rev(conc),cif=rev(cif),
		time=breaks,call=m,type="time",logor=logor[k:1,])
###	coef=lapply(res,function(x) x$all),
###    class(res) <- ""
    return(res)    
} ## }}} 


##' Fits two-stage binomial for describing depdendence in binomial data
##' using marginals that are on logistic form using the binomial.twostage funcion, but
##' call is different and easier and the data manipulation is build into the function.
##' Useful in particular for family design data. 
##'
##' If clusters contain more than two times, the algoritm uses a compososite likelihood
##' based on the pairwise bivariate models.
##'
##' The reported standard errors are based on the estimated information from the 
##' likelihood assuming that the marginals are known. This gives correct standard errors
##' in the case of the plackett distribution (OR model for dependence), but incorrect for
##' the clayton-oakes types model. The OR model is often known as the ALR model. 
##' Our fitting procedures gives correct standard errors due to the ortogonality and is 
##' fast. 
##'
##' @examples
##' data(twinstut)
##' twinstut0 <- subset(twinstut, tvparnr<2300000)
##' twinstut <- twinstut0
##' twinstut$binstut <- (twinstut$stutter=="yes")*1
##' theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
##' margbin <- glm(binstut~factor(sex)+age,data=twinstut,family=binomial())
##' bin <- binomial.twostage(margbin,data=twinstut,var.link=1,
##' 		         clusters=twinstut$tvparnr,theta.des=theta.des,detail=0,
##' 	                 score.method="fisher.scoring")
##' summary(bin)
##' lava::estimate(coef=bin$theta,vcov=bin$var.theta,f=function(p) exp(p))
##' 
##' twinstut$cage <- scale(twinstut$age)
##' theta.des <- model.matrix( ~-1+factor(zyg)+cage,data=twinstut)
##' bina <- binomial.twostage(margbin,data=twinstut,var.link=1,
##' 		         clusters=twinstut$tvparnr,theta.des=theta.des,detail=0)
##' summary(bina)
##' 
##' theta.des <- model.matrix( ~-1+factor(zyg)+factor(zyg)*cage,data=twinstut)
##' bina <- binomial.twostage(margbin,data=twinstut,var.link=1,
##' 		         clusters=twinstut$tvparnr,theta.des=theta.des)
##' summary(bina)
##' 
##' out <- easy.binomial.twostage(stutter~factor(sex)+age,data=twinstut,
##'                               response="binstut",id="tvparnr",var.link=1,
##' 			          theta.formula=~-1+factor(zyg1))
##' summary(out)
##' 
##' ## refers to zygosity of first subject in eash pair : zyg1
##' ## could also use zyg2 (since zyg2=zyg1 within twinpair's))
##' desfs <- function(x,num1="zyg1",namesdes=c("mz","dz","os"))
##'     c(x[num1]=="mz",x[num1]=="dz",x[num1]=="os")*1
##' 
##' out3 <- easy.binomial.twostage(binstut~factor(sex)+age,
##'                                data=twinstut, response="binstut",id="tvparnr",
##'                                var.link=1,theta.formula=desfs,
##'                                desnames=c("mz","dz","os"))
##' summary(out3)
##'
##' \donttest{ ## Reduce Ex.Timings
##' n <- 10000
##' set.seed(100)
##' dd <- simBinFam(n,beta=0.3) 
##' binfam <- fast.reshape(dd,varying=c("age","x","y"))
##' ## mother, father, children  (ordered)
##' head(binfam)
##'
##' ########### ########### ########### ########### ########### ###########
##' ####  simple analyses of binomial family data 
##' ########### ########### ########### ########### ########### ###########
##' desfs <- function(x,num1="num1",num2="num2")
##' {  
##'      pp <- 1*(((x[num1]=="m")*(x[num2]=="f"))|(x[num1]=="f")*(x[num2]=="m"))
##'      pc <- (x[num1]=="m" | x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1
##'      cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1
##'      c(pp,pc,cc)
##' } 
##'
##' ud <- easy.binomial.twostage(y~+1,data=binfam,
##'      response="y",id="id",
##'      theta.formula=desfs,desnames=c("pp","pc","cc"))
##' summary(ud)
##'
##' udx <- easy.binomial.twostage(y~+x,data=binfam,
##'      response="y",id="id",
##'      theta.formula=desfs,desnames=c("pp","pc","cc"))
##' summary(udx)
##'
##' ########### ########### ########### ########### ########### ###########
##' ####  now allowing parent child POR to be different for mother and father 
##' ########### ########### ########### ########### ########### ###########
##'
##' desfsi <- function(x,num1="num1",num2="num2")
##' { 
##'     pp <- (x[num1]=="m")*(x[num2]=="f")*1
##'     mc <- (x[num1]=="m")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     fc <- (x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     c(pp,mc,fc,cc)
##' }
##'
##' udi <- easy.binomial.twostage(y~+1,data=binfam,
##'      response="y",id="id",
##'      theta.formula=desfsi,desnames=c("pp","mother-child","father-child","cc"))
##' summary(udi)
##'
##' ##now looking to see if interactions with age or age influences marginal models 
##' ##converting factors to numeric to make all involved covariates numeric
##' ##to use desfai2 rather then desfai that works on binfam 
##'
##' nbinfam <- binfam
##' nbinfam$num <- as.numeric(binfam$num)
##' head(nbinfam)
##'
##' desfsai <- function(x,num1="num1",num2="num2")
##' { 
##'     pp <- (x[num1]=="m")*(x[num2]=="f")*1
##' ### av age for pp=1 i.e parent pairs
##'     agepp <- ((as.numeric(x["age1"])+as.numeric(x["age2"]))/2-30)*pp 
##'     mc <- (x[num1]=="m")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     fc <- (x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1
##'     agecc <- ((as.numeric(x["age1"])+as.numeric(x["age2"]))/2-12)*cc 
##'     c(pp,agepp,mc,fc,cc,agecc)
##' } 
##'
##' desfsai2 <- function(x,num1="num1",num2="num2")
##' { 
##'     pp <- (x[num1]==1)*(x[num2]==2)*1
##'     agepp <- (((x["age1"]+x["age2"]))/2-30)*pp ### av age for pp=1 i.e parent pairs
##'     mc <- (x[num1]==1)*(x[num2]==3 | x[num2]==4)*1
##'     fc <- (x[num1]==2)*(x[num2]==3 | x[num2]==4)*1
##'     cc <- (x[num1]==3)*(x[num2]==3 | x[num2]==4)*1
##'     agecc <- ((x["age1"]+x["age2"])/2-12)*cc ### av age for children 
##'     c(pp,agepp,mc,fc,cc,agecc)
##' } 
##'
##' udxai2 <- easy.binomial.twostage(y~+x+age,data=binfam,
##'      response="y",id="id",
##'      theta.formula=desfsai,
##'      desnames=c("pp","pp-age","mother-child","father-child","cc","cc-age"))
##' summary(udxai2)
##' }
##' @keywords binomial regression 
##' @export
##' @param margbin Marginal binomial model 
##' @param data data frame
##' @param response name of response variable in data frame
##' @param id name of cluster variable in data frame
##' @param score.method Scoring method
##' @param Nit Number of iterations
##' @param detail Detail for more output for iterations 
##' @param silent Debug information
##' @param weights Weights for log-likelihood, can be used for each type of outcome in 2x2 tables. 
##' @param control Optimization arguments
##' @param theta Starting values for variance components
##' @param theta.formula design for depedence, either formula or design function
##' @param desnames names for dependence parameters
##' @param deshelp if 1 then prints out some data sets that are used, on on which the design function operates
##' @param var.link Link function for variance 
##' @param iid Calculate i.i.d. decomposition
##' @param step Step size
##' @param model model
##' @param marginal.p vector of marginal probabilities 
##' @param strata strata for fitting 
##' @param max.clust max clusters
##' @param se.clusters clusters for iid decomposition for roubst standard errors
easy.binomial.twostage <- function(margbin=NULL,data=sys.parent(),score.method="fisher.scoring",
                                   response="response",id="id",
                                   Nit=60,detail=0, silent=1,weights=NULL,control=list(),
                                   theta=NULL,theta.formula=NULL,desnames=NULL,deshelp=0,var.link=1,iid=1,
                                   step=1.0,model="plackett",marginal.p=NULL,
				   strata=NULL,max.clust=NULL,se.clusters=NULL)
{ ## {{{
    if (class(margbin)[1]=="glm") ps <- predict(margbin,type="response") 
    else if (class(margbin)=="formula") {
        margbin <- glm(margbin,data=data,family=binomial())
        ps <- predict(margbin,type="response")
    }  else if (is.null(marginal.p)) stop("without marginal model, marginal p's must be given\n"); 

    if (!is.null(marginal.p)) ps <- marginal.p

    data <- cbind(data,ps)

### make all pairs in the families,
    fam <- familycluster.index(data[,id])
    data.fam <- data[c(fam$familypairindex),]
    data.fam$subfam <- c(fam$subfamilyindex)

### make dependency design using wide format for all pairs 
    data.fam.clust <- fast.reshape(data.fam,id="subfam")
    if (is.function(theta.formula)) {
        desfunction <- compiler::cmpfun(theta.formula)
	if (deshelp==1){
            cat("These names appear in wide version of pairs for dependence \n")
            cat("design function must be defined in terms of these: \n")
            cat(names(data.fam.clust)); cat("\n")
            cat("Here is head of wide version with pairs\n")
            print(head(data.fam.clust)); cat("\n")
	}
###	des.theta <- Reduce("rbind",lapply(seq(nrow(data.fam.clust)),function(i) unlist(desfunction(data.fam.clust[i,] ))))
        des.theta <- t(apply(data.fam.clust,1, function(x) desfunction(x)))
	colnames(des.theta) <- desnames
	desnames <- desnames
    } else {
        if (is.null(theta.formula)) theta.formula <- ~+1
        des.theta <- model.matrix(theta.formula,data=data.fam.clust)
        desnames <- colnames(des.theta); 
    }
    data.fam.clust <- cbind(data.fam.clust,des.theta)
    if (deshelp==1) {
        cat("These names appear in wide version of pairs for dependence \n")
        print(head(data.fam.clust))
    }

### back to long format keeping only needed variables
    data.fam <- fast.reshape(data.fam.clust,varying=c(response,id,"ps"))
    if (deshelp==1) {
	cat("Back to long format for binomial.twostage (head)\n"); 
        print(head(data.fam)); 
	cat("\n")
	cat(paste("binomial.twostage, called with reponse",response,"\n")); 
	cat(paste("cluster=",id,",  subcluster (pairs)=subfam \n")); 
	cat(paste("design variables =")); 
	cat(desnames)
	cat("\n  theta design (head) \n")
        print(head(data.fam[,desnames]));
	cat("\n  response (head) \n")
        print(head(data.fam[,response]));
	cat("\n")
    } 
    print(head(data.fam))

    out <- binomial.twostage(data.fam[,response],data=data.fam,
                             clusters=data.fam$subfam,
                             theta.des=data.fam[,desnames],
                             detail=detail, score.method=score.method, Nit=Nit,step=step,
                             iid=iid,theta=theta,var.link=var.link,model=model, 
                             max.clust=max.clust,
                             marginal.p=data.fam[,"ps"],se.clusters=data.fam[,id])

    return(out)
} ## }}}

##' @export
simBinPlack <- function(n,beta=0.3,theta=1,...) { ## {{{ 
    x1 <- rbinom(n,1,0.5)
    x2 <- rbinom(n,1,0.5)
###
    p1 <- exp(0.5+x1*beta)
    p2 <- exp(0.5+x2*beta)
    p1 <- p1/(1+p1)
    p2 <- p2/(1+p2)
###
    p11 <- plack.cif2(p1,p2,theta)
    p10 <- p1-p11
    p01 <- p2-p11
    p00 <- 1- p10-p01-p11
###
    y1 <- rbinom(n,1,p1)
    y2 <- (y1==1)*rbinom(n,1,p11/p1)+(y1==0)*rbinom(n,1,p01/(1-p1))
    list(x1=x1,x2=x2,y1=y1,y2=y2,id=1:n)
} ## }}} 
 
##' @export 
simBinFam <- function(n,beta=0.0,rhopp=0.1,rhomb=0.7,rhofb=0.1,rhobb=0.7) { ## {{{ 
    xc <- runif(n)*0.5
    xm <- rbinom(n,1,0.5+xc); 
    xf <- rbinom(n,1,0.5+xc); 
    xb1 <- rbinom(n,1,0.3+xc); 
    xb2 <- rbinom(n,1,0.3+xc); 
###
    rn <- matrix(rnorm(n*4),n,4)
    corm <- matrix( c(1,rhopp,rhomb,rhomb, rhopp,1,rhofb,rhofb, rhomb,rhofb,1,rhobb, rhomb,rhofb,rhobb,1),4,4)
    rnn <- t( corm %*% t(rn))
    zm <- exp(rnn[,1]); zf <- exp(rnn[,2]); zb1 <- exp(rnn[,3]); zb2 <- exp(rnn[,4]); 
    pm <- exp(0.5+xm*beta+zm)
    pf <- exp(0.5+xf*beta+zf)
    pf <- pf/(1+pf)
    pm <- pm/(1+pm)
    pb1 <- exp(0.5+xb1*beta+zb1)
    pb1 <- pb1/(1+pb1)
    pb2 <- exp(0.5+xb2*beta+zb2)
    pb2 <- pb2/(1+pb2)
    ym <- rbinom(n,1,pm)
    yf <- rbinom(n,1,pf)
    yb1 <- rbinom(n,1,pb1)
    yb2 <- rbinom(n,1,pb2)
                                        #
    agem <- 20+runif(n)*10
    ageb1 <- 5+runif(n)*10
    data.frame(agem=agem,agef=agem+3+rnorm(n)*2,
               ageb1=ageb1,ageb2=ageb1+1+runif(n)*3,xm=xm,xf=xf,xb1=xb1,xb2=xb2,ym=ym,yf=yf,yb1=yb1,yb2=yb2,id=1:n)
} ## }}} 

##' @export
simBinFam2 <- function(n,beta=0.0,alpha=0.5,lam1=1,lam2=1,...) { ## {{{ 
    x1 <- rbinom(n,1,0.5); x2 <- rbinom(n,1,0.5); 
    x3 <- rbinom(n,1,0.5); x4 <- rbinom(n,1,0.5); 
### random effects speicification of model via gamma distributions 
    zf <- rgamma(n,shape=lam1); zb <- rgamma(n,shape=lam2); 
    if (length(alpha)!=4) alpha <- rep(alpha[1],4)
    if (length(beta)!=4) beta <- rep(beta[1],4)
    pm <- exp(alpha[1]+x1*beta[1]+zf)
    pf <- exp(alpha[2]+x2*beta[2]+zf)
    pf <- pf/(1+pf)
    pm <- pm/(1+pm)
    pb1 <- exp(alpha[3]+x1*beta[3]+zf+zb)
    pb1 <- pb1/(1+pb1)
    ym <- rbinom(n,1,pm)
    yf <- rbinom(n,1,pf)
    yb1 <- rbinom(n,1,pb1)
    yb2 <- rbinom(n,1,pb1)
                                        #
    data.frame(x1=x1,x2=x2,ym=ym,yf=yf,yb1=yb1,yb2=yb2,id=1:n)
} ## }}} 

##' @export
simbinClaytonOakes.family.ace <- function(K,varg,varc,beta=NULL,alpha=NULL)  ## {{{ 
{
  ## K antal clustre (families), n=antal i clustre
###  K <- 10000
  n <- 4 # family members with ace structure
  ## total variance 1/(varg+varc)
  ## {{{  random effects 
  ###  means varg/(varg+varc) and variances varg/(varg+varc)^2
  eta <- varc+varg
  ### mother and father share environment
  ### children share half the genes with mother and father and environment 
  mother.g <-  cbind(rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta)
  father.g <-  cbind(rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta, rgamma(K,varg*0.25)/eta)
  env <- rgamma(K,varc)/eta 
  mother <- apply(mother.g,1,sum)+env
  father <- apply(father.g,1,sum)+env
  child1 <- apply(cbind(mother.g[,c(1,2)],father.g[,c(1,2)]),1,sum) + env
  child2 <- apply(cbind(mother.g[,c(1,3)],father.g[,c(1,3)]),1,sum) + env
  Gam1 <- cbind(mother,father,child1,child2)
  ## }}} 
  apply(Gam1,2,mean); apply(Gam1,2,var); cor(Gam1)

  ## {{{ marginals p's and conditional p's given random effects 
  ### marginals p's for mother, father, children
  xb1 <- rbinom(K,1,0.5)
  xb2 <- rbinom(K,1,0.5)
  xm <- rbinom(K,1,0.5)
  xf <- rbinom(K,1,0.5)
###
###
  if (is.null(beta)) beta <- rep(0.3,4)
  if (is.null(alpha)) alpha <- rep(0.5,4)
  pm <- exp(alpha[1]+xm*beta[1]); pf <- exp(alpha[2]+xf*beta[2])
  pb1 <- exp(alpha[3]+xb1*beta[3]); pb2 <- exp(alpha[4]+xb2*beta[4])
  p <- cbind(pm,pf,pb1,pb2)
  p <- p/(1+p)

  vartot <- eta
  pgivenZ <- ilap(vartot,p)
###  pgivenZ <- ilap(vartot,p)
  pgivenZ <- exp(- Gam1*pgivenZ)
  ## }}} 

  Ybin <- matrix(rbinom(n*K,1,c(pgivenZ)),K,n)
  Ybin <- t(Ybin)
  xs <- cbind(xm,xf,xb1,xb2)
  type <- rep(c("mother","father","child","child"),K)

  ud <- data.frame(ybin=c(Ybin),x=c(t(xs)),type=type,cluster=rep(1:K,each=n))

###names(ud)<-c("ybin","x","cluster","type")
return(ud)
} ## }}} 

##' @export
simbinClaytonOakes.twin.ace <- function(K,varg,varc,beta=NULL,alpha=NULL)  ## {{{ 
{
  ## K antal clustre (families), n=antal i clustre
###  K <- 10000
  n <- 2 # twins with ace structure
  ## total variance 1/(varg+varc)
  ## {{{  random effects 
  ###  means varg/(varg+varc) and variances varg/(varg+varc)^2
  eta <- varc+varg
  ### mother and father share environment
  ### children share half the genes with mother and father and environment 
  dz.g <-  cbind(rgamma(K/2,varg*0.5)/eta, rgamma(K/2,varg*0.5)/eta, rgamma(K/2,varg*0.5)/eta)
  mz.g <-  rgamma(K/2,varg)/eta
  env1 <- rgamma(K/2,varc)/eta 
  env2 <- rgamma(K/2,varc)/eta 
  mz   <- mz.g+env1
  dz1 <- apply(dz.g[,c(1,2)],1,sum) + env2
  dz2 <- apply(dz.g[,c(1,3)],1,sum) + env2
  Gam1 <- rbind(cbind(mz,mz),cbind(dz1,dz2))
  ## }}} 

  ## {{{ marginals p's and conditional p's given random effects 
  ### marginals p's for mother, father, children
  xb1 <- rbinom(K,1,0.5)
  xb2 <- rbinom(K,1,0.5)
###
###
  if (is.null(beta)) beta <- rep(0.3,2)
  if (is.null(alpha)) alpha <- rep(0.5,2)
  pb1 <- exp(alpha[1]+xb1*beta[1]); 
  pb2 <- exp(alpha[2]+xb2*beta[2])
  p <- cbind(pb1,pb2)
  p <- p/(1+p)

  vartot <- eta
  pgivenZ <- ilap(vartot,p)

###  pgivenZ <- ilap(vartot,p)
  pgivenZ <- exp(- Gam1*pgivenZ)
  ## }}} 

  Ybin <- matrix(rbinom(n*K,1,c(pgivenZ)),K,n)
  Ybin <- t(Ybin)
  xs <- cbind(xb1,xb2)
  type <- rep(c("twin1","twin2"),K)
  zygosity <- rep(c("MZ","DZ"),each=K)

  ud <- data.frame(ybin=c(Ybin),x=c(t(xs)),type=type,cluster=rep(1:K,each=n),zygosity=zygosity)

###names(ud)<-c("ybin","x","cluster","type")
return(ud)
} ## }}} 

##' @export
simbinClaytonOakes.pairs <- function(K,varc,beta=NULL,alpha=NULL)  ## {{{ 
{
  ## K antal clustre (families), n=antal i clustre
###  K <- 10000
  n <- 2 # twins with ace structure
  ## total variance 1/(varg+varc)
  ## {{{  random effects 
  ###  means varg/(varg+varc) and variances varg/(varg+varc)^2
  eta <- varc <- 1/varc
  ### mother and father share environment
  ### children share half the genes with mother and father and environment 
  rans <- rgamma(K,varc)/eta
  Gam1 <- cbind(rans,rans)
  ## }}} 
###  apply(Gam1,2,mean); apply(Gam1,2,var); cor(Gam1)

  ## {{{ marginals p's and conditional p's given random effects 
  ### marginals p's for mother, father, children
  x1 <- rbinom(K,1,0.5)
  x2 <- rbinom(K,1,0.5)
###
###
  if (is.null(beta)) beta <- rep(0.3,2)
  if (is.null(alpha)) alpha <- rep(0.5,2)
  p1 <- exp(alpha[1]+x1*beta[1]); 
  p2 <- exp(alpha[1]+x2*beta[2])
  p <- cbind(p1,p2)
  p <- p/(1+p)

  vartot <- eta
###  pgivenZ <- mets:::ilap(vartot,p)
  pgivenZ <- ilap(vartot,p)
  pgivenZ <- exp(- Gam1*pgivenZ)
  ## }}} 

  Ybin <- matrix(rbinom(n*K,1,c(pgivenZ)),K,n)
  Ybin <- t(Ybin)

  xs <- cbind(x1,x2)
  ud <- cbind(c(Ybin),c(t(xs)),rep(1:K,each=n))  

  ud <- cbind(ud)
  ud <- data.frame(ud)

names(ud)<-c("ybin","x","cluster")
return(ud)
} ## }}} 

### pairwise POR model based on case-control data
##' @export
CCbinomial.twostage <- function(margbin=NULL,data=sys.parent(),score.method="nlminb",
    response="response",id="id",num="num",case.num=0,
    Nit=60,detail=0, silent=1,weights=NULL, control=list(),
    theta=NULL,theta.formula=NULL,desnames=NULL,
    deshelp=0,var.link=1,iid=1,
    step=0.5,model="plackett",marginal.p=NULL,
    strata=NULL,max.clust=NULL,se.clusters=NULL)
{ ## {{{
  ### under construction 

    if (class(margbin)[1]=="glm") ps <- predict(margbin,type="response") 
    else if (class(margbin)=="formula") {
        margbin <- glm(margbin,data=data,family=binomial())
        ps <- predict(margbin,type="response")
    }  else if (is.null(marginal.p)) stop("without marginal model, marginal p's must be given\n"); 

    if (!is.null(marginal.p)) ps <- marginal.p

    data <- cbind(data,ps)

### make all pairs in the families,
    fam <- familycluster.index(data[,id])
    data.fam <- data[fam$familypairindex,]
    data.fam$subfam <- fam$subfamilyindex

### make dependency design using wide format for all pairs 
    data.fam.clust <- fast.reshape(data.fam,id="subfam")
    if (is.function(theta.formula)) {
        desfunction <- compiler::cmpfun(theta.formula)
	if (deshelp==1){
            cat("These names appear in wide version of pairs for dependence \n")
            cat("design function must be defined in terms of these: \n")
            cat(names(data.fam.clust)); cat("\n")
            cat("Here is head of wide version with pairs\n")
            print(head(data.fam.clust)); cat("\n")
	}
###	des.theta <- Reduce("rbind",lapply(seq(nrow(data.fam.clust)),function(i) unlist(desfunction(data.fam.clust[i,] ))))
        des.theta <- t(apply(data.fam.clust,1, function(x) desfunction(x)))
	colnames(des.theta) <- desnames
	desnames <- desnames
    } else {
        if (is.null(theta.formula)) theta.formula <- ~+1
        des.theta <- model.matrix(theta.formula,data=data.fam.clust)
        desnames <- colnames(des.theta); 
    }
    data.fam.clust <- cbind(data.fam.clust,des.theta)
    if (deshelp==1) {
        cat("These names appear in wide version of pairs for dependence \n")
        print(head(data.fam.clust))
    }

### back to long format keeping only needed variables
    data.fam <- fast.reshape(data.fam.clust,varying=c(response,id,"ps"))
    if (deshelp==1) {
	cat("Back to long format for binomial.twostage (head)\n"); 
        print(head(data.fam)); 
	cat("\n")
	cat(paste("binomial.twostage, called with reponse",response,"\n")); 
	cat(paste("cluster=",id,",  subcluster (pairs)=subfam \n")); 
	cat(paste("design variables =")); 
	cat(desnames)
	cat("\n")
    } 

    out <- binomial.twostage(data.fam[,response],data=data.fam,
                             clusters=data.fam$subfam,
                             theta.des=data.fam[,desnames],
                             detail=detail, score.method=score.method, Nit=Nit,step=step,
                             iid=iid,theta=theta, var.link=var.link,model=model, 
                             max.clust=max.clust,
                             marginal.p=data.fam[,"ps"], se.clusters=data.fam[,id])
    return(out)
} ## }}}


