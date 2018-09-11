##' @title Twostage survival model for multivariate survival data 
##'
##' @description 
##' Fits Clayton-Oakes or bivariate Plackett models for bivariate survival data 
##' using marginals that are on Cox form. The dependence can be 
##' modelled via 
##' \enumerate{
##' \item  Regression design on dependence parameter. 
##' \item  Random effects, additive gamma model. 
##' }
##'
##' If clusters contain more than two subjects, we use a composite likelihood
##' based on the pairwise bivariate models, for MLE see twostageMLE. 
##'
##' The two-stage model is constructed such that 
##' given the gamma distributed random effects it is assumed that the survival functions 
##' are indpendent, and that the marginal survival functions are on Cox form (or additive form)
##' \deqn{
##' P(T > t| x) = S(t|x)= exp( -exp(x^T \beta) A_0(t) )
##' }
##'
##' One possibility is to model the variance within clusters via a regression design, and
##' then one can specify a regression structure for the indenpendent gamma distributed 
##' random effect for each cluster, such that the variance is given by 
##' \deqn{
##'  \theta = z_j^T \alpha
##' }
##' where \eqn{z} is specified by theta.des 
##' The reported standard errors are based on the estimated information from the 
##' likelihood assuming that the marginals are known. 
##'
##' Can also fit a structured additive gamma random effects model, such
##' as the ACE, ADE model for survival data.  In this case the 
##' random.design specificies the random effects for each subject within a cluster. This is
##' a matrix of 1's and 0's with dimension n x d.  With d random effects. 
##' For a cluster with two subjects, we let the random.design rows be 
##'  \eqn{v_1} and \eqn{v_2}. 
##' Such that the random effects for subject 
##' 1 is \deqn{v_1^T (Z_1,...,Z_d)}, for d random effects. Each random effect
##' has an associated parameter \eqn{(\lambda_1,...,\lambda_d)}. 
##' By construction subjects 1's random effect are Gamma distributed with 
##' mean \eqn{\lambda_j/v_1^T \lambda}
##' and variance \eqn{\lambda_j/(v_1^T \lambda)^2}. Note that the random effect 
##' \eqn{v_1^T (Z_1,...,Z_d)} has mean 1 and variance \eqn{1/(v_1^T \lambda)}.
##' It is here asssumed that  \eqn{lamtot=v_1^T \lambda} is fixed within clusters 
##' as it would be for the ACE model below.
##'
##' Based on these parameters the relative contribution (the heritability, h) is 
##' equivalent to  the expected values of the random effects: \eqn{\lambda_j/v_1^T \lambda}
##'
##' The DEFAULT parametrization (var.par=1) uses the variances of the random effecs 
##' \deqn{
##' \theta_j  = \lambda_j/(v_1^T \lambda)^2
##' }
##' For alternative parametrizations one can specify how the parameters relate to \eqn{\lambda_j}
##' with the argument var.par=0.
##'
##' For both types of models the basic model assumptions are that 
##' given the random effects of the clusters the survival distributions within a cluster 
##' are independent and ' on the form 
##' \deqn{
##' P(T > t| x,z) = exp( -Z \cdot Laplace^{-1}(lamtot,lamtot,S(t|x)) )  
##' }
##' with the inverse laplace of the gamma distribution with mean 1 and variance 1/lamtot.
##'
##' The parameters \eqn{(\lambda_1,...,\lambda_d)} are related to the parameters of the model
##' by a regression construction \eqn{pard} (d x k), that links the \eqn{d} 
##' \eqn{\lambda} parameters
##' with the (k) underlying \eqn{\theta} parameters 
##' \deqn{
##' \lambda = theta.des  \theta 
##' }
##' here using theta.des to specify these low-dimension association. Default is a diagonal matrix. 
##' This can be used to make structural assumptions about the variances of the random-effects 
##' as is needed for the ACE model for example. 
##'
##' The case.control option that can be used with the pair specification of the pairwise parts
##' of the estimating equations. Here it is assumed that the second subject of each pair is the proband. 
##'
##' @references
##' Estimating heritability for cause specific mortality based on twins studies
##' Scheike, Holst, Hjelmborg (2014), LIDA  
##' 
##' Measuring early or late dependence for bivariate twin data
##' Scheike, Holst, Hjelmborg (2015), LIDA  
##' 
##' Twostage modelling of additive gamma frailty models for survival data. 
##' Scheike and Holst, working paper 
##' 
##' Additive Gamma frailty models for competing risks data, Biometrics (2015)
##' Eriksson and Scheike (2015), 
##' 
##' Shih and Louis (1995) Inference on the association parameter in copula models for bivariate
##' survival data, Biometrics, (1995).
##' 
##' Glidden (2000), A Two-Stage estimator of the dependence
##' parameter for the Clayton Oakes model, LIDA, (2000).
##' 
##' @examples
##' data(diabetes)
##' 
##' # Marginal Cox model  with treat as covariate
##' margph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
##' ### Clayton-Oakes, MLE 
##' fitco1<-twostageMLE(margph,data=diabetes,theta=1.0)
##' summary(fitco1)
##' 
##' ### Plackett model
##' mph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
##' fitp <- survival.twostage(mph,data=diabetes,theta=3.0,Nit=40,
##'                clusters=diabetes$id,var.link=1,model="plackett")
##' summary(fitp)
##' 
##' ### Clayton-Oakes
##' fitco2 <- survival.twostage(mph,data=diabetes,theta=0.0,detail=0,
##'                  clusters=diabetes$id,var.link=1,model="clayton.oakes")
##' summary(fitco2)
##' fitco3 <- survival.twostage(margph,data=diabetes,theta=1.0,detail=0,
##'                  clusters=diabetes$id,var.link=0,model="clayton.oakes")
##' summary(fitco3)
##'
##' ### without covariates but with stratafied 
##' marg <- phreg(Surv(time,status)~+strata(treat)+cluster(id),data=diabetes)
##' fitpa <- survival.twostage(marg,data=diabetes,theta=1.0,
##'                 clusters=diabetes$id,score.method="optimize")
##' summary(fitpa)
##' 
##' fitcoa <- survival.twostage(marg,data=diabetes,theta=1.0,clusters=diabetes$id,
##'                  model="clayton.oakes")
##' summary(fitcoa)
##' 
##' ### Piecewise constant cross hazards ratio modelling
##' ########################################################
##' 
##' d <- subset(simClaytonOakes(2000,2,0.5,0,stoptime=2,left=0),!truncated)
##' udp <- piecewise.twostage(c(0,0.5,2),data=d,score.method="optimize",
##'                           id="cluster",timevar="time",
##'                           status="status",model="clayton.oakes",silent=0)
##' summary(udp)
##' 
##' \donttest{ ## Reduce Ex.Timings
##' ### Same model using the strata option, a bit slower
##' ########################################################
##' ## makes the survival pieces for different areas in the plane 
##' ##ud1=surv.boxarea(c(0,0),c(0.5,0.5),data=d,id="cluster",timevar="time",status="status")
##' ##ud2=surv.boxarea(c(0,0.5),c(0.5,2),data=d,id="cluster",timevar="time",status="status")
##' ##ud3=surv.boxarea(c(0.5,0),c(2,0.5),data=d,id="cluster",timevar="time",status="status")
##' ##ud4=surv.boxarea(c(0.5,0.5),c(2,2),data=d,id="cluster",timevar="time",status="status")
##' 
##' ## everything done in one call 
##' ud <- piecewise.data(c(0,0.5,2),data=d,timevar="time",status="status",id="cluster")
##' ud$strata <- factor(ud$strata); 
##' ud$intstrata <- factor(ud$intstrata)
##' 
##' ## makes strata specific id variable to identify pairs within strata
##' ## se's computed based on the id variable across strata "cluster"
##' ud$idstrata <- ud$id+(as.numeric(ud$strata)-1)*2000
##' 
##' marg2 <- aalen(Surv(boxtime,status)~-1+factor(num):factor(intstrata),
##'                data=ud,n.sim=0,robust=0)
##' tdes <- model.matrix(~-1+factor(strata),data=ud)
##' fitp2 <- survival.twostage(marg2,data=ud,se.clusters=ud$cluster,clusters=ud$idstrata,
##'                 score.method="fisher.scoring",model="clayton.oakes",
##'                 theta.des=tdes,step=0.5)
##' summary(fitp2)
##' 
##' ### now fitting the model with symmetry, i.e. strata 2 and 3 same effect
##' ud$stratas <- ud$strata; 
##' ud$stratas[ud$strata=="0.5-2,0-0.5"] <- "0-0.5,0.5-2"
##' tdes2 <- model.matrix(~-1+factor(stratas),data=ud)
##' fitp3 <- survival.twostage(marg2,data=ud,clusters=ud$idstrata,se.cluster=ud$cluster,
##'                 score.method="fisher.scoring",model="clayton.oakes",
##'                 theta.des=tdes2,step=0.5)
##' summary(fitp3)
##' 
##' ### same model using strata option, a bit slower 
##' fitp4 <- survival.twostage(marg2,data=ud,clusters=ud$cluster,se.cluster=ud$cluster,
##'                 score.method="fisher.scoring",model="clayton.oakes",
##'                 theta.des=tdes2,step=0.5,strata=ud$strata)
##' summary(fitp4)
##' }
##' 
##' \donttest{ ## Reduce Ex.Timings
##' ### structured random effects model additive gamma ACE 
##' ### simulate structured two-stage additive gamma ACE model
##' data <- simClaytonOakes.twin.ace(4000,2,1,0,3)
##' out <- twin.polygen.design(data,id="cluster")
##' pardes <- out$pardes
##' pardes 
##' des.rv <- out$des.rv
##' head(des.rv)
##' aa <- phreg(Surv(time,status)~x+cluster(cluster),data=data,robust=0)
##' ts <- survival.twostage(aa,data=data,clusters=data$cluster,detail=0,
##' 	       theta=c(2,1),var.link=0,step=0.5,
##' 	       random.design=des.rv,theta.des=pardes)
##' summary(ts)
##' }
##' 
##' @keywords survival
##' @author Thomas Scheike
##' @param margsurv Marginal model 
##' @param data data frame
##' @param score.method Scoring method "fisher.scoring", "nlminb", "optimize", "nlm"
##' @param Nit Number of iterations
##' @param detail Detail
##' @param clusters Cluster variable
##' @param silent Debug information
##' @param weights Weights
##' @param control Optimization arguments
##' @param theta Starting values for variance components
##' @param theta.des design for dependence parameters, when pairs are given this is could be a
##' (pairs) x (numer of parameters)  x (max number random effects) matrix
##' @param var.link Link function for variance 
##' @param iid Calculate i.i.d. decomposition
##' @param step Step size
##' @param model model
##' @param marginal.trunc marginal left truncation probabilities
##' @param marginal.survival optional vector of marginal survival probabilities 
##' @param marginal.status related to marginal survival probabilities 
##' @param strata strata for fitting, see example
##' @param se.clusters for clusters for se calculation with iid
##' @param numDeriv to get numDeriv version of second derivative, otherwise uses sum of squared score 
##' @param random.design random effect design for additive gamma model, when pairs are given this is
##' a (pairs) x (2) x (max number random effects) matrix, see pairs.rvs below
##' @param pairs matrix with rows of indeces (two-columns) for the pairs considered in the pairwise
##' composite score, useful for case-control sampling when marginal is known.
##' @param pairs.rvs for additive gamma model and random.design and theta.des are given as arrays,
##' this specifice number of random effects for each pair. 
##' @param numDeriv.method uses simple to speed up things and second derivative not so important.
##' @param additive.gamma.sum for two.stage=0, this is specification of the lamtot in the models via
##' a matrix that is multiplied onto the parameters theta (dimensions=(number random effects x number
##' of theta parameters), when null then sums all parameters.
##' @param var.par is 1 for the default parametrization with the variances of the random effects,
##' var.par=0 specifies that the \eqn{\lambda_j}'s are used as parameters. 
##' @param cr.models competing risks models for two.stage=0, should be given as a list with models for each cause
##' @param case.control assumes case control structure for "pairs" with second column being the probands,
##' when this options is used the twostage model is profiled out via the paired estimating equations for the
##' survival model. 
##' @param ascertained if the pair are sampled only when there is an event. This is in contrast to
##' case.control sampling where a proband is given. This can be combined with control probands. Pair-call
##' of twostage is needed  and second column of pairs are the first jump time with an event for ascertained pairs,
##' or time of control proband.
##' @param shut.up to make the program more silent in the context of iterative procedures for case-control
##' and ascertained sampling
##' @aliases survival.twostage survival.twostage.fullse twostage.aalen twostage.cox.aalen twostage.coxph twostage.phreg randomDes
##' @export survival.twostage
survival.twostage <- function(margsurv,data=sys.parent(),
    score.method="fisher.scoring",Nit=60,detail=0,clusters=NULL,
    silent=1,weights=NULL, control=list(),theta=NULL,theta.des=NULL,
    var.link=1,iid=1,step=0.5,model="clayton.oakes",
    marginal.trunc=NULL,marginal.survival=NULL,marginal.status=NULL,strata=NULL,
    se.clusters=NULL,numDeriv=0,random.design=NULL,pairs=NULL,pairs.rvs=NULL,
    numDeriv.method="simple",additive.gamma.sum=NULL,var.par=1,cr.models=NULL,
    case.control=0,ascertained=0,shut.up=0)
{## {{{
## {{{ seting up design and variables
two.stage <- 1; rate.sim <- 1; sym=1; var.func <- NULL
if (model=="clayton.oakes" || model=="gamma") dep.model <- 1
else if (model=="plackett" || model=="or") dep.model <- 2
else stop("Model must by either clayton.oakes or plackett \n"); 
start.time <- NULL; ptrunc <- NULL; psurvmarg <- NULL; status <- NULL; score.iid <- NULL
fix.baseline <- 0; convergence.bp <- 1;  ### to control if baseline profiler converges

  if ((!is.null(margsurv)) | (!is.null(marginal.survival))) fix.baseline <- 1
  antpers <- nrow(data); RR <-  rep(1,antpers); 

  if (!is.null(margsurv)) {
     rrr <-  readmargsurv(margsurv,data,clusters)
     psurvmarg <- rrr$psurvmarg; ptrunc <- rrr$ptrunc; start.time <- rrr$entry;
     time2 <- rrr$exit; status <- rrr$status; clusters <- rrr$clusters
  }

  if (is.null(psurvmarg)) psurvmarg <- rep(1,antpers);
  if (!is.null(marginal.survival)) psurvmarg <- marginal.survival 
  if (!is.null(marginal.trunc)) ptrunc <- marginal.trunc 
  if (is.null(ptrunc)) ptrunc <- rep(1,length(psurvmarg))
  if (!is.null(marginal.status)) status <- marginal.status 
  if (is.null(status) & is.null(cr.models)) stop("must give status variable for survival via either margninal model (margsurv), marginal.status or as cr.models \n"); 
  if (is.null(weights)==TRUE) weights <- rep(1,antpers); 
  if (is.null(strata)==TRUE) strata<- rep(1,antpers); 
  if (length(strata)!=antpers) stop("Strata must have length equal to number of data points \n"); 

  ## {{{ cluster set up
  cluster.call <- clusters
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

###  if ((!is.null(max.clust))) if (max.clust< antiid) {
###        coarse.clust <- TRUE
###	qq <- unique(quantile(se.clusters, probs = seq(0, 1, by = 1/max.clust)))
###	qqc <- cut(se.clusters, breaks = qq, include.lowest = TRUE)    
###	se.clusters <- as.integer(qqc)-1
###	max.clusters <- length(unique(se.clusters))
###	maxclust <- max.clust    
###	antiid <- max.clusters
###  }                                                        
###  if (maxclust==1) stop("No clusters, maxclust size=1\n"); 
  ## }}} 

  ### setting design for random variables, in particular with pairs are given
  ddd <- randomDes(dep.model,random.design,theta.des,theta,antpers,additive.gamma.sum,pairs,pairs.rvs,var.link,clusterindex)
  random.design=ddd$random.design;clusterindex=ddd$clusterindex;
  antpairs=ddd$antpairs; pairs.rvs=ddd$pairs.rvs;
  theta=ddd$theta;ptheta=ddd$ptheta;theta.des=ddd$theta.des
  pair.structure=ddd$pair.structure; dep.model=ddd$dep.model
  dim.rv <- ddd$dim.rv; additive.gamma.sum=ddd$additive.gamma.sum

  theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 
  ## }}}

###  setting up arguments for Aalen baseline profile estimates
 if (fix.baseline==0)  { ## {{{  when baseline is estimated 
 if (is.null(cr.models)) stop("give hazard models for different causes, ex cr.models=list(Surv(time,status==1)~+1,Surv(time,status==2)~+1) \n")

      if (case.control==0 & ascertained==0) { ## {{{ 

       	#### organize subject specific random variables and design
        ###  for additive gamma model
	## {{{ 
	dimt <- dim(theta.des[,,1,drop=FALSE])[-3]
	dimr <- dim(random.design[,,,drop=FALSE])
	mtheta.des <- array(0,c(dimt,nrow(data)))
	mrv.des <- array(0,c(dimr[1]/2,dimr[2],nrow(data)))
	nrv.des <- rep(0,nrow(data))
	nrv.des[pairs[,1]] <- pairs.rvs
	nrv.des[pairs[,2]] <- pairs.rvs
	mtheta.des[,,pairs[,1]] <- theta.des
	mtheta.des[,,pairs[,2]] <- theta.des
	mrv.des[,,pairs[,1]] <- random.design[1:(dimr[1]/2),,,drop=FALSE]
	mrv.des[,,pairs[,2]] <- random.design[(dimr[1]/2+1):dimr[1],,,drop=FALSE]
	### array thetades to jump times (subjects)
	mtheta.des <- mtheta.des[,,ids,drop=FALSE]
	### array randomdes to jump times (subjects)
	mrv.des <- mrv.des[,,ids,drop=FALSE]
	nrv.des <- pairs.rvs[ids]
	## }}} 
      } ## }}} 

      if (case.control==1 || ascertained==1) { ## {{{ 

###      print(dim(data)); print(summary(pairs))
         data1 <-        data[pairs[,1],]
         data.proband <- data[pairs[,2],]
	 weights1 <-     weights[pairs[,1]]
###	print(summary(data.proband)); print(summary(data1))

## {{{ setting up designs for jump times 
        timestatus <- all.vars(cr.models[[1]])
        if (is.null(status)) status <- data[,timestatus[2]]
	alltimes      <- data[,timestatus[1]]
	times         <- data1[,timestatus[1]]
	lstatus       <- data1[,timestatus[2]]
	timescase     <- data.proband[,timestatus[1]]
	lstatuscase   <- data.proband[,timestatus[2]]
	### organize increments according to overall jump-times 
	jumps         <- lstatus!=0
	dtimes        <- times[jumps]
	dtimescase    <- timescase[jumps]
	st            <- order(dtimes)
	dtimesst      <- dtimes[st]
	dtimesstcase  <- dtimescase[st]
	dcauses       <- lstatus[jumps][st]
	dcausescase   <- lstatuscase[jumps][st]
	ids           <- (1:nrow(data1))[jumps][st]
	###
	### delayed entry for case because of ascertained sampling
	### controls are however control probands, and have entry=0 
	entry <- timescase*lstatuscase
	data1$entry <- entry
	cr.models2 <- list()
	if (ascertained==1) {
           for (i in 1:length(cr.models)) {
	      cr.models2[[i]] <- update(cr.models[[i]],as.formula(paste("Surv(entry,",timestatus[1],",",timestatus[2],")~.",sep="")))
	    }
	} else cr.models2 <- cr.models

        nc <- 0
	for (i in 1:length(cr.models)) {
              X <- aalen.des(as.formula(cr.models[[i]]),data=data1)$X
	      nc <- nc+ncol(X)
	}
	dBaalen    <- matrix(0,length(dtimes),nc)
	xjump      <- array(0,c(length(cr.models),nc,length(ids)))
	xjumpcase  <- array(0,c(length(cr.models),nc,length(ids)))

	## first compute marginal aalen models for all causes
	a <- list(); da <- list(); 
        ### starting values for iteration 
        Bit <- Bitcase <- c()
	for (i in 1:length(cr.models)) { ## {{{ 
	      a[[i]]  <-  aalen(as.formula(cr.models2[[i]]),data=data1,robust=0,weights=weights1)
	      da[[i]] <- apply(a[[i]]$cum[,-1,drop=FALSE],2,diff)
	      jumpsi  <- (1:length(dtimes))[dcauses==i]
              X       <- aalen.des(as.formula(cr.models[[i]]),data=data1)$X
              Xcase   <- aalen.des(as.formula(cr.models[[i]]),data=data.proband)$X
	      if (i==1) fp <- 1 
	      indexc  <- fp:(fp+ncol(X)-1)
	      dBaalen[jumpsi,indexc] <- da[[i]]
              xjump[i,indexc,]       <- t(X[ids,])
              xjumpcase[i,indexc,]   <- t(Xcase[ids,])
	      fp <- fp+ncol(X)
	      ### starting values 
	      Bit <- cbind(Bit,Cpred(a[[i]]$cum,dtimesst)[,-1,drop=FALSE])
       } ## }}} 
       Bit.ini <- Bit
       ## }}} 

####  organize subject specific random variables and design
####  already done in basic pairwise setup 
	mtheta.des <- theta.des[,,ids,drop=FALSE]
	### array randomdes to jump times (subjects)
	mrv.des <- random.design[,,ids,drop=FALSE]
	nrv.des <- pairs.rvs[ids]
      } ## }}} 

 }  else {
	 mrv.des <- array(0,c(1,1,1)); mtheta.des <- array(0,c(1,1,1)); margthetades <- array(0,c(1,1,1)); 
	 xjump <- array(0,c(1,1,1)); dBaalen <- matrix(0,1,1); nrv.des <- 3
 } ## }}} 

  loglike <- function(par) 
  { ## {{{

     if (pair.structure==0 | dep.model!=3) Xtheta <- as.matrix(theta.des) %*% matrix(c(par),nrow=ptheta,ncol=1);
     if (pair.structure==1 & dep.model==3) Xtheta <- matrix(0,antpers,1); ## not needed 
     DXtheta <- array(0,c(1,1,1));

      if (var.link==1 & dep.model==3) epar <- c(exp(par)) else epar <- c(par)
      partheta <- epar

      if (var.par==1 & dep.model==3) {
	 ## from variances to 
	 if (is.null(var.func)) {
	    sp <- sum(epar)
	    partheta <- epar/sp^2 
         } else partheta <- epar ## par.func(epar)
      } 


      if (pair.structure==0) {# {{{
	      outl<-.Call("twostageloglikeRV", ## {{{ only two stage model for this option
	      icause=status,ipmargsurv=psurvmarg, 
	      itheta=c(partheta),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),ithetades=theta.des,
	      icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
	      ivarlink=var.link,iid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
	      itrunkp=ptrunc,istrata=as.numeric(strata),iseclusters=se.clusters,iantiid=antiid,
	      irvdes=random.design,iags=additive.gamma.sum,iascertained=ascertained,
              PACKAGE="mets") ## }}}
      }# }}}
      else { ## pair-structure 
	     ## twostage model,    case.control option,   profile out baseline 
	     ## conditional model, case.control option,   profile out baseline 
         if (fix.baseline==0) ## if baseline is not given 
         {
          cum1 <- cbind(dtimesst,Bit)
          if ( (case.control==1 || ascertained==1) & (convergence.bp==1)) { ## {{{  profiles out baseline under case-control/ascertainment sampling

###	      ## initial values , only one cr.model for survival 
###           Bit <- cbind(Cpred(a[[1]]$cum,dtimesst)[,-1])
             if (detail>1) plot(dtimesst,Bit,type="l",main="Bit")

             if (ncol(Bit)==0) Bit <- Bit.ini
             Bitcase  <- Cpred(cbind(dtimesst,Bit),dtimesstcase)[,-1,drop=FALSE]
             Bitcase <- .Call("MatxCube",Bitcase,dim(xjumpcase),xjumpcase,PACKAGE="mets")$X

             for (j in 1:5) { ## {{{ profile via iteration 
                   cncc <- .Call("BhatAddGamCC",1,dBaalen,dcauses,dim(xjump),xjump,
		                 c(partheta),dim(mtheta.des),mtheta.des,additive.gamma.sum,var.link, 
				 dim(mrv.des),mrv.des,nrv.des,1,Bit,Bitcase,dcausescase,PACKAGE="mets")
		   d <- max(abs(Bit-cncc$B))
		   if (detail>1) print(d)
		   Bit <- cncc$B
###		   if (detail>1) print(c(par,epar,partheta)); 
###		   if (detail>1) print(summary(Bit)); 
		   if (detail>1) print(summary(cncc$caseweights))
		   cum1 <- cbind(dtimesst,cncc$B)
		   Bitcase  <-cbind(Cpred(cum1,dtimesstcase)[,-1])
###		   if (detail>1) print(summary(Bitcase))
		   if (detail>1) lines(dtimesst,Bit,col=j+1); 
		   if (is.na(d)) { 
			  if (shut.up==0) cat("Baseline profiler gives missing values\n");  
		          Bit <- Bit.ini; cum1 <- cbind(dtimesst,Bit); convergence.bp <<- 0; break;
		   }
		   Bitcase <- .Call("MatxCube",Bitcase,dim(xjumpcase),xjumpcase,PACKAGE="mets")$X
		   if (d<0.00001) break; 
           } ## }}} 


	   nulrow    <- rep(0,ncol(Bit)+1)
	   pbases    <- Cpred(rbind(nulrow,cbind(dtimesst,Bit)),alltimes)[,-1,drop=FALSE]
           X         <- aalen.des(as.formula(cr.models[[1]]),data=data)$X
	   psurvmarg <- exp(-apply(X*pbases,1,sum))  ## psurv given baseline 
	   if (ascertained==1) {
              Xcase     <- aalen.des(as.formula(cr.models[[1]]),data=data.proband)$X
              X         <- aalen.des(as.formula(cr.models[[1]]),data=data1)$X
	      pba.case  <- Cpred(rbind(nulrow,cbind(dtimesst,Bit)),entry)[,-1,drop=FALSE]
	      ptrunc    <- rep(0,nrow(data))
	      ### for control probands ptrunc=1, thus no adjustment
	      ptrunc[pairs[,1]]     <- exp(-apply(X*    pba.case,1,sum)*lstatuscase) ## delayed entry at time of ascertainment proband 
 	      ptrunc[pairs[,2]]     <- exp(-apply(Xcase*pba.case,1,sum)*lstatuscase)
	   }

          } ## }}} 
          }


          outl<-.Call("twostageloglikeRVpairs", ## {{{
          icause=status,ipmargsurv=psurvmarg, 
          itheta=c(partheta),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),
	  ithetades=theta.des,
          icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
          ivarlink=var.link,iiid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
          itrunkp=ptrunc,istrata=as.numeric(strata),iseclusters=se.clusters,iantiid=antiid,
          irvdes=random.design,
          idimthetades=dim(theta.des),idimrvdes=dim(random.design),irvs=pairs.rvs,
	  iags=additive.gamma.sum, 
	  iascertained=ascertained,PACKAGE="mets")
	  ## }}} 

          if (fix.baseline==0)  {
              outl$baseline <- cum1; 
	      outl$marginal.surv <- psurvmarg; 
	      outl$marginal.trunc <- ptrunc
	  }

      } ## }}} 

    if (detail==3) print(c(partheta,outl$loglike))

    ## variance parametrization, and inverse.link 
    if (dep.model==3) {# {{{
    if (var.par==1) {
         ## from variances to and with sum for all random effects 
         if (is.null(var.func)) {
	 if (var.link==0)  {
###		 print(c(sp,epar))
             mm <- matrix(-epar*2*sp,length(epar),length(epar))
	     diag(mm) <- sp^2-epar*2*sp
###	    print(mm)
	 } else {
            mm <- -c(epar) %o% c(epar)*2*sp
            diag(mm) <- epar*sp^2-epar^2*2*sp
###	    print(mm)
	 }
	    mm <- mm/sp^4
	 } else mm  <- numDeriv::hessian(var.func,par)
      } else {
	   if (var.link==0) mm <- diag(length(epar)) else 
			  mm <- diag(length(c(epar)))*c(epar)
      }
      }# }}}

    if (dep.model==3) {# {{{
       outl$score <-  t(mm) %*% outl$score
       outl$Dscore <- t(mm) %*% outl$Dscore %*% mm 
       if (iid==1) { outl$theta.iid <- t(t(mm) %*% t(outl$theta.iid))
               if (class(margsurv)=="phreg") {  
	          outl$D1thetal  <- t(t(mm) %*% t(outl$D1thetal))
	          outl$D2thetal  <- t(t(mm) %*% t(outl$D2thetal))
	       }
       }
    }# }}}

    attr(outl,"gradient") <-outl$score 
    if (oout==0) ret <- c(-1*outl$loglike) else if (oout==1) ret <- sum(outl$score^2) else if (oout==2) ret <- outl else ret <- outl$score
    return(ret)
  } ## }}}

  if (score.method=="optimize" && ptheta!=1) {
     cat("optimize only works for d==1, score.mehod set to nlminb \n"); 
     score.method <- "nlminb";
  } 

  score1 <- NULL
  theta.iid <- NULL
  logl <- NULL
  p <- theta

    if (score.method=="fisher.scoring") { ## {{{
        oout <- 2;  ### output control for obj
        if (Nit>0) 
            for (i in 1:Nit)
            {
                    out <- loglike(p)
		    ## updating starting values for cumulative baselines
		    if (fix.baseline==0) Bit <- out$baseline[,-1,drop=FALSE]
		    if (fix.baseline==1) hess <-  out$Dscore
                    ###  uses simple second derivative for computing derivative of score
		    if (numDeriv==2 || (((fix.baseline==0)) & (i==1))) {
			    oout <- 3
			    hess <- numDeriv::jacobian(loglike,p,method="simple")
			    oout <- 2
		    }
                    if (!is.na(sum(hess))) hessi <- lava::Inverse(hess) else hessi <- hess 
                    if (detail==1) {## {{{
                        cat(paste("Fisher-Scoring ===================: it=",i,"\n")); 
                        cat("theta:");print(c(theta))
                        cat("loglike:");cat(c(out$loglike),"\n"); 
                        cat("score:");cat(c(out$score),"\n"); 
                        cat("hess:\n"); cat(hess,"\n"); 
                    }## }}}
                    delta <- step*( hessi %*% out$score )
		    ### update p, but note that score and derivative in fact related to previous p 
		    ### unless Nit=0, 
	            if (Nit>0) {
                    p <- p - delta
                    theta <- p; 
		    }
                    if (is.nan(sum(out$score))) break; 
                    if (sum(abs(out$score))<0.00001) break; 
                    if (max(theta)>20 & var.link==1) { cat("theta too large lacking convergence \n"); break; }
           }
        if (!is.nan(sum(p))) { 
            if (detail==1 && iid==1) cat("iid decomposition\n"); 
            out <- loglike(p) 
            logl <- out$loglike
            score1 <- score <- out$score
            oout <- 0; 
            hess1 <- hess  <- -1*out$Dscore 
           if (iid==1) {  ## {{{
		theta.iid <- out$theta.iid
		ucc <-  unique(cluster.call) 
		if (length(ucc)== nrow(theta.iid))
		    rownames(theta.iid) <- unique(cluster.call) 
	        score.iid <-   theta.iid

                if (class(margsurv)=="phreg") {  ## {{{

                  if (is.null(margsurv$opt) | is.null(margsurv$coef)) fixbeta<- 1 else fixbeta <- 0
		  xx <- margsurv$cox.prep
                  id  <-  xx$id+1
		  cumhazt <- rrr$cum
		  rr <- rrr$RR

		  ## these refers to id given in cluster.call 
		  xx <- margsurv$cox.prep
		  D1thetal<- out$D1thetal
		  D2thetal<- out$D2thetal
		  Dlamthetal <- (D1thetal+D2thetal)
                  Fbeta <- - t(D1thetal+D2thetal) %*% (margsurv$X *cumhazt*psurvmarg*rr)

                  ### timeordered 
                  Gbase <- apply(-Dlamthetal[xx$ord+1,,drop=FALSE]*psurvmarg[xx$ord+1]*rr[xx$ord+1],2,revcumsumstrata,xx$strata,xx$nstrata)
###		  print(c(Gbase)); print(Fbeta)

                  if (fixbeta==0) { # {{{ iid after beta of marginal model 
		  Z <- xx$X
		  U <- E <- matrix(0,nrow(xx$X),margsurv$p)
		  E[xx$jumps+1,] <- margsurv$E
		  U[xx$jumps+1,] <- margsurv$U
		  invhess <- -solve(margsurv$hessian)
		  S0i <- rep(0,length(xx$strata))
		  S0i[xx$jumps+1] <- 1/margsurv$S0
		  cumhaz <- c(cumsumstrata(S0i,xx$strata,xx$nstrata))
		  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
		  rr <- c(xx$sign*exp(Z %*% coef(margsurv) + xx$offset))
		  ### Martingale  as a function of time and for all subjects to handle strata 
		  MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0)*rr*c(xx$weights)
		  UUbeta <- apply(MGt,2,sumstrata,id-1,max(id))
	          uxid <- unique(margsurv$id)
		  rownames(UUbeta) <- uxid 
		  UUbeta <- UUbeta %*% invhess
		  GbaseEdLam0 <- t(Gbase) %*% (E*S0i) 
		  Fbeta   <-  Fbeta -  GbaseEdLam0
                  Fbetaiid <- UUbeta %*% t(Fbeta)
	          }# }}}

	  ###  \int_0^\tau (GBase(s) / S_0(s)) dN_i(s) - dLamba_i(s)
	  xx <- margsurv$cox.prep
	  S0i  <- S0i2 <- rep(0,length(xx$strata))
	  S0i[xx$jumps+1] <- 1/margsurv$S0
	  S0i2[xx$jumps+1] <- 1/margsurv$S0^2

	  GbasedN <- Gbase*S0i
	  GbasedLam0 <- apply(Gbase*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)

	  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(margsurv) + xx$offset)) else rr <- c(xx$sign*exp( xx$offset))
	  MGt <- (GbasedN[,drop=FALSE]-GbasedLam0*rr)*c(xx$weights)
	  MGti <- apply(MGt,2,sumstrata,id-1,max(id))
	  if (fixbeta==1) uxid <- unique(margsurv$id)
	  rownames(MGti) <- uxid

###          if (fixbeta==0) print(cbind(theta.iid,MGti,Fbetaiid)) else print(cbind(theta.iid,MGti))

	  ## add id's after rownames that may not be exactly the same in marginal model and in pairs for fitting 
	  if (all.equal(uxid,ucc)==TRUE) {
	     theta.iid <-  theta.iid + MGti
	     if (fixbeta==0) theta.iid  <-  theta.iid + Fbetaiid
	  } else {
	     allid <-  unique(c(uxid,ucc))
	     theta.iid <- matrix(0,length(allid),ncol(theta.iid))
	     rownames(theta.iid) <- allid
	     theta.iid[rownames(score.iid),] <- score.iid 
	     theta.iid[rownames(MGti),] <- theta.iid[rownames(MGti),] + MGti
	     if (fixbeta==0) theta.iid[rownames(MGti),] <-  theta.iid[rownames(MGti),]+Fbetaiid
	  }
	  out$theta.iid <- theta.iid
	  } ## }}}
   } ## }}} 


            if (detail==1 && iid==1) cat("finished iid decomposition\n"); 
	    ### for profile solutions update second derivative at final 
	    if (numDeriv==2 || ((fix.baseline==0))) {
                    oout <- 3
                    hess <- numDeriv::jacobian(loglike,p,method="simple")
		    oout <- 2
		    }
        }
        if (numDeriv>=1) {
            if (detail==1 ) cat("starting numDeriv for second derivative \n"); 
            oout <- 0; 
            score2 <- numDeriv::jacobian(loglike,p)
	    score1 <- matrix(score2,ncol=1)
            oout <- 3
            hess <- numDeriv::jacobian(loglike,p,method="simple")
            if (detail==1 ) cat("finished numDeriv for second derivative \n"); 
        }
        if (detail==1 & Nit==0) {## {{{
            cat(paste("Fisher-Scoring ===================: final","\n")); 
            cat("theta:");print(c(p))
            cat("loglike:");cat(c(out$loglike),"\n"); 
            cat("score:");cat(c(out$score),"\n"); 
            cat("hess:\n"); cat(hess,"\n"); 
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
    hess1 <- hess <- -1* out$Dscore
    if (iid==1) { theta.iid <- out$theta.iid; score.iid <- out$theta.iid }
    if (numDeriv==1) {
    if (detail==1 ) cat("numDeriv hessian start\n"); 
      oout <- 3; ## returns score 
      hess <- numDeriv::jacobian(loglike,opt$par)
    if (detail==1 ) cat("numDeriv hessian done\n"); 
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
    hess1 <- hess <- -1* out$Dscore
    if (numDeriv==1) {
    if (detail==1 ) cat("numDeriv hessian start\n"); 
      oout <- 3;  ## to get jacobian
      hess <- numDeriv::jacobian(loglike,theta)
    if (detail==1 ) cat("numDeriv hessian done\n"); 
    }
    hessi <- lava::Inverse(hess); 
    if (iid==1) { theta.iid <- out$theta.iid; score.iid <- out$theta.iid }
  ## }}}
  } else if (score.method=="nlm") { ## {{{ nlm optimizer
    iid <- 0; oout <- 0; 
    tryCatch(opt <- nlm(loglike,theta,hessian=TRUE,print.level=detail),error=function(x) NA)
    iid <- 1; 
    hess <- opt$hessian
    score <- opt$gradient
    if (detail==1) print(opt); 
    hessi <- lava::Inverse(hess); 
    theta <- opt$estimate
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    oout <- 2
    out <- loglike(opt$estimate)
    logl <- out$loglike
    score1 <- out$score
    hess1 <- out$Dscore
    if (iid==1) { theta.iid <- out$theta.iid; score.iid <- out$theta.iid }
  ## }}}
  }  else stop("score.methods = optimize(dim=1) nlm nlminb fisher.scoring\n"); 

## {{{ handling output
  loglikeiid <- NULL
  robvar.theta <- NULL
  likepairs <- NULL
  if (fix.baseline==1)  { marginal.surv <- psurvmarg; marginal.trunc <- ptrunc; } else { marginal.surv <- out$marginal.surv; marginal.trunc <- out$marginal.trunc;} 

  if (iid==1) {
  if (dep.model==3 & pair.structure==1) likepairs <- out$likepairs
  if (dep.model==3 & two.stage==0) {# {{{
	  hessi <- -1*hessi
	  all.likepairs <- out$all.likepairs
	  colnames(all.likepairs) <- c("surv","dt","ds","dtds","cause1","cause2")
  }# }}}
     theta.iid <- out$theta.iid %*% hessi
### rownames not set to make more robust 
###     if (is.null(call.secluster)) rownames(theta.iid) <- unique(cluster.call) else rownames(theta.iid) <- unique(se.clusters)
     robvar.theta  <- crossprod(theta.iid) 
     loglikeiid <- out$loglikeiid
    } else { all.likepairs <- NULL}
  var.theta <- robvar.theta
  if (is.null(robvar.theta)) var.theta <- hessi

  if (!is.null(colnames(theta.des))) thetanames <- colnames(theta.des) else thetanames <- paste("dependence",1:length(theta),sep="")
  theta <- matrix(theta,length(c(theta)),1)
  if (length(thetanames)==nrow(theta)) { rownames(theta) <- thetanames; rownames(var.theta) <- colnames(var.theta) <- thetanames; }
  if (!is.null(logl)) logl <- -1*logl

  if (convergence.bp==0) theta <- rep(NA,length(theta))
  ud <- list(theta=theta,score=score,hess=hess,hessi=hessi,var.theta=var.theta,
	     model=model,robvar.theta=robvar.theta,
             theta.iid=theta.iid,loglikeiid=loglikeiid,likepairs=likepairs,
	     thetanames=thetanames,loglike=logl,score1=score1,Dscore=out$Dscore,
	     marginal.surv=marginal.surv,marginal.trunc=marginal.trunc,
	     baseline=out$baseline,se=diag(robvar.theta)^.5,
	     score.iid=score.iid,theta.des=theta.des,random.design=random.design)
  class(ud) <- "mets.twostage" 
  attr(ud,"response") <- "survival"
  attr(ud,"Formula") <- formula
  attr(ud,"clusters") <- clusters
  attr(ud,"cluster.call") <- cluster.call
  attr(ud,"secluster") <- c(se.clusters)
  attr(ud,"sym")<-sym; 
  attr(ud,"var.link")<-var.link; 
  attr(ud,"var.par")<-var.par; 
  attr(ud,"var.func")<-var.func; 
  attr(ud,"ptheta")<-ptheta
  attr(ud,"antpers")<-antpers; 
  attr(ud,"antclust")<-antclust; 
  attr(ud,"Type") <- model
  attr(ud,"twostage") <- two.stage
  attr(ud,"additive-gamma") <- (dep.model==3)*1
  if (!is.null(marginal.trunc)) attr(ud,"trunclikeiid")<- out$trunclikeiid
  if (dep.model==3 & two.stage==0) attr(ud,"all.likepairs")<- all.likepairs
  if (dep.model==3 ) attr(ud,"additive.gamma.sum") <- additive.gamma.sum
  #likepairs=likepairs,##  if (dep.model==3 & pair.structure==1) attr(ud, "likepairs") <- c(out$likepairs)
  if (dep.model==3 & pair.structure==0) attr(ud, "pardes") <-    theta.des
  if (dep.model==3 & pair.structure==0) attr(ud, "theta.des") <- theta.des
  if (dep.model==3 & pair.structure==1) attr(ud, "pardes") <-    theta.des[,,1]
  if (dep.model==3 & pair.structure==1) attr(ud, "theta.des") <- theta.des[,,1]
  if (dep.model==3 & pair.structure==0) attr(ud, "rv1") <- random.design[1,,drop=FALSE]
  if (dep.model==3 & pair.structure==1) attr(ud, "rv1") <- random.design[,,1]
  return(ud);
  ## }}}

} ## }}}


##' @export
randomDes <- function(dep.model,random.design,theta.des,theta,antpers,ags,pairs,pairs.rvs,var.link,clusterindex)
{# {{{
   additive.gamma.sum <- ags

   if (!is.null(random.design)) { ### different parameters for Additive random effects # {{{
     dep.model <- 3
     dim.rv <- ncol(random.design); 
     if (is.null(theta.des)) theta.des<-diag(dim.rv);
 } else { random.design <- matrix(0,1,1);  dim.rv <- 1; 
           additive.gamma.sum <- matrix(1,1,1); 
   }

  if (is.null(theta.des)) ptheta<-1; 
  if (is.null(theta.des)) theta.des<-matrix(1,antpers,ptheta); ###  else theta.des<-as.matrix(theta.des); 

  if (length(dim(theta.des))==3) ptheta<-dim(theta.des)[2] else if (length(dim(theta.des))==2) ptheta<-ncol(theta.des)
###  if (nrow(theta.des)!=antpers & dep.model!=3 ) stop("Theta design does not have correct dim");

  if (length(dim(theta.des))!=3) theta.des <- as.matrix(theta.des)

  if (is.null(theta)==TRUE) {
         if (var.link==1) theta<- rep(-0.7,ptheta);  
         if (var.link==0) theta<- rep(exp(-0.7),ptheta);   
  }       

  if (length(theta)!=ptheta) {
###	 warning("dimensions of theta.des and theta do not match\n"); 
         theta<-rep(theta[1],ptheta); 
  }
  theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 
  antpairs <- 1; ### to define 

  if (is.null(additive.gamma.sum)) additive.gamma.sum <- matrix(1,dim.rv,ptheta)
  if (!is.null(pairs)) { pair.structure <- 1;} else  pair.structure <- 0;  

  if (pair.structure==1 & dep.model==3) { ## {{{ 
### something with dimensions of rv.des 
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
       clusterindex <- pairs-1; 
  } ## }}} 

  if (pair.structure==1 & dep.model!=3) {
       clusterindex <- pairs-1; 
       antpairs <- nrow(pairs); 
       pairs.rvs <- 1
  }# }}}
  if (is.null(pairs.rvs)) pairs.rvs <- 1

  return(list(random.design=random.design,clusterindex=clusterindex,
	 antpairs=antpairs,pair.structure=pair.structure,
	 dep.model=dep.model,dim.rv=dim.rv,
	 additive.gamma.sum=additive.gamma.sum,
	 pairs.rvs=pairs.rvs,theta=theta,ptheta=ptheta,theta.des=theta.des)) 

}# }}}

readmargsurv <- function(margsurv,data,clusters)
{# {{{
start.time <- 0
if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen") {
	 formula<-attr(margsurv,"Formula");
	 beta.fixed <- attr(margsurv,"beta.fixed")
	 if (is.null(beta.fixed)) beta.fixed <- 1; 
	 ldata<-aalen.des(formula,data=data,model="cox.aalen");
	 id <- attr(margsurv,"id"); 
	 mclusters <- attr(margsurv,"cluster.call")
	 X<-ldata$X; 
	 time<-ldata$time2; 
	 Z<-ldata$Z;  
	 status<-ldata$status;
	 time2 <- attr(margsurv,"stop"); 
	 start.time <- attr(margsurv,"start")
	 antpers<-nrow(X);
         if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else { 
		     Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
         if (npar==TRUE) {Z<-matrix(0,antpers,1); pz<-1; fixed<-0;} else {fixed<-1;pz<-ncol(Z);}
	 px<-ncol(X);

         if (is.null(clusters) && is.null(mclusters)) stop("No cluster variabel specified in marginal or twostage call\n"); 
         if (is.null(clusters)) clusters <- mclusters 
	 if (nrow(X)!=length(clusters)) stop("Length of Marginal survival data not consistent with cluster length\n"); 
	 } else if (class(margsurv)=="phreg") { 
         antpers <- length(margsurv$id)
	start.time <- margsurv$entry 
} else { ### coxph 
	  antpers <- margsurv$n
	  id <- 0:(antpers-1)
	  mt <- model.frame(margsurv)
	  Y <- model.extract(mt, "response")
	  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
	  if (attr(Y, "type") == "right") {
	      time2 <- Y[, "time"]; 
	      status <- Y[, "status"]
		start.time <- rep(0,antpers);
		} else {
		 start.time <- Y[, 1]; 
		 time2 <- Y[, 2];
		 status <- Y[, 3];
		}
###	   Z <- na.omit(model.matrix(margsurv)[,-1]) ## Discard intercept
	   Z <- matrix(1,antpers,length(coef(margsurv)));

	   if (is.null(clusters)) stop("must give clusters for coxph\n");
	   cluster.call <- clusters; 
	   X <- matrix(1,antpers,1); ### Z <- matrix(0,antpers,1); ### no use for these
	   px <- 1; pz <- ncol(Z); 
	   start <- rep(0,antpers);
	   beta.fixed <- 0
	   semi <- 1
}

  if (any(abs(start.time)>0)) lefttrunk <- 1 else lefttrunk <- 0
  if (!is.null(start.time)) {
      if (any(abs(start.time)>0)) lefttrunk <- 1  else lefttrunk <- 0;  
  } else lefttrunk <- 0

if (!is.null(margsurv))  {
  if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen")  { 
         resi <- residualsTimereg(margsurv,data=data) 
         RR  <- resi$RR
	 cum <- resi$cumhaz/RR
	 psurvmarg <- exp(-resi$cumhaz); 
         ptrunc <- rep(1,length(psurvmarg)); 
	 if (lefttrunk==1) ptrunc <- exp(-resi$cumhazleft); 
  } 
  else if (class(margsurv)=="coxph") {  
    ### some problems here when data is different from data used in margsurv
       residuals <- residuals(margsurv)
       cum <- cumhaz <- status-residuals
       psurvmarg <- exp(-cumhaz); 
       cumhazleft <- rep(0,antpers)
       ptrunc <- rep(1,length(psurvmarg)); 
       RR<- exp(margsurv$linear.predictors+sum(margsurv$means*coef(margsurv)))
       cum <- cum/RR
        if ((lefttrunk==1)) { 
         stop("Use cox.aalen function for truncation case \n"); 
         baseout <- survival::basehaz(margsurv,centered=FALSE); 
         cumt <- cbind(baseout$time,baseout$hazard)
	 cumt <- Cpred(cumt,start.time)[,2]
	 ptrunc <- exp(-cumt * RR)
	}
  } else if (class(margsurv)=="phreg") {  
	xx <- margsurv$cox.prep
	S0i2 <- S0i <- rep(0,length(xx$strata))
	S0i[xx$jumps+1] <-  1/margsurv$S0
	if (!is.null(margsurv$coef)) 
		RR <- c(exp(margsurv$X  %*% margsurv$coef))
        else RR <- rep(1,antpers) 
	cumhazt <- cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum
###	cumhazt <- cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum
	### back to original order
	orig.o <- (1:nrow(xx$X))[xx$ord+1]
	bto <- order(orig.o)
	cum <- cumhazt <- cumhazt[bto]
###	cumhazt <- Cpred(cum,margsurv$exit)[,2]
	psurvmarg <- c(exp(-cumhazt*RR))
        ptrunc <- rep(1,length(psurvmarg)); 
	start.time <- c(margsurv$entry)
	status <- c(margsurv$status)
	time2 <- margsurv$exit
  } 
} 

if (is.null(clusters) &  class(margsurv)!="phreg") stop("must give clusters")

return(list(psurvmarg=psurvmarg,ptrunc=ptrunc,entry=start.time,exit=time2,
	    status=status,clusters=clusters,cum=cum,RR=RR))

} ## }}}


survival.twostageG <- function(margsurv,data=sys.parent(),score.method="fisher.scoring",Nit=60,detail=0,clusters=NULL,
             silent=1,weights=NULL, control=list(),theta=NULL,theta.des=NULL,
             var.link=1,iid=1,step=0.5,notaylor=0,model="clayton.oakes",
             marginal.trunc=NULL,marginal.survival=NULL,marginal.status=NULL,strata=NULL,
             se.clusters=NULL,numDeriv=0,random.design=NULL,pairs=NULL,pairs.rvs=NULL,numDeriv.method="simple",
             additive.gamma.sum=NULL,var.par=1,cr.models=NULL,case.control=0,ascertained=0,shut.up=0)
{## {{{
## {{{ seting up design and variables
score.iid <- NULL; two.stage <- 1; rate.sim <- 1; sym=1; var.func <- NULL
if (model=="clayton.oakes" || model=="gamma") dep.model <- 1
else if (model=="plackett" || model=="or") dep.model <- 2
else stop("Model must by either clayton.oakes or plackett \n"); 
start.time <- NULL; ptrunc <- NULL; psurvmarg <- NULL; status <- NULL
fix.baseline <- 0; 
convergence.bp <- 1;  ### to control if baseline profiler converges
if ((!is.null(margsurv)) | (!is.null(marginal.survival))) fix.baseline <- 1


if (!is.null(margsurv)) {
if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen") { ## {{{
	 formula<-attr(margsurv,"Formula");
	 beta.fixed <- attr(margsurv,"beta.fixed")
	 if (is.null(beta.fixed)) beta.fixed <- 1; 
	 ldata<-aalen.des(formula,data=data,model="cox.aalen");
	 id <- attr(margsurv,"id"); 
	 mclusters <- attr(margsurv,"cluster.call")
	 X<-ldata$X; 
	 time<-ldata$time2; 
	 Z<-ldata$Z;  
	 status<-ldata$status;
	 time2 <- attr(margsurv,"stop"); 
	 start.time <- attr(margsurv,"start")
	 antpers<-nrow(X);
         if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else { Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
         if (npar==TRUE) {Z<-matrix(0,antpers,1); pz<-1; fixed<-0;} else {fixed<-1;pz<-ncol(Z);}
	 px<-ncol(X);

         if (is.null(clusters) && is.null(mclusters)) stop("No cluster variabel specified in marginal or twostage call\n"); 
         if (is.null(clusters)) clusters <- mclusters 
	 if (nrow(X)!=length(clusters)) stop("Length of Marginal survival data not consistent with cluster length\n"); 
## }}}
	 } else if (class(margsurv)=="phreg") { ## {{{
	        antpers <- length(margsurv$id)
		start.time <- margsurv$entry
} else { ### coxph ## {{{
	  notaylor <- 1
	  antpers <- margsurv$n
	  id <- 0:(antpers-1)
	  mt <- model.frame(margsurv)
	  Y <- model.extract(mt, "response")
	  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
	  if (attr(Y, "type") == "right") {
	      time2 <- Y[, "time"]; 
	      status <- Y[, "status"]
		start.time <- rep(0,antpers);
		} else {
		 start.time <- Y[, 1]; 
		 time2 <- Y[, 2];
		 status <- Y[, 3];
		}
###	   Z <- na.omit(model.matrix(margsurv)[,-1]) ## Discard intercept
	   Z <- matrix(1,antpers,length(coef(margsurv)));

	   if (is.null(clusters)) stop("must give clusters for coxph\n");
	   cluster.call <- clusters; 
	   X <- matrix(1,antpers,1); ### Z <- matrix(0,antpers,1); ### no use for these
	   px <- 1; pz <- ncol(Z); 
	   start <- rep(0,antpers);
	   beta.fixed <- 0
	   semi <- 1
###	   start.time <- 0
} ## }}}
} else { start.time<- 0}

###print("00"); print(summary(psurvmarg)); print(summary(ptrunc))

  if (any(abs(start.time)>0)) lefttrunk <- 1 else lefttrunk <- 0
  if (!is.null(start.time)) {
      if (any(abs(start.time)>0)) lefttrunk <- 1  else lefttrunk <- 0;  
  } else lefttrunk <- 0

if (!is.null(margsurv))  {
  if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen")  { ## {{{
         resi <- residualsTimereg(margsurv,data=data) 
         RR  <- resi$RR
	 psurvmarg <- exp(-resi$cumhaz); 
         ptrunc <- rep(1,length(psurvmarg)); 
	 if (lefttrunk==1) ptrunc <- exp(-resi$cumhazleft); 
  } ## }}}
  else if (class(margsurv)=="coxph") {  ## {{{
	### some problems here when data is different from data used in margsurv
       notaylor <- 1
       residuals <- residuals(margsurv)
       cumhaz <- status-residuals
       psurvmarg <- exp(-cumhaz); 
       cumhazleft <- rep(0,antpers)
       ptrunc <- rep(1,length(psurvmarg)); 
       RR<- exp(margsurv$linear.predictors-sum(margsurv$means*coef(margsurv)))
        if ((lefttrunk==1)) { 
         stop("Use cox.aalen function for truncation case \n"); 
         baseout <- survival::basehaz(margsurv,centered=FALSE); 
         cum <- cbind(baseout$time,baseout$hazard)
	 cum <- Cpred(cum,start.time)[,2]
	 ptrunc <- exp(-cum * RR)
	}
  }
  else if (class(margsurv)=="phreg") {  ## {{{
	xx <- margsurv$cox.prep
	S0i2 <- S0i <- rep(0,length(xx$strata))
	S0i[xx$jumps+1] <-  1/margsurv$S0
	if (!is.null(margsurv$coef)) 
		rr <- c(exp(margsurv$X  %*% margsurv$coef))
        else rr <- rep(1,antpers) 
###	### back to original order
	cumhazt <- cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum
	orig.o <- (1:nrow(xx$X))[xx$ord+1]
	bto <- order(orig.o)
	cumhazt <- cumhazt[bto]
	psurvmarg <- c(exp(-cumhazt*rr))
        ptrunc <- rep(1,length(psurvmarg)); 
	start.time <- c(margsurv$entry)
	status <- c(margsurv$status)
  } ## }}} 
} ## }}}

###     print(head(cbind(psurvmarg,status)))

  antpers <- nrow(data); ## mydim(marginal.survival)[1]
  RR <-  rep(1,antpers); 

  if (is.null(psurvmarg)) psurvmarg <- rep(1,antpers);
  if (!is.null(marginal.survival)) psurvmarg <- marginal.survival 
  if (!is.null(marginal.trunc)) ptrunc <- marginal.trunc 
  if (is.null(ptrunc)) ptrunc <- rep(1,length(psurvmarg))

###  print(summary(psurvmarg)); print(summary(ptrunc))

  if (!is.null(marginal.status)) status <- marginal.status 
  if (is.null(status) & is.null(cr.models)) stop("must give status variable for survival via either margninal model (margsurv), marginal.status or as cr.models \n"); 

  if (is.null(weights)==TRUE) weights <- rep(1,antpers); 
  if (is.null(strata)==TRUE) strata<- rep(1,antpers); 
  if (length(strata)!=antpers) stop("Strata must have length equal to number of data points \n"); 

  ## {{{ cluster set up
  cluster.call <- clusters
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

###  if ((!is.null(max.clust))) if (max.clust< antiid) {
###        coarse.clust <- TRUE
###	qq <- unique(quantile(se.clusters, probs = seq(0, 1, by = 1/max.clust)))
###	qqc <- cut(se.clusters, breaks = qq, include.lowest = TRUE)    
###	se.clusters <- as.integer(qqc)-1
###	max.clusters <- length(unique(se.clusters))
###	maxclust <- max.clust    
###	antiid <- max.clusters
###  }                                                        
  ## }}} 

###  pxz <- px + pz;

   if (!is.null(random.design)) { ### different parameters for Additive random effects # {{{
     dep.model <- 3

###  if (is.null(random.design)) random.design <- matrix(1,antpers,1); 
     dim.rv <- ncol(random.design); 
     if (is.null(theta.des)) theta.des<-diag(dim.rv);


###     ptheta <- dimpar <- ncol(theta.des); 
###   if (dim(theta.des)[2]!=ncol(random.design)) 
###   stop("nrow(theta.des)!= ncol(random.design),\nspecifies restrictions on paramters, if theta.des not given =diag (free)\n"); 
 } else { random.design <- matrix(0,1,1);  dim.rv <- 1; 
           additive.gamma.sum <- matrix(1,1,1); 
   }

  if (is.null(theta.des)) ptheta<-1; 
  if (is.null(theta.des)) theta.des<-matrix(1,antpers,ptheta); ###  else theta.des<-as.matrix(theta.des); 
###    ptheta<-ncol(theta.des); 
###    if (nrow(theta.des)!=antpers) stop("Theta design does not have correct dim");

  if (length(dim(theta.des))==3) ptheta<-dim(theta.des)[2] else if (length(dim(theta.des))==2) ptheta<-ncol(theta.des)
  if (nrow(theta.des)!=antpers & dep.model!=3 ) stop("Theta design does not have correct dim");

  if (length(dim(theta.des))!=3) theta.des <- as.matrix(theta.des)

  if (is.null(theta)==TRUE) {
         if (var.link==1) theta<- rep(-0.7,ptheta);  
         if (var.link==0) theta<- rep(exp(-0.7),ptheta);   
  }       

  if (length(theta)!=ptheta) {
###	 warning("dimensions of theta.des and theta do not match\n"); 
###         print(theta); 
         theta<-rep(theta[1],ptheta); 
  }
  theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 

  if (maxclust==1) stop("No clusters, maxclust size=1\n"); 

  antpairs <- 1; ### to define 
  if (is.null(additive.gamma.sum)) additive.gamma.sum <- matrix(1,dim.rv,ptheta)

  if (!is.null(pairs)) { pair.structure <- 1;} else  pair.structure <- 0;  

## ppprint 
###  print(c(pair.structure,dep.model,fix.baseline))
###  print(head(theta.des))
###  print(c(case.control,ascertained))

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
###	       print("laver 3-dim design "); 
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
  }# }}}
  
  ## }}}

###  setting up arguments for Aalen baseline profile estimates
 if (fix.baseline==0)  { ## {{{  when baseline is estimated when baseline is estimated 
 if (is.null(cr.models)) stop("give hazard models for different causes, ex cr.models=list(Surv(time,status==1)~+1,Surv(time,status==2)~+1) \n")

      if (case.control==0 & ascertained==0) { ## {{{ 
## {{{ setting up random effects and covariates for marginal modelling
        timestatus <- all.vars(cr.models[[1]])
	times <- data[,timestatus[1]]
	if (is.null(status)) status <- data[,timestatus[2]]
	lstatus <- data[,timestatus[2]]
	### organize increments according to overall jump-times 
	jumps <- lstatus!=0
	dtimes <- times[jumps]
	st <- order(dtimes)
	dtimesst <- dtimes[st]
	dcauses <- lstatus[jumps][st]
	ids <- (1:nrow(data))[jumps][st]

        nc <- 0
	for (i in 1:length(cr.models)) {
              a2 <- aalen.des(as.formula(cr.models[[i]]),data=data)
	      X <- a2$X
	      nc <- nc+ncol(X)
	}
	dBaalen <- matrix(0,length(dtimes),nc)
	xjump <- array(0,c(length(cr.models),nc,length(ids)))

	## first compute marginal aalen models for all causes
	a <- list(); da <- list(); 
	for (i in 1:length(cr.models)) {
	      a[[i]] <-  aalen(as.formula(cr.models[[i]]),data=data,robust=0,weights=weights)
              a2 <- aalen.des(as.formula(cr.models[[i]]),data=data)
	      X <- a2$X
	      da[[i]] <- apply(a[[i]]$cum[,-1,drop=FALSE],2,diff)
	      jumpsi <- (1:length(dtimes))[dcauses==i]
	      if (i==1) fp <- 1 
	      indexc <- fp:(fp+ncol(X)-1)
	      dBaalen[jumpsi,indexc] <- da[[i]]
              xjump[i,indexc,] <- t(X[ids,])
	      fp <- ncol(X)+1
        }

	## }}} 

       	#### organize subject specific random variables and design
        ###  for additive gamma model
	## {{{ 
	dimt <- dim(theta.des[,,1,drop=FALSE])[-3]
	dimr <- dim(random.design[,,,drop=FALSE])
	mtheta.des <- array(0,c(dimt,nrow(data)))
	mrv.des <- array(0,c(dimr[1]/2,dimr[2],nrow(data)))
	nrv.des <- rep(0,nrow(data))
	nrv.des[pairs[,1]] <- pairs.rvs
	nrv.des[pairs[,2]] <- pairs.rvs
	mtheta.des[,,pairs[,1]] <- theta.des
	mtheta.des[,,pairs[,2]] <- theta.des
	mrv.des[,,pairs[,1]] <- random.design[1:(dimr[1]/2),,,drop=FALSE]
	mrv.des[,,pairs[,2]] <- random.design[(dimr[1]/2+1):dimr[1],,,drop=FALSE]
	### array thetades to jump times (subjects)
	mtheta.des <- mtheta.des[,,ids,drop=FALSE]
	### array randomdes to jump times (subjects)
	mrv.des <- mrv.des[,,ids,drop=FALSE]
	nrv.des <- pairs.rvs[ids]
	## }}} 

###       	#### organize subject specific random variables and design
###        ###  for additive gamma model
###	## {{{ 
###	dimt <- dim(theta.des[,,1])
###	dimr <- dim(random.design[,,])
###	mtheta.des <- array(0,c(dimt,nrow(data)))
###	mrv.des <- array(0,c(dimr[1]/2,dimr[2],nrow(data)))
###	mtheta.des[,,pairs[,1]] <- theta.des
###	mtheta.des[,,pairs[,2]] <- theta.des
###	mrv.des[,,pairs[,1]] <- random.design[1:(dimr[1]/2),,,drop=FALSE]
###	mrv.des[,,pairs[,2]] <- random.design[(dimr[1]/2+1):dimr[1],,,drop=FALSE]
###	nrv.des <- rep(0,nrow(data))
###	nrv.des[pairs[,1]] <- pairs.rvs
###	nrv.des[pairs[,2]] <- pairs.rvs
###	### array thetades to jump times (subjects)
###	mtheta.des <- mtheta.des[,,ids,drop=FALSE]
###	### array randomdes to jump times (subjects)
###	mrv.des <- mrv.des[,,ids,drop=FALSE]
###	nrv.des <- pairs.rvs[ids]
###	## }}} 

      } ## }}} 

      if (case.control==1 || ascertained==1) { ## {{{ 

###      print(dim(data)); print(summary(pairs))
         data1 <-        data[pairs[,1],]
         data.proband <- data[pairs[,2],]
	 weights1 <-     weights[pairs[,1]]
###	print(summary(data.proband)); print(summary(data1))

## {{{ setting up designs for jump times 
        timestatus <- all.vars(cr.models[[1]])
        if (is.null(status)) status <- data[,timestatus[2]]
	alltimes      <- data[,timestatus[1]]
	times         <- data1[,timestatus[1]]
	lstatus       <- data1[,timestatus[2]]
	timescase     <- data.proband[,timestatus[1]]
	lstatuscase   <- data.proband[,timestatus[2]]
	### organize increments according to overall jump-times 
	jumps         <- lstatus!=0
	dtimes        <- times[jumps]
	dtimescase    <- timescase[jumps]
	st            <- order(dtimes)
	dtimesst      <- dtimes[st]
	dtimesstcase  <- dtimescase[st]
	dcauses       <- lstatus[jumps][st]
	dcausescase   <- lstatuscase[jumps][st]
	ids           <- (1:nrow(data1))[jumps][st]
	###
	### delayed entry for case because of ascertained sampling
	### controls are however control probands, and have entry=0 
	entry <- timescase*lstatuscase
	data1$entry <- entry
	cr.models2 <- list()
	if (ascertained==1) {
           for (i in 1:length(cr.models)) {
	      cr.models2[[i]] <- update(cr.models[[i]],as.formula(paste("Surv(entry,",timestatus[1],",",timestatus[2],")~.",sep="")))
	    }
	} else cr.models2 <- cr.models

        nc <- 0
	for (i in 1:length(cr.models)) {
              X <- aalen.des(as.formula(cr.models[[i]]),data=data1)$X
	      nc <- nc+ncol(X)
	}
	dBaalen    <- matrix(0,length(dtimes),nc)
	xjump      <- array(0,c(length(cr.models),nc,length(ids)))
	xjumpcase  <- array(0,c(length(cr.models),nc,length(ids)))

	## first compute marginal aalen models for all causes
	a <- list(); da <- list(); 
        ### starting values for iteration 
        Bit <- Bitcase <- c()
	for (i in 1:length(cr.models)) { ## {{{ 
	      a[[i]]  <-  aalen(as.formula(cr.models2[[i]]),data=data1,robust=0,weights=weights1)
	      da[[i]] <- apply(a[[i]]$cum[,-1,drop=FALSE],2,diff)
	      jumpsi  <- (1:length(dtimes))[dcauses==i]
              X       <- aalen.des(as.formula(cr.models[[i]]),data=data1)$X
              Xcase   <- aalen.des(as.formula(cr.models[[i]]),data=data.proband)$X
	      if (i==1) fp <- 1 
	      indexc  <- fp:(fp+ncol(X)-1)
	      dBaalen[jumpsi,indexc] <- da[[i]]
              xjump[i,indexc,]       <- t(X[ids,])
              xjumpcase[i,indexc,]   <- t(Xcase[ids,])
	      fp <- fp+ncol(X)
	      ### starting values 
	      Bit <- cbind(Bit,Cpred(a[[i]]$cum,dtimesst)[,-1,drop=FALSE])
       } ## }}} 
       Bit.ini <- Bit
       ## }}} 

####  organize subject specific random variables and design
####  already done in basic pairwise setup 
	mtheta.des <- theta.des[,,ids,drop=FALSE]
	### array randomdes to jump times (subjects)
	mrv.des <- random.design[,,ids,drop=FALSE]
	nrv.des <- pairs.rvs[ids]
      } ## }}} 

 }  else {
	 mrv.des <- array(0,c(1,1,1)); mtheta.des <- array(0,c(1,1,1)); margthetades <- array(0,c(1,1,1)); 
	 xjump <- array(0,c(1,1,1)); dBaalen <- matrix(0,1,1); nrv.des <- 3
 } ## }}} 

  loglike <- function(par) 
  { ## {{{

     if (pair.structure==0 | dep.model!=3) Xtheta <- as.matrix(theta.des) %*% matrix(c(par),nrow=ptheta,ncol=1);
     if (pair.structure==1 & dep.model==3) Xtheta <- matrix(0,antpers,1); ## not needed 
     DXtheta <- array(0,c(1,1,1));

      if (var.link==1 & dep.model==3) epar <- c(exp(par)) else epar <- c(par)
      partheta <- epar

      if (var.par==1 & dep.model==3) {
	 ## from variances to 
	 if (is.null(var.func)) {
	    sp <- sum(epar)
	    partheta <- epar/sp^2 
         } else partheta <- epar ## par.func(epar)
      } 


      if (pair.structure==0) {
	      outl<-.Call("twostageloglikeRV", ## {{{ only two stage model for this option
	      icause=status,ipmargsurv=psurvmarg, 
	      itheta=c(partheta),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),ithetades=theta.des,
	      icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
	      ivarlink=var.link,iid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
	      itrunkp=ptrunc,istrata=as.numeric(strata),iseclusters=se.clusters,iantiid=antiid,
	      irvdes=random.design,iags=additive.gamma.sum,iascertained=ascertained,
              PACKAGE="mets") ## }}}
      }
      else { ## pair-structure 
	     ## twostage model,    case.control option,   profile out baseline 
	     ## conditional model, case.control option,   profile out baseline 
         if (fix.baseline==0) ## if baseline is not given 
         {
          cum1 <- cbind(dtimesst,Bit)
          if ( (case.control==1 || ascertained==1) & (convergence.bp==1)) { ## {{{  profiles out baseline under case-control/ascertainment sampling

###	      ## initial values , only one cr.model for survival 
###           Bit <- cbind(Cpred(a[[1]]$cum,dtimesst)[,-1])
             if (detail>1) plot(dtimesst,Bit,type="l",main="Bit")

             if (ncol(Bit)==0) Bit <- Bit.ini
             Bitcase  <- Cpred(cbind(dtimesst,Bit),dtimesstcase)[,-1,drop=FALSE]
             Bitcase <- .Call("MatxCube",Bitcase,dim(xjumpcase),xjumpcase,PACKAGE="mets")$X

             for (j in 1:5) { ## {{{ profile via iteration 
                   cncc <- .Call("BhatAddGamCC",1,dBaalen,dcauses,dim(xjump),xjump,
		                 c(partheta),dim(mtheta.des),mtheta.des,additive.gamma.sum,var.link, 
				 dim(mrv.des),mrv.des,nrv.des,1,Bit,Bitcase,dcausescase,PACKAGE="mets")
		   d <- max(abs(Bit-cncc$B))
		   if (detail>1) print(d)
		   Bit <- cncc$B
###		   if (detail>1) print(c(par,epar,partheta)); 
###		   if (detail>1) print(summary(Bit)); 
		   if (detail>1) print(summary(cncc$caseweights))
		   cum1 <- cbind(dtimesst,cncc$B)
		   Bitcase  <-cbind(Cpred(cum1,dtimesstcase)[,-1])
###		   if (detail>1) print(summary(Bitcase))
		   if (detail>1) lines(dtimesst,Bit,col=j+1); 
		   if (is.na(d)) { 
			  if (shut.up==0) cat("Baseline profiler gives missing values\n");  
		          Bit <- Bit.ini; cum1 <- cbind(dtimesst,Bit); convergence.bp <<- 0; break;
		   }
		   Bitcase <- .Call("MatxCube",Bitcase,dim(xjumpcase),xjumpcase,PACKAGE="mets")$X
		   if (d<0.00001) break; 
           } ## }}} 


	   nulrow    <- rep(0,ncol(Bit)+1)
	   pbases    <- Cpred(rbind(nulrow,cbind(dtimesst,Bit)),alltimes)[,-1,drop=FALSE]
           X         <- aalen.des(as.formula(cr.models[[1]]),data=data)$X
	   psurvmarg <- exp(-apply(X*pbases,1,sum))  ## psurv given baseline 
	   if (ascertained==1) {
              Xcase     <- aalen.des(as.formula(cr.models[[1]]),data=data.proband)$X
              X         <- aalen.des(as.formula(cr.models[[1]]),data=data1)$X
	      pba.case  <- Cpred(rbind(nulrow,cbind(dtimesst,Bit)),entry)[,-1,drop=FALSE]
	      ptrunc    <- rep(0,nrow(data))
	      ### for control probands ptrunc=1, thus no adjustment
	      ptrunc[pairs[,1]]     <- exp(-apply(X*    pba.case,1,sum)*lstatuscase) ## delayed entry at time of ascertainment proband 
 	      ptrunc[pairs[,2]]     <- exp(-apply(Xcase*pba.case,1,sum)*lstatuscase)
	   }
###	   print(head(cbind(psurvmarg,ptrunc)))
###	   print(summary(psurvmarg))
###	   print(summary(ptrunc))
###	   print(dim(pbases)); 

          } ## }}} 
          }

###      browser()
###      print(dim(random.design))
###      print(summary(status)); print(summary(psurvmarg)); print(summary(clusters)); print(summary(random.design));
###      print(random.design)
###      print(theta.des)

          outl<-.Call("twostageloglikeRVpairs", ## {{{
          icause=status,ipmargsurv=psurvmarg, 
          itheta=c(partheta),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),
	  ithetades=theta.des,
          icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
          ivarlink=var.link,iiid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
          itrunkp=ptrunc,istrata=as.numeric(strata),iseclusters=se.clusters,iantiid=antiid,
          irvdes=random.design,
          idimthetades=dim(theta.des),idimrvdes=dim(random.design),irvs=pairs.rvs,iags=additive.gamma.sum, 
	  iascertained=ascertained,PACKAGE="mets")
	  ## }}} 


          if (fix.baseline==0)  {
              outl$baseline <- cum1; 
	      outl$marginal.surv <- psurvmarg; 
	      outl$marginal.trunc <- ptrunc
	  }

      } ## }}} 

    if (detail==3) print(c(partheta,outl$loglike))

    ## variance parametrization, and inverse.link 
    if (dep.model==3) {# {{{
    if (var.par==1) {
         ## from variances to and with sum for all random effects 
         if (is.null(var.func)) {
	 if (var.link==0)  {
###		 print(c(sp,epar))
             mm <- matrix(-epar*2*sp,length(epar),length(epar))
	     diag(mm) <- sp^2-epar*2*sp
###	    print(mm)
	 } else {
            mm <- -c(epar) %o% c(epar)*2*sp
            diag(mm) <- epar*sp^2-epar^2*2*sp
###	    print(mm)
	 }
	    mm <- mm/sp^4
	 } else mm  <- numDeriv::hessian(var.func,par)
      } else {
	   if (var.link==0) mm <- diag(length(epar)) else 
			  mm <- diag(length(c(epar)))*c(epar)
      }
      }# }}}

###    print(c(var.link,dep.model,var.par))
###    print("hh"); print(mm); print(outl$score)

    if (dep.model==3) {# {{{
       outl$score <-  t(mm) %*% outl$score
       outl$Dscore <- t(mm) %*% outl$Dscore %*% mm 
       if (iid==1) { outl$theta.iid <- t(t(mm) %*% t(outl$theta.iid))
               if (class(margsurv)=="phreg") {  
###	    print(dim(mm))
	          outl$D1thetal  <- t(t(mm) %*% t(outl$D1thetal))
	          outl$D2thetal  <- t(t(mm) %*% t(outl$D2thetal))
	       }
       } 

###       print(crossprod(outl$theta.iid)); print(outl$Dscore)
###       print(c(outl$score))
###       print(apply(outl$theta.iid,2,sum))
    }# }}}

    attr(outl,"gradient") <-outl$score 
    if (oout==0) ret <- c(-1*outl$loglike) else if (oout==1) ret <- sum(outl$score^2) else if (oout==2) ret <- outl else ret <- outl$score
    return(ret)
  } ## }}}

  if (score.method=="optimize" && ptheta!=1) {
     cat("optimize only works for d==1, score.mehod set to nlminb \n"); 
     score.method <- "nlminb";
  } 

  score1 <- NULL
  theta.iid <- NULL
  logl <- NULL
  p <- theta

    if (score.method=="fisher.scoring") { ## {{{
        oout <- 2;  ### output control for obj
        if (Nit>0) 
            for (i in 1:Nit)
            {
                    out <- loglike(p)
		    ## updating starting values for cumulative baselines
		    if (fix.baseline==0) Bit <- out$baseline[,-1,drop=FALSE]
		    if (fix.baseline==1) hess <-  out$Dscore
                    ###  uses simple second derivative for computing derivative of score
		    if (numDeriv==2 || (((fix.baseline==0)) & (i==1))) {
			    oout <- 3
			    hess <- numDeriv::jacobian(loglike,p,method="simple")
			    oout <- 2
		    }
                    if (!is.na(sum(hess))) hessi <- lava::Inverse(hess) else hessi <- hess 
                    if (detail==1) {## {{{
                        cat(paste("Fisher-Scoring ===================: it=",i,"\n")); 
                        cat("theta:");print(c(theta))
                        cat("loglike:");cat(c(out$loglike),"\n"); 
                        cat("score:");cat(c(out$score),"\n"); 
                        cat("hess:\n"); cat(hess,"\n"); 
                    }## }}}
                    delta <- step*( hessi %*% out$score )
		    ### update p, but note that score and derivative in fact related to previous p 
		    ### unless Nit=0, 
	            if (Nit>0) {
                    p <- p - delta
                    theta <- p; 
		    }
                    if (is.nan(sum(out$score))) break; 
                    if (sum(abs(out$score))<0.00001) break; 
                    if (max(theta)>20 & var.link==1) { cat("theta too large lacking convergence \n"); break; }
           }
        if (!is.nan(sum(p))) { 
            if (detail==1 && iid==1) cat("iid decomposition\n"); 
            out <- loglike(p) 
            logl <- out$loglike
            score1 <- score <- out$score
            oout <- 0; 
            hess1 <- hess  <- -1*out$Dscore 

           if (iid==1) {  ## {{{
		theta.iid <- out$theta.iid
	        score.iid <-   theta.iid

                if (class(margsurv)=="phreg") {  ## {{{

                  if (is.null(margsurv$opt) | is.null(margsurv$coef)) fixbeta<- 1 else fixbeta <- 0
                  id  <-  xx$id+1

		  xx <- margsurv$cox.prep
		  D1thetal<- out$D1thetal
		  D2thetal<- out$D2thetal
		  Dlamthetal <- (D1thetal+D2thetal)
###		  print(mydim(Dlamthetal))

                  Fbeta <- - t(D1thetal+D2thetal) %*% (margsurv$X *cumhazt*psurvmarg*rr)

                  ### time-ordered 
                  Gbase <- apply(-Dlamthetal[xx$ord+1,,drop=FALSE]*psurvmarg[xx$ord+1]*rr[xx$ord+1],2,revcumsumstrata,xx$strata,xx$nstrata)

###		  print(c(Gbase)) print(Fbeta)

                  if (fixbeta==0) {
			  Z <- xx$X
			  U <- E <- matrix(0,nrow(xx$X),margsurv$p)
			  E[xx$jumps+1,] <- margsurv$E
			  U[xx$jumps+1,] <- margsurv$U
			  invhess <- -solve(margsurv$hessian)
			  S0i <- rep(0,length(xx$strata))
			  S0i[xx$jumps+1] <- 1/margsurv$S0
			  cumhaz <- c(cumsumstrata(S0i,xx$strata,xx$nstrata))
			  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
			  rr <- c(xx$sign*exp(Z %*% coef(margsurv) + xx$offset))
			  ### Martingale  as a function of time and for all subjects to handle strata 
			  MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0)*rr*c(xx$weights)
			  UUbeta <- apply(MGt,2,sumstrata,id-1,max(id))
			  UUbeta <- UUbeta %*% invhess
			  GbaseEdLam0 <- t(Gbase) %*% (E*S0i) 
			  Fbeta   <-  Fbeta -  GbaseEdLam0
			  Fbetaiid <- UUbeta %*% t(Fbeta)
	          }

		  ###  \int_0^\tau (GBase(s) / S_0(s)) dN_i(s) - dLamba_i(s)
		  xx <- margsurv$cox.prep
		  S0i <- rep(0,length(xx$strata))
		  S0i[xx$jumps+1] <- 1/margsurv$S0
		  S0i2[xx$jumps+1] <- 1/margsurv$S0^2

		  GbasedN <- Gbase*S0i
		  GbasedLam0 <- apply(Gbase*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)

		  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(margsurv) + xx$offset)) else rr <- c(xx$sign*exp( xx$offset))
		  MGt <- (GbasedN[,drop=FALSE]-GbasedLam0*rr)*c(xx$weights)
		  MGti <- apply(MGt,2,sumstrata,id-1,max(id))

	###	  print(cbind(theta.iid,MGti,Fbetaiid))

		  theta.iid <-  theta.iid + MGti
		  if (fixbeta==0) theta.iid  <-  theta.iid + Fbetaiid
		  out$theta.iid <- theta.iid
	  } ## }}}
   } ## }}} 


            if (detail==1 && iid==1) cat("finished iid decomposition\n"); 
	    ### for profile solutions update second derivative at final 
	    if (numDeriv==2 || ((fix.baseline==0))) {
                    oout <- 3
                    hess <- numDeriv::jacobian(loglike,p,method="simple")
		    oout <- 2
		    }
        }
        if (numDeriv>=1) {
            if (detail==1 ) cat("starting numDeriv for second derivative \n"); 
            oout <- 0; 
            score2 <- numDeriv::jacobian(loglike,p)
	    score1 <- matrix(score2,ncol=1)
            oout <- 3
            hess <- numDeriv::jacobian(loglike,p,method="simple")
            if (detail==1 ) cat("finished numDeriv for second derivative \n"); 
        }
        if (detail==1 & Nit==0) {## {{{
            cat(paste("Fisher-Scoring ===================: final","\n")); 
            cat("theta:");print(c(p))
            cat("loglike:");cat(c(out$loglike),"\n"); 
            cat("score:");cat(c(out$score),"\n"); 
            cat("hess:\n"); cat(hess,"\n"); 
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
    hess1 <- hess <- -1* out$Dscore
    if (iid==1) { theta.iid <- out$theta.iid; score.iid <- out$theta.iid }
    if (numDeriv==1) {
    if (detail==1 ) cat("numDeriv hessian start\n"); 
      oout <- 3; ## returns score 
      hess <- numDeriv::jacobian(loglike,opt$par)
    if (detail==1 ) cat("numDeriv hessian done\n"); 
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
    hess1 <- hess <- -1* out$Dscore
    if (numDeriv==1) {
    if (detail==1 ) cat("numDeriv hessian start\n"); 
      oout <- 3;  ## to get jacobian
      hess <- numDeriv::jacobian(loglike,theta)
    if (detail==1 ) cat("numDeriv hessian done\n"); 
    }
    hessi <- lava::Inverse(hess); 
    if (iid==1) { theta.iid <- out$theta.iid; score.iid <- out$theta.iid }
  ## }}}
  } else if (score.method=="nlm") { ## {{{ nlm optimizer
    iid <- 0; oout <- 0; 
    tryCatch(opt <- nlm(loglike,theta,hessian=TRUE,print.level=detail),error=function(x) NA)
    iid <- 1; 
    hess <- opt$hessian
    score <- opt$gradient
    if (detail==1) print(opt); 
    hessi <- lava::Inverse(hess); 
    theta <- opt$estimate
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    oout <- 2
    out <- loglike(opt$estimate)
    logl <- out$loglike
    score1 <- out$score
    hess1 <- out$Dscore
    if (iid==1) { theta.iid <- out$theta.iid; score.iid <- out$theta.iid }
  ## }}}
  }  else stop("score.methods = optimize(dim=1) nlm nlminb fisher.scoring\n"); 

## {{{ handling output
  loglikeiid <- NULL
  robvar.theta <- NULL
  likepairs <- NULL
  if (fix.baseline==1)  { marginal.surv <- psurvmarg; marginal.trunc <- ptrunc; } else { marginal.surv <- out$marginal.surv; marginal.trunc <- out$marginal.trunc;} 

  if (iid==1) {
  if (dep.model==3 & pair.structure==1) likepairs <- out$likepairs
  if (dep.model==3 & two.stage==0) {# {{{
	  hessi <- -1*hessi
	  all.likepairs <- out$all.likepairs
	  colnames(all.likepairs) <- c("surv","dt","ds","dtds","cause1","cause2")
  }# }}}
###  print(crossprod(out$theta.iid) %*% hessi)
     theta.iid <- out$theta.iid %*% hessi
     if (is.null(call.secluster)) rownames(theta.iid) <- unique(cluster.call) else rownames(theta.iid) <- unique(se.clusters)
     robvar.theta  <- crossprod(theta.iid) 
     loglikeiid <- out$loglikeiid
    } else { all.likepairs <- NULL}
  var.theta <- robvar.theta
  if (is.null(robvar.theta)) var.theta <- hessi

  if (!is.null(colnames(theta.des))) thetanames <- colnames(theta.des) else thetanames <- paste("dependence",1:length(theta),sep="")
  theta <- matrix(theta,length(c(theta)),1)
  if (length(thetanames)==nrow(theta)) { rownames(theta) <- thetanames; rownames(var.theta) <- colnames(var.theta) <- thetanames; }
  if (!is.null(logl)) logl <- -1*logl

  if (convergence.bp==0) theta <- rep(NA,length(theta))
  ud <- list(theta=theta,score=score,hess=hess,hessi=hessi,var.theta=var.theta,
	     model=model,robvar.theta=robvar.theta,
             theta.iid=theta.iid,loglikeiid=loglikeiid,likepairs=likepairs,
	     thetanames=thetanames,loglike=logl,score1=score1,Dscore=out$Dscore,
	     marginal.surv=marginal.surv,marginal.trunc=marginal.trunc,
	     baseline=out$baseline,se=diag(robvar.theta)^.5,
	     score.iid=score.iid)
  class(ud) <- "mets.twostage" 
  attr(ud,"response") <- "survival"
  attr(ud,"Formula") <- formula
  attr(ud,"clusters") <- clusters
  attr(ud,"cluster.call") <- cluster.call
  attr(ud,"secluster") <- c(se.clusters)
  attr(ud,"sym")<-sym; 
  attr(ud,"var.link")<-var.link; 
  attr(ud,"var.par")<-var.par; 
  attr(ud,"var.func")<-var.func; 
  attr(ud,"ptheta")<-ptheta
  attr(ud,"antpers")<-antpers; 
  attr(ud,"antclust")<-antclust; 
  attr(ud,"Type") <- model
  attr(ud,"twostage") <- two.stage
  attr(ud,"additive-gamma") <- (dep.model==3)*1
  if (!is.null(marginal.trunc)) attr(ud,"trunclikeiid")<- out$trunclikeiid
  if (dep.model==3 & two.stage==0) attr(ud,"all.likepairs")<- all.likepairs
  if (dep.model==3 ) attr(ud,"additive.gamma.sum") <- additive.gamma.sum
  #likepairs=likepairs,##  if (dep.model==3 & pair.structure==1) attr(ud, "likepairs") <- c(out$likepairs)
  if (dep.model==3 & pair.structure==0) attr(ud, "pardes") <-    theta.des
  if (dep.model==3 & pair.structure==0) attr(ud, "theta.des") <- theta.des
  if (dep.model==3 & pair.structure==1) attr(ud, "pardes") <-    theta.des[,,1]
  if (dep.model==3 & pair.structure==1) attr(ud, "theta.des") <- theta.des[,,1]
  if (dep.model==3 & pair.structure==0) attr(ud, "rv1") <- random.design[1,]
  if (dep.model==3 & pair.structure==1) attr(ud, "rv1") <- random.design[,,1]
  return(ud);
  ## }}}

} ## }}}


##' @title Twostage survival model fitted by pseudo MLE 
##'
##' @description 
##' Fits Clayton-Oakes clustered  survival data 
##' using marginals that are on Cox form in the likelihood for the dependence parameter 
##' as in Glidden (2000). The dependence can be modelled via  a 
##' \enumerate{
##' \item  Regression design on dependence parameter. 
##' }
##'
##' We allow a regression structure for the indenpendent gamma distributed 
##' random effects  and their variances that may depend on cluster covariates. So
##' \deqn{
##'  \theta = h( z_j^T \alpha)
##' }
##' where \eqn{z} is specified by theta.des . The link function can be the exp when var.link=1
##' @references
##' 
##' Measuring early or late dependence for bivariate twin data
##' Scheike, Holst, Hjelmborg (2015), LIDA  
##' 
##' Twostage modelling of additive gamma frailty models for survival data. 
##' Scheike and Holst, working paper 
##' 
##' Shih and Louis (1995) Inference on the association parameter in copula models for bivariate
##' survival data, Biometrics, (1995).
##' 
##' Glidden (2000), A Two-Stage estimator of the dependence
##' parameter for the Clayton Oakes model, LIDA, (2000).
##' 
##' @examples
##' data(diabetes)
##' dd <- phreg(Surv(time,status==1)~treat+cluster(id),diabetes)
##' oo <- twostageMLE(dd,data=diabetes)
##' summary(oo)
##' 
##' theta.des <- model.matrix(~-1+factor(adult),diabetes)
##' 
##' oo <-twostageMLE(dd,data=diabetes,theta.des=theta.des)
##' summary(oo)
##' @keywords survival
##' @author Thomas Scheike
##' @param margsurv Marginal model from phreg 
##' @param data data frame
##' @param theta Starting values for variance components
##' @param theta.des design for dependence parameters, when pairs are given this is could be a
##' (pairs) x (numer of parameters)  x (max number random effects) matrix
##' @param var.link Link function for variance  if 1 then uses exp link 
##' @param method type of opitmizer, default is Newton-Raphson "NR"
##' @param no.opt to not optimize, for example to get score and iid for specific theta 
##' @param weights cluster specific weights, but given with length equivalent to data-set, weights for score equations 
##' @param ... arguments to be passed to  optimizer
##' @export 
twostageMLE <-function(margsurv,data=sys.parent(),
theta=NULL,theta.des=NULL,var.link=0,method="NR",no.opt=FALSE,weights=NULL,...)
{# {{{
 if (class(margsurv)!="phreg")  stop("Must use phreg for this \n"); 

 clusters <- margsurv$cox.prep$id
 n <- nrow(margsurv$cox.prep$X)
 secluster <- NULL

 if (is.null(theta.des)==TRUE) ptheta<-1; 
 if (is.null(theta.des)==TRUE) theta.des<-matrix(1,n,ptheta) else theta.des<-as.matrix(theta.des); 
  ptheta<-ncol(theta.des); 
  if (nrow(theta.des)!=n) stop("Theta design does not have correct dim");

  if (is.null(theta)==TRUE) {
      if (var.link==1) theta<- rep(0.0,ptheta); 
      if (var.link==0) theta<- rep(1.0,ptheta); 
  }
  if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta); 
  theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 

  max.clust <- length(unique(clusters))
  theta.iid <- matrix(0,max.clust,ptheta)

  xx <- margsurv$cox.prep
  nn <- length(xx$strata)

  if (is.null(weights)) weights <-  rep(1,nn)
  if (length(weights)!=nn) stop("Weights do not have right length")

  statusxx <- rep(0,length(xx$strata))
  statusxx[xx$jumps+1] <- 1
  xx$status <- statusxx
  mid <- max(xx$id)+1

  ## Ni.(t-), Ni.(t)
  Nsum <- cumsumstratasum(statusxx,xx$id,mid,type="all")
  ## Ni.(tau)
  Ni.tau <- sumstrata(statusxx,xx$id,mid)
  ## cumahz(T_i-)
  S0i2 <- S0i <- rep(0, length(xx$strata))
  S0i[xx$jumps + 1] <- 1/margsurv$S0
  cumhazD <- cumsumstratasum(S0i, xx$strata, xx$nstrata)$lagsum
  if (!is.null(margsurv$coef)) RR <- exp(xx$X %*% margsurv$coef) else RR  <-  rep(1,nn)
  H <- c(cumhazD * RR)

  cc <- cluster.index(xx$id)
  firstid <- cc$firstclustid+1
  if (max(cc$cluster.size)==1) stop("No clusters !, maxclust size=1\n"); 

  ### order after time-sorted data
  theta.des <- theta.des[xx$ord+1,,drop=FALSE]
  weightsid <- weights <- weights[xx$ord+1]
  ### clusterspecific weights
  weights <- weights[firstid]

###  par <- c(2,2,0.1) 
###  par <- 2
obj <- function(par,all=FALSE)
{# {{{

if (var.link==1) epar <- c(exp(par)) else epar <- c(par)

thetaX<- as.matrix(theta.des[firstid,,drop=FALSE])
thetav <- c(as.matrix(theta.des) %*%  c(epar))
thetai <-thetav[firstid]

Hs <- sumstrata(H*exp(thetav*H),xx$id,mid)
R <- sumstrata((exp(thetav*H)-1),xx$id,mid) + 1
H2 <- sumstrata(H^2*exp(thetav*H),xx$id,mid) 

l1 <- sumstrata(log(1+thetav*Nsum$lagsum)*statusxx,xx$id,mid)
l2 <- sumstrata(statusxx*H,xx$id,mid)
l3 <- -((thetai)^{-1}+Ni.tau) * log(R)
logliid <- (l1 + thetai* l2 + l3)*c(weights)
###cbind(logliid,fit2$loglikeiid,Ni.tau)
logl <- sum(logliid)
ploglik <- logl

## first derivative 
l1s <- sumstrata(Nsum$lagsum/(1+thetav*Nsum$lagsum)*statusxx,xx$id,mid)
l2s <- (thetai^{-2}) * log(R) 
l3s <- -(thetai^{-1}+Ni.tau) * Hs / R
Dltheta <- (l1s+l2s+l3s+l2)*c(weights)
scoreiid <- thetaX* c(Dltheta)

## second derivative 
D2N <- -sumstrata(Nsum$lagsum^2/(1+thetav*Nsum$lagsum)^2*statusxx,xx$id,mid)
Dhes  <- c(D2N+(2/thetai^2)*Hs/R-(2/thetai^3)*log(R)-(1/thetai+Ni.tau)*(H2*R-Hs*Hs)/R^2)
Dhes  <-  Dhes * c(weights)

if (var.link==1) { 
	scoreiid  <-  scoreiid * c(exp(par));
	Dhes<- Dhes* thetai^2 + thetai *  Dltheta
} 

gradient <- apply(scoreiid,2,sum)
hessian <-  crossprod(thetaX,thetaX*c(Dhes))
hess2 <-    crossprod(scoreiid)

val <- list(id=xx$id,score.iid=scoreiid,logl.iid=logliid,ploglik=ploglik,
            gradient=-gradient,hessian=hessian,hess2=hess2)

if (all) return(val); 

with(val, structure(-ploglik,gradient=-gradient,hessian=hessian))
}# }}}


  opt <- NULL
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
          opt <- lava::NR(theta,obj,...)
###          opt <- lava::NR(theta,obj)
          opt$estimate <- opt$par
      } else {
          opt <- nlm(obj,theta,...)
          opt$method <- "nlm"
      }
      cc <- opt$estimate;  
###   names(cc) <- colnames(theta.des)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
  } else val <- c(list(coef=theta),obj(theta,all=TRUE))
  val$score <- val$gradient


  theta <- matrix(c(val$coef),length(c(val$coef)),1)
  if (!is.null(colnames(theta.des))) thetanames <- colnames(theta.des) else thetanames <- paste("dependence",1:length(c(theta)),sep="")
  if (length(thetanames)==length(c(theta))) { rownames(theta) <- thetanames; rownames(var.theta) <- colnames(var.theta) <- thetanames; }

  hessianI <- solve(val$hessian)
  val$theta.iid.naive  <-  val$score.iid %*% hessianI

  ### iid due to Marginal model 
  biid <- 1
  if (biid==1) {
# {{{
  if (is.null(margsurv$opt) | is.null(margsurv$coef)) fixbeta<- 1 else fixbeta <- 0

   id  <-  xx$id+1
   if (var.link==1) epar <- c(exp(theta)) else epar <- c(theta)

   thetaX<- as.matrix(theta.des[firstid,,drop=FALSE])
   thetav <- c(as.matrix(theta.des) %*%  c(epar))

   Hs <- c(sumstrata(H*exp(thetav*H),xx$id,mid))
   R <- c(sumstrata((exp(thetav*H)-1),xx$id,mid)) + 1
   H2 <- c(sumstrata(H^2*exp(thetav*H),xx$id,mid) )

Ft <- ((1/(thetav*R[id]))*exp(thetav*H)-(1/thetav+Ni.tau[id])*(1+thetav*H)*exp(thetav*H)/R[id]+statusxx+(1+thetav*Ni.tau[id])*exp(thetav*H)*Hs[id]/R[id]^2)

   if (var.link==1) { 
	Ft  <-  Ft * c(exp(par));
   } 
   Gt <- c(RR *Ft * xx$sign * weightsid)
   Ft  <- c( Ft * H * weightsid )
###  oi <- order(id)
###print(round(cbind(id,Ft,thetav,H,R[id],Hs[id],statusxx,Ni.tau[id],xx$X)[oi,],4))
###print(round(cbind(id,Ft,xx$X,Ft*xx$X)[oi,],4))

###   pit <- revcumsumstrata(Ft* RR* theta.des,xx$strata,xx$nstrata)
  Gbase <- apply(Gt* theta.des,2,revcumsumstrata,xx$strata,xx$nstrata)
  Fbeta <-  t(theta.des) %*% ( xx$X * Ft )

###  print(c(Gbase)); print(Fbeta)

  ### beta part 
  if (fixbeta==0) {
	  Z <- xx$X
	  U <- E <- matrix(0,nrow(xx$X),margsurv$p)
	  E[xx$jumps+1,] <- margsurv$E
	  U[xx$jumps+1,] <- margsurv$U
          invhess <- -solve(margsurv$hessian)
	  S0i <- rep(0,length(xx$strata))
	  S0i[xx$jumps+1] <- 1/margsurv$S0
	  cumhaz <- c(cumsumstrata(S0i,xx$strata,xx$nstrata))
	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
	  rr <- c(xx$sign*exp(Z %*% coef(margsurv) + xx$offset))
	  ### Martingale  as a function of time and for all subjects to handle strata 
	  MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0)*rr*c(xx$weights)
          UUbeta <- apply(MGt,2,sumstrata,id-1,max(id))
	  UUbeta <- UUbeta %*% invhess
	  GbaseEdLam0 <- t(Gbase) %*% (E*S0i) 
	  Fbeta   <-  Fbeta -  GbaseEdLam0
          Fbetaiid <- UUbeta %*% t(Fbeta)
  }

  ###  \int_0^\tau (GBase(s) / S_0(s)) dN_i(s) - dLamba_i(s)
  xx <- margsurv$cox.prep
  S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <- 1/margsurv$S0
  S0i2[xx$jumps+1] <- 1/margsurv$S0^2

  GbasedN <- Gbase*S0i
  GbasedLam0 <- apply(Gbase*S0i2,2,cumsumstrata,xx$strata,xx$nstrata)

  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(margsurv) + xx$offset)) else rr <- c(xx$sign*exp( xx$offset))
  MGt <- (GbasedN[,drop=FALSE]-GbasedLam0*rr)*c(xx$weights)
  MGti <- apply(MGt,2,sumstrata,id-1,max(id))

###  if (fixbeta==0) print(cbind(val$score.iid,MGti,Fbetaiid)) else print(cbind(val$score.iid,MGti))

  theta.iid <-  val$score.iid + MGti
  if (fixbeta==0) theta.iid  <-  theta.iid + Fbetaiid
  theta.iid  <-  theta.iid %*% hessianI 
  ## take names from phreg call 
  if (!is.null(margsurv$id)) rownames(theta.iid) <- unique(margsurv$id) 
  val$theta.iid <- theta.iid
# }}}
  }

  var <- robvar.theta  <-  var.theta  <-  crossprod(val$theta.iid)
  naive.var  <-   crossprod(val$theta.iid.naive)

  val  <- c(val,list(theta=theta,var.theta=var.theta,robvar.theta=robvar.theta,
		     var=var,thetanames=thetanames,model="clayton.oakes",
	             se=diag(robvar.theta)^.5),
	             var.naive=naive.var)
  class(val) <- "mets.twostage" 
  attr(val,"clusters") <- clusters
###  attr(val,"secluster") <- c(se.clusters)
  attr(val,"var.link")<-var.link; 
  attr(val,"ptheta")<-ptheta
  attr(val,"n")<-n ; 
  attr(val,"response")<- "survival"
  attr(val,"additive-gamma")<-0 
  attr(val,"twostage") <- "two.stage"

  return(val)  
}# }}}

##' @title Survival model for multivariate survival data 
##'
##' @description 
##' Fits additive gamma frailty model 
##' with additive hazard condtional on the random effects
##' \deqn{
##' \lambda_{ij} = (V_{ij}^T Z) (X_{ij}^T \alpha(t)) 
##' }
##' The baseline \eqn{\alpha(t)} is profiled out using
##' marginal modelling adjusted for the random effects structure as in Eriksson and Scheike (2015).
##' One advantage of the standard frailty model is that one can deal with competing risks 
##' for this model. 
##'
##' For all models the 
##' standard errors do not reflect this uncertainty of the baseline estimates, and might therefore be a bit to small.
##' To remedy this one can do bootstrapping or use survival.twostage.fullse function when possible.
##'
##' If clusters contain more than two times, we use a composite likelihood
##' based on the pairwise bivariate models. Can also fit a additive gamma random
##' effects model described in detail below.
##'
##' We allow a regression structure for the indenpendent gamma distributed 
##' random effects  and their variances that may depend on cluster covariates. So
##' \deqn{
##'  \theta = z_j^T \alpha
##' }
##' where \eqn{z} is specified by theta.des 
##' The reported standard errors are based on the estimated information from the 
##' likelihood assuming that the marginals are known. 
##'
##' Can also fit a structured additive gamma random effects model, such
##' as the ACE, ADE model for survival data. 
##'
##' Now random.design specificies the random effects for each subject within a cluster. This is
##' a matrix of 1's and 0's with dimension n x d.  With d random effects. 
##' For a cluster with two subjects, we let the random.design rows be 
##'  \eqn{v_1} and \eqn{v_2}. 
##' Such that the random effects for subject 
##' 1 is \deqn{v_1^T (Z_1,...,Z_d)}, for d random effects. Each random effect
##' has an associated parameter \eqn{(\lambda_1,...,\lambda_d)}. 
##' By construction subjects 1's random effect are Gamma distributed with 
##' mean \eqn{\lambda_j/v_1^T \lambda}
##' and variance \eqn{\lambda_j/(v_1^T \lambda)^2}. Note that the random effect 
##' \eqn{v_1^T (Z_1,...,Z_d)} has mean 1 and variance \eqn{1/(v_1^T \lambda)}.
##' It is here asssumed that  \eqn{lamtot=v_1^T \lambda} is fixed over all clusters
##' as it would be for the ACE model below.
##' The lamtot parameter may be specified separately for some sets of the parameter
##' is the additive.gamma.sum (ags) matrix is specified and then lamtot for the 
##' j'th random effect is \eqn{ags_j^T \lambda}.
##'
##' Based on these parameters the relative contribution (the heritability, h) is 
##' equivalent to  the expected values of the random effects  \eqn{\lambda_j/v_1^T \lambda}
##'
##' The DEFAULT parametrization uses the variances of the random effecs 
##' \deqn{
##' \theta_j  = \lambda_j/(v_1^T \lambda)^2
##' }
##' For alternative parametrizations one can specify how the parameters relate to \eqn{\lambda_j}
##' with the function 
##'
##' Given the random effects the survival distributions with a cluster are independent and
##' on the form 
##' \deqn{
##' P(T > t| x,z) = exp( -Z  A(t) \exp( Z^t beta))
##' }
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
##' The case.control option that can be used with the pair specification of the pairwise parts
##' of the estimating equations. Here it is assumed that the second subject of each pair is the
##' proband. 
##' @keywords survival
##' @author Thomas Scheike
##' @param margsurv Marginal model 
##' @param data data frame
##' @param score.method Scoring method "fisher.scoring", "nlminb", "optimize", "nlm"
##' @param Nit Number of iterations
##' @param detail Detail
##' @param clusters Cluster variable
##' @param silent Debug information
##' @param weights Weights
##' @param control Optimization arguments
##' @param theta Starting values for variance components
##' @param theta.des design for dependence parameters, when pairs are given this is could be a
##' (pairs) x (numer of parameters)  x (max number random effects) matrix
##' @param var.link Link function for variance 
##' @param iid Calculate i.i.d. decomposition
##' @param step Step size
##' @param model model
##' @param marginal.trunc marginal left truncation probabilities
##' @param marginal.survival optional vector of marginal survival probabilities 
##' @param marginal.status related to marginal survival probabilities 
##' @param strata strata for fitting, see example
##' @param se.clusters for clusters for se calculation with iid
##' @param max.clust max se.clusters for se calculation with iid
##' @param numDeriv to get numDeriv version of second derivative, otherwise uses sum of squared score 
##' @param random.design random effect design for additive gamma model, when pairs are given this is
##' a (pairs) x (2) x (max number random effects) matrix, see pairs.rvs below
##' @param pairs matrix with rows of indeces (two-columns) for the pairs considered in the pairwise
##' composite score, useful for case-control sampling when marginal is known.
##' @param pairs.rvs for additive gamma model and random.design and theta.des are given as arrays,
##' this specifice number of random effects for each pair. 
##' @param numDeriv.method uses simple to speed up things and second derivative not so important.
##' @param additive.gamma.sum for two.stage=0, this is specification of the lamtot in the models via
##' a matrix that is multiplied onto the parameters theta (dimensions=(number random effects x number
##' of theta parameters), when null then sums all parameters.
##' @param var.par is 1 for the default parametrization with the variances of the random effects,
##' var.par=0 specifies that the \eqn{\lambda_j}'s are used as parameters. 
##' @param cr.models competing risks models for two.stage=0, should be given as a list with models for each cause
##' @param case.control assumes case control structure for "pairs" with second column being the probands,
##' when this options is used the twostage model is profiled out via the paired estimating equations for the
##' survival model. 
##' @param ascertained if the pair are sampled only when there is an event. This is in contrast to
##' case.control sampling where a proband is given. This can be combined with control probands. Pair-call
##' of twostage is needed  and second column of pairs are the first jump time with an event for ascertained pairs,
##' or time of control proband.
##' @param shut.up to make the program more silent in the context of iterative procedures for case-control
##' and ascertained sampling
##' @export
survival.iterative <- function(margsurv,data=sys.parent(),score.method="fisher.scoring",Nit=60,detail=0,clusters=NULL,
             silent=1,weights=NULL, control=list(),theta=NULL,theta.des=NULL,
             var.link=1,iid=1,step=0.5,model="clayton.oakes",
             marginal.trunc=NULL,marginal.survival=NULL,marginal.status=NULL,strata=NULL,
             se.clusters=NULL,max.clust=NULL,numDeriv=0,random.design=NULL,pairs=NULL,pairs.rvs=NULL,numDeriv.method="simple",
             additive.gamma.sum=NULL,var.par=1,cr.models=NULL,case.control=0,ascertained=0,shut.up=0)
{ ## {{{
## {{{ seting up design and variables
two.stage <- 0; 
rate.sim <- 1; sym=1; var.func <- NULL
if (model=="clayton.oakes" || model=="gamma") dep.model <- 1
else if (model=="plackett" || model=="or") dep.model <- 2
else stop("Model must by either clayton.oakes or plackett \n"); 
start.time <- NULL; ptrunc <- NULL; psurvmarg <- NULL; status <- NULL
fix.baseline <- 0; 
convergence.bp <- 1;  ### to control if baseline profiler converges
if ((!is.null(margsurv)) | (!is.null(marginal.survival))) fix.baseline <- 1


if (!is.null(margsurv)) {
if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen") { ## {{{
	 formula<-attr(margsurv,"Formula");
	 beta.fixed <- attr(margsurv,"beta.fixed")
	 if (is.null(beta.fixed)) beta.fixed <- 1; 
	 ldata<-aalen.des(formula,data=data,model="cox.aalen");
	 id <- attr(margsurv,"id"); 
	 mclusters <- attr(margsurv,"cluster.call")
	 X<-ldata$X; 
	 time<-ldata$time2; 
	 Z<-ldata$Z;  
	 status<-ldata$status;
	 time2 <- attr(margsurv,"stop"); 
	 start.time <- attr(margsurv,"start")
	 antpers<-nrow(X);
         if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else { Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
         if (npar==TRUE) {Z<-matrix(0,antpers,1); pz<-1; fixed<-0;} else {fixed<-1;pz<-ncol(Z);}
	 px<-ncol(X);

         if (is.null(clusters) && is.null(mclusters)) stop("No cluster variabel specified in marginal or twostage call\n"); 
         if (is.null(clusters)) clusters <- mclusters 
###	 else if (sum(abs(clusters-mclusters))>0) 
###         cat("Warning: Clusters for marginal model different than those specified for two.stage\n"); 

###         if (!is.null(attr(margsurv,"max.clust")))
###         if ((attr(margsurv,"max.clust")< attr(margsurv,"orig.max.clust")) && (!is.null(mclusters))) 
###		  cat("Warning: Probably want to estimate marginal model with max.clust=NULL\n"); 

	 if (nrow(X)!=length(clusters)) stop("Length of Marginal survival data not consistent with cluster length\n"); 
## }}}
} else { ### coxph ## {{{
	  notaylor <- 1
	  antpers <- margsurv$n
	  id <- 0:(antpers-1)
	  mt <- model.frame(margsurv)
	  Y <- model.extract(mt, "response")
	  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
	  if (attr(Y, "type") == "right") {
	      time2 <- Y[, "time"]; 
	      status <- Y[, "status"]
		start.time <- rep(0,antpers);
		} else {
		 start.time <- Y[, 1]; 
		 time2 <- Y[, 2];
		 status <- Y[, 3];
		}
###	   Z <- na.omit(model.matrix(margsurv)[,-1]) ## Discard intercept
	   Z <- matrix(1,antpers,length(coef(margsurv)));

	   if (is.null(clusters)) stop("must give clusters for coxph\n");
	   cluster.call <- clusters; 
	   X <- matrix(1,antpers,1); ### Z <- matrix(0,antpers,1); ### no use for these
	   px <- 1; pz <- ncol(Z); 
	   start <- rep(0,antpers);
	   beta.fixed <- 0
	   semi <- 1
###	   start.time <- 0
} ## }}}
} else { start.time<- 0}

###print("00"); print(summary(psurvmarg)); print(summary(ptrunc))

  if (any(abs(start.time)>0)) lefttrunk <- 1 else lefttrunk <- 0
  if (!is.null(start.time)) {
      if (any(abs(start.time)>0)) lefttrunk <- 1  else lefttrunk <- 0;  
  } else lefttrunk <- 0

if (!is.null(margsurv))  {
  if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen")  { ## {{{
         resi <- residualsTimereg(margsurv,data=data) 
         RR  <- resi$RR
	 psurvmarg <- exp(-resi$cumhaz); 
         ptrunc <- rep(1,length(psurvmarg)); 
	 if (lefttrunk==1) ptrunc <- exp(-resi$cumhazleft); 
  } ## }}}
  else if (class(margsurv)=="coxph") {  ## {{{
	### some problems here when data is different from data used in margsurv
       notaylor <- 1
       residuals <- residuals(margsurv)
       cumhaz <- status-residuals
       psurvmarg <- exp(-cumhaz); 
       cumhazleft <- rep(0,antpers)
       ptrunc <- rep(1,length(psurvmarg)); 
       RR<- exp(margsurv$linear.predictors-sum(margsurv$means*coef(margsurv)))
        if ((lefttrunk==1)) { 
         stop("Use cox.aalen function for truncation case \n"); 
         baseout <- survival::basehaz(margsurv,centered=FALSE); 
         cum <- cbind(baseout$time,baseout$hazard)
	 cum <- Cpred(cum,start.time)[,2]
	 ptrunc <- exp(-cum * RR)
	}
  } ## }}}
}

  antpers <- nrow(data); ## mydim(marginal.survival)[1]
  RR <-  rep(1,antpers); 

  if (is.null(psurvmarg)) psurvmarg <- rep(1,antpers);
  if (!is.null(marginal.survival)) psurvmarg <- marginal.survival 
  if (!is.null(marginal.trunc)) ptrunc <- marginal.trunc 
  if (is.null(ptrunc)) ptrunc <- rep(1,length(psurvmarg))

###  print(summary(psurvmarg)); print(summary(ptrunc))

  if (!is.null(marginal.status)) status <- marginal.status 
  if (is.null(status) & is.null(cr.models)) stop("must give status variable for survival via either margninal model (margsurv), marginal.status or as cr.models \n"); 

  if (is.null(weights)==TRUE) weights <- rep(1,antpers); 
  if (is.null(strata)==TRUE) strata<- rep(1,antpers); 
  if (length(strata)!=antpers) stop("Strata must have length equal to number of data points \n"); 

  ## {{{ cluster set up
  cluster.call <- clusters
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
  ## }}} 

###  pxz <- px + pz;

   if (!is.null(random.design)) { ### different parameters for Additive random effects # {{{
     dep.model <- 3

###  if (is.null(random.design)) random.design <- matrix(1,antpers,1); 
     dim.rv <- ncol(random.design); 
     if (is.null(theta.des)) theta.des<-diag(dim.rv);


###     ptheta <- dimpar <- ncol(theta.des); 
###   if (dim(theta.des)[2]!=ncol(random.design)) 
###   stop("nrow(theta.des)!= ncol(random.design),\nspecifies restrictions on paramters, if theta.des not given =diag (free)\n"); 
 } else { random.design <- matrix(0,1,1);  dim.rv <- 1; 
           additive.gamma.sum <- matrix(1,1,1); 
   }

  if (is.null(theta.des)) ptheta<-1; 
  if (is.null(theta.des)) theta.des<-matrix(1,antpers,ptheta); ###  else theta.des<-as.matrix(theta.des); 
###    ptheta<-ncol(theta.des); 
###    if (nrow(theta.des)!=antpers) stop("Theta design does not have correct dim");

  if (length(dim(theta.des))==3) ptheta<-dim(theta.des)[2] else if (length(dim(theta.des))==2) ptheta<-ncol(theta.des)
  if (nrow(theta.des)!=antpers & dep.model!=3 ) stop("Theta design does not have correct dim");

  if (length(dim(theta.des))!=3) theta.des <- as.matrix(theta.des)

  if (is.null(theta)==TRUE) {
         if (var.link==1) theta<- rep(-0.7,ptheta);  
         if (var.link==0) theta<- rep(exp(-0.7),ptheta);   
  }       

  if (length(theta)!=ptheta) {
###	 warning("dimensions of theta.des and theta do not match\n"); 
###         print(theta); 
         theta<-rep(theta[1],ptheta); 
  }
  theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 

  if (maxclust==1) stop("No clusters, maxclust size=1\n"); 

  antpairs <- 1; ### to define 
  if (is.null(additive.gamma.sum)) additive.gamma.sum <- matrix(1,dim.rv,ptheta)

  if (!is.null(pairs)) { pair.structure <- 1;} else  pair.structure <- 0;  

## ppprint 
###  print(c(pair.structure,dep.model,fix.baseline))
###  print(head(theta.des))
###  print(c(case.control,ascertained))

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
###	       print("laver 3-dim design "); 
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
  }# }}}
  
  ## }}}

###  setting up arguments for Aalen baseline profile estimates
 if (fix.baseline==0)  { ## {{{  when baseline is estimated when baseline is estimated 
 if (is.null(cr.models)) stop("give hazard models for different causes, ex cr.models=list(Surv(time,status==1)~+1,Surv(time,status==2)~+1) \n")

      if (case.control==0 & ascertained==0) { ## {{{ 
## {{{ setting up random effects and covariates for marginal modelling
        timestatus <- all.vars(cr.models[[1]])
	times <- data[,timestatus[1]]
	if (is.null(status)) status <- data[,timestatus[2]]
	lstatus <- data[,timestatus[2]]
	### organize increments according to overall jump-times 
	jumps <- lstatus!=0
	dtimes <- times[jumps]
	st <- order(dtimes)
	dtimesst <- dtimes[st]
	dcauses <- lstatus[jumps][st]
	ids <- (1:nrow(data))[jumps][st]

        nc <- 0
	for (i in 1:length(cr.models)) {
              a2 <- aalen.des(as.formula(cr.models[[i]]),data=data)
	      X <- a2$X
	      nc <- nc+ncol(X)
	}
	dBaalen <- matrix(0,length(dtimes),nc)
	xjump <- array(0,c(length(cr.models),nc,length(ids)))

	## first compute marginal aalen models for all causes
	a <- list(); da <- list(); 
	for (i in 1:length(cr.models)) {
	      a[[i]] <-  aalen(as.formula(cr.models[[i]]),data=data,robust=0,weights=weights)
              a2 <- aalen.des(as.formula(cr.models[[i]]),data=data)
	      X <- a2$X
	      da[[i]] <- apply(a[[i]]$cum[,-1,drop=FALSE],2,diff)
	      jumpsi <- (1:length(dtimes))[dcauses==i]
	      if (i==1) fp <- 1 
	      indexc <- fp:(fp+ncol(X)-1)
	      dBaalen[jumpsi,indexc] <- da[[i]]
              xjump[i,indexc,] <- t(X[ids,])
	      fp <- ncol(X)+1
        }

	## }}} 

       	#### organize subject specific random variables and design
        ###  for additive gamma model
	## {{{ 
	dimt <- dim(theta.des[,,1,drop=FALSE])[-3]
	dimr <- dim(random.design[,,,drop=FALSE])
	mtheta.des <- array(0,c(dimt,nrow(data)))
	mrv.des <- array(0,c(dimr[1]/2,dimr[2],nrow(data)))
	nrv.des <- rep(0,nrow(data))
	nrv.des[pairs[,1]] <- pairs.rvs
	nrv.des[pairs[,2]] <- pairs.rvs
	mtheta.des[,,pairs[,1]] <- theta.des
	mtheta.des[,,pairs[,2]] <- theta.des
	mrv.des[,,pairs[,1]] <- random.design[1:(dimr[1]/2),,,drop=FALSE]
	mrv.des[,,pairs[,2]] <- random.design[(dimr[1]/2+1):dimr[1],,,drop=FALSE]
	### array thetades to jump times (subjects)
	mtheta.des <- mtheta.des[,,ids,drop=FALSE]
	### array randomdes to jump times (subjects)
	mrv.des <- mrv.des[,,ids,drop=FALSE]
	nrv.des <- pairs.rvs[ids]
	## }}} 

###       	#### organize subject specific random variables and design
###        ###  for additive gamma model
###	## {{{ 
###	dimt <- dim(theta.des[,,1])
###	dimr <- dim(random.design[,,])
###	mtheta.des <- array(0,c(dimt,nrow(data)))
###	mrv.des <- array(0,c(dimr[1]/2,dimr[2],nrow(data)))
###	mtheta.des[,,pairs[,1]] <- theta.des
###	mtheta.des[,,pairs[,2]] <- theta.des
###	mrv.des[,,pairs[,1]] <- random.design[1:(dimr[1]/2),,,drop=FALSE]
###	mrv.des[,,pairs[,2]] <- random.design[(dimr[1]/2+1):dimr[1],,,drop=FALSE]
###	nrv.des <- rep(0,nrow(data))
###	nrv.des[pairs[,1]] <- pairs.rvs
###	nrv.des[pairs[,2]] <- pairs.rvs
###	### array thetades to jump times (subjects)
###	mtheta.des <- mtheta.des[,,ids,drop=FALSE]
###	### array randomdes to jump times (subjects)
###	mrv.des <- mrv.des[,,ids,drop=FALSE]
###	nrv.des <- pairs.rvs[ids]
###	## }}} 

      } ## }}} 

      if (case.control==1 || ascertained==1) { ## {{{ 

###      print(dim(data)); print(summary(pairs))
         data1 <-        data[pairs[,1],]
         data.proband <- data[pairs[,2],]
	 weights1 <-     weights[pairs[,1]]
###	print(summary(data.proband)); print(summary(data1))

## {{{ setting up designs for jump times 
        timestatus <- all.vars(cr.models[[1]])
        if (is.null(status)) status <- data[,timestatus[2]]
	alltimes      <- data[,timestatus[1]]
	times         <- data1[,timestatus[1]]
	lstatus       <- data1[,timestatus[2]]
	timescase     <- data.proband[,timestatus[1]]
	lstatuscase   <- data.proband[,timestatus[2]]
	### organize increments according to overall jump-times 
	jumps         <- lstatus!=0
	dtimes        <- times[jumps]
	dtimescase    <- timescase[jumps]
	st            <- order(dtimes)
	dtimesst      <- dtimes[st]
	dtimesstcase  <- dtimescase[st]
	dcauses       <- lstatus[jumps][st]
	dcausescase   <- lstatuscase[jumps][st]
	ids           <- (1:nrow(data1))[jumps][st]
	###
	### delayed entry for case because of ascertained sampling
	### controls are however control probands, and have entry=0 
	entry <- timescase*lstatuscase
	data1$entry <- entry
	cr.models2 <- list()
	if (ascertained==1) {
           for (i in 1:length(cr.models)) {
	      cr.models2[[i]] <- update(cr.models[[i]],as.formula(paste("Surv(entry,",timestatus[1],",",timestatus[2],")~.",sep="")))
	    }
	} else cr.models2 <- cr.models

        nc <- 0
	for (i in 1:length(cr.models)) {
              X <- aalen.des(as.formula(cr.models[[i]]),data=data1)$X
	      nc <- nc+ncol(X)
	}
	dBaalen    <- matrix(0,length(dtimes),nc)
	xjump      <- array(0,c(length(cr.models),nc,length(ids)))
	xjumpcase  <- array(0,c(length(cr.models),nc,length(ids)))

	## first compute marginal aalen models for all causes
	a <- list(); da <- list(); 
        ### starting values for iteration 
        Bit <- Bitcase <- c()
	for (i in 1:length(cr.models)) { ## {{{ 
	      a[[i]]  <-  aalen(as.formula(cr.models2[[i]]),data=data1,robust=0,weights=weights1)
	      da[[i]] <- apply(a[[i]]$cum[,-1,drop=FALSE],2,diff)
	      jumpsi  <- (1:length(dtimes))[dcauses==i]
              X       <- aalen.des(as.formula(cr.models[[i]]),data=data1)$X
              Xcase   <- aalen.des(as.formula(cr.models[[i]]),data=data.proband)$X
	      if (i==1) fp <- 1 
	      indexc  <- fp:(fp+ncol(X)-1)
	      dBaalen[jumpsi,indexc] <- da[[i]]
              xjump[i,indexc,]       <- t(X[ids,])
              xjumpcase[i,indexc,]   <- t(Xcase[ids,])
	      fp <- fp+ncol(X)
	      ### starting values 
	      Bit <- cbind(Bit,Cpred(a[[i]]$cum,dtimesst)[,-1,drop=FALSE])
       } ## }}} 
       Bit.ini <- Bit
       ## }}} 

####  organize subject specific random variables and design
####  already done in basic pairwise setup 
	mtheta.des <- theta.des[,,ids,drop=FALSE]
	### array randomdes to jump times (subjects)
	mrv.des <- random.design[,,ids,drop=FALSE]
	nrv.des <- pairs.rvs[ids]
      } ## }}} 

 }  else {
	 mrv.des <- array(0,c(1,1,1)); mtheta.des <- array(0,c(1,1,1)); margthetades <- array(0,c(1,1,1)); 
	 xjump <- array(0,c(1,1,1)); dBaalen <- matrix(0,1,1); nrv.des <- 3
 } ## }}} 

  loglike <- function(par) 
  { ## {{{

     if (pair.structure==0 | dep.model!=3) Xtheta <- as.matrix(theta.des) %*% matrix(c(par),nrow=ptheta,ncol=1);
     if (pair.structure==1 & dep.model==3) Xtheta <- matrix(0,antpers,1); ## not needed 
     DXtheta <- array(0,c(1,1,1));

      if (var.link==1 & dep.model==3) epar <- c(exp(par)) else epar <- c(par)
      partheta <- epar

      if (var.par==1 & dep.model==3) {
	 ## from variances to 
	 if (is.null(var.func)) {
	    sp <- sum(epar)
	    partheta <- epar/sp^2 
         } else partheta <- epar ## par.func(epar)
      } 


      if (pair.structure==0) {
	      outl<-.Call("twostageloglikeRV", ## {{{ only two stage model for this option
	      icause=status,ipmargsurv=psurvmarg, 
	      itheta=c(partheta),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),ithetades=theta.des,
	      icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
	      ivarlink=var.link,iid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
	      itrunkp=ptrunc,istrata=as.numeric(strata),iseclusters=se.clusters,iantiid=antiid,
	      irvdes=random.design,iags=additive.gamma.sum,iascertained=ascertained,
              PACKAGE="mets") ## }}}
      }
      else { ## pair-structure 
	     ## twostage model,    case.control option,   profile out baseline 
	     ## conditional model, case.control option,   profile out baseline 
      if (two.stage==1) { ## {{{ two-stage model 

         if (fix.baseline==0) ## if baseline is not given 
         {
          cum1 <- cbind(dtimesst,Bit)
          if ( (case.control==1 || ascertained==1) & (convergence.bp==1)) { ## {{{  profiles out baseline under case-control/ascertainment sampling

###	      ## initial values , only one cr.model for survival 
###           Bit <- cbind(Cpred(a[[1]]$cum,dtimesst)[,-1])
             if (detail>1) plot(dtimesst,Bit,type="l",main="Bit")

             if (ncol(Bit)==0) Bit <- Bit.ini
             Bitcase  <- Cpred(cbind(dtimesst,Bit),dtimesstcase)[,-1,drop=FALSE]
             Bitcase <- .Call("MatxCube",Bitcase,dim(xjumpcase),xjumpcase,PACKAGE="mets")$X

             for (j in 1:5) { ## {{{ profile via iteration 
                   cncc <- .Call("BhatAddGamCC",1,dBaalen,dcauses,dim(xjump),xjump,
		                 c(partheta),dim(mtheta.des),mtheta.des,additive.gamma.sum,var.link, 
				 dim(mrv.des),mrv.des,nrv.des,1,Bit,Bitcase,dcausescase,PACKAGE="mets")
		   d <- max(abs(Bit-cncc$B))
		   if (detail>1) print(d)
		   Bit <- cncc$B
###		   if (detail>1) print(c(par,epar,partheta)); 
###		   if (detail>1) print(summary(Bit)); 
		   if (detail>1) print(summary(cncc$caseweights))
		   cum1 <- cbind(dtimesst,cncc$B)
		   Bitcase  <-cbind(Cpred(cum1,dtimesstcase)[,-1])
###		   if (detail>1) print(summary(Bitcase))
		   if (detail>1) lines(dtimesst,Bit,col=j+1); 
		   if (is.na(d)) { 
			  if (shut.up==0) cat("Baseline profiler gives missing values\n");  
		          Bit <- Bit.ini; cum1 <- cbind(dtimesst,Bit); convergence.bp <<- 0; break;
		   }
		   Bitcase <- .Call("MatxCube",Bitcase,dim(xjumpcase),xjumpcase,PACKAGE="mets")$X
		   if (d<0.00001) break; 
           } ## }}} 


	   nulrow    <- rep(0,ncol(Bit)+1)
	   pbases    <- Cpred(rbind(nulrow,cbind(dtimesst,Bit)),alltimes)[,-1,drop=FALSE]
           X         <- aalen.des(as.formula(cr.models[[1]]),data=data)$X
	   psurvmarg <- exp(-apply(X*pbases,1,sum))  ## psurv given baseline 
	   if (ascertained==1) {
              Xcase     <- aalen.des(as.formula(cr.models[[1]]),data=data.proband)$X
              X         <- aalen.des(as.formula(cr.models[[1]]),data=data1)$X
	      pba.case  <- Cpred(rbind(nulrow,cbind(dtimesst,Bit)),entry)[,-1,drop=FALSE]
	      ptrunc    <- rep(0,nrow(data))
	      ### for control probands ptrunc=1, thus no adjustment
	      ptrunc[pairs[,1]]     <- exp(-apply(X*    pba.case,1,sum)*lstatuscase) ## delayed entry at time of ascertainment proband 
 	      ptrunc[pairs[,2]]     <- exp(-apply(Xcase*pba.case,1,sum)*lstatuscase)
	   }
###	   print(head(cbind(psurvmarg,ptrunc)))
###	   print(summary(psurvmarg))
###	   print(summary(ptrunc))
###	   print(dim(pbases)); 

          } ## }}} 
          }

###      browser()
###      print(dim(random.design))


          outl<-.Call("twostageloglikeRVpairs", ## {{{
          icause=status,ipmargsurv=psurvmarg, 
          itheta=c(partheta),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),ithetades=theta.des,
          icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
          ivarlink=var.link,iiid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
          itrunkp=ptrunc,istrata=as.numeric(strata),iseclusters=se.clusters,iantiid=antiid,
          irvdes=random.design,
          idimthetades=dim(theta.des),idimrvdes=dim(random.design),irvs=pairs.rvs,iags=additive.gamma.sum, 
	  iascertained=ascertained,PACKAGE="mets")
	  ## }}} 

          if (fix.baseline==0)  {
              outl$baseline <- cum1; 
	      outl$marginal.surv <- psurvmarg; 
	      outl$marginal.trunc <- ptrunc
	  }

      } ## }}} 
      else { ## {{{  survival model two.stage==0

###      cum1 <- cbind(dtimesst,Bit)
	 entry.cause  <- rep(0,nrow(data))
###	 proband.time <- rep(0,nrow(data))

         ### update aalen type baseline  
	 if (fix.baseline==0)  { ## {{{ 

         if ((case.control==1 || ascertained==1) & (convergence.bp==1)) { ## {{{  profiles out baseline under case-control/ascertainment sampling

	     if (detail>1) print(summary(Bit))
             if (detail>1) matplot(dtimesst,Bit,type="l",main="Bit",ylim=c(0,2))
             if (ncol(Bit)==0) Bit <- Bit.ini
             Bitcase  <- Cpred(cbind(dtimesst,Bit),dtimesstcase)[,-1,drop=FALSE]
             Bitcase <- .Call("MatxCube",Bitcase,dim(xjumpcase),xjumpcase,PACKAGE="mets")$X

            for (j in 1:10) { ## {{{ profile via iteration 
                   profile.baseline <- .Call("BhatAddGamCC",0,dBaalen,dcauses,dim(xjump),xjump,
	                               c(partheta), dim(mtheta.des),mtheta.des, additive.gamma.sum,var.link, 
	                               dim(mrv.des),mrv.des,nrv.des,1,Bit,Bitcase,dcausescase,PACKAGE="mets")
                   d <- max(abs(Bit-profile.baseline$B))
		   Bit <- profile.baseline$B
		   cum1 <- cbind(dtimesst,Bit)
		   Bitcase  <-cbind(Cpred(cum1,dtimesstcase)[,-1])
###		   matlines(dtimesst,Bit,type="l",col=j+1)

		   if (detail>1) matlines(dtimesst,Bit,col=j+1); 
		   if (is.na(d)) { 
			  if (shut.up==0) cat("Baseline profiler gives missing values\n");  
		          Bit <- Bit.ini; cum1 <- cbind(dtimesst,Bit); convergence.bp <<- 0; break;
		   }
		   Bitcase <- .Call("MatxCube",Bitcase,dim(xjumpcase),xjumpcase,PACKAGE="mets")$X
                   if (d<0.00001) break; 
           } ## }}} 

###	     print("profile slut"); 
###	   plot(cum1); abline(c(0,1))

           ### makes cumulative hazard for all subjects
	   nulrow <- rep(0,ncol(Bit)+1)
	   pbases <- Cpred(rbind(nulrow,cbind(dtimesst,Bit)),alltimes)[,-1,drop=FALSE]
	   psurvmarg <- c()     ### update psurvmarg
	   if (ascertained==1 || case.control==1)  { ### sets up truncation probabilities to match situation
		   pbase.case  <- Cpred(rbind(nulrow,cbind(dtimesst,Bit)),timescase)[,-1,drop=FALSE]
		   ptrunc <- c() ### update ptrunc
	   }
           for (i in 1:length(cr.models)) {
                  if (i==1) fp <- 1 
                  a2 <- aalen.des(as.formula(cr.models[[i]]),data=data)
	          X <- a2$X
		  indexc <- fp:(fp+ncol(X)-1)
	          psurvmarg <- cbind(psurvmarg,apply(X*pbases[,indexc],1,sum))
	          if (ascertained==1 || case.control==1) {
                     Xcase     <- aalen.des(as.formula(cr.models[[i]]),data=data.proband)$X
                     X         <- aalen.des(as.formula(cr.models[[i]]),data=data1)$X
###		     print(dim(Xcase)); print(dim(pbase.case[,indexc,drop=FALSE]))
		     ptrunc.new <- rep(0,length(alltimes))
		     ## delayed entry at time of ascertainment proband, for case control no adjustment for first 
		     if (ascertained==1) ptrunc.new[pairs[,1]] <- apply(X*pbase.case[,indexc,drop=FALSE],1,sum)*lstatuscase 
	             else   ptrunc.new[pairs[,1]] <- 0
		     ## for second component adjustment for marginal or ascertainment
		     ptrunc.new[pairs[,2]] <- apply(Xcase*pbase.case[,indexc,drop=FALSE],1,sum)
		     ptrunc <- cbind(ptrunc,ptrunc.new)
		  }
	          fp <- fp+ncol(X)
	   }
	      
          } ## }}} 
	  else { ## {{{ profile out baseline, recursive estimator


          profile.baseline  <- .Call("BhatAddGam",recursive=1,
          dBaalen,dcauses,dim(xjump),xjump,c(partheta),
	  dim(mtheta.des),mtheta.des, 
	  additive.gamma.sum,0,dim(mrv.des),mrv.des,0,matrix(0,1,1),PACKAGE="mets")

###	  print(summary(cbind(dtimesst,profile.baseline$B)))
###	  matplot(dtimesst,profile.baseline$B,type="l")
###	  print(head(profile.baseline$B))
###	  abline(c(0,0.2),col=3)
###	  abline(c(0,0.4),col=3)

	  marg.model <- "no-cox"
          if (marg.model=="cox") {# {{{
	      Bit <- profile.baseline$B
              Bit <- .Call("MatxCube",Bit,dim(xjump),xjump,PACKAGE="mets")$X
	      caseweights <- profile.baseline$caseweights

	      pp <- timereg.formula(as.formula(cr.models[[i]]))
	      profile.cox <- cox.aalen(pp,data=data,robust=0,caseweight=caseweights)
	      print(summary(profile.cox))
	      plot(profile.cox)
	  }# }}}

	   nulrow <- rep(0,ncol(dBaalen)+1)
	   pbases <- Cpred(rbind(nulrow,cbind(dtimesst,profile.baseline$B)),times)[,-1,drop=FALSE]
	   psurvmarg <- c()
           for (i in 1:length(cr.models)) {
	          if (i==1) fp <- 1 
                  a2 <- aalen.des(as.formula(cr.models[[i]]),data=data)
	          X <- a2$X
	          indexc <- fp:(fp+ncol(X)-1)
	          psurvmarg <- cbind(psurvmarg,apply(X*pbases[,indexc],1,sum))
	          fp <- fp+ncol(X)
	  }
          ### no truncation in this case ? 
	  ptrunc <- 0*psurvmarg

	  } ## }}} 

	 } ## }}}  

###	 if (fix.baseline==1) if (is.null(psurvmarg)) stop("must provide baselines or set fix.baseline=0\n"); 

###	 print("er vi her surv "); 
###	 print(fix.baseline)
###         print(dim(as.matrix(psurvmarg))); print(dim(as.matrix(ptrunc)))
###	 print(dim(pairs))
###	 print(head(pairs))
###         print(summary(psurvmarg[pairs,])); print(summary(ptrunc[pairs,]))

      ### cumulative hazard for this model  when fix.baseline==1
      if (fix.baseline==1 ) {
	      psurvmarg <- -log(psurvmarg); 
	      ptrunc <- -log(ptrunc); 
      }


      outl<-.Call("survivalloglikeRVpairs",icause=status,ipmargsurv=as.matrix(psurvmarg), 
      itheta=c(partheta),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),ithetades=theta.des,
      icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
      iiid=iid,iweights=weights,isilent=silent,idepmodel=dep.model,
      itrunkp=as.matrix(ptrunc),istrata=as.numeric(strata),iseclusters=se.clusters,iantiid=antiid,
      irvdes=random.design,
      idimthetades=dim(theta.des),idimrvdes=dim(random.design),
      irvs=pairs.rvs,iags=additive.gamma.sum,ientry.cause=entry.cause,iascertained=(ascertained+case.control>0)*1,
      PACKAGE="mets")

      if (fix.baseline==0)  {
	      outl$baseline <- cbind(dtimesst,profile.baseline$B); 
	      outl$marginal.surv <- psurvmarg; 
	      outl$marginal.trunc <- ptrunc
      }

      } ## }}} 
      }  

    if (detail==3) print(c(partheta,outl$loglike))

    ## variance parametrization, and inverse.link 
    if (dep.model==3) {# {{{
    if (var.par==1) {
         ## from variances to and with sum for all random effects 
         if (is.null(var.func)) {
	 if (var.link==0)  {
###		 print(c(sp,epar))
             mm <- matrix(-epar*2*sp,length(epar),length(epar))
	     diag(mm) <- sp^2-epar*2*sp
###	    print(mm)
	 } else {
            mm <- -c(epar) %o% c(epar)*2*sp
            diag(mm) <- epar*sp^2-epar^2*2*sp
###	    print(mm)
	 }
	    mm <- mm/sp^4
	 } else mm  <- numDeriv::hessian(var.func,par)
      } else {
	   if (var.link==0) mm <- diag(length(epar)) else 
			  mm <- diag(length(c(epar)))*c(epar)
      }
      }# }}}

###    print(c(var.link,dep.model,var.par))
###    print("hh"); print(mm); print(outl$score)

    if (dep.model==3) {# {{{
###	    print(dim(mm))
       outl$score <-  t(mm) %*% outl$score
       outl$Dscore <- t(mm) %*% outl$Dscore %*% mm
       if (iid==1) outl$theta.iid <- t(t(mm) %*% t(outl$theta.iid))

###       print(crossprod(outl$theta.iid)); print(outl$Dscore)
###       print(c(outl$score))
###       print(apply(outl$theta.iid,2,sum))
    }# }}}

    attr(outl,"gradient") <-outl$score 
    if (oout==0) ret <- c(-1*outl$loglike) else if (oout==1) ret <- sum(outl$score^2) else if (oout==2) ret <- outl else ret <- outl$score
    return(ret)
  } ## }}}

  if (score.method=="optimize" && ptheta!=1) {
     cat("optimize only works for d==1, score.mehod set to nlminb \n"); 
     score.method <- "nlminb";
  } 

  score1 <- NULL
  theta.iid <- NULL
  logl <- NULL
  p <- theta

    if (score.method=="fisher.scoring") { ## {{{
        oout <- 2;  ### output control for obj
        if (Nit>0) 
            for (i in 1:Nit)
            {
                    out <- loglike(p)
		    ## updating starting values for cumulative baselines
		    if (fix.baseline==0) Bit <- out$baseline[,-1,drop=FALSE]
		    if (fix.baseline==1) hess <-  out$Dscore
                    ###  uses simple second derivative for computing derivative of score
		    if (numDeriv==2 || (((fix.baseline==0)) & (i==1))) {
			    oout <- 3
			    hess <- numDeriv::jacobian(loglike,p,method="simple")
			    oout <- 2
		    }
                    if (!is.na(sum(hess))) hessi <- lava::Inverse(hess) else hessi <- hess 
                    if (detail==1) {## {{{
                        cat(paste("Fisher-Scoring ===================: it=",i,"\n")); 
                        cat("theta:");print(c(theta))
                        cat("loglike:");cat(c(out$loglike),"\n"); 
                        cat("score:");cat(c(out$score),"\n"); 
                        cat("hess:\n"); cat(hess,"\n"); 
                    }## }}}
                    delta <- step*( hessi %*% out$score )
		    ### update p, but note that score and derivative in fact related to previous p 
		    ### unless Nit=0, 
	            if (Nit>0) {
                    p <- p - delta
                    theta <- p; 
		    }
                    if (is.nan(sum(out$score))) break; 
                    if (sum(abs(out$score))<0.00001) break; 
                    if (max(theta)>20 & var.link==1) { cat("theta too large lacking convergence \n"); break; }
           }
        if (!is.nan(sum(p))) { 
            if (detail==1 && iid==1) cat("iid decomposition\n"); 
            out <- loglike(p) 
            logl <- out$loglike
            score1 <- score <- out$score
            oout <- 0; 
            hess1 <- hess  <- -1*out$Dscore 
###            if (numDeriv==2 || (baseline.fix==0 && two.stage==0)) {
###                    oout <- 3
###                    hess <- numDeriv::jacobian(loglike,p,method="simple")
###		    oout <- 2
###	    }
            if (iid==1) { theta.iid <- out$theta.iid
                if (class(margsurv)=="phreg") {  ## {{{
		       browser()
		       ## order after time
		       D1dltheta1 <- out$D1dltheta1[xx$order+1,]
		       D2dltheta1 <- out$D1dltheta1[xx$order+1,]

		       ### baseline iid 
			xx <- margsurv$cox.prep
			S0i2 <- S0i <- rep(0,length(xx$strata))
			S0i[xx$jumps+1] <-  1/margsurv$S0
			rr <- exp(margsurv$cox.prep$X  %*% margsurv$coef)
			cumhazt <- cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum
			psurvmarg <- exp(-cumhazt*rr)
	        } ## }}}
	   }

            if (detail==1 && iid==1) cat("finished iid decomposition\n"); 
	    ### for profile solutions update second derivative at final 
	    if (numDeriv==2 || ((fix.baseline==0))) {
                    oout <- 3
                    hess <- numDeriv::jacobian(loglike,p,method="simple")
		    oout <- 2
		    }
        }
        if (numDeriv>=1) {
            if (detail==1 ) cat("starting numDeriv for second derivative \n"); 
            oout <- 0; 
            score2 <- numDeriv::jacobian(loglike,p)
	    score1 <- matrix(score2,ncol=1)
            oout <- 3
            hess <- numDeriv::jacobian(loglike,p,method="simple")
            if (detail==1 ) cat("finished numDeriv for second derivative \n"); 
        }
        if (detail==1 & Nit==0) {## {{{
            cat(paste("Fisher-Scoring ===================: final","\n")); 
            cat("theta:");print(c(p))
            cat("loglike:");cat(c(out$loglike),"\n"); 
            cat("score:");cat(c(out$score),"\n"); 
            cat("hess:\n"); cat(hess,"\n"); 
        }## }}}
        if (!is.na(sum(hess))) hessi <- lava::Inverse(hess) else hessi <- diag(nrow(hess))
        ## }}}
  } else if (score.method=="nlminb") { ## {{{ nlminb optimizer
    oout <- 0; 
    if (two.stage==0) oout <- 1 ## score 
    tryCatch(opt <- nlminb(theta,loglike,control=control),error=function(x) NA)
    if (detail==1) print(opt); 
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    oout <- 2
    theta <- opt$par
    out <- loglike(opt$par)
    logl <- out$loglike
    score1 <- score <- out$score
    hess1 <- hess <- -1* out$Dscore
    if (iid==1) theta.iid <- out$theta.iid
    if (numDeriv==1) {
    if (detail==1 ) cat("numDeriv hessian start\n"); 
      oout <- 3; ## returns score 
      hess <- numDeriv::jacobian(loglike,opt$par)
    if (detail==1 ) cat("numDeriv hessian done\n"); 
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
    hess1 <- hess <- -1* out$Dscore
    if (numDeriv==1) {
    if (detail==1 ) cat("numDeriv hessian start\n"); 
      oout <- 3;  ## to get jacobian
      hess <- numDeriv::jacobian(loglike,theta)
    if (detail==1 ) cat("numDeriv hessian done\n"); 
    }
    hessi <- lava::Inverse(hess); 
    if (iid==1) theta.iid <- out$theta.iid
  ## }}}
  } else if (score.method=="nlm") { ## {{{ nlm optimizer
    iid <- 0; oout <- 0; 
    tryCatch(opt <- nlm(loglike,theta,hessian=TRUE,print.level=detail),error=function(x) NA)
    iid <- 1; 
    hess <- opt$hessian
    score <- opt$gradient
    if (detail==1) print(opt); 
    hessi <- lava::Inverse(hess); 
    theta <- opt$estimate
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    oout <- 2
    out <- loglike(opt$estimate)
    logl <- out$loglike
    score1 <- out$score
    hess1 <- out$Dscore
    if (iid==1) theta.iid <- out$theta.iid
  ## }}}
  }  else stop("score.methods = optimize(dim=1) nlm nlminb fisher.scoring\n"); 

## {{{ handling output
  loglikeiid <- NULL
  robvar.theta <- NULL
  likepairs <- NULL
  if (fix.baseline==1)  { marginal.surv <- psurvmarg; marginal.trunc <- ptrunc; } else { marginal.surv <- out$marginal.surv; marginal.trunc <- out$marginal.trunc;} 

  if (iid==1) {
  if (dep.model==3 & pair.structure==1) likepairs <- out$likepairs
  if (dep.model==3 & two.stage==0) {
	  hessi <- -1*hessi
	  all.likepairs <- out$all.likepairs
	  colnames(all.likepairs) <- c("surv","dt","ds","dtds","cause1","cause2")
  }
###  print(crossprod(out$theta.iid) %*% hessi)
     theta.iid <- out$theta.iid %*% hessi
     if (is.null(call.secluster) & is.null(max.clust)) rownames(theta.iid) <- unique(cluster.call) else rownames(theta.iid) <- unique(se.clusters)
     robvar.theta  <- crossprod(theta.iid) 
     loglikeiid <- out$loglikeiid
    } else { all.likepairs <- NULL}
  var.theta <- robvar.theta
  if (is.null(robvar.theta)) var.theta <- hessi

  if (!is.null(colnames(theta.des))) thetanames <- colnames(theta.des) else thetanames <- paste("dependence",1:length(theta),sep="")
  theta <- matrix(theta,length(c(theta)),1)
  if (length(thetanames)==nrow(theta)) { rownames(theta) <- thetanames; rownames(var.theta) <- colnames(var.theta) <- thetanames; }
  if (!is.null(logl)) logl <- -1*logl

  if (convergence.bp==0) theta <- rep(NA,length(theta))
  ud <- list(theta=theta,score=score,hess=hess,hessi=hessi,var.theta=var.theta,model=model,robvar.theta=robvar.theta,
             theta.iid=theta.iid,loglikeiid=loglikeiid,likepairs=likepairs,
	     thetanames=thetanames,loglike=logl,score1=score1,Dscore=out$Dscore,
	     marginal.surv=marginal.surv,marginal.trunc=marginal.trunc,baseline=out$baseline,
	     se=diag(robvar.theta)^.5)
  class(ud) <- "mets.twostage" 
  attr(ud,"response") <- "survival"
  attr(ud,"Formula") <- formula
  attr(ud,"clusters") <- clusters
  attr(ud,"cluster.call") <- cluster.call
  attr(ud,"secluster") <- c(se.clusters)
  attr(ud,"sym")<-sym; 
  attr(ud,"var.link")<-var.link; 
  attr(ud,"var.par")<-var.par; 
  attr(ud,"var.func")<-var.func; 
  attr(ud,"ptheta")<-ptheta
  attr(ud,"antpers")<-antpers; 
  attr(ud,"antclust")<-antclust; 
  attr(ud,"Type") <- model
  attr(ud,"twostage") <- two.stage
  attr(ud,"additive-gamma") <- (dep.model==3)*1
  if (!is.null(marginal.trunc)) attr(ud,"trunclikeiid")<- out$trunclikeiid
  if (dep.model==3 & two.stage==0) attr(ud,"all.likepairs")<- all.likepairs
  if (dep.model==3 ) attr(ud,"additive.gamma.sum") <- additive.gamma.sum
  #likepairs=likepairs,##  if (dep.model==3 & pair.structure==1) attr(ud, "likepairs") <- c(out$likepairs)
  if (dep.model==3 & pair.structure==0) attr(ud, "pardes") <-    theta.des
  if (dep.model==3 & pair.structure==0) attr(ud, "theta.des") <- theta.des
  if (dep.model==3 & pair.structure==1) attr(ud, "pardes") <-    theta.des[,,1]
  if (dep.model==3 & pair.structure==1) attr(ud, "theta.des") <- theta.des[,,1]
  if (dep.model==3 & pair.structure==0) attr(ud, "rv1") <- random.design[1,]
  if (dep.model==3 & pair.structure==1) attr(ud, "rv1") <- random.design[,,1]
  return(ud);
  ## }}}

} ## }}}


##' @export
summary.mets.twostage <- function(object,digits = 3,silent=0,...) 
{ ## {{{
  if (!(inherits(object,"mets.twostage"))) stop("Must be a Two-Stage object")
  
  var.link<-attr(object,"var.link");
  var.par  <- attr(object,"var.par"); 
  model <- object$model
  if ((model=="plackett" || model=="or") )         model <- "or"
  if ((model=="clayton.oakes" || model=="gamma") ) model <- "gamma"
  if ((attr(object,"additive-gamma")==1) & (silent==0)) addgam <- TRUE else addgam <- FALSE

  if ((model=="or") && (silent==0)) cat("Dependence parameter for Odds-Ratio (Plackett) model \n"); 
  if (attr(object,"response")=="binomial") response <- "binomial" else response <- "survival"
  if ((model=="gamma") && (silent==0)) {
	  cat("Dependence parameter for Clayton-Oakes model\n")
	  if (var.par==1 || !addgam) 
	  cat("Variance of Gamma distributed random effects \n"); 
	  if (var.par==0 && addgam)
	  cat("Inverse of variance of Gamma distributed random effects \n"); 
  }
  if (var.link==1 && silent==0) cat("With log-link \n")
  if ((sum(abs(object$score))>0.0001) & (silent==0))  {
	  cat("    Variance parameters did not converge, allow more iterations.\n"); 
	  cat(paste("    Score:",object$score,"  \n")); 
  }
  theta <- object$theta
  if (is.null(rownames(theta)))  rownames(theta) <- paste("dependence",1:nrow(theta),sep="")

###  print(theta)

  coefs <- coef.mets.twostage(object,response=response,...);

  if (attr(object,"additive-gamma")==1 & (!is.null(object$robvar.theta))  ) {
      var.link <- attr(object,"var.link"); 
      var.par  <- attr(object,"var.par"); 
      rv1      <- attr(object,"rv1"); 
      if (is.matrix(rv1)) rv1 <- rv1[1,]
      theta.des <- attr(object,"pardes"); 
      ags <- attr(object,"additive.gamma.sum"); 
      ptheta <- attr(object,"ptheta")
      npar <- nrow(object$theta)
      theta <- object$theta[seq(1,ptheta),1,drop=FALSE]
###   print(theta.des); print(theta); print(theta.des %*% theta); print(rv1); 
      robvar.theta <- object$robvar.theta[seq(1,ptheta),seq(1,ptheta)]
      if (var.link==1) par <- theta.des %*% exp(theta) else  par <- theta.des %*% theta

      if (model=="or" || model=="plackett") var.par<-1 ## same as var.par=0 for this model

      if (attr(object,"twostage")==0) {
         ###         cat("MLE estimates of marginal parameters\n"); 
###         var.gamma <- object$robvar.theta[seq(ptheta+1,npar),seq(ptheta+1,npar)]
###         obj.marginal <- list(gamma=object$theta[seq(ptheta+1,npar),1],var.gamma=var.gamma,robvar.gamma=var.gamma)
###	 marginal <- coefBase(obj.marginal)
      }

      if (var.par==0) { ## {{{ 
	      if (var.link==1) { 
		     fp <- function(p,d,t){  res <- exp(p*t)/(sum(rv1* c(theta.des %*% matrix(exp(p),ncol=1))))^d; 
					     if (t==0) res <- res[1]; return(res); }
             pare <-   lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) exp(p),labels=rownames(theta))
	      } else {
		      fp <- function(p,d,t) {  res <- (p^t)/(sum(rv1* c(theta.des %*% matrix(p,ncol=1))))^d;
					     if (t==0) res <- res[1]; return(res); }
	      pare <- NULL
	      }

             e      <-   lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,1,1),labels=rownames(theta))
             vare <-   lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,2,1),labels=rownames(theta))
             vartot <- lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,1,0))
      } ## }}} 
      if (var.par==1) { ## {{{ 

	   if (var.link==1) { ## {{{ 
	     fp <- function(p,d,t){  res <- exp(p*t)/(sum(rv1* c(theta.des %*% matrix(exp(p),ncol=1))))^d; 
                                     if (t==0) res <- res[1]; return(res); }
             e      <-   lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,1,1),labels=rownames(theta))
             vare   <-   lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) exp(p),labels=rownames(theta))
             vartot <-   lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) sum(exp(p)))
###	     names(e) <- names(vare) <- colnames(coefs)
      } else {
              fp <- function(p,d,t) {  res <- (p^t)/(sum(rv1* c(theta.des %*% matrix(p,ncol=1))))^d;
                                     if (t==0) res <- res[1]; return(res); }
              e     <-  lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,1,1),labels=rownames(theta))
              pare  <-  lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) fp(p,2,1),labels=rownames(theta))
              vartot <- lava::estimate(coef=theta,vcov=robvar.theta,f=function(p) sum(p))
	      vare <- NULL
###	      names(e) <- names(pare) <- colnames(coefs)
	       } ## }}} 
      } ## }}} 
      res <- list(estimates=coefs, type=attr(object,"Type"),h=e,vare=vare,vartot=vartot)

  } else {
	 if (var.link==1) { ## {{{ 
	     if (model=="or") {
		or <-  lava::estimate(coef=object$theta,vcov=object$var.theta,f=function(p) exp(p),labels=rownames(theta))
	        res <- list(estimates=coefs,or=or,type=attr(object,"Type"))
	     } else {
	        vargam <-  lava::estimate(coef=object$theta,vcov=object$var.theta,f=function(p) exp(p),labels=rownames(theta))
###             names(vargam) <- colnames(coefs) 
                res <- list(estimates=coefs,vargam=vargam, type=attr(object,"Type"))
	     }
      } else {
	     if (model=="or") {
             lor <- lava::estimate(coef=object$theta,vcov=object$var.theta,f=function(p) log(p),labels=rownames(theta))
###          names(lor) <- colnames(coefs) 
	     res <- list(estimates=coefs,log.or=lor,type=attr(object,"Type"))
	     } else  {
             res <- list(estimates=coefs,type=attr(object,"Type"))
	     }
      } ## }}} 
  }

###  if (attr(object,"twostage")==0)  res <- c(res,list(marginal=marginal))
  class(res) <- "summary.mets.twostage"
  res
} ## }}}

##' @export
coef.mets.twostage <- function(object,var.link=NULL,response="survival",...)
{ ## {{{
  pt <- attr(object,"ptheta")
  theta <- object$theta[seq(pt),1]

  if (is.null(object$robvar.theta)) 
	  var.theta  <-  object$var.theta[seq(1,pt),seq(1,pt),drop=FALSE] else 
	  var.theta  <-   object$robvar.theta[seq(1,pt),seq(1,pt),drop=FALSE]
  se <- diag(var.theta)^.5

  if (is.null(var.link))
     if (attr(object,"var.link")==1) vlink <- 1 else vlink <- 0
     else vlink <- var.link
  res <- cbind(theta, se )
  wald <- theta/se
  waldp <- (1 - pnorm(abs(wald))) * 2
  if (response=="survival") { 
       if (object$model=="plackett") {
       spearman <- alpha2spear(theta,link=vlink)
       Dspear <- numDeriv::jacobian(alpha2spear,theta,link=vlink) 
       var.spearman <- Dspear %*% var.theta %*%  Dspear
       se.spearman <- diag(var.spearman)^.5
       res <- as.matrix(cbind(res, wald, waldp,spearman,se.spearman))
       if (vlink==1) colnames(res) <- c("log-Coef.", "SE","z", "P-val","Spearman Corr.","SE")
	  else colnames(res) <- c("Coef.", "SE","z", "P-val","Spearman Corr.","SE")
	  if ((!is.null(object$thetanames)) & (nrow(res)==length(object$thetanames))) rownames(res)<-object$thetanames
       }
       if (object$model=="clayton.oakes") {
       kendall <- alpha2kendall(theta,link=vlink)
       Dken <- numDeriv::jacobian(alpha2kendall,theta,link=vlink) 
       var.kendall<- Dken %*% var.theta %*%  Dken
       se.kendall <- diag(var.kendall)^.5
       res <- as.matrix(cbind(res, wald, waldp,kendall,se.kendall))
       if (vlink==1) colnames(res) <- c("log-Coef.", "SE","z", "P-val","Kendall tau","SE")
       else colnames(res) <- c("Coef.", "SE","z", "P-val","Kendall tau","SE")
       if ((!is.null(object$thetanames)) & (nrow(res)==length(object$thetanames))) rownames(res)<-object$thetanames
       }
  }
  return(res)
} ## }}}

##' @export
print.mets.twostage<-function(x,digits=3,...)
{ ## {{{
###  print(x$call); 
  cat("\n")
  print(summary(x,silent=0)); 
} ## }}} 

##' @export
plot.mets.twostage<-function(x,pointwise.ci=1,robust=0,specific.comps=FALSE,
		level=0.05, 
		start.time=0,stop.time=0,add.to.plot=FALSE,mains=TRUE,
                xlab="Time",ylab ="Cumulative regression function",...) 
{ ## {{{
  if (!(inherits(x, 'two.stage'))) stop("Must be a Two-Stage object")
  object <- x; rm(x);  
 
  B<-object$cum; V<-object$var.cum; p<-dim(B)[[2]]; 
  if (robust>=1) V<-object$robvar.cum; 

  if (sum(specific.comps)==FALSE) comp<-2:p else comp<-specific.comps+1
  if (stop.time==0) stop.time<-max(B[,1]);

  med<-B[,1]<=stop.time & B[,1]>=start.time
  B<-B[med,]; Bs<-B[1,];  B<-t(t(B)-Bs); B[,1]<-B[,1]+Bs[1];
  V<-V[med,]; Vs<-V[1,]; V<-t( t(V)-Vs); 
  Vrob<-object$robvar.cum; 
  Vrob<-Vrob[med,]; Vrobs<-Vrob[1,]; Vrob<-t( t(Vrob)-Vrobs); 

  c.alpha<- qnorm(1-level/2)
  for (v in comp) { 
    c.alpha<- qnorm(1-level/2)
    est<-B[,v];ul<-B[,v]+c.alpha*V[,v]^.5;nl<-B[,v]-c.alpha*V[,v]^.5;
    if (add.to.plot==FALSE) 
      {
        plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlab,ylab=ylab) 
        if (mains==TRUE) title(main=colnames(B)[v]); }
    else lines(B[,1],est,type="s"); 
    if (pointwise.ci>=1) {
      lines(B[,1],ul,lty=pointwise.ci,type="s");
      lines(B[,1],nl,lty=pointwise.ci,type="s"); }
    if (robust>=1) {
      lines(B[,1],ul,lty=robust,type="s"); 
      lines(B[,1],nl,lty=robust,type="s"); }
    abline(h=0); 
  }
}  ## }}}

##' @export
matplot.mets.twostage <- function(object,...)
{ ## {{{ 
B <- object$baseline
matplot(B[,1],B[,-1],type="s",...)
} ## }}} 

##' @export
predict.mets.twostage <- function(object,X=NULL,Z=NULL,times=NULL,times2=NULL,theta.des=NULL,diag=TRUE,...)
{ ## {{{
time.coef <- data.frame(object$cum)
if (!is.null(times)) {
cum <- Cpred(object$cum,times);
cum2 <- Cpred(object$cum,times);
} else { cum <- object$cum; cum2 <- object$cum }
if (!is.null(times2)) cum2 <- Cpred(object$cum,times2);

if (is.null(X)) X <- 1;
if (is.null(X) & (!is.null(Z))) { Z <- as.matrix(Z);  X <- matrix(1,nrow(Z),1)}
if (is.null(Z) & (!is.null(X)))  {X <- as.matrix(X);  Z <- matrix(0,nrow(X),1); gamma <- 0}

if (diag==FALSE) {
   time.part <-  X %*% t(cum[,-1]) 
   time.part2 <-  X %*% t(cum2[,-1]) 
   if (!is.null(object$gamma)) { RR <- exp( Z %*% gamma ); 
       cumhaz <- t( t(time.part) * RR ); cumhaz2 <- t( t(time.part2) * RR )}
	    else { cumhaz <- time.part;  cumhaz2 <- time.part2;   }
} else { 
	time.part <-  apply(as.matrix(X*cum[,-1]),1,sum) 
	time.part2 <-  apply(as.matrix(X*cum2[,-1]),1,sum) 
}

if (!is.null(object$gamma)) {
	RR<- exp(Z%*%gamma); 
	cumhaz <- t( t(time.part) * RR );  
	cumhaz2 <- t( t(time.part2) * RR )} else {
		cumhaz <- time.part;  cumhaz2 <- time.part2; 
} 
S1 <- exp(-cumhaz); S2 <- exp(-cumhaz2)

if (attr(object,"var.link")==1) theta  <- exp(object$theta) else theta <- object$theta
if (!is.null(theta.des)) theta <- c(theta.des %*% object$theta)

if (diag==FALSE) St1t2<- (outer(c(S1)^{-(theta)},c(S2)^{-(theta)},FUN="+") - 1)^(-(1/theta)) else 
St1t2<- ((S1^{-(theta)}+S2^{-(theta)})-1)^(-(1/theta))

out=list(St1t2=St1t2,S1=S1,S2=S2,times=times,times2=times2,theta=theta)
return(out)
} ## }}}

##' @export ascertained.pairs
ascertained.pairs <-function (pairs,data,cr.models,bin=FALSE) 
{# {{{
      timestatus <- all.vars(cr.models)
      ### let first event by second column and only
      ### use pairs where first is event 
      apairs <- pairs
      if (bin==TRUE) fj <- ifelse(data[pairs[,1],timestatus[1]] > data[pairs[,2],timestatus[1]],1,2)
      else fj <- ifelse(data[pairs[,1],timestatus[1]] < data[pairs[,2],timestatus[1]],1,2)
      ### change order when 1st comes first 
      apairs[fj==1,1] <- pairs[fj==1,2]
      apairs[fj==1,2] <- pairs[fj==1,1]
      ### only take pairs where first is a jump
      if (bin==FALSE) { 
	      jmpf <- (data[apairs[,2],timestatus[2]]==1)
	      apairs <- apairs[data[apairs[,2],timestatus[2]]==1,]
              attr(pairs,"jump-first") <- jmpf
      }
      pairs <- apairs
      return(pairs)
} # }}}

##' @export
alpha2spear <- function(theta,link=1) { ## {{{ 
   if (link==1) theta <- exp(theta)
if (length(theta)>1) {
   out <- c()
   for (thet in theta) {
   if (thet!=1) out <- c(out,( (thet+1)/(thet-1) -2* thet* log(thet)/ (thet-1)^2))
   else out <- c(out,0)
   }
} else { if (theta!=1) out <- ( (theta+1)/(theta-1) -2* theta* log(theta)/ (theta-1)^2) }

return(out)
} ## }}} 

##' @export
alpha2kendall <- function(theta,link=0) {  ## {{{ 
   if (link==1) theta <- exp(theta)
   return(1/(1+2/theta)) 
} ## }}} 

##' @export piecewise.twostage
piecewise.twostage <- function(cut1,cut2,data=sys.parent(),timevar="time",status="status",id="id",covars=NULL,covars.pairs=NULL,num=NULL,
            score.method="optimize",Nit=100,detail=0,silent=1,weights=NULL,
            control=list(),theta=NULL,theta.des=NULL,var.link=1,iid=1,step=0.5,model="plackett",data.return=0)
{ ## {{{
ud <- list()
if (missing(cut2)) cut2 <- cut1; 
nc1 <- length(cut1); nc2 <- length(cut2)
names1 <- names2 <- c()
theta.mat <- se.theta.mat <- cor.mat <- score.mat <- se.cor.mat <- matrix(0,nc1-1,nc2-1); 
clusters <- data[,id]
cluster.call <- clusters
idi <- unique(data[,id]); 
###print(head(idi))

## {{{ 
###   se.clusters=NULL,max.clust=1000,
###  evt saette cluster se max.clust paa 
###  if (is.null(se.clusters)) { se.clusters <- clusters; antiid <- nrow(clusterindex);} else  {
###      iids <-  unique(seclusters); 
###      antiid <- length(iids); 
###      if (is.numeric(seclusters)) se.clusters <-  fast.approx(iids,se.clusters)-1
###       else se.clusters <- as.integer(factor(se.clusters, labels = seq(antiid)))-1
###  }
###  if (length(se.clusters)!=length(clusters)) stop("Length of seclusters and clusters must be same\n"); 
###
###  if ((!is.null(max.clust))) if (max.clust< antiid) {
###        coarse.clust <- TRUE
###	qq <- unique(quantile(se.clusters, probs = seq(0, 1, by = 1/max.clust)))
###	qqc <- cut(se.clusters, breaks = qq, include.lowest = TRUE)    
###	se.clusters <- as.integer(qqc)-1
###	max.clusters <- length(unique(se.clusters))
###	maxclust <- max.clust    
###	antiid <- max.clusters
###  }                                                         
## }}} 

if (iid==1) { theta.iid <- matrix(0,length(idi),(nc1-1)*(nc2-1));
              rownames(theta.iid) <- idi
            } else theta.iid <- NULL

thetal <- c()
k <- 0; 
for (i1 in 2:nc1)
for (i2 in 2:nc2)
{
k <-(i1-2)*(nc2-1)+(i2-1)
if (silent<=0) cat(paste("Data-set ",k,"out of ",(nc1-1)*(nc2-1)),"\n");
datalr <- surv.boxarea(c(cut1[i1-1],cut2[i2-1]),c(cut1[i1],cut2[i2]),data,timevar=timevar,
			status=status,id=id,covars=covars,covars.pairs=covars.pairs,num=num,silent=silent) 
if (silent<=-1) print("back in piecewise.twostage"); 
if (silent<=-1) print(summary(datalr)); 
if (silent<=-1) print(head(datalr)); 
if (silent<=-1) print(summary(datalr[,id])); 
 boxlr <- list(left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]))
datalr$tstime <- datalr[,timevar]
datalr$tsstatus <- datalr[,status]
datalr$tsid <- datalr[,id]
###

if (is.null(covars)) 
f <- as.formula(with(attributes(datalr),paste("Surv(",time,",",status,")~-1+factor(",num,")")))
else f <- as.formula(with(attributes(datalr),paste("Surv(",time,",",status,")~-1+factor(",num,"):",covars)))
marg1 <- aalen(f,data=datalr,n.sim=0,robust=0)

fitlr<-  survival.twostage(marg1,data=datalr,clusters=datalr$tsid,model=model,score.method=score.method,
             Nit=Nit,detail=detail,silent=silent,weights=weights,
             control=control,theta=theta,theta.des=theta.des,var.link=var.link,iid=iid,step=step)
####
coef <- coef(fitlr)
theta.mat[i1-1,i2-1] <- fitlr$theta
se.theta.mat[i1-1,i2-1] <- fitlr$var.theta^.5
cor.mat[i1-1,i2-1] <- coef[1,5]
se.cor.mat[i1-1,i2-1] <- coef[1,6]
score.mat[i1-1,i2-1] <- fitlr$score
if (data.return==0) 
ud[[k]] <- list(index=c(i1,i2),left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]),fitlr=fitlr)
if (data.return==1) 
ud[[k]] <- list(index=c(i1,i2),left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]),fitlr=fitlr,data=datalr)
if (i2==2) names1 <- c(names1, paste(cut1[i1-1],"-",cut1[i1]))
if (i1==2) names2 <- c(names2, paste(cut2[i2-1],"-",cut2[i2]))
thetal <- c(thetal,fitlr$theta)

if ((silent<=-1) & (iid==1)) print(head(fitlr$theta.iid)); 
if ((silent<=-1) & (iid==1)) {
print(idi) ; print(datalr$tsid)
print(dim(fitlr$theta.iid))
print(head(fitlr$theta.iid))
print(dim(theta.iid))
print(length( idi %in% unique(datalr$tsid)))
}
if (iid==1) theta.iid[idi %in% unique(datalr$tsid),k] <-c(fitlr$theta.iid) 
###if (iid==1) theta.iid[rownames(fitlr$theta.iid),k] <-  fitlr$theta.iid 
}

var.thetal <- NULL
if (iid==1)  var.thetal <- t(theta.iid) %*% theta.iid

colnames(score.mat) <- colnames(cor.mat) <-  colnames(se.cor.mat)  <- colnames(se.theta.mat) <- colnames(theta.mat) <- names1; 
rownames(score.mat) <- rownames(cor.mat) <-  rownames(se.cor.mat) <-  rownames(se.theta.mat) <- rownames(theta.mat) <- names2; 

ud <- list(model.fits=ud,theta=theta.mat,var.theta=se.theta.mat^2,
	   se.theta=se.theta.mat,thetal=thetal,thetal.iid=theta.iid,var.thetal=var.thetal,model=model,
	   cor=cor.mat,se.cor=se.cor.mat,score=score.mat); 
class(ud)<-"pc.twostage" 
attr(ud,"var.link")<-var.link; 
attr(ud, "Type") <- model
return(ud);
} ## }}}

##' @export piecewise.data
piecewise.data <- function(cut1,cut2,data=sys.parent(),timevar="time",status="status",id="id",covars=NULL,covars.pairs=NULL,num=NULL,silent=1)
{ ## {{{
ud <- list()
if (missing(cut2)) cut2 <- cut1; 
nc1 <- length(cut1); nc2 <- length(cut2)
dataud <- c()

k <- 0; 
for (i1 in 2:nc1)
for (i2 in 2:nc2)
{
k <-(i1-2)*(nc2-1)+(i2-1)
if (silent<=0) cat(paste("Data-set ",k,"out of ",(nc1-1)*(nc2-1)),"\n"); 
 datalr <- surv.boxarea(c(cut1[i1-1],cut2[i2-1]),c(cut1[i1],cut2[i2]),data,timevar=timevar,
			status=status,id=id,covars=covars,covars.pairs=covars.pairs,num=num,silent=silent) 
if (silent<=-1) print(summary(datalr)); 
if (silent<=-1) print(head(datalr)); 
datalr$tstime <- datalr[,timevar]
datalr$tsstatus <- datalr[,status]
datalr$tsid <- datalr[,id]
###
datalr$strata <- paste( c(cut1[i1-1],cut2[i2-1]),c(cut1[i1],cut2[i2]),collapse=",",sep="-")
datalr$intstrata <- 
c(paste(c(cut1[i1-1],cut1[i1]),collapse=",",sep="-"),paste( c(cut2[i2-1],cut2[i2]),collapse=",",sep="-"))

if (silent<=-1) print(head(datalr)); 
dataud <- rbind(dataud,datalr)
}

return(data.frame(dataud))
} ## }}}

##' @export
summary.pc.twostage <- function(object,var.link=NULL,...)
{ ## {{{
  if (!(inherits(object,"pc.twostage"))) stop("Must be a Piecewise constant two-Stage object")
  
  res <- list(estimates=object$theta,se=object$se.theta,cor=object$cor,se.cor=object$se.cor,
	      model=object$model,score=object$score)
  class(res) <- "summary.pc.twostage"
  attr(res,"var.link")<-attr(object,"var.link"); 
  attr(res, "Type") <- object$model
  res
} ## }}}

##' @export
print.pc.twostage <- function(x,var.link=NULL,...)
{ ## {{{
   if (!(inherits(x,"pc.twostage"))) stop("Must be a Piecewise constant two-Stage object")
   print( summary(x,var.link=var.link,...))
} ## }}}

##' @export
print.summary.pc.twostage <- function(x,var.link=NULL, digits=3,...)
{ ## {{{
  
  if (is.null(var.link)) { if (attr(x,"var.link")==1) vlink <- 1 else vlink <- 0; } else vlink <- var.link
  print(vlink)

  if (x$model=="plackett") cat("Dependence parameter for Plackett model \n"); 
  if (x$model=="clayton.oakes") cat("Dependence parameter for Clayton-Oakes model \n"); 
 
  if (max(x$score)>0.001) { cat("Score of log-likelihood for parameter estimates (too large?)\n"); print(x$score);cat("\n\n");}

  if (vlink==1) cat("log-coefficient for dependence parameter (SE) \n")  else cat("Dependence parameter (SE) \n");
  print(coefmat(x$estimate,x$se,digits=digits,...))
  cat("\n") 

  if (x$model=="plackett") {cat("Spearman Correlation (SE) \n");cor.type <- "Spearman Correlation"; }
  if (x$model=="clayton.oakes") {cat("Kendall's tau (SE) \n"); cor.type <- "Kendall's tau";}

  print(coefmat(x$cor,x$se.cor,digits,...))
  cat("\n") 
} ## }}}

##' @export
coefmat <- function(est,stderr,digits=3,...) { ## {{{
  myest <- round(10^digits*(est))/10^digits;
  myest <- paste(ifelse(myest<0,""," "),myest,sep="")
  mysd <- round(10^digits*(stderr))/10^digits;  
  res <- matrix(paste(format(myest)," (",format(mysd),")",sep=""),ncol=ncol(est))
  dimnames(res) <- dimnames(est)
  colnames(res) <- paste("",colnames(res))
  noquote(res)
} ## }}}

##' Wrapper for easy fitting of Clayton-Oakes or bivariate Plackett models for bivariate survival data 
##'
##' Fits two-stage model for describing depdendence in survival data
##' using marginals that are on cox or aalen form using the twostage funcion, but
##' call is different and easier and the data manipulation  build into the function.
##' Useful in particular for family design data. 
##'
##' If clusters contain more than two times, the algoritm uses a composite likelihood
##' based on the pairwise bivariate models.
##'
##' The reported standard errors are based on the estimated information from the 
##' likelihood assuming that the marginals are known. 
##'
##' @examples
##' library(mets)
##' data("prt",package="mets")
##' margp <- coxph(Surv(time,status==1)~factor(country),data=prt)
##' fitco <- survival.twostage(margp,data=prt,clusters=prt$id)
##' summary(fitco)
##' 
##' des <- model.matrix(~-1+factor(zyg),data=prt); 
##' fitco <- survival.twostage(margp,data=prt,theta.des=des,clusters=prt$id)
##' summary(fitco)
##' 
##' dfam <- simSurvFam(1000)
##' dfam <- fast.reshape(dfam,var=c("x","time","status"))
##' 
##' desfs <- function(x,num1="num1",num2="num2")
##' { 
##' pp <- (x[num1]=="m")*(x[num2]=="f")*1   ## mother-father 
##' pc <- (x[num1]=="m" | x[num1]=="f")*(x[num2]=="b1" | x[num2]=="b2")*1 ## mother-child
##' cc <- (x[num1]=="b1")*(x[num2]=="b1" | x[num2]=="b2")*1               ## child-child
##' c(pp,pc,cc)
##' } 
##' 
##' marg <- coxph(Surv(time,status)~factor(num),data=dfam)
##' out3 <- easy.survival.twostage(marg,data=dfam,time="time",status="status",id="id",
##'              deshelp=0,
##'              score.method="fisher.scoring",theta.formula=desfs,
##'              model="plackett",
##'              desnames=c("parent-parent","parent-child","child-cild"),iid=1)
##' summary(out3)
##' 
##' @keywords survival twostage 
##' @export easy.survival.twostage
##' @param margsurv model 
##' @param data data frame
##' @param score.method Scoring method
##' @param status Status at exit time
##' @param time Exit time
##' @param entry Entry time
##' @param id name of cluster variable in data frame
##' @param Nit Number of iterations
##' @param detail Detail for more output for iterations 
##' @param silent Debug information
##' @param weights Weights for log-likelihood, can be used for each type of outcome in 2x2 tables. 
##' @param control Optimization arguments
##' @param theta Starting values for variance components
##' @param theta.formula design for depedence, either formula or design function
##' @param desnames names for dependence parameters
##' @param deshelp if 1 then prints out some data sets that are used, on on which the design function operates
##' @param var.link Link function for variance (exp link)
##' @param iid Calculate i.i.d. decomposition
##' @param step Step size for newton-raphson
##' @param model plackett or clayton-oakes model
##' @param marginal.surv vector of marginal survival probabilities 
##' @param strata strata for fitting 
##' @param se.clusters clusters for iid decomposition for roubst standard errors
easy.survival.twostage <- function(margsurv=NULL,data=sys.parent(),score.method="nlminb",
status="status",time="time",entry=NULL,id="id", Nit=60,detail=0, silent=1,weights=NULL, control=list(),
theta=NULL,theta.formula=NULL,desnames=NULL,deshelp=0,var.link=1,iid=1,
step=0.5,model="plackett",marginal.surv=NULL,strata=NULL,se.clusters=NULL)
{ ## {{{
### marginal trunction probabilty, to be computed from model 
pentry <- NULL

if (is.null(marginal.surv))
if (class(margsurv)[1]=="coxph")
{ ## {{{
###    ps <- survfit(margsurv)$surv
    coxformula <- margsurv$formula
    X <- model.matrix(coxformula,data=data)[,-1]; 
    baseout <- survival::basehaz(margsurv,centered=FALSE); 
    baseout <- cbind(baseout$time,baseout$hazard)
    cumh <-  Cpred(baseout,data[,time])[,2]
    RR<-exp(X %*% coef(margsurv))
    ps<-exp(-cumh*RR)
    ## }}}
  } else if (class(margsurv)[1]=="phreg")
  {  ## {{{
	if (!is.null(margsurv$coef)) 
		rr <- c(exp(margsurv$X  %*% margsurv$coef))
	else rr <- rep(1,nrow(margsurv$X))
        ps <- exp(-rr * Cpred(margsurv$cumhaz,margsurv$exit)[,2])
  } ## }}} 
  else stop("marginal survival probabilities must be given as marginal.sur or margsurv \n"); 

  data <- cbind(data,ps)
  if (!is.null(pentry)) data <- cbind(data,pentry)

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
    des.theta  <- t( apply(data.fam.clust,1,desfunction)) 
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
     if (is.null(pentry))
    data.fam <- fast.reshape(data.fam.clust,varying=c(id,"ps",status))
    else data.fam <- fast.reshape(data.fam.clust,varying=c(id,"ps",status,"pentry"))
    if (deshelp==1) {
	cat("Back to long format for twostage (head)\n"); 
        print(head(data.fam)); 
	cat("\n")
###	cat(paste("twostage, called with reponse",response,"\n")); 
	cat(paste("cluster=",id,",  subcluster (pairs)=subfam \n")); 
	cat(paste("design variables =")); 
	cat(desnames)
	cat("\n")
    } 

    if (is.null(pentry)) ptrunc <- NULL else ptrunc <- data.fam[,pentry]

    out <- survival.twostage(NULL,data=data.fam,
                    clusters=data.fam$subfam,
		    theta.des=as.matrix(data.fam[,desnames]),
                    detail=detail, score.method=score.method, Nit=Nit,step=step,
                    iid=iid,theta=theta, var.link=var.link,model=model, 
                    marginal.survival=data.fam[,"ps"],
                    marginal.status=data.fam[,status],
		    marginal.trunc=ptrunc,
		    se.clusters=data.fam[,id])

   return(out)
} ## }}}

##' @export
simSurvFam <- function(n,beta=0.0,theta=1,lam0=0.5,lam1=1,lam2=1,ctime=10,...) { ## {{{ 
###	n=10; beta=0; theta=1; lam1=1;lam2=1; ctime=10; lam0=0.5
xm <- rbinom(n,1,0.5); xf <- rbinom(n,1,0.5); 
xb1 <- rbinom(n,1,0.5); xb2 <- rbinom(n,1,0.5); 
###
zf <- rgamma(n,shape=lam1); zb <- rgamma(n,shape=lam2); 
tm <- rexp(n)/(zf*exp(xm*beta)*lam0)
tf <- rexp(n)/(zf*exp(xf*beta)*lam0)
tb1 <- rexp(n)/((zf+zb)*exp(xb1*beta)*2*lam0)
tb2 <- rexp(n)/((zf+zb)*exp(xb2*beta)*2*lam0)
cm <- ifelse(tm<ctime,1,0); cf <- ifelse(tf<ctime,1,0); 
cb1 <- ifelse(tb1<ctime,1,0); cb2 <- ifelse(tb2<ctime,1,0); 
tm <- ifelse(tm<ctime,tm,ctime); tf <- ifelse(tf<ctime,tf,ctime)
tb1 <- ifelse(tb1<ctime,tb1,ctime); tb2 <- ifelse(tb2<ctime,tb2,ctime)
#
data.frame(xm=xm,xf=xf,xb1=xb1,xb2=xb2,timem=tm,timef=tf,timeb1=tb1,timeb2=tb2,statusm=cm,statusf=cf,
	   statusb1=cb1,statusb2=cb2,id=1:n)
} ## }}} 

##' @export
object.defined <- function(object)
{ 
   exists(as.character(substitute(object)))
}

##' @export survival.twostage.fullse
survival.twostage.fullse <- function(margsurv,data=sys.parent(),
   score.method="fisher.scoring",Nit=60,detail=0,clusters=NULL,
   silent=1,weights=NULL, control=list(),theta=NULL,
   theta.des=NULL,var.link=1,iid=1,
   step=0.5,model="clayton.oakes",
   marginal.trunc=NULL, marginal.survival=NULL,
   marginal.status=NULL,strata=NULL,
   se.clusters=NULL,numDeriv=0,
   random.design=NULL,fdetail=1,...)
{ ## {{{ 
  if (is.null(margsurv$gamma.iid)) stop("Call marginal model with resample.iid=1, only Cox model via cox.aalen \n"); 
  beta.iid <- margsurv$gamma.iid
  base.iid <- margsurv$B.iid

  ### takes cluster grouping of marginal models for se
  se.clusters <- attr(margsurv,"cluster")

  udtwo <- survival.twostage(margsurv,data=data,
  score.method=score.method,Nit=Nit,detail=detail,
  clusters=clusters,
  silent=silent,weights=weights,control=control,
  theta=theta,theta.des=theta.des,
  var.link=var.link,iid=iid,step=step,
  model=model,marginal.trunc=marginal.trunc,
  marginal.survival=marginal.survival,
  marginal.status=marginal.status,strata=strata,
  se.clusters=se.clusters,numDeriv=numDeriv,random.design=random.design,...)
  par <- margsurv$gamma
  theta <- udtwo$theta
###
 
 if (object.defined(margsurv$time.sim.resolution))
 {
 parbase <- Cpred(margsurv$cum,margsurv$time.sim.resolution) 
 margsurv$cum <- parbase
 parbase <- parbase[,2]
 } else parbase <- margsurv$cum[,2]

twobeta  <- function(par,beta=1)
{ ## {{{ 
   if (beta==1) margsurv$gamma <- par else margsurv$cum[,2] <- par

   udl <- survival.twostage(margsurv,data=data,
                    score.method=score.method,Nit=0,detail=detail,clusters=clusters,
		    silent=silent,weights=weights, control=control,theta=theta,theta.des=theta.des,
		    var.link=var.link,iid=0,
                    step=step,model=model,
		    marginal.trunc=marginal.trunc,marginal.survival=marginal.survival,
		    marginal.status=marginal.status,strata=strata,
		    se.clusters=se.clusters,numDeriv=0,random.design=random.design,...)
###  udl <- two.stage(margsurv,data=data,Nit=1,clusters=clusters,theta=theta) ###  udl$theta.score
  udl$score
} ## }}} 

if (fdetail==1) cat("Ready for numDeriv wrt beta and baseline\n")
DUbeta <-  numDeriv::jacobian(twobeta,par,beta=1) #,method="complex")
DUbase <-  numDeriv::jacobian(twobeta,parbase,beta=0) #,method="complex")
if (fdetail==1) cat("Finished numDeriv wrt beta and baseline\n")

biid <- c()
for (i in 1:length(base.iid)) biid <- cbind(biid,base.iid[[i]])
					 
IDUbeta <-  udtwo$hessi %*% DUbeta
IDUbase <-  udtwo$hessi %*% DUbase
betaiid <-t( IDUbeta %*% t(beta.iid))
baseiid <- t( IDUbase %*% biid)
###
iidfull <- udtwo$theta.iid
iidfull <- iidfull+betaiid+baseiid
var2 <- t(iidfull) %*% iidfull
se <- cbind(diag(var2)^.5); colnames(se) <- "se"

se.naive=coef(udtwo)[,2,drop=FALSE]; 
colnames(se.naive) <- "se.naive"
###
res <- list(theta.iid=udtwo$theta.iid, iid=iidfull,coef=udtwo$theta,var=var2,se=se, se.naive=se.naive)
attr(res,"DUbeta") <- DUbeta; 
attr(res,"DUbase") <- DUbase; 
attr(res,"DUthetainv") <- udtwo$hessi

class(res) <- "survival.twostage.fullse"
return(res)
} ## }}} 

##' @export
summary.survival.twostage.fullse <- function(object,digits=3,...)
{ ## {{{ 
	tval <- object$coef/object$se
	pval <- 2*(1-pnorm(abs(tval)))
       	res <- cbind(object$coef,object$se,object$se.naive,pval)
	colnames(res) <- c("coef","se","se.naive","pval")
return(res)
} ## }}} 

##' @export
coef.survival.twostage.fullse <- function(object,digits=3,...)
{ ## {{{ 
summary.survival.twostage.fullse(object)
} ## }}} 

##' @export
print.survival.twostage.fullse  <-  function(x,...)
{  ## {{{ 
summary.survival.twostage.fullse(x)
} ## }}} 

##' @export
twin.polygen.design <-function (data,id="id",zyg="DZ",zygname="zyg",type="ace",tv=NULL,...) { ## {{{
  ### twin case 
  nid <- table(data[,id])
  id <- data[,id]
  tv <- diff(c(NA,id))
  tv[tv!=0 | is.na(tv)] <- 1
  tv[tv==0] <- 2

  zygbin <- (data[,zygname]==zyg)*1
  zygdes=model.matrix(~-1+factor(zygbin),data)
  n <- length(zygbin)

  if (type=="ace") { ### ace ## {{{ 
  ### random effects for each cluster
  DZns <- cbind((zygbin==1)*(tv==1)*cbind(rep(1,n),rep(0,n))+
		(zygbin==1)*(tv==2)*cbind(rep(0,n),rep(1,n)))
  des.rv <- cbind(zygdes,DZns,1)
  colnames(des.rv) <- c("MZ","DZ","DZns1","DZns2","env")
  pard <- rbind(c(1,0), c(0.5,0),c(0.5,0), c(0.5,0), c(0,1))
  } ## }}} 

  if (type=="ae") { ### ae ## {{{ 
  ###AE model 
  DZns <- cbind((zygbin==1)*(tv==1)*cbind(rep(1,n),rep(0,n))+
		(zygbin==1)*(tv==2)*cbind(rep(0,n),rep(1,n)))
  des.rv <- cbind(zygdes,DZns)
  colnames(des.rv) <- c("MZ","DZ","DZns1","DZns2")
  pard <- rbind(c(1,0), c(0.5,0),c(0.5,0), c(0.5,0))[,1,drop=FALSE]
  } ## }}} 

  if (type=="dce") { ### dce ## {{{ 
  ### DCE  
  ### random effects for each cluster
  DZns <- cbind((zygbin==1)*(tv==1)*cbind(rep(1,n),rep(0,n))+(zygbin==1)*(tv==2)*cbind(rep(0,n),rep(1,n)))
  des.rv <- cbind(zygdes,DZns,1)
  colnames(des.rv) <- c("MZ","DZ","DZns1","DZns2","env")
  pard <- rbind(c(1,0), c(0.25,0),c(0.75,0), c(0.75,0), c(0,1))
  } ## }}} 

  if (type=="ade") { ### ade ## {{{ 
  #ADE
  DZns <- cbind((zygbin==1)*(tv==1)*cbind(rep(1,n),rep(0,n))+(zygbin==1)*(tv==2)*cbind(rep(0,n),rep(1,n)))
  des.rv <- cbind(zygdes,DZns,DZns,1)
  pard <- rbind(c(1,0),c(0.25,0),c(0.75,0),c(0.75,0),c(0,1),c(0,0.5),c(0,0.5),c(0,0.5) )
  pardes <- matrix(pard,n,16,byrow=TRUE)
  des.rv <- NULL
  } ## }}} 

  if (type=="adce") { ### adce ## {{{ 
###zygdes=model.matrix(~-1+factor(zygbin),prtwomen)
###n <- nrow(prtwomen)
###des.rv <- cbind( zygdes[,c(2,1)], (prtwomen$zygbin==0)*(prtwomen$tv==1), 
###	(prtwomen$zygbin==0)*(prtwomen$tv==2),
###	zygdes[,c(2,1)], (prtwomen$zygbin==0)*(prtwomen$tv==1), 
###	(prtwomen$zygbin==0)*(prtwomen$tv==2),1
###	)
###pard <- rbind(c(1,0,0), c(0.25,0,0),c(0.75,0,0), c(0.75,0,0), 
###	      c(0,1,0), c(0,0.5,0),c(0,0.5,0), c(0,0.5,0) ,c(0,0,1))
###pardes <- matrix(pard,n,27,byrow=TRUE)
  } ## }}} 

  if (type=="de") { ### ae ## {{{ 
  DZns <- cbind((zygbin==1)*(tv==1)*cbind(rep(1,n),rep(0,n))+
		(zygbin==1)*(tv==2)*cbind(rep(0,n),rep(1,n)))
  des.rv <- cbind(zygdes,DZns)
  colnames(des.rv) <- c("MZ","DZ","DZns1","DZns2")
  pard <- rbind(c(1,0), c(0.25,0),c(0.75,0), c(0.75,0))[,1,drop=FALSE]
  } ## }}} 

  if (type=="un") { ### ae ## {{{ 
  des.rv <- cbind(zygdes)
  colnames(des.rv) <- c("MZ","DZ")
  pard <- rbind(c(1,0),c(0,1))
  } ## }}} 

res <- list(pardes=pard,des.rv=des.rv)
return(res)
} ## }}}

##' @export
make.pairwise.design  <- function(pairs,kinship,type="ace")
{ ## {{{ 
### makes pairwise random effects design for shared and non-shared random effects
### kinship gives shared genes for each pair

if (type=="ace" | type=="ade") {
theta.des  <- array(0,c(4,2,nrow(pairs)))
random.des <- array(0,c(2,4,nrow(pairs)))
rvs <- c()
}
if (type=="ae") {
theta.des  <- array(0,c(3,1,nrow(pairs)))
random.des <- array(0,c(2,3,nrow(pairs)))
rvs <- c()
}
if (type=="simple") {
theta.des  <- array(1,c(1,1,nrow(pairs)))
random.des <- array(1,c(2,1,nrow(pairs)))
rvs <- rep(1,nrow(pairs))
}


for (i in 1:nrow(pairs))
{ 
	if (type=="ace" | type=="ade") {
         ### only 3 random variables for ace 
         ### (gene, shared, non-shared, environment
         ### kinship gives amount of genes shared 
	if (type=="ace") 
	 theta.des[,,i] <- rbind(c(kinship[i],0),
				 c(1-kinship[i],0),
				 c(1-kinship[i],0),
				 c(0,1))
	if (type=="ade") 
	theta.des[,,i] <- rbind(c((kinship[i]==1)+(kinship[i]!=1)*kinship[i]*0.5,0),
				c(1-(kinship[i]==1)+(kinship[i]!=1)*kinship[i]*0.5,0),
				c(1-(kinship[i]==1)+(kinship[i]!=1)*kinship[i]*0.5,0),
				c(0,1))
       	 random.des[,,i] <- rbind(c(1,1,0,1),c(1,0,1,1))
	 rvs <- c(rvs,4)
	} 
        if (type=="ae" | type=="ad") {
         ### only 2 random variables for ace 
         ### (gene, shared, not-shared  
         ### kinship gives amount of genes shared 
	if (type=="ae") 
	 theta.des[,,i] <- matrix(c(kinship[i],1-kinship[i],1-kinship[i]),nrow=3,ncol=1)
        if (type=="ad") 
	theta.des[,,i] <- matrix(c((kinship[i]==1)+(kinship[i]!=1)*kinship[i]*0.5),
				 c(1-(kinship[i]==1)+(kinship[i]!=1)*kinship[i]*0.5),
				 c(1-(kinship[i]==1)+(kinship[i]!=1)*kinship[i]*0.5),3,1)
       	 random.des[,,i] <- rbind(c(1,1,0),c(1,0,1))
	 rvs <- c(rvs,3)
	}
} 

return(list(random.design=random.des,theta.des=theta.des,ant.rvs=rvs))
} ## }}} 

##' @export
ace.family.design <-function (data,id="id",member="type",mother="mother",father="father",child="child",child1="child",type="ace",...) { 
## {{{
  ### standard family case 
  nid <- table(data[,id])
  id <- data[,id]
###  tv <- diff(c(NA,id))
###  tv[tv!=0 | is.na(tv)] <- 1
###  tv[tv==0] <- 2

###  zygbin <- (data[,zygname]==zyg)*1
###  zygdes=model.matrix(~-1+factor(zygbin),data)
###  n <- length(zygbin)

  if (type=="ace") { ### ace ## {{{ 
  ### random effects for each cluster
	  mom <- 1*(data[,member]==mother)
	  mo <- cbind(mom*1,1,1,1)*(mom)
	  fa <- (data[,member]==father)
	  fad <- cbind(fa,1,1,1)*fa
	  ch1 <- (data[,member]==child)*(data[,child1]==1)
	  cc1 <- cbind(ch1,1,0,0,1,1,0,0)*ch1
	  ch2 <- (data[,member]==child)*(data[,child1]==0)
	  cc2 <- cbind(ch2,0,1,0,1,0,1,0)*ch2
	  des.rv <- cbind(cbind(mo,fad)+cc1+cc2,1)

	  colnames(des.rv) <- c("m1","m2","m3","m4","f1","f2","f3","f4","env")
	  pard <- rbind( c(0.25,0), c(0.25,0),c(0.25,0), c(0.25,0), 
			 c(0.25,0), c(0.25,0),c(0.25,0), c(0.25,0), c(0,1))
  } ## }}} 

  if (type=="ae") { ### ae ## {{{ 
  ###AE model 
          mom <- 1*(data[,member]==mother)
	  mo <- cbind(mom*1,1,1,1)*(mom)
	  fa <- (data[,member]==father)
	  fad <- cbind(fa,1,1,1)*fa
	  ch1 <- (data[,member]==child)*(data[,child1]==1)
	  cc1 <- cbind(ch1,1,0,0,1,1,0,0)*ch1
	  ch2 <- (data[,member]==child)*(data[,child1]==0)
	  cc2 <- cbind(ch2,0,1,0,1,0,1,0)*ch2
	  des.rv <- cbind(cbind(mo,fad)+cc1+cc2)
	  colnames(des.rv) <- c("m1","m2","m3","m4","f1","f2","f3","f4")
	  pard <- rbind( c(0.25), c(0.25),c(0.25), c(0.25), 
			 c(0.25), c(0.25),c(0.25), c(0.25))
  } ## }}} 

  if (type=="dce") { ### dce ## {{{ 
  ### DCE  
  ### random effects for each cluster
	  stop("not done yet"); 
          mom <- 1*(data[,member]==mother)
	  mo <- cbind(mom*1,1,1,1)*(mom)
	  fa <- (data[,member]==father)
	  fad <- cbind(fa,1,1,1)*fa
	  ch1 <- (data[,member]==child)*(data[,child1]==1)
	  cc1 <- cbind(ch1,1,0,0,1,1,0,0)*ch1
	  ch2 <- (data[,member]==child)*(data[,child1]==0)
	  cc2 <- cbind(ch2,0,1,0,1,0,1,0)*ch2
	  des.rv <- cbind(cbind(mo,fad)+cc1+cc2)
	  colnames(des.rv) <- c("m1","m2","m3","m4","f1","f2","f3","f4")
	  pard <- rbind( c(0.25), c(0.25),c(0.25), c(0.25), 
			 c(0.25), c(0.25),c(0.25), c(0.25))

  } ## }}} 

  if (type=="ade") { ### ade ## {{{ 
  #ADE
	  stop("not done yet"); 
###  pard <- rbind(c(1,0),c(0.25,0),c(0.75,0),c(0.75,0),c(0,1),c(0,0.5),c(0,0.5),c(0,0.5) )
###  pardes <- matrix(pard,n,16,byrow=TRUE)
  } ## }}} 

  if (type=="adce") { ### adce ## {{{ 
	  stop("not done yet"); 
  } ## }}} 

  if (type=="de") { ### ae ## {{{ 
	  stop("not done yet"); 
###  pard <- rbind(c(1,0), c(0.25,0),c(0.75,0), c(0.75,0))[,1,drop=FALSE]
  } ## }}} 

  if (type=="un") { ### ae ## {{{ 
	 stop("not done yet"); 
         pard <- diag(4)
  } ## }}} 

res <- list(pardes=pard,des.rv=des.rv)
return(res)
} ## }}}

##' @export
make.pairwise.design.competing <- function(pairs,kinship,type="ace",compete=2,overall=1)
 ## {{{ 
 {
 ### makes pairwise random effects design for shared and non-shared random effects
 ### kinship gives shared genes for each pair
 ### overall ace + 1 ace , 2 ace (6 pars) 
 
 if (type=="ace") {
    theta.des  <- array(0,c((2*compete+overall*1)*4+overall*2,(compete)*2,nrow(pairs))) 
    random.des <- array(0,c(2*compete,(2*compete+overall*1)*4+overall*2,nrow(pairs))) 
 }
 
 if (type=="simple") {
    theta.des <- array(0,c((3*compete+overall),(compete+overall),nrow(pairs)))
    random.des <- array(0,c(2*compete,(3*compete+overall),nrow(pairs)))
 }
 
 if (type=="ae") {
    theta.des  <- array(0,c(6,1,nrow(pairs)))
    random.des <- array(0,c(2*compete,6,nrow(pairs)))
 }
 
 print(dim(theta.des))
 print(dim(random.des))

 rvs <- c()
 for (i in 1:nrow(pairs))
 { 
 	if (type=="ace") {  ## {{{ 
          ### only 3 random variables for ace 
          ### (gene, shared, non-shared, environment
          ### kinship gives amount of genes shared 
 	if (compete==3) { ## {{{ 
 		stop("to do ")
 	} ## }}} 
 	if (compete==2) { ## {{{ 
 		stop("to do ")
 		k <- kinship
 		if (overall==1)  { ## {{{ 
### 		 random.des[,,i] <- rbind( c(1,1,0,1,rep(0,4),c(1,1,0,0,0,1)),
 		 rd           <- rbind( c(1,1,0,1,rep(0,4),c(1,1,0,0,0,1)),
 					   c(rep(0,4),1,1,0,1,c(1,0,1,0,0,1)),
 					   c(1,0,1,1,rep(0,4),c(1,0,0,1,0,1)),
 					   c(rep(0,4),1,0,1,1,c(1,0,0,0,1,1)))
 			theta.des[,,i] <- matrix(0,14,6)
 			mini <- rbind( c(k[i],  0), c(1-k[i],0), c(1-k[i],0), c(0,     1))
 			minis <- rbind( c(k[i],0), 
 				        c((1-k[i]),0), c((1-k[i]),0),
 				        c((1-k[i]),0), c((1-k[i]),0),
### 				        c(0.5*(1-k[i]),0), c(0.5*(1-k[i]),0),
### 				        c(0.5*(1-k[i]),0), c(0.5*(1-k[i]),0),
 				        c(0,1))
 			theta.des[1:4,1:2,i] <- mini
 			theta.des[5:8,3:4,i] <- mini
 			theta.des[9:14,5:6,i] <- minis
### 		 random.des[,,i] <- rbind( c(1,1,0,1,rep(0,4),c(1,1,0,0,0,1,0,0,0,1)),
### 					   c(rep(0,4),1,1,0,1,c(1,0,1,0,0,0,1,0,0,1)),
### 					   c(1,0,1,1,rep(0,4),c(1,0,0,1,0,0,0,1,0,1)),
### 					   c(rep(0,4),1,0,1,1,c(1,0,0,0,1,0,0,0,1,1)))
 		 rvs <- c(rvs,14)
 		} ## }}} 
 		if (overall==0)  { ## {{{ 
 		 theta.des[,,i] <- rbind( 
 				 c(kinship[i],0,0,0),
 				 c(1-kinship[i],0,0,0),
 				 c(1-kinship[i],0,0,0),
 				 c(0,1,0,0),
 				 c(0,0,kinship[i],0),
 				 c(0,0,1-kinship[i],0),
 				 c(0,0,1-kinship[i],0),
 				 c(0,0,0,1))
 		 random.des[,,i] <- rbind( c(1,1,0,1,rep(0,4)),
 					   c(rep(0,4),1,1,0,1),
 					   c(1,0,1,1,rep(0,4)),
 					   c(rep(0,4),1,0,1,1))
 		 rvs <- c(rvs,8)
 		} ## }}} 
 	} ## }}} 
 	} 
 
        if (type=="simple") {  ## {{{ 
          ### only 3 random variables 
 	if (compete==3) { ## {{{ 
        	if (overall==1)  { ## {{{ 
 		 theta.des[,,i] <- rbind( c(1,0,0),
 		                          c(1,0,0),
 		                          c(1,0,0),
 		                          c(1,0,0),
 		                          c(0,1,0),
 		                          c(0,1,0),
 		                          c(0,1,0),
 		                          c(0,1,0),
 		                          c(0,0,1))
 		 random.des[,,i] <- rbind(c(1,0,0,1,0,0,1),
 					  c(0,1,0,0,1,0,1),
 					  c(0,0,1,0,0,1,1),
					  c(1,0,0,1,0,0,1),
 					  c(0,1,0,0,1,0,1),
 					  c(0,0,1,0,0,1,1))
 		 rvs <- c(rvs,7)
 		} ## }}} 
 		if (overall==0)  { ## {{{ 
 			stop("not done"); 
 		 rvs <- c(rvs,2)
 		} ## }}} 

 	} ## }}} 
 	if (compete==2) { ## {{{ 
 		if (overall==1)  { ## {{{ 
 		 theta.des[,,i] <- rbind( c(1,0,0),
 		                          c(1,0,0),
 		                          c(1,0,0),
 		                          c(0,1,0),
 		                          c(0,1,0),
 		                          c(0,1,0),
 		                          c(0,0,1))
 		 random.des[,,i] <- rbind(c(1,0,0,0,1,0,1),
 					  c(0,1,0,1,0,0,1),
					  c(1,0,0,0,0,1,1),
					  c(0,0,1,1,0,0,1))
 		 rvs <- c(rvs,7)
 		} ## }}} 
 		if (overall==0)  { ## {{{ 
 			stop("not done"); 
 		 rvs <- c(rvs,2)
 		} ## }}} 
 	} ## }}} 
 	}
         if (type=="ae") {
          ### only 2 random variables for ace 
          ### (gene, shared, not-shared  
          ### kinship gives amount of genes shared 
 	 theta.des[,,i] <- matrix(c(kinship[i],1-kinship[i],1-kinship[i]),nrow=3,ncol=1)
        	 random.des[,,i] <- rbind(c(1,1,0),c(1,0,1))
 	 rvs <- c(rvs,3)
 	}
 } 
 
 return(list(random.design=random.des,theta.des=theta.des,ant.rvs=rvs))
 } 
## }}} 
## }}}
## }}}

##' Relative risk for additive gamma model
##'
##' Computes the relative risk for additive gamma model at time 0
##' 
##' @references
##' 
##' Eriksson and Scheike (2015), Additive Gamma frailty models for competing risks data, Biometrics (2015)
##' 
##' @examples
##' lam0 <- c(0.5,0.3)
##' pars <- c(1,1,1,1,0,1)
##' ## genetic random effects, cause1, cause2 and overall 
##' parg <- pars[c(1,3,5)]
##' ## environmental random effects, cause1, cause2 and overall 
##' parc <- pars[c(2,4,6)]
##' 
##' ## simulate competing risks with two causes with hazards 0.5 and 0.3
##' ## ace for each cause, and overall ace 
##' out <- simCompete.twin.ace(10000,parg,parc,0,2,lam0=lam0,overall=1,all.sum=1)
##' 
##' ## setting up design for running the model 
##' mm <- familycluster.index(out$cluster)
##' head(mm$familypairindex,n=10)
##' pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
##' tail(pairs,n=12)
##' #
##' kinship <- (out[pairs[,1],"zyg"]=="MZ")+ (out[pairs[,1],"zyg"]=="DZ")*0.5
##' 
##' # dout <- make.pairwise.design.competing(pairs,kinship,
##' #          type="ace",compete=length(lam0),overall=1)
##' # head(dout$ant.rvs)
##' ## MZ
##' # dim(dout$theta.des)
##' # dout$random.design[,,1]
##' ## DZ
##' # dout$theta.des[,,nrow(pairs)]
##' # dout$random.design[,,nrow(pairs)]
##' #
##' # thetades <- dout$theta.des[,,1]
##' # x <- dout$random.design[,,1]
##' # x
##' ##EVaddGam(rep(1,6),x[1,],x[3,],thetades,matrix(1,18,6))
##' 
##' # thetades <- dout$theta.des[,,nrow(out)/2]
##' # x <- dout$random.design[,,nrow(out)/2]
##' ##EVaddGam(rep(1,6),x[1,],x[4,],thetades,matrix(1,18,6))
##' @author Thomas Scheike 
##' @export
##' @param theta theta 
##' @param x1 x1 
##' @param x2 x2
##' @param thetades thetades
##' @param ags ags
EVaddGam <- function(theta,x1,x2,thetades,ags)
{ ## {{{ 
	pars <- thetades %*% theta
	lamtot <- ags %*% theta

	mvar <- pars/lamtot
	vvar <- pars/lamtot^2
###	print(pars)

	x1mvar <- sum(x1 * mvar)
	x2mvar <- sum(x2 * mvar)
	x1x2vvar <- sum(x1*x2*vvar)
###	print(x1x2vvar)

	list(x1m=x1mvar,mx2=x2mvar,
	     dN=x1x2vvar/x2mvar)
###	     pars=pars,mvar=mvar,vvar=vvar)
} ## }}} 


##' @export
twostage.aalen <- function(object,...) survival.twostage(object,...)

##' @export
twostage.cox.aalen <- function(object,...) survival.twostage(object,...)

##' @export
twostage.coxph <- function(object,...) survival.twostage(object,...)

##' @export
twostage.phreg <- function(object,...) survival.twostage(object,...)
