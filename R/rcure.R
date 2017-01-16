#' Robust cure model
#'
#' @description Fits robust cure model by incorporating a weakly informative prior distribution for uncure probability part in cure models
#' @param formula a formula object for the survival part in cure model. left must be a survival object as returned by the Surv function
#' @param cureform specifies the variables in the uncure probability part in cure model
#' @param offset.s variable(s) with coefficient 1 in PH model or AFT model
#' @param data a a data.frame
#' @param na.action a missing-data filter function. By default na.action = na.omit
#' @param model specifies survival part in cure model, "ph" or "aft"
#' @param link specifies the link in incidence part. The "logit", "probit" or complementary loglog ("cloglog") links are available. By default link = "logit".
#' @param Var If it is TRUE, the program returns Std.Error by bootstrap method. If set to False, the program only returns estimators of coefficients. By default, Var = TRUE
#' @param emmax specifies the maximum iteration number. If the convergence criterion is not met, the EM iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default emmax = 100.
#' @param eps sets the convergence criterion. The default is eps = 1e-7. The iterations are considered to be converged when the maximum relative change in the parameters and likelihood estimates between iterations is less than the value specified.
#' @param nboot specifies the number of bootstrap sampling. The default nboot = 100
#' @param family a description of the error distribution and link function to be used in the model. Default is binomial(link="logit")
#' @param method the method to be used in fitting the glmbayes model. The default method "glm.fit" uses iteratively reweighted least squares (IWLS). The only current alternative is "model.frame" which returns the model frame and does no fitting
#' @param prior.mean prior mean for the coefficients: default is 0. Can be a vector of length equal to the number of predictors (not counting the intercept, if any). If it is a scalar, it is expanded to the length of this vector.
#' @param prior.scale prior scale for the coefficients: default is NULL; if is NULL, for a logit model, prior.scale is 2.5; for a probit model, prior scale is 2.5*1.6. Can be a vector of length equal to the number of predictors (not counting the intercept, if any). If it is a scalar, it is expanded to the length of this vector
#' @param prior.df prior degrees of freedom for the coefficients. For t distribution: default is 1 (Cauchy). Set to Inf to get normal prior distributions. Can be a vector of length equal to the number of predictors (not counting the intercept, if any). If it is a scalar, it is expanded to the length of this vector
#' @param prior.mean.for.intercept prior mean for the intercept: default is 0. 
#' @param prior.scale.for.intercept prior scale for the intercept: default is NULL; for a logit model, prior scale for intercept is 10; for probit model, prior scale for intercept is rescaled as 10*1.6.
#' @param prior.df.for.intercept prior degrees of freedom for the intercept: default is 1.
#' @param min.prior.scale Minimum prior scale for the coefficients: default is 1e-12.
#' @param scaled scaled=TRUE, the scales for the prior distributions of the coefficients are determined as follows: For a predictor with only one value, we just use prior.scale. For a predictor with two values, we use prior.scale/range(x). For a predictor with more than two values, we use prior.scale/(2*sd(x)). If the response is Gaussian, prior.scale is also multiplied by 2 * sd(y). Default is TRUE
#' @param Warning default is TRUE, which will show the error messages of not convergence and separation
#' @param n.iter integer giving the maximal number of bayesglm IWLS iterations, default is 100.
#' @param ... further arguments passed to or from other methods
#' @author Xiaoxia Han
#' @references Cai, C., Zou, Y., Peng, Y., & Zhang, J. (2012). smcure: An R-Package for estimating semiparametric mixture cure models. Computer methods and programs in biomedicine, 108(3), 1255-1260.
#' @references Gelman, A., Jakulin, A., Pittau, M. G., & Su, Y. S. (2008). A weakly informative default prior distribution for logistic and other regression models. The Annals of Applied Statistics, 1360-1383.
#' @examples
#' library(survival)
#' library(smcure)
#' library(arm)
#' data(e1684)
#' 
#' # fit PH robust cure model
#' pd <- rcure(Surv(FAILTIME,FAILCENS)~TRT+SEX+AGE,cureform=~TRT+SEX+AGE,
#' data=e1684,model="ph",Var = FALSE,
#' method = "glm.fit", prior.mean = 0, prior.scale = NULL, prior.df = 1,
#' prior.mean.for.intercept = 0, prior.scale.for.intercept = NULL,
#' prior.df.for.intercept = 1, min.prior.scale = 1e-12,
#' scaled = TRUE, n.iter = 100, Warning=F)
#' printrcure(pd,Var = FALSE)
#' # plot predicted survival curves for male with median centered age by treatment groups
#' predm=predictrcure(pd,newX=cbind(c(1,0),c(0,0),c(0.579,0.579)),
#' newZ=cbind(c(1,0),c(0,0),c(0.579,0.579)),model="ph")
#' plotpredictrcure(predm,model="ph")
#' 
#' # just a test:  this should be identical to classical cure model
#' pd2<-smcure(Surv(FAILTIME,FAILCENS)~TRT+SEX+AGE,cureform=~TRT+SEX+AGE,
#' data=e1684,model="ph",Var = FALSE,
#' method = "glm.fit", prior.mean = 0, prior.scale = Inf, prior.df = Inf,
#' prior.mean.for.intercept = 0, prior.scale.for.intercept = Inf,
#' prior.df.for.intercept = Inf, Warning=F)
#' pd3 <- smcure(Surv(FAILTIME,FAILCENS)~TRT+SEX+AGE,cureform=~TRT+SEX+AGE,
#' data=e1684,model="ph",Var = FALSE)
#' 
#' data(bmt)
#' # fit AFT robust cure model
#' bmtfit <- rcure(formula = Surv(Time, Status) ~ TRT, cureform = ~TRT,
#' data = bmt, model = "aft", Var = FALSE,
#' method = "glm.fit", prior.mean = 0, prior.scale = NULL, prior.df = 1,
#' prior.mean.for.intercept = 0, prior.scale.for.intercept = NULL,
#' prior.df.for.intercept = 1, min.prior.scale = 1e-12,
#' scaled = TRUE, n.iter = 100, Warning=F)
#' printrcure(bmtfit,Var = FALSE)
#' # plot predicted Survival curves by treatment groups
#' predbmt=predictrcure(bmtfit,newX=c(0,1),newZ=c(0,1),model="aft")
#' plotpredictrcure(predbmt,model="aft")
#' @importFrom stats binomial coef model.extract model.frame model.matrix na.omit optim pnorm var
#' @export

rcure <-function(formula, cureform, offset.s=NULL, data, na.action=na.omit, model= c("aft", "ph"), link="logit",
                 Var=TRUE,emmax=50,eps=1e-7,nboot=100,
                 family = binomial(link="logit"),
                 method = "glm.fit", prior.mean = 0, prior.scale = NULL, prior.df = 1,
                 prior.mean.for.intercept = 0,
                 prior.scale.for.intercept = NULL, prior.df.for.intercept = 1,
                 min.prior.scale = 1e-12, scaled = TRUE, n.iter = 100, Warning=TRUE,...)

{
	call <- match.call()

	model <- match.arg(model)
	cat("Program is running..be patient...")
	## prepare data
	data <- na.action(data)
	n <- dim(data)[1]
	mf <- model.frame(formula,data)
	cvars <- all.vars(cureform)
	Z <- as.matrix(cbind(rep(1,n),data[,cvars]))
	colnames(Z) <- c("(Intercept)",cvars)
	if(!is.null(offset.s)) {
	offsetvar <- all.vars(offset.s)
	offsetvar<-data[,offsetvar]} else
    offsetvar <- NULL
	Y <- model.extract(mf,"response")
	X <- model.matrix(attr(mf,"terms"), mf)
	if (!inherits(Y, "Surv")) stop("Response must be a survival object")
   	Time <- Y[,1]
   	Status <- Y[,2]
	bnm <- colnames(Z)
	nb <- ncol(Z)
	if(model == "ph") {
			betanm <- colnames(X)[-1]
			nbeta <- ncol(X)-1}
	if(model == "aft"){
			betanm <- colnames(X)
			nbeta <- ncol(X)}
	## initial value
	w <- Status

  bayesglm.fit0<-bayesglm(w ~ Z[,-1], family=binomial(link=link),
                          prior.mean = prior.mean,
                          prior.scale = prior.scale,
                          prior.df = prior.df,
                          prior.mean.for.intercept = prior.mean.for.intercept,
                          prior.scale.for.intercept = prior.scale.for.intercept,
                          prior.df.for.intercept = prior.df.for.intercept,
                          min.prior.scale=min.prior.scale,
                          scaled = scaled, maxit = n.iter, Warning=Warning)
  b<-coef(bayesglm.fit0)
  if(model=="ph") beta <- coxph(Surv(Time, Status)~X[,-1]+offset(log(w)), subset=w!=0, method="breslow")$coef
	if(model=="aft") beta <- survreg(Surv(Time,Status)~X[,-1])$coef
     	## do EM algo
	emfit <- rem(Time,Status,X,Z,offsetvar,b,beta,model,link,emmax,eps,
                    prior.mean, prior.scale, prior.df, prior.mean.for.intercept, prior.scale.for.intercept,
	                  prior.df.for.intercept, min.prior.scale, scaled, n.iter, Warning)
		b <- emfit$b
		beta <- emfit$latencyfit
		s <- emfit$Survival
		logistfit <- emfit$logistfit
    nn.iter<-emfit$nn.iter

  ## bootstrap to get the sd, they bootstrap Status=1 and Status=0 seperately
		  if(Var){
		   if(model=="ph") {b_boot<-matrix(rep(0,nboot*nb), nrow=nboot)
					beta_boot<-matrix(rep(0,nboot*(nbeta)), nrow=nboot)
					nn.iter_boot <- matrix(rep(0,nboot),ncol=1)
		   }

		   if(model=="aft") {b_boot<-matrix(rep(0,nboot*nb), nrow=nboot)
					beta_boot<-matrix(rep(0,nboot*(nbeta)), nrow=nboot)
					nn.iter_boot <- matrix(rep(0,nboot),ncol=1)
		   }

		  ## stratified bootstrap to calculate variance of parameter
		  tempdata <- cbind(Time,Status,X,Z)
		  data1<-subset(tempdata,Status==1);data0<-subset(tempdata,Status==0)
		  n1<-nrow(data1);n0<-nrow(data0)
		  i<-1

		while (i<=nboot){
             id1<-sample(1:n1,n1,replace=TRUE);id0<-sample(1:n0,n0,replace=TRUE)
			bootdata<-rbind(data1[id1,],data0[id0,])
			bootZ <- bootdata[,bnm]
			if(model=="ph") bootX <- as.matrix(cbind(rep(1,n),bootdata[,betanm]))
			if(model=="aft") bootX <- bootdata[,betanm]
			bootfit <- rem(bootdata[,1],bootdata[,2],bootX,bootZ,offsetvar,b,beta,model,link,emmax,eps,
			                    prior.mean, prior.scale, prior.df, prior.mean.for.intercept, prior.scale.for.intercept,
			                    prior.df.for.intercept, min.prior.scale, scaled, n.iter, Warning)
			b_boot[i,] <- bootfit$b
		  beta_boot[i,] <- bootfit$latencyfit
			nn.iter_boot[i,]<-bootfit$nn.iter
   		if (bootfit$tau<eps) i<-i+1
		}

		b_var <- apply(b_boot, 2, var)
		beta_var <- apply(beta_boot, 2, var)
		b_sd <- sqrt(b_var)
		beta_sd <- sqrt(beta_var)

		nn.iter_mean<-apply(nn.iter_boot, 2, mean)
		nn.iter_var<-apply(nn.iter_boot, 2, var)
		nn.iter_sd<-sqrt(nn.iter_var)
	}

	fit<-list()
	class(fit) <- c("rcure")
	fit$logistfit <- logistfit
	fit$b <- b
	fit$beta <- beta
	fit$nn.iter<-nn.iter

	if(Var){
	fit$b_var <- b_var
	fit$b_sd <- b_sd
	fit$b_zvalue <- fit$b/b_sd
	fit$b_pvalue <- (1-pnorm(abs(fit$b_zvalue)))*2
	fit$beta_var <- beta_var
	fit$beta_sd <- beta_sd
	fit$beta_zvalue <- fit$beta/beta_sd
	fit$beta_pvalue <- (1-pnorm(abs(fit$beta_zvalue)))*2

	fit$nn.iter_boot<-nn.iter_boot
	fit$nn.iter_mean<-nn.iter_mean
	fit$nn.iter_sd<-nn.iter_sd
	}
	cat(" done.\n")
	fit$call <- call
	fit$bnm <- bnm
	fit$betanm <- betanm
	fit$s <- s
	fit$Time <- Time
	if(model=="aft"){
	error <- drop(log(Time)-beta%*%t(X))
	fit$error <- error}
	fit
	printrcure(fit,Var)
}


