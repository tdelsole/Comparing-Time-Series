compare.anncyc = function(ts1,ts2,order.pic,nharm.tot,period=12,first.step=c(1,1),alpha=0.05) {
## PERFORMS DIFFERENCE IN ARX(P) MODELS, WERE 'X' IS THE ANNUAL CYCLE
## INPUT:
##   TS1[N1]: 1ST TIME SERIES: N1=TIME 
##   TS2[N2]: 2ND TIME SERIES: N2=TIME
##   ORDER.PIC: ORDER OF THE AR PROCESS
##	 NHARM.TOT: NUMBER OF HARMONICS (SINE AND COSINE COUNT AS "1" HARMONIC)
##   PERIOD: PERIOD OF THE ANNUAL CYCLE (DEFAULT = 12 FOR MONTHLY DATA; 365.24 FOR DAILY DATA)
##   FIRST.STEP[2]: POSITION WITHIN THE ANNUAL CYCLE OF THE FIRST TIME STEP OF EACH TIME SERIES
##   ALPHA: SIGNIFICANCE LEVEL

##############################################################
######## BEGIN FUNCTION
##############################################################

	n1.all = dim(as.matrix(ts1))[1]
	n2.all = dim(as.matrix(ts2))[1]
	npic1  = (1+order.pic):n1.all
	npic2  = (1+order.pic):n2.all
	ntot1  = length(npic1)
	ntot2  = length(npic2)
	j1     = rep(1,ntot1)
	j2     = rep(1,ntot2)
	
	if (nharm.tot * 2 >= period) stop('number of harmonics cannot exceed Nyquist frequency')
	alpha.step = 1-(1-alpha)^(1/3)
	
	no1 = npic1 + first.step[1] - 1
	no2 = npic2 + first.step[2] - 1
	
	e1  = NULL
	e2  = NULL
	for (nh in 1:nharm.tot) e1  = cbind(e1,cos(2*pi*no1*nh/period),sin(2*pi*no1*nh/period))
	for (nh in 1:nharm.tot) e2  = cbind(e2,cos(2*pi*no2*nh/period),sin(2*pi*no2*nh/period))
	edim = 2*nharm.tot

	estar      = NULL
	tstar      = 1:ceiling(period) 
	for (nh in 1:nharm.tot) estar = cbind(estar,cos(2*pi*tstar*nh/period),sin(2*pi*tstar*nh/period))
	ete.diag   = diag(t(estar) %*% estar)
		
	y1     = ts1[npic1]
	y2     = ts2[npic2]
	x1.lag = NULL
	x2.lag = NULL
	if (order.pic > 0) for (lag in 1:order.pic) x1.lag = cbind(x1.lag,ts1[npic1-lag])
	if (order.pic > 0) for (lag in 1:order.pic) x2.lag = cbind(x2.lag,ts2[npic2-lag])
	
	colnames(e1    ) = paste('E',1:(2*nharm.tot),sep='')
	colnames(e2    ) = paste('E',1:(2*nharm.tot),sep='')
	colnames(x1.lag) = paste('Lag',1:order.pic  ,sep='')
	colnames(x2.lag) = paste('Lag',1:order.pic  ,sep='')

	x1        = cbind(e1,x1.lag)
	x2        = cbind(e2,x2.lag)	
	lm.1      = lm(y1 ~ x1)
	lm.2      = lm(y2 ~ x2)
	sse1      = sum(residuals(lm.1)^2)
	sse2      = sum(residuals(lm.2)^2)
	loglik.H3 = ntot1 * log(sse1/ntot1) + ntot2 * log(sse2/ntot1)
	loglik.H2 = (ntot1 + ntot2) * log ( (sse1 + sse2)/(ntot1 + ntot2) )
	parm.H3   = dim(x1)[2] + dim(x2)[2] + 4  ## 4 = 2 means and 2 variance
	parm.H2   = dim(x1)[2] + dim(x2)[2] + 3  ## 3 = 2 means and 1 variance
	x1.save   = cbind(rep(1,ntot1),x1)
	x2.save   = cbind(rep(1,ntot2),x2)
	varn1     = sse1 / lm.1$df.residual
	varn2     = sse2 / lm.2$df.residual
		
	y         = c(y1,y2)
	x1        = cbind(e1,x1.lag,array(0,dim=c(ntot1,order.pic)),j1,1-j1)
	x2        = cbind(e2,array(0,dim=c(ntot2,order.pic)),x2.lag,1-j2,j2)
	x         = rbind(x1,x2)
	lm.H1     = lm(y~x-1)
	sse.H1    = sum(residuals(lm.H1)^2)
	loglik.H1 = (ntot1 + ntot2) * log ( sse.H1 / (ntot1 + ntot2) )
	parm.H1   = dim(x)[2] + 1
	
	x1        = cbind(e1,x1.lag,j1,1-j1)
	x2        = cbind(e2,x2.lag,1-j2,j2)
	x         = rbind(x1,x2)
	lm.H0     = lm(y~x-1)
	sse.H0    = sum(residuals(lm.H0)^2)
	loglik.H0 = (ntot1 + ntot2) * log ( sse.H0 / (ntot1 + ntot2) )
	parm.H0   = dim(x)[2] + 1
	
	############################################
	##### CALCULATE DEVIANCES
	############################################
	dev.01    = loglik.H0 - loglik.H1
	dev.12    = loglik.H1 - loglik.H2
	dev.23    = loglik.H2 - loglik.H3
	dev.13    = loglik.H1 - loglik.H3
	dev.03    = loglik.H0 - loglik.H3

	############################################
	##### DIAGNOSE DIFFERENCES IN E-COEFFICIENTS
	############################################
	beta1.H2   = coef(lm.1)
	beta2.H2   = coef(lm.2)
	beta1.H1.l = c(rep(TRUE ,edim),rep(TRUE ,order.pic),rep(FALSE,order.pic),FALSE,FALSE)
	beta2.H1.l = c(rep(TRUE ,edim),rep(FALSE,order.pic),rep(TRUE ,order.pic),FALSE,FALSE)
	intr1.H1.l = c(rep(FALSE,edim),rep(FALSE,order.pic),rep(FALSE,order.pic),TRUE ,FALSE)
	intr2.H1.l = c(rep(FALSE,edim),rep(FALSE,order.pic),rep(FALSE,order.pic),FALSE,TRUE )
	beta1.H1   = c(coef(lm.H1)[intr1.H1.l],coef(lm.H1)[beta1.H1.l])
	beta2.H1   = c(coef(lm.H1)[intr2.H1.l],coef(lm.H1)[beta2.H1.l])
	sse.diff   = sse.H1 - sse1 - sse2
	pred1      = x1.save %*% (beta1.H2 - beta1.H1)
	pred2      = x2.save %*% (beta2.H2 - beta2.H1)
	sse.diff2  = sum(pred1^2) + sum(pred2^2)
	if (!all.equal(sse.diff,sse.diff2)) stop('inconsistent errors')
	
	a.mat      = cbind(rep(0,edim),diag( 1,nrow=edim),array(0,dim=c(edim,order.pic)),
	                   rep(0,edim),diag(-1,nrow=edim),array(0,dim=c(edim,order.pic)))
	b.mat      = c(beta1.H2,beta2.H2)
	x1.dum     = cbind(x1.save,array(0,dim=dim(x1.save)))
	x2.dum     = cbind(array(0,dim=dim(x2.save)),x2.save)
	x.mat      = rbind(x1.dum,x2.dum)
	xtx.inv    = chol2inv(chol(t(x.mat) %*% x.mat))	
	cov.ab     = chol2inv(chol( a.mat %*% xtx.inv %*% t(a.mat) ))
	ab.mat     = a.mat %*% b.mat
	sse.diff3  = as.numeric(t(ab.mat) %*% cov.ab %*% ab.mat)
	
	
	l.i        = c(TRUE ,rep(FALSE,edim),rep(FALSE,order.pic))
	l.e        = c(FALSE,rep(TRUE ,edim),rep(FALSE,order.pic))
	l.c        = c(FALSE,rep(FALSE,edim),rep(TRUE ,order.pic))

	intercept  = c    (beta1.H2[l.i],beta2.H2[l.i])
	e.hat      = cbind(beta1.H2[l.e],beta2.H2[l.e])
	c.hat      = cbind(beta1.H2[l.c],beta2.H2[l.c])

	anncyc1    = estar %*% e.hat[,1]
	anncyc2    = estar %*% e.hat[,2]
	covm.star  = estar %*% cov.ab %*% t(estar) / ete.diag[1]^2
	sse.diff3.check = as.numeric(t(anncyc1-anncyc2) %*% covm.star %*% (anncyc1-anncyc2))
	if (!all.equal(sse.diff3,sse.diff3.check)) stop('diff3 not consistent')
	
	
	############################################
	##### ANNUAL CYCLE RESPONSE
	############################################
	anncyc.response = array(0,dim=c(12,2))
	for (n in 1:2) {
		amps = e.hat[,n]
		phis = c.hat[,n]
		for (nh in 1:nharm.tot) {
			ak = amps[2*nh-1] - 1i * amps[2*nh]
			ck = ak / (1 - sum(phis * exp(-1i * 2 * pi * nh * (1:order.pic) / 12)))
			anncyc.response[,n] = anncyc.response[,n] + Re(ck * exp(1i * 2 * pi * nh * (1:12)/12))
		}
	}

			
	############################################
	##### DIAGNOSE DIFFERENCES IN NOISES
	############################################
	fval = (sse1/lm.1$df.residual) / (sse2/lm.2$df.residual)
	fval.upper = qf(alpha.step/2,lm.1$df.residual,lm.2$df.residual,lower.tail=FALSE)
	fval.lower = qf(alpha.step/2,lm.1$df.residual,lm.2$df.residual,lower.tail=TRUE)
	fscale     = fval * lm.1$df.residual / lm.2$df.residual
	dev.23.check = ntot1 * log(ntot1/(ntot1+ntot2)) + ntot1 * log( 1 + 1/fscale ) +
	               ntot2 * log(ntot2/(ntot1+ntot2)) + ntot2 * log( 1 + fscale   )
	if (abs(dev.23-dev.23.check) > 1.e-10) stop('inconsistent noise ratio')
	# if (!all.equal(dev.23,dev.23.check)) stop('inconsistent noise ratio')
	
	############################################
	##### RESIDUALS
	############################################
	residuals.data = cbind(residuals(lm.1),residuals(lm.2)) 
	
	############################################
	##### WRITE OUT RESULTS
	############################################
	dev.01.crit = qchisq(alpha.step,parm.H1-parm.H0,lower.tail=FALSE)
	dev.12.crit = qchisq(alpha.step,parm.H2-parm.H1,lower.tail=FALSE)
	dev.23.crit = qchisq(alpha.step,parm.H3-parm.H2,lower.tail=FALSE)
	dev.03.crit = qchisq(alpha     ,parm.H3-parm.H0,lower.tail=FALSE)

	
	list.dev = list(dev.01 = dev.01, dev.01.crit = dev.01.crit,
	     			dev.03 = dev.03, dev.03.crit = dev.03.crit,
	     			dev.12 = dev.12, dev.12.crit = dev.12.crit,
	     			dev.23 = dev.23, dev.23.crit = dev.23.crit,
	     			parm.H0 = parm.H0, parm.H1   = parm.H1,
	     			parm.H2 = parm.H2, parm.H3   = parm.H3,
	     			loglik.H0 = loglik.H0, loglik.H1 = loglik.H1,
	     			loglik.H2 = loglik.H2, loglik.H3 = loglik.H3,
	     			fval    = fval, fval.upper = fval.upper, fval.lower = fval.lower,
	     			anncyc1 = anncyc1, anncyc2 = anncyc2,
	     			beta1.H2 = beta1.H2, beta2.H2 = beta2.H2,
	     			varn1   = varn1, varn2 = varn2, anncyc.response = anncyc.response,
	     			alpha.step = alpha.step, intercept = intercept,
	     			e.hat = e.hat, c.hat = c.hat, residuals.data = residuals.data)

##############################################################
######## END FUNCTION
##############################################################



}