diff.var.test = function(ts1,ts2,order.pic,alpha=0.05,monte.carlo.test=FALSE,ntrials=10000,lag.lbtest=10) {
## PERFORMS DIFFERENCE IN VAR(P) PROCESSES TEST
## INPUT:
##   TS1[N1,SDIM]: 1ST TIME SERIES: N1=TIME AND SDIM=VARIABLES
##   TS2[N2,SDIM]: 2ND TIME SERIES: N2=TIME AND SDIM=VARIABLES
##   ORDER.PIC: ORDER OF THE AR PROCESS
##   ALPHA: SIGNIFICANCE LEVEL
##   MONTE.CARLO.TEST[LOGICAL]: DO MONTE CARLO TO FIND CRITICAL THESHOLD? 
## COMMENTS:
##   1) IF DIMENSIONS ARE CONSTANT FOR EACH CALL, SET MONTE.CARLO.TEST=TRUE *ONCE* 
##      TO COMPUTE SIGNIFCANCE LEVELS, THEN SET = FALSE FOR ALL OTHER CALCULATIONS
##	 2) THIS FUNCTION CALLS TEST.EQUALITY.COV.HIERARCHICAL 

ts1     = as.matrix(ts1)
ts2     = as.matrix(ts2)
n1.all  = dim(ts1)[1]
n2.all  = dim(ts2)[1]
sdim    = dim(ts1)[2]
if (sdim != dim(ts2)[2]) stop('ts1 and ts2 not dimensioned correctly')

npic1   = (1+order.pic):n1.all
npic2   = (1+order.pic):n2.all
ntot1   = length(npic1)
ntot2   = length(npic2)

one     = c(rep(1,ntot1),rep(0,ntot2))
one     = cbind(one,1-one)
y1      = as.matrix(ts1[npic1,])
y2      = as.matrix(ts2[npic2,])
x1      = NULL
x2      = NULL
for (lag in 1:order.pic) x1 = cbind(x1,ts1[npic1-lag,])
for (lag in 1:order.pic) x2 = cbind(x2,ts2[npic2-lag,])
x1.svd  = svd(cbind(rep(1,ntot1),x1))
x2.svd  = svd(cbind(rep(1,ntot2),x2))
x0.svd  = svd(cbind(one,rbind(x1,x2)))
err1    = y1 - x1.svd$u %*% ( t(x1.svd$u) %*% y1)
err2    = y2 - x2.svd$u %*% ( t(x2.svd$u) %*% y2)
err0    = rbind(y1,y2) - x0.svd$u %*% ( t(x0.svd$u) %*% rbind(y1,y2))
sse1    = t(err1) %*% err1
sse2    = t(err2) %*% err2
sse0    = t(err0) %*% err0
nu1     = dim(x1.svd$u)[1] - dim(x1.svd$v)[1]
nu2     = dim(x2.svd$u)[1] - dim(x2.svd$v)[1]
nu0     = dim(x0.svd$u)[1] - dim(x0.svd$v)[1]
cova    = (sse1 + sse2)/(nu1+nu2)
d.noise = (nu1+nu2)*log(det(cova)) - nu1*log(det(sse1/nu1)) - nu2*log(det(sse2/nu2))
d.parms = (nu1+nu2)* ( log(det(sse0/(nu1+nu2))) - log(det(cova)))
d.total = d.noise + d.parms

### DISCRIMINANTS OF D.PARMS
cda.noise = gev(sse1/nu1,sse2/nu2)
cda.parms = gev((sse0-sse1-sse2)/(nu1+nu2),(sse1+sse2)/(nu1+nu2))

cda.noise.deviance = nu1 * log(nu1 + nu2 / cda.noise$lambda) + nu2 * log(nu1 * cda.noise$lambda + nu2) - (nu1 + nu2) * log(nu1+nu2)
cda.parms.deviance = (nu1 + nu2) * log(1 + cda.parms$lambda)

#### SVD OF NORMALIZED DIFFERENCE IN BETAS
x1    = t(t(x1)-colMeans(x1))
x2    = t(t(x2)-colMeans(x2))
covmh = chol2inv(chol(chol2inv(chol(t(x1) %*% x1)) + chol2inv(chol(t(x2) %*% x2))))/(nu1+nu2)
cova.eigen  = eigen(cova )
covmh.eigen = eigen(covmh)
covmh.sqrt  = covmh.eigen$vectors %*% (sqrt(covmh.eigen$values) * t(covmh.eigen$vectors))
covai.sqrt  = cova.eigen$vectors %*% (t(cova.eigen$vectors)/sqrt(cova.eigen$values))

beta1 = (x1.svd$v %*% ((t(x1.svd$u) / x1.svd$d) %*% y1))
beta2 = (x2.svd$v %*% ((t(x2.svd$u) / x2.svd$d) %*% y2))
int1  = beta1[ 1,]
int2  = beta2[ 1,]
beta1 = beta1[-1,,drop=FALSE]
beta2 = beta2[-1,,drop=FALSE]

delta.norm = covmh.sqrt %*% (beta2 - beta1) %*% covai.sqrt
delta.svd  = svd(delta.norm)
sval       = delta.svd$d
px         = covmh.sqrt %*% delta.svd$u
qx         = chol2inv(chol(covmh)) %*% px
qy         = covai.sqrt %*% delta.svd$v
py         = cova %*% qy

#### COMPUTE VARIATES
r1.noise = err1 %*% cda.noise$q
r2.noise = err2 %*% cda.noise$q

#### COMPUTE LOADING VECTORS
p.noise = (sse2/nu2) %*% cda.noise$q
p.parms = ((sse1+sse2)/(nu1+nu2)) %*% cda.parms$q

##### PROPORTIONALITY TEST
if (sdim > 1) {
	smat      = array(NA,dim=c(sdim,sdim,2))
	smat[,,1] = sse1/nu1
	smat[,,2] = sse2/nu2
	dof.all   = c(nu1,nu2)
	prop.list = test.equality.cov.hierarchical(smat,dof.all,alpha=alpha)	
} else prop.list = list(d.prop=NA,d.equal.var=NA,d.equal.cor=NA,crit.d.prop=NA,crit.d.equal.var=NA,crit.d.equal.cor=NA)

###### PORTMANTAU TEST
lbtest1 = LjungBox(err1,lags=lag.lbtest,order=order.pic)
lbtest2 = LjungBox(err2,lags=lag.lbtest,order=order.pic)

#### DECOMPOSE RESIDUAL VARIANCE = VARIANCE ( 1- R.SQUARE )
# y1      = t(t(y1)-colMeans(y1))
# y2      = t(t(y2)-colMeans(y2))
# r1.y    = y1 %*% cda.noise$q
# r2.y    = y2 %*% cda.noise$q
# var.tot1 = colMeans(r1.y^2)
# var.tot2 = colMeans(r2.y^2)
# var.err1 = colSums(r1.noise^2)/nu1
# var.err2 = colSums(r2.noise^2)/nu2
# rsqr1    = 1 - var.err1/var.tot1
# rsqr2    = 1 - var.err2/var.tot2

### DECOMPOSE NOISE MATRICES INTO VARIANCE * ( 1- CCA)
### ABANDON: D.NOISE.RSQ OFTEN IS NEGATIVE, SHOWING THIS IS NOT A MEANINGFUL DECOMPOSITION.
# cov1 = t(y1) %*% y1 / nu1
# cov2 = t(y2) %*% y2 / nu2
# cca1 = chol2inv(chol(cov1)) %*% (sse1/nu1)
# cca2 = chol2inv(chol(cov2)) %*% (sse2/nu2)
# ccaa = chol2inv(chol((nu1*cov1 + nu2*cov2)/(nu1+nu2))) %*% cova

# d.noise.rsq = (nu1+nu2)*log(det(ccaa)) - nu1*log(det(cca1)) - nu2*log(det(cca2))
# d.noise.var = (nu1+nu2)*log(det((nu1*cov1 + nu2*cov2)/(nu1+nu2))) - nu1*log(det(cov1)) - nu2*log(det(cov2))


## consistency checks
# print(all.equal(t(qy) %*% cova %*% qy, diag(1,sdim)))
# print(all.equal(beta2-beta1,qx %*% diag(sval) %*% t(py)))
# print(all.equal(t(px) %*% ( beta2 - beta1),diag(sval) %*% t(py)))
# print(all.equal(t(px) %*% ( beta2 - beta1) %*% qy, diag(sval)))
# print(all.equal(abs(p.parms),abs(py)))


if (monte.carlo.test) {
	print(paste('doing montecarlo',nu1,nu2))
	cov.identity = diag(1,nrow=sdim,ncol=sdim)
	d.noise.mc = as.numeric(rep(NA,ntrials))
	d.parms.mc = as.numeric(rep(NA,ntrials))
	d.lambd.mc = array(NA,dim=c(sdim,ntrials))
	lambda.parms.mc = array(NA,dim=c(sdim,ntrials))
	for (nt in 1:ntrials) {
		cov1 = as.matrix(rWishart(1,nu1,cov.identity)[,,1])
		cov2 = as.matrix(rWishart(1,nu2,cov.identity)[,,1])
		cov3 = as.matrix(rWishart(1,nu0-nu1-nu2,cov.identity)[,,1])
		cova = (cov1 + cov2)/(nu1+nu2)
		d.noise.mc[nt] = (nu1+nu2)*log(det(cova)) - nu1*log(det(cov1/nu1)) - nu2*log(det(cov2/nu2))
		d.parms.mc[nt] = (nu1+nu2)* ( log(det((cov1+cov2+cov3)/(nu1+nu2))) - log(det(cova)))
		d.lambd.mc[,nt] = gev(cov1/nu1,cov2/nu2)$lambda
		lambda.parms.mc[,nt] = gev(cov3/(nu1+nu2),cova)$lambda
	}
	d.noise.crit = quantile(d.noise.mc,probs=1-alpha)
	d.parms.crit = quantile(d.parms.mc,probs=1-alpha)
	d.total.crit = quantile(d.noise.mc+d.parms.mc,probs=1-alpha)
	lambda.noise.crit = array(NA,dim=c(sdim,2*length(alpha)))
	lambda.parms.crit = array(NA,dim=c(sdim,2*length(alpha)))
	for (ns in 1:sdim) lambda.noise.crit[ns,] = quantile(     d.lambd.mc[ns,],probs=c(rev(alpha),1-alpha))
	for (ns in 1:sdim) lambda.parms.crit[ns,] = quantile(lambda.parms.mc[ns,],probs=c(rev(alpha),1-alpha))
	colnames(lambda.noise.crit) = paste(100*c(rev(alpha),1-alpha),'%',sep='')
	colnames(lambda.parms.crit) = paste(100*c(rev(alpha),1-alpha),'%',sep='')
} else {
	d.noise.crit = 'must set monte.carlo.test=TRUE'
	d.parms.crit = 'must set monte.carlo.test=TRUE'
	d.total.crit = 'must set monte.carlo.test=TRUE'
	lambda.noise.crit = 'must set monte.carlo.test=TRUE'
	lambda.parms.crit = 'must set monte.carlo.test=TRUE'
}

nu.parms = sdim*(nu0-nu1-nu2)
nu.noise = sdim*(sdim+1)/2
nu.total = nu.parms + nu.noise
d.parms.crit.asym = qchisq(alpha,df=nu.parms,lower.tail=FALSE)
d.noise.crit.asym = qchisq(alpha,df=nu.noise,lower.tail=FALSE)
d.total.crit.asym = qchisq(alpha,df=nu.total,lower.tail=FALSE)

# SANITY CHECK: VARIANCE RATIO = DISCRIMINANT RATIOS
# print(all.equal(colSums(r1.noise^2)/colSums(r2.noise^2)*nu2/nu1,cda.noise$lambda))
# print(colSums(r2.noise^2)/nu2) ## should be 1
# print(sum(cda.noise.deviance))
# print(d.noise)
# if (abs(sum(cda.parms.deviance)-d.parms) > 1.e-10) stop('parms components do not add up correctly')
# if (abs(sum(cda.noise.deviance)-d.noise) > 1.e-10) stop('noise components do not add up correctly')

### FILL TIME SERIES WITH NAs TO PRESERVE DIMENSION
na.fill  = array(NA,dim=c(order.pic,sdim))
r1.noise = rbind(na.fill,r1.noise)
r2.noise = rbind(na.fill,r2.noise)

#### OPTIMIZE INITIAL CONDITIONS
x1    = t(t(x1)-colMeans(x1))
x2    = t(t(x2)-colMeans(x2))
beta1 = (x1.svd$v %*% ((t(x1.svd$u) / x1.svd$d) %*% y1))
beta2 = (x2.svd$v %*% ((t(x2.svd$u) / x2.svd$d) %*% y2))
int1  = beta1[ 1,]
int2  = beta2[ 1,]
beta1 = beta1[-1,,drop=FALSE]
beta2 = beta2[-1,,drop=FALSE]

covmh = chol2inv(chol(chol2inv(chol(t(x1) %*% x1)) + chol2inv(chol(t(x2) %*% x2))))
cda.initc.vec  = covmh %*% (beta1 - beta2) %*% cda.parms$q %*% diag(cda.parms$lambda)


##############################################
#### SOLVE SEPARATE EIGENVALUE PROBLEM FOR SANITY CHECK
#### (there are 2 ways to get optimal initial conditions)
##############################################
# # CHECK DELTA IS THE SAME
# delta.check = t(beta1 - beta2) %*% covmh %*% (beta1 - beta2)
# delta       = sse0 - sse1 - sse2
# print(all.equal(delta,delta.check))

# # SEPARATE EIGENVALUE PROBLEM
# imat1       = (beta1 - beta2) %*% chol2inv(chol(cova)) %*% t(beta1 - beta2) / (nu1+nu2)
# imat2       = chol2inv(chol(covmh))
# cda.initc   = gev(imat1,imat2)
# cda.initc.deviance = (nu1 + nu2) * log(1 + cda.initc$lambda)

# #### CHECK THAT PARMS IS THE SAME USING DELTA AND COVA
# cda.parms.check = gev(delta.check/(nu1+nu2),cova)
# print(all.equal(cda.parms.check$lambda,cda.parms$lambda))

# ## CHECK THAT EIGENVALUES FROM PARMS AND INITC ARE THE SAME
# cda.initc$lambda = ifelse(abs(cda.initc$lambda) < 1.e-10,0,cda.initc$lambda)
# print(cda.parms$lambda)
# print(cda.initc$lambda)

# ## CHECK THAT EIGENVECTORS ARE THE SAME, UP TO A CONSTANT FACTOR
# print(cda.initc.vec/cda.initc$q[,1:sdim])


list(d.noise=d.noise,d.parms=d.parms,d.total=d.total,cda.noise=cda.noise,cda.parms=cda.parms,
     d.noise.crit=d.noise.crit,d.parms.crit=d.parms.crit,d.total.crit=d.total.crit,
     lambda.noise.crit=lambda.noise.crit,lambda.parms.crit=lambda.parms.crit,
     d.parms.crit.asym = d.parms.crit.asym , d.noise.crit.asym = d.noise.crit.asym, 
     d.total.crit.asym = d.total.crit.asym,
     sdim=sdim,order.pic=order.pic,
     cda.noise.deviance=cda.noise.deviance,
     cda.parms.deviance=cda.parms.deviance,
     cda.initc.vec=cda.initc.vec,
     r1.noise=r1.noise,r2.noise=r2.noise,
     err1=err1,err2=err2,
     beta1=beta1,beta2=beta2,int1=int1,int2=int2,
     p.noise=p.noise,p.parms=p.parms,
     px=px,qx=qx,py=py,qy=qy,sval=sval,
     prop.list=prop.list,lbtest1=lbtest1,lbtest2=lbtest2)
	
}






