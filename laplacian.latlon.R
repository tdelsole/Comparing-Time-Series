laplacian.latlon = function(lon,lat,field=NULL,keeper=NULL,neof=50,luniform=TRUE,idate=NULL,dname=NULL,rm.clim=TRUE) {
### THIS FUNCTION COMPUTES AN ORTHOGONAL BASIS SET IN AN
### ARBITRARY DOMAIN DEFINED BY KEEPER=TRUE, 
### AND PROJECTS THEM ONTO 'FIELD' IF IT IS DEFINED.  
### 
### THE BASIS SET IS COMPUTED FROM THE ORTHOGONAL EIGENVECTORS OF THE
### GREENS FUNCTION ON A SPHERE WITHIN THE DOMAIN DEFINED BY KEEPER = TRUE.  
##  INPUT:
##    LON[NLON]: LONGITUDE POINTS (ASSUMED TO BE EQUALLY SPACED)
##    LAT[NLAT]: LATITUDE POINTS (ASSUMED TO BE EQUALLY SPACED)
##    FIELD[NLON,NLAT,NTIME]: FIELD FOR PROJECTING ORTHOGONAL EIGENVECTORS (CAN BE OMITTED)
##    KEEPER[NLON,NLAT]: A LOGICAL MASK, GRIDPOINTS SET TO 'FALSE' WILL BE OMITTED.
##    NEOF: NUMBER OF ORTHOGONAL BASIS VECTORS TO RETAIN (RESET TO A SMALLER VALUE IF DIMENSION < NEOF)
##    LUNIFORM: A LOGICAL VARIABLE.  = TRUE TO INCLUDE ZERO-EIGENVALUE AS AN EIGENVECTOR.
##    IDATE = A VECTOR INDICATING THE STARTING TIME AND INCREMENT (E.G., IDATA = C(1, "JAN", 1979, "1MO") )
##    DNAME = NAME OF DATA (E.G., HADSST, ERSSTV, NCEP/NCAR REAN TEMP, ETC)
# OUTPUT: LIST
#	$EOF[NLON*NLAT,NEOF]: THE FIRST NEOF SCALED ORTHOGONAL EOFS OF GREENS FUNCTION, NORMALIZED TO UNIT AREA WEIGHTED VARIANCE
#	$PC[NTIM,NEOF]: THE PCS, NORMALIZED TO UNIT VARIANCE [NTOT X NEOF]
#	$EVAL: THE EIGENVALUES
#	$FEXPVAR: FRACTION OF EXPLAINED VARIANCE FOR EACH EOF.
#	$EOFI: PSEUDO INVERSE OF EOF (I.E.,T(EOFI) %*% EOF = I)
#	$BAD[NLON*NLAT]: LOGICAL ARRAY INDICATING WHICH POINTS WERE DROPPED FROM THE EOF ANALYSIS
#   $LON: THE LONGITUDES OF THE SPATIAL FIELD
#   $LAT: THE LATITUDES OF THE SPATIAL FIELD
#   $IDATE: REPEATS INPUT "IDATE"
#   $DNAME: REPEATS INPUT "DNAME"
#   $TOTWAVNUM: APPROXIMATE TOTAL WAVENUMBER (CAN BE NEGATIVE DUE TO NUMERICAL APPROXIMATION)
#   $AREA: AREA FACTOR FOR NORMALIZING EOFS AND EXPLAINED VARIANCES
### COMMENT:
#    1) ORTHOGONALITY T(EOF) %*% DIAG(AREA) %*% EOF = IDENTITY
#    2) EXPLAINED VARIANCE IS RELATIVE TO TRACE(FIELD^T %*% DIAG(AREA) %*% FIELD)/N

#####################################
#######  DEFINE GRID
#####################################
nlon    = length(lon)
nlat    = length(lat)
lat.all = rep(lat,each=nlon) * pi / 180
lon.all = rep(lon,     nlat) * pi / 180

if ( is.null(keeper)) keeper=rep(TRUE,length(lon)*length(lat))

#####################################
#######  CHECK FOR UNEQUAL GRID
#####################################
if ( var(diff(lat)) > 1.e-6) print('latitudes not equally space')
if ( var(diff(lon)) > 1.e-6) stop('longitudes not equally space')

### IDENTIFY MISSING DATA AND MASK THOSE TOO
if (!is.null(field)) {
	if (length(field) %% (nlon*nlat) != 0) stop('field not consistent with lat-lon grid')
	ntime = length(field)/nlon/nlat
	dim(field) = c(nlon*nlat,ntime)
	keeper     = keeper & !is.na(rowSums(field))
}


ngood    = sum(keeper)
lon.good = lon.all[keeper]
lat.good = lat.all[keeper]

weight   = sqrt(cos(lat.good))
dphi     = abs(lon[2]-lon[1]) * pi / 180 
dtheta   = mean(diff(lat)) * pi / 180

print(paste('Number of non-missing grid points=',ngood))
#####################################
#######  RE-SET NEOF TO THE MINIMUM OF (NEOF,NGOOD)
#####################################
neof     = min(neof,ngood)

#####################################
#######  COMPUTE GREAT CIRCLE DISTANCES
#####################################
sinsqr.gcdist.over2 = array(NA,dim=c(ngood,ngood))
for ( j in 1:ngood) {
	i = 1:j
	delphi  = lat.good[i]-lat.good[j]
	dlambda = lon.good[i]-lon.good[j]
	sinsqr.gcdist.over2[i,j] = sin(delphi/2)^2 + cos(lat.good[i]) * cos(lat.good[j]) * sin(dlambda/2)^2
	sinsqr.gcdist.over2[j,i] = sinsqr.gcdist.over2[i,j]
}

#####################################
#### CONSTRUCT OMEGA MATRIX, REMOVE INFINITE DIAGONAL ELEMENTS
#####################################
omega       = -log(2*sinsqr.gcdist.over2)/4/pi * dphi * dtheta

#####################################
#### SYMMETERIZE THE GREENS FUNCTION
#####################################
omega  = omega  * rep(weight,ngood) * rep(weight,each=ngood)

#####################################
#### CORRECT DIAGONAL ELEMENTS
#####################################
rho0        = sqrt(dphi * dtheta * cos(lat.good) /pi)
diag(omega) = rho0^2/4*( 1 - 2*log(rho0/sqrt(2)) ) 


#####################################
#### COMPUTE ORTHOGONAL COMPLEMENT TO CONSTANT VECTOR
#####################################
if (luniform) {
	vone  = rep(1,ngood)
	uvec  = diag(weight,nrow=ngood) %*% vone
	uvec  = uvec / sqrt(sum(uvec^2))
	ou    = omega %*% uvec
	omega = omega - ou %*% t(uvec) - uvec %*% t(ou) + uvec %*% t(uvec) * sum(ou*uvec)
}


#####################################
#### COMPUTE EIGENVECTORS #####
#####################################
omega.eigen = eigen(omega,symmetric=TRUE)

#####################################
#### IDENTIFY THE ZERO EIGENVECTOR
#####################################
npic.zero = which(abs(omega.eigen$values) < 1.e-10)
if ( length(npic.zero) != 1 ) print(omega.eigen$values)
if ( length(npic.zero) != 1 ) stop('difficulty finding zero eigenvalue')
omega.eigen$values = c(omega.eigen$values[npic.zero],omega.eigen$values[-npic.zero])
omega.eigen$vectors = cbind(omega.eigen$vectors[,npic.zero],omega.eigen$vectors[,-npic.zero])

#####################################
#### FORCE CONSTANT PATTERN TO BE POSITIVE
#####################################
if (mean(omega.eigen$vectors[,1]) < 0) omega.eigen$vectors[,1] = -omega.eigen$vectors[,1]

#####################################
#### FORCE FIRST UNDEFINED POINT TO BE POSITIVE
#####################################
first.defined.point = min(which(!is.na(omega.eigen$vectors[,1])))
for ( n in 1:neof) if (omega.eigen$vectors[first.defined.point,n] < 0) omega.eigen$vectors[,n] = - omega.eigen$vectors[,n]

#####################################
#### CONSTRUCT BASIS VECTORS AND PSEUDOINVERSE
#####################################
bvec          =  weight/sqrt(sum(weight^2))
eof           =  array(NA,dim=c(nlon*nlat,neof))
eofi          =  array(NA,dim=c(nlon*nlat,neof))
eof [keeper,] =  omega.eigen$vectors[,1:neof] / bvec 
eofi[keeper,] =  omega.eigen$vectors[,1:neof] * bvec 
eval          =  omega.eigen$values 

totwavnum     = 1/(eval-eval[ngood])
totwavnum[1]  = 0
eval[1]       = Inf


#####################################
#### NORMALIZE EIGENVECTORS TO HAVE UNIT AREA WEIGHTED VARIANCE
#####################################
area = bvec^2

#####################################
#### LENGTH SCALE
#####################################
ae               = 6.37e3
length.scale     = ae*pi/sqrt(totwavnum[1:max.num.lapl])

#####################################
#### PROJECT ONTO DATA
#####################################
if ( !is.null(field)) {
	if (rm.clim) clim = rowMeans(field) else clim = rep(0,nlon*nlat)
	field      = field - clim
	pc         = t(field[keeper,]) %*% eofi[keeper,]
	fexpvar    = colSums(pc^2)/sum(field[keeper]^2*area)
	
	list(eof=eof,eofi=eofi,eval=eval,ngood=ngood,bad=!keeper,neof=neof,
	   pc=pc,fexpvar=fexpvar,lon=lon,lat=lat,idate=idate,dname=dname,
	   totwavnum=totwavnum,area=area,length.scale)
} else {
	list(eof=eof,eofi=eofi,eval=eval,ngood=ngood,bad=!keeper,neof=neof,
	   lon=lon,lat=lat,idate=idate,dname=dname,totwavnum=totwavnum,area=area,
	   length.scale)
}

}