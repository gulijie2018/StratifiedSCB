######## Coverage frequencies of SCBs under proportional allocation based on bootstrap method ######
## Two strata (Table 5)

timestart=Sys.time()
library(stats)
library(MASS)

kintegral=function(u)  ## distrution of kernal function
{
  Ktilde = (3/16*u^5-5/8*u^3+15/16*u+1/2)*(abs(u)<=1)+(u>1)
}

kernelcdf = function(h,z,zfixed) ##KDE
{
  n = nrow(z)
  m = ncol(zfixed)
  x = matrix(1,n,1)%*%zfixed-z%*%matrix(1,1,m)
  result = kintegral(x/h)
  F = colMeans(result)
}


bpopn=function(y,N)  ## generate bootstrap population from sample
{
n=length(y)
re=floor(N/n)
by1=kronecker(rep(1,re),y)
by2=sample(y,N-nrow(by1),replace=FALSE,prob=NULL)
popn=c(by1,by2)
}


alpha=c(2,3)  ## parameters of Gamma distribution
beta=c(0.5,1)
Ns=c(2000,1000)   ## stratum sizes
N=sum(Ns);N       ## total population size
 
S=2    ##there are two strata

n=150  ## total sample size n=150,300,600,900
r=400  ## (r+1) is the number of grid points
m=1000 ## replications

D = matrix(1,m,2)
v1=rep(0,4)
v2=rep(0,4)

p1=rep(0,4)
p2=rep(0,4)

Ws=as.matrix(Ns/N) ## stratum weights
lamda=(1/n-1/N)^(-1/2) ## corrected scale factor

#### stratum sample sizes under proportional allocation ####
ns=round(n*Ns/N);ns;sum(ns)
indmax=which(ns== max(ns), arr.ind = TRUE)
ns[indmax]=n-sum(ns[-indmax]);ns;sum(ns) ##propotional

set.seed(123456)
############### population #####################
subp1=matrix(rgamma(Ns[1]*m,shape=alpha[1],scale=beta[1]),Ns[1],m)
subp2=matrix(rgamma(Ns[2]*m,shape=alpha[2],scale=beta[2]),Ns[2],m)
popn=rbind(subp1,subp2)

for (j in 1:m)
{
#####################generate sample########################
y1=sample(subp1[,j],ns[1],replace=FALSE,prob=NULL)
y2=sample(subp2[,j],ns[2],replace=FALSE,prob=NULL)

#### grid points ########
a=min(y1,y2);a
b=max(y1,y2);b
leftpoint=a-(b-a)/(n^2)
int=(b-leftpoint)/r
grid=leftpoint+int*c(0:r)

#### sample EDF ########
Fn1=ecdf(y1)(grid)
Fn2=ecdf(y2)(grid)
Fn=Ws[1]*Fn1+Ws[2]*Fn2  ##EDF

############ sample KDE with bandwidth h1/h2 ################
#h1= (quantile(z)[4]-quantile(z)[2])*((1/n-1/N)^(-1/2))^(-2)
#h2= (quantile(z)[4]-quantile(z)[2])*((1/n-1/N)^(-1/2))^(-2/3)

h11= (quantile(y1)[4]-quantile(y1)[2])*((1/ns[1]-1/Ns[1])^(-1/2))^(-2)
h12= (quantile(y2)[4]-quantile(y2)[2])*((1/ns[2]-1/Ns[2])^(-1/2))^(-2)

Fhat1=kernelcdf(h11,data.matrix(y1),t(as.matrix(grid))) 
Fhat2=kernelcdf(h12,data.matrix(y2),t(as.matrix(grid))) 
Fhat=Ws[1]*Fhat1+Ws[2]*Fhat2   #KDE

######################### population distribution ################################
FN1=ecdf(subp1[,j])(grid)
FN2=ecdf(subp2[,j])(grid)
FN=Ws[1]*FN1+Ws[2]*FN2

################ critical value ################

##### generate bootstrap population #####
bpopn1=bpopn(y1,Ns[1])
bpopn2=bpopn(y2,Ns[2])

Lstar1=matrix(0,m,1)
Lstar2=matrix(0,m,1)  
Bs1=matrix(0,(r+1),1)
Bs2=matrix(0,(r+1),1)

for (i in 1:m)
{
### resample from bootstrap population ########
by1=sample(bpopn1,ns[1],replace=FALSE,prob=NULL)
by2=sample(bpopn2,ns[2],replace=FALSE,prob=NULL)

## EDF
bFn1=ecdf(by1)(grid)
bFn2=ecdf(by2)(grid)
bFn=Ws[1]*bFn1+Ws[2]*bFn2

## KDE
bFhat1=kernelcdf(h11,data.matrix(by1),t(as.matrix(grid))) 
bFhat2=kernelcdf(h12,data.matrix(by2),t(as.matrix(grid))) 
bFhat=Ws[1]*bFhat1+Ws[2]*bFhat2

#### taking maximal absolute values ######
Bs1=(abs(bFn-Fn))*lamda
Bs2=(abs(bFhat-Fhat))*lamda
Lstar1[i]=max(Bs1)
Lstar2[i]=max(Bs2)
}  
Lstar1=sort(Lstar1)
Lstar2=sort(Lstar2)
qvalue1=Lstar1[c(0.99*m,0.95*m,0.9*m,0.8*m)]  ## nonsmooth critical values
qvalue2=Lstar2[c(0.99*m,0.95*m,0.9*m,0.8*m)]  ## smooth critical values


###################### SCB coverage ###########################
D[j,1] = max(abs(Fn-FN))
D[j,2] = max(abs(Fhat-FN))

 for (k in 1:4)
{
if (D[j,1]<=qvalue1[k]/lamda)
 v1[k]=v1[k]+1
if (D[j,2]<=qvalue2[k]/lamda)
 v2[k]=v2[k]+1
 }
}
p1=v1/m;p2=v2/m
p1;p2  ##results listed in Table 5

timeend=Sys.time()
runtime=timeend-timestart;runtime

