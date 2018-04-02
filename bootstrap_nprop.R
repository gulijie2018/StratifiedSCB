######## Coverage frequencies of SCBs under proportional allocation based on bootstrap method ######
## Four strata (Table 4)

timestart=Sys.time()
library(stats)
library(MASS)

kintegral=function(u)  ## distrution of kernal function
{
  Ktilde = (3/16*u^5-5/8*u^3+15/16*u+1/2)*(abs(u)<=1)+(u>1)
}

kernelcdf = function(h,z,zfixed)  ## KDE
{
  n = nrow(z)
  m = ncol(zfixed)
  x = matrix(1,n,1)%*%zfixed-z%*%matrix(1,1,m)
  result = kintegral(x/h)
  F = colMeans(result)
}

bpopn=function(y,N) ## generate bootstrap population from sample
{
n=length(y)
re=floor(N/n)
by1=kronecker(rep(1,re),y)
by2=sample(y,N-nrow(by1),replace=FALSE,prob=NULL)
popn=c(by1,by2)
}

alpha=c(1.41,1.45,0.67,0.76) ## parameters of Gamma distribution
beta=c(0.66,2.25,2.98,9.56)
Ns=c(200,1000,1400,400)  ## stratum sizes
N=sum(Ns);N   ## total population size
 
S=4    ##there are four strata

n=150  ## total sample size n=150,300,600,900
r=400  ## (r+1) is the number of grid points
m=1000 ## replications

D = matrix(1,m,8)
v1=rep(0,4)
v2=rep(0,4)

p1=rep(0,4)
p2=rep(0,4)

Fn1=matrix(0,(r+1),m)
Fn2=matrix(0,(r+1),m)
Fn3=matrix(0,(r+1),m)
Fn4=matrix(0,(r+1),m)
Fn=matrix(0,(r+1),m)

Fhat1=matrix(0,(r+1),m)
Fhat2=matrix(0,(r+1),m)
Fhat3=matrix(0,(r+1),m)
Fhat4=matrix(0,(r+1),m)
Fhat=matrix(0,(r+1),m)

FN1=matrix(0,(r+1),m)
FN2=matrix(0,(r+1),m)
FN3=matrix(0,(r+1),m)
FN4=matrix(0,(r+1),m)
FN=matrix(0,(r+1),m)

grid=matrix(0,(r+1),m)
qvalue1=matrix(0,m,4)
qvalue2=matrix(0,m,4)

Ws=as.matrix(Ns/N)     ## stratum weights
lamda=(1/n-1/N)^(-1/2) ## corrected scale factor

#### stratum sample sizes under proportional allocation ####
ns=round(n*Ns/N);ns;sum(ns)
indmax=which(ns== max(ns), arr.ind = TRUE)
ns[indmax]=n-sum(ns[-indmax]);ns;sum(ns) ##propotional


set.seed(123456)
############### population #######################
subp1=matrix(rgamma(Ns[1]*m,shape=alpha[1],scale=beta[1]),Ns[1],m)
subp2=matrix(rgamma(Ns[2]*m,shape=alpha[2],scale=beta[2]),Ns[2],m)
subp3=matrix(rgamma(Ns[3]*m,shape=alpha[3],scale=beta[3]),Ns[3],m)
subp4=matrix(rgamma(Ns[4]*m,shape=alpha[4],scale=beta[4]),Ns[4],m)
popn=rbind(subp1,subp2,subp3,subp4)


for (j in 1:m)
{
#####################generate sample########################
y1=sample(subp1[,j],ns[1],replace=FALSE,prob=NULL)
y2=sample(subp2[,j],ns[2],replace=FALSE,prob=NULL)
y3=sample(subp3[,j],ns[3],replace=FALSE,prob=NULL)
y4=sample(subp4[,j],ns[4],replace=FALSE,prob=NULL)

#### grid points ########
a=min(y1,y2,y3,y4);a
b=max(y1,y2,y3,y4);b
leftpoint=a-(b-a)/(n^2)
int=(b-leftpoint)/r
grid[,j]=leftpoint+int*c(0:r)

#### sample EDF ########
Fn1[,j]=ecdf(y1)(grid[,j])
Fn2[,j]=ecdf(y2)(grid[,j])
Fn3[,j]=ecdf(y3)(grid[,j])
Fn4[,j]=ecdf(y4)(grid[,j])
Fn[,j]=Ws[1]*Fn1[,j]+Ws[2]*Fn2[,j]+Ws[3]*Fn3[,j]+Ws[4]*Fn4[,j]  ##EDF

#### sample KDE with bandwidth h1/h2 ########
#h1= (quantile(z)[4]-quantile(z)[2])*((1/n-1/N)^(-1/2))^(-2)
#h2= (quantile(z)[4]-quantile(z)[2])*((1/n-1/N)^(-1/2))^(-2/3)

h11= (quantile(y1)[4]-quantile(y1)[2])*((1/ns[1]-1/Ns[1])^(-1/2))^(-2)
h12= (quantile(y2)[4]-quantile(y2)[2])*((1/ns[2]-1/Ns[2])^(-1/2))^(-2)
h13= (quantile(y3)[4]-quantile(y3)[2])*((1/ns[3]-1/Ns[3])^(-1/2))^(-2)
h14= (quantile(y4)[4]-quantile(y4)[2])*((1/ns[4]-1/Ns[4])^(-1/2))^(-2)

Fhat1[,j]=kernelcdf(h11,data.matrix(y1),t(as.matrix(grid[,j]))) 
Fhat2[,j]=kernelcdf(h12,data.matrix(y2),t(as.matrix(grid[,j]))) 
Fhat3[,j]=kernelcdf(h13,data.matrix(y3),t(as.matrix(grid[,j]))) 
Fhat4[,j]=kernelcdf(h14,data.matrix(y4),t(as.matrix(grid[,j]))) 
Fhat[,j]=Ws[1]*Fhat1[,j]+Ws[2]*Fhat2[,j]+Ws[3]*Fhat3[,j]+Ws[4]*Fhat4[,j]   #KDE

######################### population distribution ################################
FN1[,j]=ecdf(subp1[,j])(grid[,j])
FN2[,j]=ecdf(subp2[,j])(grid[,j])
FN3[,j]=ecdf(subp3[,j])(grid[,j])
FN4[,j]=ecdf(subp4[,j])(grid[,j])
FN[,j]=Ws[1]*FN1[,j]+Ws[2]*FN2[,j]+Ws[3]*FN3[,j]+Ws[4]*FN4[,j]


################ critical value ################

##### generate bootstrap population #####
bpopn1=bpopn(y1,Ns[1])
bpopn2=bpopn(y2,Ns[2])
bpopn3=bpopn(y3,Ns[3])
bpopn4=bpopn(y4,Ns[4])

Lstar1=matrix(0,m,1)
Lstar2=matrix(0,m,1)  
Bs1=matrix(0,(r+1),1)
Bs2=matrix(0,(r+1),1)

for (i in 1:m)
{
### resample from bootstrap population ########
by1=sample(bpopn1,ns[1],replace=FALSE,prob=NULL)
by2=sample(bpopn2,ns[2],replace=FALSE,prob=NULL)
by3=sample(bpopn3,ns[3],replace=FALSE,prob=NULL)
by4=sample(bpopn4,ns[4],replace=FALSE,prob=NULL)

## EDF
bFn1=ecdf(by1)(grid[,j])
bFn2=ecdf(by2)(grid[,j])
bFn3=ecdf(by3)(grid[,j])
bFn4=ecdf(by4)(grid[,j])
bFn=Ws[1]*bFn1+Ws[2]*bFn2+Ws[3]*bFn3+Ws[4]*bFn4

## KDE
bFhat1=kernelcdf(h11,data.matrix(by1),t(as.matrix(grid[,j]))) 
bFhat2=kernelcdf(h12,data.matrix(by2),t(as.matrix(grid[,j]))) 
bFhat3=kernelcdf(h13,data.matrix(by3),t(as.matrix(grid[,j]))) 
bFhat4=kernelcdf(h14,data.matrix(by4),t(as.matrix(grid[,j]))) 
bFhat=Ws[1]*bFhat1+Ws[2]*bFhat2+Ws[3]*bFhat3+Ws[4]*bFhat4

#### taking maximal absolute values ######
Bs1=(abs(bFn-Fn[,j]))*lamda
Bs2=(abs(bFhat-Fhat[,j]))*lamda
Lstar1[i]=max(Bs1)
Lstar2[i]=max(Bs2)
}  
Lstar1=sort(Lstar1)
Lstar2=sort(Lstar2)
qvalue1[j,]=Lstar1[c(0.99*m,0.95*m,0.9*m,0.8*m)]  ## nonsmooth critical values
qvalue2[j,]=Lstar2[c(0.99*m,0.95*m,0.9*m,0.8*m)]  ## smooth critical values

###################### SCB coverage ###########################
D[j,1] = max(abs(Fn[,j]-FN[,j]))
D[j,2] = max(abs(Fhat[,j]-FN[,j]))

for (k in 1:4)
{
if (D[j,1]<=qvalue1[j,k]/lamda)
 v1[k]=v1[k]+1
if (D[j,2]<=qvalue2[j,k]/lamda)
 v2[k]=v2[k]+1
 }
}
p1=v1/m;p2=v2/m
p1;p2   ##results listed in Table 4

timeend=Sys.time()
runtime=timeend-timestart;runtime
