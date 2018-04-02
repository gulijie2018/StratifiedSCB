###### SCBs under proportional allocation based on limiting method ######
## Four strata (Table 4)

timestart=Sys.time()
library(stats)
library(MASS)

kintegral=function(u) ## the distribution of kernal funtion
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

rbridge=function(x,m)  ## generate Brownian Bridge
{
r=length(x)
x=as.matrix(x,r,1)
var=c(x[1],diff(x))
w=matrix(0,r,1)
w=sqrt(var)*matrix(rnorm(r*m,0,1),r,m)
B=apply(w,2,cumsum)-x%*%colSums(w)
return(B)
}


alpha=c(1.41,1.45,0.67,0.76)  ## parameters of Gamma distribution
beta=c(0.66,2.25,2.98,9.56)
Ns=c(200,1000,1400,400)       ## each stratum size
N=sum(Ns);N                   ## total population size
 
S=4    ##there are four strata

n=150  ## total sample size n=150,300,600,900
r=400  ## (r+1) is the number of grid points
m=1000 ## replication for coverage frequencies
L=1000 ## replication for critical values

D = matrix(1,m,8)
v1=rep(0,4)
v2=rep(0,4)
v3=rep(0,4)
v4=rep(0,4)
v5=rep(0,4)
v6=rep(0,4)
v7=rep(0,4)
v8=rep(0,4)

p1=rep(0,4)
p2=rep(0,4)
p3=rep(0,4)
p4=rep(0,4)
p5=rep(0,4)
p6=rep(0,4)
p7=rep(0,4)
p8=rep(0,4)

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

F1=matrix(0,(r+1),m)
F2=matrix(0,(r+1),m)
F3=matrix(0,(r+1),m)
F4=matrix(0,(r+1),m)
F=matrix(0,(r+1),m)


grid=matrix(0,(r+1),m)
qvalue=matrix(0,m,4)

SRSFhat=matrix(0,(r+1),m)
SRSFn=matrix(0,(r+1),m)


Ws=as.matrix(Ns/N)                ## the stratum weights
lamda=(1/n-1/N)^(-1/2)            ## corrected scale factor
SRSqvalue=c(1.63,1.36,1.22,1.07)  ##quantiles of Kolmogorov distribution 

#### stratum sample sizes under proportional allocation #######
ns=round(n*Ns/N);ns;sum(ns)
indmax=which(ns== max(ns), arr.ind = TRUE)
ns[indmax]=n-sum(ns[-indmax]);ns;sum(ns) 


set.seed(123456)
############### popupation #################
subp1=matrix(rgamma(Ns[1]*m,shape=alpha[1],scale=beta[1]),Ns[1],m)
subp2=matrix(rgamma(Ns[2]*m,shape=alpha[2],scale=beta[2]),Ns[2],m)
subp3=matrix(rgamma(Ns[3]*m,shape=alpha[3],scale=beta[3]),Ns[3],m)
subp4=matrix(rgamma(Ns[4]*m,shape=alpha[4],scale=beta[4]),Ns[4],m)
popn=rbind(subp1,subp2,subp3,subp4)

################## SRS ####################
SRSyn=popn[sample(1:N,n,replace=FALSE,prob=NULL),]

for (j in 1:m)
{
##################### generate stratified sample ########################
y1=sample(subp1[,j],ns[1],replace=FALSE,prob=NULL)
y2=sample(subp2[,j],ns[2],replace=FALSE,prob=NULL)
y3=sample(subp3[,j],ns[3],replace=FALSE,prob=NULL)
y4=sample(subp4[,j],ns[4],replace=FALSE,prob=NULL)
strsample=c(y1,y2,y3,y4)

#### grid points ######
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

#### sample KDE with bandwidth h1 or h2 ########
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

######################### superpopulation distribution ############################
F1[,j]=pgamma(grid[,j],shape=alpha[1],scale=beta[1])
F2[,j]=pgamma(grid[,j],shape=alpha[2],scale=beta[2])
F3[,j]=pgamma(grid[,j],shape=alpha[3],scale=beta[3])
F4[,j]=pgamma(grid[,j],shape=alpha[4],scale=beta[4])
F[,j]=Ws[1]*F1[,j]+Ws[2]*F2[,j]+Ws[3]*F3[,j]+Ws[4]*F4[,j]


################ critical value ################
weight=sqrt(Ws) ## under proportional allocation
B1=rbridge(Fn1[,j],L)
B2=rbridge(Fn2[,j],L)
B3=rbridge(Fn3[,j],L)
B4=rbridge(Fn4[,j],L)
Bstar1=weight[1]*B1+weight[2]*B2+weight[3]*B3+weight[4]*B4  ## weighted sum of Brownian bridges
Lstar=apply(abs(Bstar1),2,max) ## taking maximal absolute values 
Lstar=sort(Lstar)
qvalue[j,]=Lstar[c(0.99*L,0.95*L,0.9*L,0.8*L)]  ## critical values of SCB

################## treat the stratified sample as SRS ####################
SRSh1= (quantile(strsample)[4]-quantile(strsample)[2])*((1/n-1/N)^(-1/2))^(-2)
#SRSh2= (quantile(strsample)[4]-quantile(strsample)[2])*((1/n-1/N)^(-1/2))^(-2/3)
SRSFhat[,j]=kernelcdf(SRSh1,data.matrix(strsample),t(as.matrix(grid[,j]))) ##KDE
SRSFn[,j]=ecdf(strsample)(grid[,j])   ##EDF

###################### SCB coverage ###########################
D[j,1] = max(abs(Fn[,j]-FN[,j]))
D[j,2] = max(abs(Fhat[,j]-FN[,j]))
D[j,3] = max(abs(SRSFn[,j]-FN[,j]))
D[j,4] = max(abs(SRSFhat[,j]-FN[,j]))

D[j,5] = max(abs(Fn[,j]-F[,j]))
D[j,6] = max(abs(Fhat[,j]-F[,j]))
D[j,7] = max(abs(SRSFn[,j]-F[,j]))
D[j,8] = max(abs(SRSFhat[,j]-F[,j]))

 for (k in 1:4)
{
if (D[j,1]<=qvalue[j,k]/lamda)
 v1[k]=v1[k]+1
if (D[j,2]<=qvalue[j,k]/lamda)
 v2[k]=v2[k]+1
if (D[j,3]<=SRSqvalue[k]/lamda)
 v3[k]=v3[k]+1
if (D[j,4]<=SRSqvalue[k]/lamda)
 v4[k]=v4[k]+1

if (D[j,5]<=qvalue[j,k]/lamda)
 v5[k]=v5[k]+1
if (D[j,6]<=qvalue[j,k]/lamda)
 v6[k]=v6[k]+1
if (D[j,7]<=SRSqvalue[k]/lamda)
 v7[k]=v7[k]+1
if (D[j,8]<=SRSqvalue[k]/lamda)
 v8[k]=v8[k]+1

 }
}
p1=v1/m;p2=v2/m;p3=v3/m;p4=v4/m
p5=v5/m;p6=v6/m;p7=v7/m;p8=v8/m

p1;p2;p3;p4;p5;p6;p7;p8  ## results listed in Table 4

timeend=Sys.time()
runtime=timeend-timestart;runtime

######################## 95% SCB ########################
j=996  
lamda=(1/n-1/N)^(-1/2)
up1 = matrix(1,(r+1),1)
lw1 = matrix(1,(r+1),1)
for (i in 1:(r+1))
{
up1[i] = min(1, Fn[i,j]+qvalue[2]/lamda) ## upper bound of SCB based on EDF
lw1[i] = max(0, Fn[i,j]-qvalue[2]/lamda) ## lower bound of SCB based on EDF
}

up2 = matrix(1,(r+1),1)
lw2 = matrix(1,(r+1),1)
for (i in 1:(r+1))
{
up2[i] = min(1, Fhat[i,j]+qvalue[2]/lamda) ## upper bound of SCB based on KDE
lw2[i] = max(0, Fhat[i,j]-qvalue[2]/lamda) ## lower bound of SCB based on KDE
}

###################### plot SCB shown in Figure 7(b)#############################
plot(grid[,j],FN[,j],type="l",col=1,lty=1,lwd=3,,xlab="",
ylab="",xlim=c(grid[1,j],grid[(r+1),j]))
lines(grid[,j],Fhat[,j],col=3,lty=1,lwd=2)
lines(grid[,j],up2,col=5,lty=5,lwd=2)
lines(grid[,j],lw2,col=5,lty=5,lwd=2)
lines(grid[,j],Fn[,j],col=2,lty=2,lwd=2)
lines(grid[,j],up1,col=4,lty=3,lwd=2)
lines(grid[,j],lw1,col=4,lty=3,lwd=2)

legend("bottomright",c("population CDF","stratified KDE","smooth SCB",
"stratified EDF","nonsmooth SCB"),
lty=c(1,1,5,2,3),lwd=c(3,2,2,2,2),col=c(1,3,5,2,4),title="")





