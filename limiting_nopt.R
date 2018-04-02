###### SCBs under optimal allocation based on limiting method ######
## Four strata (Table 3)

timestart=Sys.time()
library(stats)
library(MASS)

kintegral=function(u)  ## the distribution of kernal funtion
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

rbridge=function(x,m) ##generate Brownian Bridge
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
Ns=c(200,1000,1400,400)    ## each stratum size
N=sum(Ns);N   #total population size
 
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
qvalue1=matrix(0,m,4)
qvalue2=matrix(0,m,4)

SRSFhat=matrix(0,(r+1),m)
SRSFn=matrix(0,(r+1),m)


Ws=as.matrix(Ns/N) ## the stratum weights
lamda=(1/n-1/N)^(-1/2)  ## corrected scale factor
SRSqvalue=c(1.63,1.36,1.22,1.07)  ##quantiles of Kolmogorov distribution 

pns=round(n*Ns/N);pns;sum(pns)
indmax=which(pns== max(pns), arr.ind = TRUE)
pns[indmax]=n-sum(pns[-indmax]);pns;sum(pns) ##propotional


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
################################ optimal allocation #########################
## pilot sample
py1=sample(subp1[,j],pns[1],replace=FALSE,prob=NULL)
py2=sample(subp2[,j],pns[2],replace=FALSE,prob=NULL)
py3=sample(subp3[,j],pns[3],replace=FALSE,prob=NULL)
py4=sample(subp4[,j],pns[4],replace=FALSE,prob=NULL)

a1=min(py1);a1
b1=max(py1);b1
c1=(b1-a1)/r
grid1=a1+c(0:r)*c1
Fn_1=ecdf(py1)(grid1)
S1=sqrt(sum(Fn_1*(1-Fn_1)*c1)*(Ns[1]-1)/Ns[1])

a2=min(py2);a2
b2=max(py2);b2
c2=(b2-a2)/r
grid2=a2+c(0:r)*c2
Fn_2=ecdf(py2)(grid2)
S2=sqrt(sum(Fn_2*(1-Fn_2)*c2)*(Ns[2]-1)/Ns[2])

a3=min(py3);a3
b3=max(py3);b3
c3=(b3-a3)/r
grid3=a3+c(0:r)*c3
Fn_3=ecdf(py3)(grid3)
S3=sqrt(sum(Fn_3*(1-Fn_3)*c3)*(Ns[3]-1)/Ns[3])

a4=min(py4);a4
b4=max(py4);b4
c4=(b4-a4)/r
grid4=a4+c(0:r)*c4
Fn_4=ecdf(py4)(grid4)
S4=sqrt(sum(Fn_4*(1-Fn_4)*c4)*(Ns[4]-1)/Ns[4])

#### stratum sample sizes under optimal allocation ####
Sh=c(S1,S2,S3,S4)
ns=n*(Sh*Ns)/sum(Sh*Ns);ns  
ns=round(ns);ns;sum(ns)
ind_max=which(ns== max(ns), arr.ind = TRUE)
ns[ind_max]=n-sum(ns[-ind_max]);ns;sum(ns) 


##################### generate stratified sample ########################
y1=sample(subp1[,j],ns[1],replace=FALSE,prob=NULL)
y2=sample(subp2[,j],ns[2],replace=FALSE,prob=NULL)
y3=sample(subp3[,j],ns[3],replace=FALSE,prob=NULL)
y4=sample(subp4[,j],ns[4],replace=FALSE,prob=NULL)
strsample=c(y1,y2,y3,y4)

#### grid points #######
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

######################### superpopulation distribution ################################
F1[,j]=pgamma(grid[,j],shape=alpha[1],scale=beta[1])
F2[,j]=pgamma(grid[,j],shape=alpha[2],scale=beta[2])
F3[,j]=pgamma(grid[,j],shape=alpha[3],scale=beta[3])
F4[,j]=pgamma(grid[,j],shape=alpha[4],scale=beta[4])
F[,j]=Ws[1]*F1[,j]+Ws[2]*F2[,j]+Ws[3]*F3[,j]+Ws[4]*F4[,j]


################ critical value ################

#Ws=as.matrix(Ns/N)
#lamda=(1/n-1/N)^(-1/2)
invlamdas=(1/ns-1/Ns)^(1/2)
weight1=(Ws*invlamdas)*lamda ## weights for FN
weight2=Ws*sqrt(n/ns)  ## weights for F
B1=rbridge(Fn1[,j],L)
B2=rbridge(Fn2[,j],L)
B3=rbridge(Fn3[,j],L)
B4=rbridge(Fn4[,j],L)
Bstar1=weight1[1]*B1+weight1[2]*B2+weight1[3]*B3+weight1[4]*B4
Bstar2=weight2[1]*B1+weight2[2]*B2+weight2[3]*B3+weight2[4]*B4
Lstar1=apply(abs(Bstar1),2,max)
Lstar2=apply(abs(Bstar2),2,max)
Lstar1=sort(Lstar1)
Lstar2=sort(Lstar2)
qvalue1[j,]=Lstar1[c(0.99*L,0.95*L,0.9*L,0.8*L)] ##critical values of SCB for FN
qvalue2[j,]=Lstar2[c(0.99*L,0.95*L,0.9*L,0.8*L)] ##critical values of SCB for F

################## treat the stratified sample as SRS ####################
SRSh1= (quantile(strsample)[4]-quantile(strsample)[2])*((1/n-1/N)^(-1/2))^(-2)
#SRSh2= (quantile(strsample)[4]-quantile(strsample)[2])*((1/n-1/N)^(-1/2))^(-2/3)
SRSFhat[,j]=kernelcdf(SRSh1,data.matrix(strsample),t(as.matrix(grid[,j])))  ## KDE
SRSFn[,j]=ecdf(strsample)(grid[,j])   ## EDF

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
if (D[j,1]<=qvalue1[j,k]/lamda)
 v1[k]=v1[k]+1
if (D[j,2]<=qvalue1[j,k]/lamda)
 v2[k]=v2[k]+1
if (D[j,3]<=SRSqvalue[k]/lamda)
 v3[k]=v3[k]+1
if (D[j,4]<=SRSqvalue[k]/lamda)
 v4[k]=v4[k]+1

if (D[j,5]<=qvalue2[j,k]/lamda)
 v5[k]=v5[k]+1
if (D[j,6]<=qvalue2[j,k]/lamda)
 v6[k]=v6[k]+1
if (D[j,7]<=SRSqvalue[k]/lamda)
 v7[k]=v7[k]+1
if (D[j,8]<=SRSqvalue[k]/lamda)
 v8[k]=v8[k]+1

 }
}
p1=v1/m;p2=v2/m;p3=v3/m;p4=v4/m
p5=v5/m;p6=v6/m;p7=v7/m;p8=v8/m

p1;p2;p3;p4;p5;p6;p7;p8  ## results listed in Table 3

timeend=Sys.time()
runtime=timeend-timestart;runtime

######################## 95% SCB ########################
j=990 
lamda=(1/n-1/N)^(-1/2)
up1 = matrix(1,(r+1),1)
lw1 = matrix(1,(r+1),1)
for (i in 1:(r+1))
{
up1[i] = min(1, Fn[i,j]+qvalue1[2]/lamda)  ## upper bound of SCB based on EDF
lw1[i] = max(0, Fn[i,j]-qvalue1[2]/lamda)  ## lower bound of SCB based on EDF
}

up2 = matrix(1,(r+1),1)
lw2 = matrix(1,(r+1),1)
for (i in 1:(r+1))
{
up2[i] = min(1, Fhat[i,j]+qvalue1[2]/lamda) ## upper bound of SCB based on KDE
lw2[i] = max(0, Fhat[i,j]-qvalue1[2]/lamda) ## lower bound of SCB based on KDE
} 


###################### plot SCB shown in Figure 7(a) #############################
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






