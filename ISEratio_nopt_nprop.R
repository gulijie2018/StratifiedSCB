####### Codes for Table 1 and Figures 1-6 ##########
timestart=Sys.time()
library(stats)
library(MASS)

kintegral=function(u) ## the distribution of kernal funtion
{
  Ktilde = (3/16*u^5-5/8*u^3+15/16*u+1/2)*(abs(u)<=1)+(u>1)
}

kernelcdf = function(N,z,zfixed)  ## KDE
{
  n = nrow(z)
  m = ncol(zfixed)
  h= (quantile(z)[4]-quantile(z)[2])*((1/n-1/N)^(-1/2))^(-2)    ##h1
  #h= (quantile(z)[4]-quantile(z)[2])*((1/n-1/N)^(-1/2))^(-2/3) ##h2
  x = matrix(1,n,1)%*%zfixed-z%*%matrix(1,1,m)
  result = kintegral(x/h)
  F = colMeans(result)
}


alpha=c(1.41,1.45,0.67,0.76) ## parameters of Gamma distribution
beta=c(0.66,2.25,2.98,9.56)
Ns=c(200,1000,1400,400)  ## each stratum size
N=sum(Ns);N   ## population size
Ws=as.matrix(Ns/N) ## the stratum weights

S=4    ## four strata
r=400  ## (r+1) is the number of grid points
m=1000 ## replications

ISE1 = matrix(1,m,4)
ISE2 = matrix(1,m,4)
ISE3 = matrix(1,m,4)
ISE4 = matrix(1,m,4)
ISE5 = matrix(1,m,4)
ISE6 = matrix(1,m,4)
ISE7 = matrix(1,m,4)
ISE8 = matrix(1,m,4)

ratio1 = matrix(1,m,4)
ratio2 = matrix(1,m,4)
ratio3 = matrix(1,m,4)
ratio4 = matrix(1,m,4)
ratio5 = matrix(1,m,4)
ratio6 = matrix(1,m,4)
ratio7 = matrix(1,m,4)
ratio8 = matrix(1,m,4)

MISE=matrix(1,4,8)


set.seed(123456)
############### generate population
subp1=matrix(rgamma(Ns[1]*m,shape=alpha[1],scale=beta[1]),Ns[1],m)
subp2=matrix(rgamma(Ns[2]*m,shape=alpha[2],scale=beta[2]),Ns[2],m)
subp3=matrix(rgamma(Ns[3]*m,shape=alpha[3],scale=beta[3]),Ns[3],m)
subp4=matrix(rgamma(Ns[4]*m,shape=alpha[4],scale=beta[4]),Ns[4],m)
popn=rbind(subp1,subp2,subp3,subp4)

n=c(150,300,600,900)  ## total sample size n has four options
lamda=(1/n-1/N)^(-1/2) ## corrected scale factor

for (i in 1:4)
{ 
## stratum sample size under proportional allocation
pns=round(n[i]*Ns/N);pns;sum(pns)
indmax=which(pns== max(pns), arr.ind = TRUE)
pns[indmax]=n[i]-sum(pns[-indmax]);pns;sum(pns) 


for (j in 1:m)
{
################# generate sample under proportional allocation ###################
py1=sample(subp1[,j],pns[1],replace=FALSE,prob=NULL)
py2=sample(subp2[,j],pns[2],replace=FALSE,prob=NULL)
py3=sample(subp3[,j],pns[3],replace=FALSE,prob=NULL)
py4=sample(subp4[,j],pns[4],replace=FALSE,prob=NULL)

#### sample ecdf ########
pa=min(py1,py2,py3,py4);pa
pb=max(py1,py2,py3,py4);pb
pleftpoint=pa-(pb-pa)/(n[i]^2)
pint=(pb-pleftpoint)/r
pgrid=pleftpoint+pint*c(0:r)

pFn1=ecdf(py1)(pgrid)
pFn2=ecdf(py2)(pgrid)
pFn3=ecdf(py3)(pgrid)
pFn4=ecdf(py4)(pgrid)
pFn=Ws[1]*pFn1+Ws[2]*pFn2+Ws[3]*pFn3+Ws[4]*pFn4  ##ecdf

#### sample KDE ########
pFhat1=kernelcdf(Ns[1],data.matrix(py1),t(as.matrix(pgrid))) 
pFhat2=kernelcdf(Ns[2],data.matrix(py2),t(as.matrix(pgrid)))
pFhat3=kernelcdf(Ns[3],data.matrix(py3),t(as.matrix(pgrid))) 
pFhat4=kernelcdf(Ns[4],data.matrix(py4),t(as.matrix(pgrid))) 
pFhat=Ws[1]*pFhat1+Ws[2]*pFhat2+Ws[3]*pFhat3+Ws[4]*pFhat4    ##KDE

######################### population ecdf ################################
pFN1=ecdf(subp1[,j])(pgrid)
pFN2=ecdf(subp2[,j])(pgrid)
pFN3=ecdf(subp3[,j])(pgrid)
pFN4=ecdf(subp4[,j])(pgrid)
pFN=Ws[1]*pFN1+Ws[2]*pFN2+Ws[3]*pFN3+Ws[4]*pFN4

######################### population cdf ################################
pF1=pgamma(pgrid,shape=alpha[1],scale=beta[1])
pF2=pgamma(pgrid,shape=alpha[2],scale=beta[2])
pF3=pgamma(pgrid,shape=alpha[3],scale=beta[3])
pF4=pgamma(pgrid,shape=alpha[4],scale=beta[4])
pF=Ws[1]*pF1+Ws[2]*pF2+Ws[3]*pF3+Ws[4]*pF4


############## optimal allocation ###############################
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

Sh=c(S1,S2,S3,S4)
ns=n[i]*(Sh*Ns)/sum(Sh*Ns);ns  ##optimal allocation
ns=round(ns);ns;sum(ns)
ind_max=which(ns== max(ns), arr.ind = TRUE)
ns[ind_max]=n[i]-sum(ns[-ind_max]);ns;sum(ns) 

################### generate sample under optimal allocation ########################
y1=sample(subp1[,j],ns[1],replace=FALSE,prob=NULL)
y2=sample(subp2[,j],ns[2],replace=FALSE,prob=NULL)
y3=sample(subp3[,j],ns[3],replace=FALSE,prob=NULL)
y4=sample(subp4[,j],ns[4],replace=FALSE,prob=NULL)

#### sample ecdf ########
a=min(y1,y2,y3,y4);a
b=max(y1,y2,y3,y4);b
leftpoint=a-(b-a)/(n[i]^2)
int=(b-leftpoint)/r
grid=leftpoint+int*c(0:r)

Fn1=ecdf(y1)(grid)
Fn2=ecdf(y2)(grid)
Fn3=ecdf(y3)(grid)
Fn4=ecdf(y4)(grid)
Fn=Ws[1]*Fn1+Ws[2]*Fn2+Ws[3]*Fn3+Ws[4]*Fn4  ##ecdf

#### sample KDE ########
Fhat1=kernelcdf(Ns[1],data.matrix(y1),t(as.matrix(grid))) 
Fhat2=kernelcdf(Ns[2],data.matrix(y2),t(as.matrix(grid))) 
Fhat3=kernelcdf(Ns[3],data.matrix(y3),t(as.matrix(grid))) 
Fhat4=kernelcdf(Ns[4],data.matrix(y4),t(as.matrix(grid))) 
Fhat=Ws[1]*Fhat1+Ws[2]*Fhat2+Ws[3]*Fhat3+Ws[4]*Fhat4   #KDE

######################### population ecdf ################################
FN1=ecdf(subp1[,j])(grid)
FN2=ecdf(subp2[,j])(grid)
FN3=ecdf(subp3[,j])(grid)
FN4=ecdf(subp4[,j])(grid)
FN=Ws[1]*FN1+Ws[2]*FN2+Ws[3]*FN3+Ws[4]*FN4

######################### population cdf ################################
F1=pgamma(grid,shape=alpha[1],scale=beta[1])
F2=pgamma(grid,shape=alpha[2],scale=beta[2])
F3=pgamma(grid,shape=alpha[3],scale=beta[3])
F4=pgamma(grid,shape=alpha[4],scale=beta[4])
F=Ws[1]*F1+Ws[2]*F2+Ws[3]*F3+Ws[4]*F4

## optimal ISE ####
ISE1[j,i] = sum((Fn-FN)^2*int)
ISE2[j,i] = sum((Fhat-FN)^2*int)
ISE3[j,i] = sum((Fn-F)^2*int)
ISE4[j,i] = sum((Fhat-F)^2*int)

## proportional ISE
ISE5[j,i] = sum((pFn-pFN)^2*pint)
ISE6[j,i] = sum((pFhat-pFN)^2*pint)
ISE7[j,i] = sum((pFn-pF)^2*pint)
ISE8[j,i] = sum((pFhat-pF)^2*pint)
}

## compare between optimal and proportional allocations
ratio1[,i]=ISE1[,i]/ISE5[,i]
ratio2[,i]=ISE2[,i]/ISE6[,i]
ratio3[,i]=ISE3[,i]/ISE7[,i]
ratio4[,i]=ISE4[,i]/ISE8[,i]

## compare between ecdf and KDE 
ratio5[,i]=ISE2[,i]/ISE1[,i]
ratio6[,i]=ISE6[,i]/ISE5[,i]
ratio7[,i]=ISE4[,i]/ISE3[,i]
ratio8[,i]=ISE8[,i]/ISE7[,i]

}
timeend=Sys.time()
runtime=timeend-timestart;runtime
apply(ISE2,2,mean)
apply(ISE6,2,mean)
apply(ISE2,2,mean)/apply(ISE6,2,mean)
apply(ISE2,2,mean)/apply(ISE1,2,mean)
apply(ISE6,2,mean)/apply(ISE5,2,mean)

apply(ISE4,2,mean)
apply(ISE8,2,mean)
apply(ISE4,2,mean)/apply(ISE8,2,mean)
apply(ISE4,2,mean)/apply(ISE3,2,mean)
apply(ISE8,2,mean)/apply(ISE7,2,mean)

boxplot(list("150"=ratio1[,1],"300"=ratio1[,2],"600"=ratio1[,3],
"900"=ratio1[,4]),xlab="Sample Size",ylab="ISE_ratio",ylim=c(0,15),
col=c("red","blue","green","purple"))
abline(1,0,lty="dotdash")

boxplot(list("150"=ratio2[,1],"300"=ratio2[,2],"600"=ratio2[,3],
"900"=ratio2[,4]),xlab="Sample Size",ylab="ISE_ratio",ylim=c(0,15),
col=c("red","blue","green","purple"))
abline(1,0,lty="dotdash")

boxplot(list("150"=ratio3[,1],"300"=ratio3[,2],"600"=ratio3[,3],
"900"=ratio3[,4]),xlab="Sample Size",ylab="ISE_ratio",ylim=c(0,15),
col=c("red","blue","green","purple"))
abline(1,0,lty="dotdash")

boxplot(list("150"=ratio4[,1],"300"=ratio4[,2],"600"=ratio4[,3],
"900"=ratio4[,4]),xlab="Sample Size",ylab="ISE_ratio",ylim=c(0,15),
col=c("red","blue","green","purple"))
abline(1,0,lty="dotdash")

boxplot(list("150"=ratio5[,1],"300"=ratio5[,2],"600"=ratio5[,3],
"900"=ratio5[,4]),xlab="Sample Size",ylab="ISE_ratio",ylim=c(0.84,1.08),
col=c("red","blue","green","purple"))
abline(1,0,lty="dotdash")

boxplot(list("150"=ratio6[,1],"300"=ratio6[,2],"600"=ratio6[,3],
"900"=ratio6[,4]),xlab="Sample Size",ylab="ISE_ratio",ylim=c(0.84,1.08),
col=c("red","blue","green","purple"))
abline(1,0,lty="dotdash")

boxplot(list("150"=ratio7[,1],"300"=ratio7[,2],"600"=ratio7[,3],
"900"=ratio7[,4]),xlab="Sample Size",ylab="ISE_ratio",ylim=c(0.84,1.08),
col=c("red","blue","green","purple"))
abline(1,0,lty="dotdash")

boxplot(list("150"=ratio8[,1],"300"=ratio8[,2],"600"=ratio8[,3],
"900"=ratio8[,4]),xlab="Sample Size",ylab="ISE_ratio",ylim=c(0.84,1.08),
col=c("red","blue","green","purple"))
abline(1,0,lty="dotdash")

boxplot(list("150"=ISE2[,1],"300"=ISE2[,2],"600"=ISE2[,3],
"900"=ISE2[,4]),xlab="Sample Size",ylab="ISE",ylim=c(0,0.08),
col=c("red","blue","green","purple"))
abline(0,0,lty="dotdash")

boxplot(list("150"=ISE6[,1],"300"=ISE6[,2],"600"=ISE6[,3],
"900"=ISE6[,4]),xlab="Sample Size",ylab="ISE",ylim=c(0,0.08),
col=c("red","blue","green","purple"))
abline(0,0,lty="dotdash")

boxplot(list("150"=ISE4[,1],"300"=ISE4[,2],"600"=ISE4[,3],
"900"=ISE4[,4]),xlab="Sample Size",ylab="ISE",ylim=c(0,0.095),
col=c("red","blue","green","purple"))
abline(0,0,lty="dotdash")

boxplot(list("150"=ISE8[,1],"300"=ISE8[,2],"600"=ISE8[,3],
"900"=ISE8[,4]),xlab="Sample Size",ylab="ISE",ylim=c(0,0.095),
col=c("red","blue","green","purple"))
abline(0,0,lty="dotdash")


