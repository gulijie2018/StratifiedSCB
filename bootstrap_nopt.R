######## Coverage frequencies of SCBs under optimal allocation based on bootstrap method ######
## Four strata (Table 3)

timestart=Sys.time() 
library(stats)
library(MASS)

kintegral=function(u)  ## distrution of kernal function
{
  Ktilde = (3/16*u^5-5/8*u^3+15/16*u+1/2)*(abs(u)<=1)+(u>1)
}

kernelcdf = function(h,z,zfixed) ## KDE
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


alpha=c(1.41,1.45,0.67,0.76) ## parameters of Gamma distribution
beta=c(0.66,2.25,2.98,9.56)
Ns=c(200,1000,1400,400)  ## stratum sizes
N=sum(Ns);N   ## total population size
 
S=4    ##there are four strata

n=150  ## total sample size n=150,300,600,900
r=400  ## (r+1) is the number of grid points
m=1000 ## replications

D = matrix(1,m,2)
v1=rep(0,4)
v2=rep(0,4)

p1=rep(0,4)
p2=rep(0,4)


Ws=as.matrix(Ns/N)         ## stratum weights
lamda=(1/n-1/N)^(-1/2)     ## corrected scale factor

pns=round(n*Ns/N);pns;sum(pns)
indmax=which(pns== max(pns), arr.ind = TRUE)
pns[indmax]=n-sum(pns[-indmax]);pns;sum(pns) ##propotional

set.seed(123456)
############### population #####################
subp1=matrix(rgamma(Ns[1]*m,shape=alpha[1],scale=beta[1]),Ns[1],m)
subp2=matrix(rgamma(Ns[2]*m,shape=alpha[2],scale=beta[2]),Ns[2],m)
subp3=matrix(rgamma(Ns[3]*m,shape=alpha[3],scale=beta[3]),Ns[3],m)
subp4=matrix(rgamma(Ns[4]*m,shape=alpha[4],scale=beta[4]),Ns[4],m)
popn=rbind(subp1,subp2,subp3,subp4)

for (j in 1:m)
{
################################ allocation #########################
py1=sample(subp1[,j],pns[1],replace=FALSE,prob=NULL)
py2=sample(subp2[,j],pns[2],replace=FALSE,prob=NULL)
py3=sample(subp3[,j],pns[3],replace=FALSE,prob=NULL)
py4=sample(subp4[,j],pns[4],replace=FALSE,prob=NULL)

a1=min(py1);a1
b1=max(py1);b1
c1=(b1-a1)/r
grid1=a1+c(0:r)*c1
Fn_1=ecdf(py1)(grid1)
S1=sqrt(sum(Fn_1*(1-Fn_1)*c1)*(Ns[1]-1)/Ns[1]) ##population variance of stratum 

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


#####################generate sample########################
y1=sample(subp1[,j],ns[1],replace=FALSE,prob=NULL)
y2=sample(subp2[,j],ns[2],replace=FALSE,prob=NULL)
y3=sample(subp3[,j],ns[3],replace=FALSE,prob=NULL)
y4=sample(subp4[,j],ns[4],replace=FALSE,prob=NULL)

#### grid points ######
a=min(y1,y2,y3,y4);a
b=max(y1,y2,y3,y4);b
leftpoint=a-(b-a)/(n^2)
int=(b-leftpoint)/r
grid=leftpoint+int*c(0:r)

#### sample EDF ########
Fn1=ecdf(y1)(grid)
Fn2=ecdf(y2)(grid)
Fn3=ecdf(y3)(grid)
Fn4=ecdf(y4)(grid)
Fn=Ws[1]*Fn1+Ws[2]*Fn2+Ws[3]*Fn3+Ws[4]*Fn4  ##EDF

#### sample KDE with bandwidth h1/h2 ########
#h1= (quantile(z)[4]-quantile(z)[2])*((1/n-1/N)^(-1/2))^(-2)
#h2= (quantile(z)[4]-quantile(z)[2])*((1/n-1/N)^(-1/2))^(-2/3)

h11= (quantile(y1)[4]-quantile(y1)[2])*((1/ns[1]-1/Ns[1])^(-1/2))^(-2)
h12= (quantile(y2)[4]-quantile(y2)[2])*((1/ns[2]-1/Ns[2])^(-1/2))^(-2)
h13= (quantile(y3)[4]-quantile(y3)[2])*((1/ns[3]-1/Ns[3])^(-1/2))^(-2)
h14= (quantile(y4)[4]-quantile(y4)[2])*((1/ns[4]-1/Ns[4])^(-1/2))^(-2)

Fhat1=kernelcdf(h11,data.matrix(y1),t(as.matrix(grid))) 
Fhat2=kernelcdf(h12,data.matrix(y2),t(as.matrix(grid))) 
Fhat3=kernelcdf(h13,data.matrix(y3),t(as.matrix(grid))) 
Fhat4=kernelcdf(h14,data.matrix(y4),t(as.matrix(grid))) 
Fhat=Ws[1]*Fhat1+Ws[2]*Fhat2+Ws[3]*Fhat3+Ws[4]*Fhat4   #KDE

######################### population distribution ################################
FN1=ecdf(subp1[,j])(grid)
FN2=ecdf(subp2[,j])(grid)
FN3=ecdf(subp3[,j])(grid)
FN4=ecdf(subp4[,j])(grid)
FN=Ws[1]*FN1+Ws[2]*FN2+Ws[3]*FN3+Ws[4]*FN4

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
bFn1=ecdf(by1)(grid)
bFn2=ecdf(by2)(grid)
bFn3=ecdf(by3)(grid)
bFn4=ecdf(by4)(grid)
bFn=Ws[1]*bFn1+Ws[2]*bFn2+Ws[3]*bFn3+Ws[4]*bFn4
 
## KDE
bFhat1=kernelcdf(h11,data.matrix(by1),t(as.matrix(grid))) 
bFhat2=kernelcdf(h12,data.matrix(by2),t(as.matrix(grid))) 
bFhat3=kernelcdf(h13,data.matrix(by3),t(as.matrix(grid))) 
bFhat4=kernelcdf(h14,data.matrix(by4),t(as.matrix(grid))) 
bFhat=Ws[1]*bFhat1+Ws[2]*bFhat2+Ws[3]*bFhat3+Ws[4]*bFhat4

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
p1;p2  ##results listed in Table 3

timeend=Sys.time()
runtime=timeend-timestart;runtime

