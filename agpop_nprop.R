######## Coverage frequencies of SCBs for agpop CDF under proportional allocation
######## Table 6, Figures 8,9

library(stats)
library(MASS)

kintegral=function(u) ## distrution of kernal function
{
  Ktilde = (3/16*u^5-5/8*u^3+15/16*u+1/2)*(abs(u)<=1)+(u>1)
}

kernelcdf = function(h,z,zfixed)   ## KDE
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

bpopn=function(y,N)   ## generate bootstrap population from sample
{
n=length(y)
re=floor(N/n)
by1=kronecker(rep(1,re),y)
by2=sample(y,N-nrow(by1),replace=FALSE,prob=NULL)
popn=c(by1,by2)
}


### reading data
data=read.table("C:/mydocuments/papers/stratifiedSCB/agpop.txt", header=TRUE,sep = ",")

head(data)

data_NE<-data[data$REGION=='NE',]
data_NC<-data[data$REGION=='NC',]
data_S<-data[data$REGION=='S',]
data_W<-data[data$REGION=='W',]

y_NE=data_NE[,3]
y_NC=data_NC[,3]
y_S=data_S[,3]
y_W=data_W[,3]


##################################### delete the missing data #############################
yNE=as.matrix(y_NE[y_NE>-99])
yNC=as.matrix(y_NC[y_NC>-99])
yS=as.matrix(y_S[y_S>-99])
yW=as.matrix(y_W[y_W>-99])
y=as.matrix(data[,3][data[,3]>-99])
ylarge=y[y>3.5*10^6];ylarge

N_NE=nrow(yNE);N_NE
N_NC=nrow(yNC);N_NC
N_S=nrow(yS);N_S
N_W=nrow(yW);N_W
Ns=c(N_NE,N_NC,N_S,N_W);Ns   ## each stratum size
N=sum(Ns);N                  ## total population size
Nsorigin=c(nrow(data_NE),nrow(data_NC),nrow(data_S),nrow(data_W));Nsorigin
TN=nrow(data);TN
nmiss=TN-N;nmiss


################ parameters of Gamma distribution ###################
## shape=alpha, scale=beta
beta1=var(yNE)*(1-1/Ns[1])/mean(yNE);beta1
alpha1=mean(yNE)/beta1;alpha1

beta2=var(yNC)*(1-1/Ns[2])/mean(yNC);beta2
alpha2=mean(yNC)/beta2;alpha2

beta3=var(yS)*(1-1/Ns[3])/mean(yS);beta3
alpha3=mean(yS)/beta3;alpha3

beta4=var(yW)*(1-1/Ns[4])/mean(yW);beta4
alpha4=mean(yW)/beta4;alpha4

alpha=c(1.41,1.45,0.67,0.76)
beta=c(0.66,2.25,2.98,9.56)

################################ allocation #########################
timestart=Sys.time()
n=150  ## total sample size n=150,300,600,900
r=400  ## (r+1) is the number of grid points
m=1000  #replications

#### stratum sample sizes under proportional allocation #######
ns=round(n*Ns/N);ns;sum(ns)
indmax=which(ns== max(ns), arr.ind = TRUE)
ns[indmax]=n-sum(ns[-indmax]);ns;sum(ns) ##propotional

Ws=as.matrix(Ns/N)      ## the stratum weights
lamda=(1/n-1/N)^(-1/2)  ## corrected scale factor

#####################generate sample########################

set.seed(123456)
y1=sample(yNE,ns[1],replace=FALSE,prob=NULL)
y2=sample(yNC,ns[2],replace=FALSE,prob=NULL)
y3=sample(yS,ns[3],replace=FALSE,prob=NULL)
y4=sample(yW,ns[4],replace=FALSE,prob=NULL)

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

#### sample KDE with bandwidth h1 or h2 ########
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
FN1=ecdf(yNE)(grid)
FN2=ecdf(yNC)(grid)
FN3=ecdf(yS)(grid)
FN4=ecdf(yW)(grid)
Ws=as.matrix(Ns/N)
FN=Ws[1]*FN1+Ws[2]*FN2+Ws[3]*FN3+Ws[4]*FN4


##################### critical value ########################
##### bootstrap method
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
by1=sample(bpopn1,ns[1],replace=FALSE,prob=NULL)
by2=sample(bpopn2,ns[2],replace=FALSE,prob=NULL)
by3=sample(bpopn3,ns[3],replace=FALSE,prob=NULL)
by4=sample(bpopn4,ns[4],replace=FALSE,prob=NULL)

bFn1=ecdf(by1)(grid)
bFn2=ecdf(by2)(grid)
bFn3=ecdf(by3)(grid)
bFn4=ecdf(by4)(grid)
bFn=Ws[1]*bFn1+Ws[2]*bFn2+Ws[3]*bFn3+Ws[4]*bFn4

bFhat1=kernelcdf(h11,data.matrix(by1),t(as.matrix(grid))) 
bFhat2=kernelcdf(h12,data.matrix(by2),t(as.matrix(grid))) 
bFhat3=kernelcdf(h13,data.matrix(by3),t(as.matrix(grid))) 
bFhat4=kernelcdf(h14,data.matrix(by4),t(as.matrix(grid))) 
bFhat=Ws[1]*bFhat1+Ws[2]*bFhat2+Ws[3]*bFhat3+Ws[4]*bFhat4

Bs1=lamda*(abs(bFn-Fn))
Bs2=lamda*(abs(bFhat-Fhat))
Lstar1[i]=max(Bs1)
Lstar2[i]=max(Bs2)
}  
Lstar1=sort(Lstar1)
Lstar2=sort(Lstar2)
level=c(0.01,0.05,0.1,0.2)
qvalue1=Lstar1[c(0.99*m,0.95*m,0.9*m,0.8*m)];qvalue1  ##nonsmooth critical value
qvalue2=Lstar2[c(0.99*m,0.95*m,0.9*m,0.8*m)];qvalue2  ##smooth critical value

#### limiting method ####################
#invlamdas=(1/ns-1/Ns)^(1/2)
#weight=(Ws*invlamdas)*lamda
weight=sqrt(Ws)
B1=rbridge(Fn1,m)
B2=rbridge(Fn2,m)
B3=rbridge(Fn3,m)
B4=rbridge(Fn4,m)
Bstar=weight[1]*B1+weight[2]*B2+weight[3]*B3+weight[4]*B4
Lstar=apply(abs(Bstar),2,max)
Lstar=sort(Lstar)
qvalue=Lstar[c(0.99*m,0.95*m,0.9*m,0.8*m)];qvalue

SRSqvalue=c(1.63,1.36,1.22,1.07) ##quantiles of Kolmogorov distribution

###################### SCB coverage ###########################
#set.seed(123456)
D = matrix(1,m,4)
v1=rep(0,4)
v2=rep(0,4)
v3=rep(0,4)
v4=rep(0,4)
v5=rep(0,4)
v6=rep(0,4)

p1=rep(0,4)
p2=rep(0,4)
p3=rep(0,4)
p4=rep(0,4)
p5=rep(0,4)
p6=rep(0,4)

for (j in 1:m)
{
y1=sample(yNE,ns[1],replace=FALSE,prob=NULL)
y2=sample(yNC,ns[2],replace=FALSE,prob=NULL)
y3=sample(yS,ns[3],replace=FALSE,prob=NULL)
y4=sample(yW,ns[4],replace=FALSE,prob=NULL)
strsample=c(y1,y2,y3,y4)

Fn1=ecdf(y1)(grid)
Fn2=ecdf(y2)(grid)
Fn3=ecdf(y3)(grid)
Fn4=ecdf(y4)(grid)
Fn=Ws[1]*Fn1+Ws[2]*Fn2+Ws[3]*Fn3+Ws[4]*Fn4

h11= (quantile(y1)[4]-quantile(y1)[2])*((1/ns[1]-1/Ns[1])^(-1/2))^(-2)
h12= (quantile(y2)[4]-quantile(y2)[2])*((1/ns[2]-1/Ns[2])^(-1/2))^(-2)
h13= (quantile(y3)[4]-quantile(y3)[2])*((1/ns[3]-1/Ns[3])^(-1/2))^(-2)
h14= (quantile(y4)[4]-quantile(y4)[2])*((1/ns[4]-1/Ns[4])^(-1/2))^(-2)
Fhat1=kernelcdf(h11,data.matrix(y1),t(as.matrix(grid))) 
Fhat2=kernelcdf(h12,data.matrix(y2),t(as.matrix(grid))) 
Fhat3=kernelcdf(h13,data.matrix(y3),t(as.matrix(grid))) 
Fhat4=kernelcdf(h14,data.matrix(y4),t(as.matrix(grid))) 
Fhat=Ws[1]*Fhat1+Ws[2]*Fhat2+Ws[3]*Fhat3+Ws[4]*Fhat4

################## treat the stratified sample as SRS ####################
SRSyn=sample(y,n,replace=FALSE,prob=NULL)  ##
SRSh1= (quantile(strsample)[4]-quantile(strsample)[2])*((1/n-1/N)^(-1/2))^(-2)
#SRSh2= (quantile(strsample)[4]-quantile(strsample)[2])*((1/n-1/N)^(-1/2))^(-2/3)
SRSFhat=kernelcdf(SRSh1,data.matrix(strsample),t(as.matrix(grid)))
SRSFn=ecdf(strsample)(grid)

D[j,1] = max(abs(Fn-FN))
D[j,2] = max(abs(Fhat-FN))
D[j,3] = max(abs(SRSFn-FN))
D[j,4] = max(abs(SRSFhat-FN))

 for (k in 1:4)
{
if (D[j,1]<=qvalue1[k]/lamda)
 v1[k]=v1[k]+1
if (D[j,2]<=qvalue2[k]/lamda)
 v2[k]=v2[k]+1
if (D[j,1]<=qvalue[k]/lamda)
 v3[k]=v3[k]+1
if (D[j,2]<=qvalue[k]/lamda)
 v4[k]=v4[k]+1
if (D[j,3]<=SRSqvalue[k]/lamda)
 v5[k]=v5[k]+1
if (D[j,4]<=SRSqvalue[k]/lamda)
 v6[k]=v6[k]+1
 }
}
p1=v1/m;p2=v2/m;p3=v3/m;p4=v4/m;p5=v5/m;p6=v6/m
p1;p2;p3;p4;p5;p6  ## results listed in Table 6

timeend=Sys.time()
runtime=timeend-timestart;runtime

######################## 95% SCB ########################
lamda=(1/n-1/N)^(-1/2)
up1 = matrix(1,(r+1),1)
lw1 = matrix(1,(r+1),1)
for (i in 1:(r+1))
{
up1[i] = min(1, Fn[i]+qvalue[2]/lamda)  ## upper bound of SCB based on Fn
lw1[i] = max(0, Fn[i]-qvalue[2]/lamda)  ## lower bound of SCB based on Fn
}

up2 = matrix(1,(r+1),1)
lw2 = matrix(1,(r+1),1)
for (i in 1:(r+1))
{
up2[i] = min(1, Fhat[i]+qvalue[2]/lamda)  ## upper bound of SCB based on KDE
lw2[i] = max(0, Fhat[i]-qvalue[2]/lamda)  ## lower bound of SCB based on KDE
}

###################### plot SCB plot SCB shown in Figure 9 #############################
plot(grid,FN,type="l",col=1,lty=1,lwd=3,,xlab="",ylab="",xlim=c(leftpoint,1*10^6),main="agpop SCB")
lines(grid,Fhat,col=3,lty=1,lwd=2)
lines(grid,up2,col=5,lty=5,lwd=2)
lines(grid,lw2,col=5,lty=5,lwd=2)
lines(grid,Fn,col=2,lty=2,lwd=2)
lines(grid,up1,col=4,lty=3,lwd=2)
lines(grid,lw1,col=4,lty=3,lwd=2)

legend("bottomright",c("population CDF","stratified KDE","smooth SCB",
"stratified EDF","nonsmooth SCB"),
lty=c(1,1,5,2,3),lwd=c(3,2,2,2,2),col=c(1,3,5,2,4),title="")

################### plot population and stratum CDFs shown in Figure 8(a) ##################
plot(grid,FN1,type="l",lty=5,lwd=2,col=1,xlab="",ylab="Empirical CDF",xlim=c(leftpoint,b))
lines(grid,FN2,col=2,lty=4,lwd=2)
lines(grid,FN3,col=3,lty=3,lwd=2)
lines(grid,FN4,col=4,lty=2,lwd=2)
lines(grid,FN,col=5,lty=1,lwd=3)
legend("bottomright",c("Northeast","North Central","South","West","population CDF"),
lty=c(5,4,3,2,1),lwd=c(2,2,2,2,3),col=c(1,2,3,4,5),title="")


#### plot density of each stratum shown in Figure 8(b) #####
fN1=density(yNE)
fN2=density(yNC)
fN3=density(yS)
fN4=density(yW)
plot(fN1,xlim=c(leftpoint,b),col=1,main="",xlab="",ylab="smoothed relative frequency",lwd=2)
lines(fN2,col=2,lty=2,lwd=2)
lines(fN3,col=3,lty=3,lwd=2)
lines(fN4,col=4,lty=4,lwd=2)
legend("topright",c("Northeast","North Central","South","West"),
lty=c(1,2,3,4),lwd=c(2,2,2,2),col=c(1,2,3,4),title="")






