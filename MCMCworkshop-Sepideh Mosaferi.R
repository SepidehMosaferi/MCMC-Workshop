###########################################################
###     Introducing Monte Carlo Methods with R          ###
###    Author: Sepideh Mosaferi                         ###
###########################################################

library(spuRs)
library(mcsm)
library(mnormt)

## generating random variable

Nsim <- 10^4
x <- runif(Nsim)
x1 <- x[-Nsim]
x2 <- x[-1]
par(mfrow=c(1,3))
hist(x)
plot(x1,x2)
acf(x)

set.seed(1)
runif(5)
set.seed(1)
runif(5)
set.seed(2)
runif(5)

## Model based simulation
empiricalMSE=rep(0,50)
for(i in 1:50)
{
  x=rnorm(1000,5,3) #var=9
  sai=rgamma(1000,4.5,scale=2)
  v=rnorm(1000,0,2) #var=4
  e=rnorm(1000,0,sqrt(sai))
  Y=1+3*x+v
  fit=lm(Y~x)
  X=cbind(rep(1,1000),x)
  beta_hat=c(coef(fit)[1],coef(fit)[2])
  Y_hat=X%*%beta_hat
  empiricalMSE[i]=sum((Y_hat-Y)^2)/1000
}
empiricalMSE
mean(empiricalMSE)


## Design based simulation
one.rep <- function(){
  ## Sample selection
  SMP.IDs <- strata(data=Employment,stratanames="STATE",
                    size=nh,method="srswor")    
  SAMPLE <- getdata(Employment,SMP.IDs)  #Output of selected sample
  SAMPLEMEAN <- tapply(SAMPLE$TOTEMP12rep,INDEX=SAMPLE$STATE,FUN=mean)
  SAMPLESD <- tapply(SAMPLE$TOTEMP12rep,INDEX=SAMPLE$STATE,FUN=sd)
  SAMPLEVAR <- (1-(nh/Nh))*(SAMPLESD^2/nh)
}
## Replicate 1000 sample:
R <- 1000
many.reps <- replicate(n=R,one.rep()); dim(many.reps)

## Model-Design based simulation
require(PracTools); require(sampling)
require(pscl); require(MCMCpack)
# Population:
x <- rgamma(2e4,shape=2,scale=5)
b <- 1.25*(x^(3/2))*((8+5*x)^(-1))
c <- 0.04*(x^(-3/2))*(8+5*x)^2
y <- rgamma(2e4,shape=c,scale=b) #Our realistic population
range(y);summary(y); summary(x)
cor(x,y)

FRAME <- data.frame(cbind(y,x))
head(FRAME)
H <- 10 #total number of strata
N <- nrow(FRAME) #population size
# Create strata
STRATUM <- cut(x,breaks=c(-Inf,quantile(x,seq(0.1,0.9,by=0.1)),Inf))
UPDATEFRAME <- data.frame(FRAME,STRATUM)
UPDATEFRAME <- UPDATEFRAME[order(STRATUM),]; head(UPDATEFRAME)
one.rep <- function(){
  #Sample selection (selecting 1 PSU)
  SMP.IDs <- strata(data=UPDATEFRAME,stratanames="STRATUM",size=ng,method="srswor")
  SAMPLE <- getdata(UPDATEFRAME,SMP.IDs)  # Output of selected sample
}
## Replicate 1000 sample:
R <- 1000
many.reps <- replicate(n=R,one.rep())



## probability integral transformation
u <- runif(Nsim,0,1)
x <- -log(1-u)
system.time(test1())
y <- rexp(Nsim)
system.time(test2())
par(mfrow=c(1,2))
hist(x,freq=F,main="Exp from Uniform")
hist(y,freq=F,main="Exp from R")

U <- runif(3*10^4)
U <- matrix(U,nrow=3)
X <- -log(U)
X <- 2*apply(X,2,sum)
system.time(test1())
X2 <- rchisq(10^4,df=6)
system.time(test2())
par(mfrow=c(1,2))
hist(X); hist(X2)

## Box-Muller Alg.
size = 100000
U1 <- runif(size)
U2 <- runif(size)
X1 <- rep(0,size) 
X2 <- rep(0,size)
for (i in 1:size){
  X1[i] <- sqrt(-2*log(U1[i]))*cos(2*pi*U2[i])
  X2[i] <- sqrt(-2*log(U1[i]))*sin(2*pi*U2[i])
}

#Y1 <- Y2 <- rnorm(Nsim)
par(mfrow=c(2,1))
plot(density(X1),main="Normal Density from X1"); plot(density(X2),main="Normal Density from X2")

## Poisson R code
Nsim <- 10^4
lambda <- 100
spread <- 3*sqrt(lambda)
t <- round(seq(max(0,lambda-spread),lambda+spread,1))
prob <- ppois(t,lambda)
X <- rep(0,Nsim)
for(i in 1:Nsim){
  u <- runif(1)
  X[i] <- t[1]+sum(prob<u)
}

## Student's t density (df=10): Continuous
df <- 10
y <- rchisq(Nsim,df,ncp=0)
x <- rnorm(Nsim,0,df/y)
par(mfrow=c(2,2))
plot(density(x),main="Density from Mixture",xlim=c(-6,6))
plot(density(rt(Nsim,df,ncp=0)),main="Density from Student's t",xlim=c(-6,6))

## Negative Binomial: Discrete
n <- 6
p <- 0.3
y <- rgamma(Nsim,n,rate=p/(1-p))
x <- rpois(Nsim,y)
hist(x,freq=FALSE,breaks = 40,main="Negative Binomial")
lines(1:50,dnbinom(1:50,n,p),col="red")

## Accept-Reject Algorithm (Betas from Uniforms)
M <- 2.669744 #Maximum of U[0,M]
Nsim <- 2500
a <- 2.7;b <- 6.3
par(mfrow=c(1,2))
y <- runif(Nsim) #generated from g (U[0,1])
u <- runif(Nsim,max=M) #generated from U[0,M]
x <- y[u<dbeta(y,a,b)] #accepted subsample
plot(y,u,pch=15,col="grey",main="Y~U[0,1] and U~U[0,M]")
polygon(density(x,bw=0.04), col="skyblue", border="red",lwd=4)
#the acceptance area as shown in the skyblue is 1/M=1/2.669=%36
rb=rbeta(Nsim,2.7,6.3)
d<-density(rb)
plot(y,u,pch=15,col="grey",main="from R beta")
polygon(d,col="skyblue", border="red",lwd=4)

## Toy function of cosine and sine
par(mfrow=c(1,2))
h <- function(x){(cos(50*x)+sin(20*x))^2}
curve(h,xlab="Function",ylab="",lwd=2)
integrate(h,0,1)
x <- h(runif(Nsim))
estint <- cumsum(x)/(1:Nsim)
esterr <- sqrt(cumsum((x-estint)^2))/(1:Nsim)
plot(estint,xlab="Mean and Error Range",type="l",lwd=2,
     ylim=mean(x)+20*c(-esterr[Nsim],esterr[Nsim]),ylab="")
lines(estint+2*esterr,col="gold",lwd=2)
lines(estint-2*esterr,col="gold",lwd=2)

## Marginal distribution related to the Beta distribution
#we consider the exp(log(f))=f.
#By looking at the image of f, we can figure out a better function for 
#simulation rather than the true posterior which is difficult for simulation.
#lambda=1;x0=y0=0.5;x=0.6
f=function(a,b){
  exp(2*(lgamma(a+b)-lgamma(a)-lgamma(b))+
        a*log(.3)+b*log(.2))}
aa=1:150 #alpha grid for image
bb=1:100 #beta grid for image
post=outer(aa,bb,f) #outer provides an array which conducts function f on aa & bb
image(aa,bb,post,xlab=expression(alpha),ylab=expression(beta),main="Posterior on Beta") #image for 3-dimensional/spatial
#data image "It basically makes a frame for outr plot"
contour(aa,bb,post,add=T)  #create a counter plot

#From the plot,we can consider t-student distribution.
x=matrix(rt(2*10^4,3),ncol=2)#T sample (df=3)
E=matrix(c(220,190,190,180),ncol=2)
image(aa,bb,post,xlab=expression(alpha),ylab=expression(beta),main="Student's t") 
##interested simulated variable:
y=t(t(chol(E))%*%t(x)+c(50,45)) #Note:x~t(3);Now, we want the distribution of
#y=(x-mu)/sigma.  x=sigam*y+mu f(y)=fx(y*sigma+mu)
points(y,cex=.6,pch=19)

#Finding out the marginal likelihood of (m(x))while g(a,b)=t-student
x=matrix(rt(2*10^4,3),ncol=2)
E=matrix(c(220,190,190,180),ncol=2)
y=t(t(chol(E))%*%t(x)+c(50,45)) 
#Hint:As a and b must essentially be positive but t distribution brings negative 
#values for them, we have to eliminate them by using y=y[ine>0]
ine=apply(y,1,min);ine
y=y[ine>0,] ;y
x=x[ine>0,];x #This just trims x for the sake of having the same number of
#variables like y and eliminated rows corresponding omitted rows'y
normx=sqrt(x[,1]^2+x[,2]^2)#here by this code we make positive x 
f=function(y){
  exp(2*(lgamma(y[,1]+y[,2])-lgamma(y[,1])-lgamma(y[,2]))+
        y[,1]*log(.3)+y[,2]*log(.2))}
h=function(y){
  exp(lgamma(y[,1]+y[,2])-lgamma(y[,1])-lgamma(y[,2])+y[,1]*log(.5)+y[,2]*(log(.5)))}
den=dt(normx,3) #density of t-distribution as g
estimator=mean(f(y)/den)/mean(h(y)/den);estimator #answer of integral

#finding the posterior expectations of the parameters 
f=function(y){
  exp(2*(lgamma(y[,1]+y[,2])-lgamma(y[,1])-lgamma(y[,2]))+
        y[,1]*log(.3)+y[,2]*log(.2))}
mean(y[,1]*f(y)/den)/mean(h(y)/den) #mean of a
mean(y[,2]*f(y)/den)/mean(h(y)/den) #mean of b

## optimization for cos+sin
#Using optimize function:
h=function(x){
  (cos(50*x)+sin(20*x))^2}
optimize(h,lower=0,upper=1,maximum=TRUE)
#According to the above code, xmax=0.379 & h(xmax)=3.8325.

#Using Uniform sampler (proposed by myself):
set.seed(1)
Nsim=10^4
x=runif(Nsim)
h=function(x){
  (cos(50*x)+sin(20*x))^2}
curve(h)
xnew=sort(x)
optimize(h,xnew,maximum=TRUE)

#Using Uniform sampler (proposed method in book):
rangom=h(matrix(runif(10^6),ncol=10^3))
monitor=t(apply(rangom,1,cummax))
plot(monitor[1,],type="l",col="white",ylab=expression(paste(beta)))
polygon(c(1:10^3,10^3:1),c(apply(monitor,2,max),
                           rev(apply(monitor,2,min))),col="grey")
abline(h=optimise(h,int=c(0,1),maximum=TRUE)$ob,col="sienna",lwd=2)

## Artificical Mini
h=function(x,y){(x*sin(20*y)+y*sin(20*x))^2*cosh(sin(10*x)
                                                 *x)+(x*cos(10*y)-y*sin(10*x))^2*cosh(cos(20*y)*y)}
x=y=seq(-3,3,le=435) #defines a grid for persp
z=outer(x,y,h)
par(bg="wheat",mar=c(1,1,1,1)) #bg stands for background
persp(x,y,z,theta=155,phi=30,col="green4",
      ltheta=-120,shade=.75,border=NA,box=FALSE)
min(z);max(z)


## Metropolis-Hastings Beta distribution
a=2.7;b=6.3;c=2.669 #initial values
Nsim=5000
X=rep(runif(1),Nsim)
for (i in 2:Nsim){
  Y=runif(1)
  rho=dbeta(Y,a,b)/dbeta(X[i-1],a,b)
  X[i]=X[i-1]+(Y-X[i-1])*(runif(1)<rho)
}
print(X)
par(mfrow=c(1,2))
hist(X,col="black",freq=F,main="Metropolis-Hastings",xlim=c(0,0.8))
curve(dbeta(x,a,b),col="red",lwd=2,add=TRUE)
hist(rbeta(5000,a,b),nclass=20,col="gray",main="Direct Generation",freq=F,xlab="X")
curve(dbeta(x,a,b),col="red",lwd=2,add=TRUE)
#Checking by Kolmogrove-Smirnov test:
ks.test(jitter(X),rbeta(5000,a,b))
#Checking by moments:
a/(a+b);a*b/(((a+b)^2)*(a+b+1)) #For Beta distribution
xbar=mean(X);xvar=sd(X)^2;xbar;xvar #For Metropolis-Hasting

#Checking for different iteration:
Nsim=4500
X1=rep(runif(1),Nsim)
for (i in 2:Nsim){
  Y=runif(1)
  rho=dbeta(Y,a,b)/dbeta(X1[i-1],a,b)
  X1[i]=X1[i-1]+(Y-X1[i-1])*(runif(1)<rho)
}
print(X1)
Nsim=4550
X2=rep(runif(1),Nsim)
for (i in 2:Nsim){
  Y=runif(1)
  rho=dbeta(Y,a,b)/dbeta(X2[i-1],a,b)
  X2[i]=X2[i-1]+(Y-X2[i-1])*(runif(1)<rho)
}
print(X2)
Nsim=4600
X3=rep(runif(1),Nsim)
for (i in 2:Nsim){
  Y=runif(1)
  rho=dbeta(Y,a,b)/dbeta(X3[i-1],a,b)
  X3[i]=X3[i-1]+(Y-X3[i-1])*(runif(1)<rho)
}
print(X3)
Nsim=4650
X4=rep(runif(1),Nsim)
for (i in 2:Nsim){
  Y=runif(1)
  rho=dbeta(Y,a,b)/dbeta(X4[i-1],a,b)
  X4[i]=X4[i-1]+(Y-X4[i-1])*(runif(1)<rho)
}
print(X4)
Nsim=4700
X5=rep(runif(1),Nsim)
for (i in 2:Nsim){
  Y=runif(1)
  rho=dbeta(Y,a,b)/dbeta(X5[i-1],a,b)
  X5[i]=X5[i-1]+(Y-X5[i-1])*(runif(1)<rho)
}
print(X5)
Nsim=4750
X6=rep(runif(1),Nsim)
for (i in 2:Nsim){
  Y=runif(1)
  rho=dbeta(Y,a,b)/dbeta(X6[i-1],a,b)
  X6[i]=X6[i-1]+(Y-X6[i-1])*(runif(1)<rho)
}
print(X6)
Nsim=4800
X7=rep(runif(1),Nsim)
for (i in 2:Nsim){
  Y=runif(1)
  rho=dbeta(Y,a,b)/dbeta(X7[i-1],a,b)
  X7[i]=X7[i-1]+(Y-X7[i-1])*(runif(1)<rho)
}
print(X7)
X=c(X1,X2,X3,X4,X5,X6,X7)
#Iteration=c(4500,4550,4600,4650,4700,4750,4800)
par(mfrow=c(1,1))
plot(X,xlim=c(4500,4800),type="l",xlab="Iterations",ylim=c(0,0.7))


#Comparison of Accept-Reject and Metropolis-Hastings algorithms 
4/4.85=0.8247423

Nsim=5000
par(mfrow=c(1,2))
hist(rgamma(Nsim,4.85,1))
hist(rgamma(Nsim,4,0.8247423))

f=function(x){dgamma(x,4.85,1)}
g=function(x){dgamma(x,4,0.8247423)}
par(mfrow=c(1,2))
plot(f, from=0, to=15, col="blue", ylab="")
plot(g, from=0, to=15, col="red",add=T)
#Trial and Error method to get a good M:
g=function(x) {1*dgamma(x,4,0.8247423)} #M=1 suits well.
f=function(x) {dgamma(x,4.85,1)}
plot(f, from=0, to=15, ylim=c(0,1),col="blue", ylab="")
plot(g, from=0, to=15, col="red",add=T)
#Accept-Reject algorithm:
Nsim=5000
x=rgamma(Nsim,4,0.8247423)
u=runif(Nsim,0,1)
ratio=dgamma(x,4.85,1)/dgamma(x,4,0.8247423)
ind=I(u<ratio)
sim<-x[ind==1]
length(sim)
par(mfrow=c(1,2))
hist(sim,freq=F,nclass=20,xlab="X",ylab="Density",main="Accept-Reject Method")
lines(density(sim,adjust=2),col="red")
#Metropolis-Hastings Algorithm
Nsim=5000
X=rep(rgamma(1,4,0.8247423),Nsim) #initiate value from target
for (i in 2:Nsim){
  Y=rgamma(1,4,0.8247423) #Candidate
  rho=dgamma(Y,4.85,1)*dgamma(X[i-1],4,0.8247423)/(dgamma(X[i-1],4.85,1)*dgamma(Y,4,0.8247423))
  X[i]=X[i-1]+(Y-X[i-1])*(runif(1)<rho)
  print(rho)
}
print(X);length(X)
hist(X,freq=F,nclass=20,xlab="X",ylab="Density",main="Metropollis-Hastings Method")
lines(density(X,adjust=2),col="red")

mean(sim);sd(sim)^2
mean(X);sd(X)^2
par(mfrow=c(1,2))
#To assess correlation in both samples:
acf(sim,ylim=c(-0.05,0.25),main="Accept-Reject Method",col="red")
acf(X,ylim=c(-0.05,0.25),main="Metropollis-Hastings Method",col="blue")

## Gibbs Samplers
library(MASS) #Required packages for dataset Energy
library(coda)
library(lattice)
library(mcsm)
Energy
x=c(91,504,557,609,693,727,764,803,857,929,970,1043,1089,1195,1384,1713,
    457,645,736,790,899,991,1104,1154,1203,1320,1417,1569,1685,1843,2296,2710)
x1=c(91,504,557,609,693,727,764,803,857,929,970,1043,1089,1195,1384,1713) #For Girls
x2=c(457,645,736,790,899,991,1104,1154,1203,1320,1417,1569,1685,1843,2296,2710) #For Boys
xbar1=mean(x1);xbar2=mean(x2) #1 for girls and 2 for boys
Nsim=5000
a1=a2=a3=b1=b2=b3=3
n1=n2=16;n=32;k=2;mu0=5
sh1=(n/2)+a1;sh2=(k/2)+a2;sh3=(1/2)+a3
theta1=theta2=mu=sigma2=tau2=sigma2mu=rep(0.5,Nsim) #init arrays
for(i in 2:Nsim){
  B1=sigma2[i-1]/(sigma2[i-1]+(n1*tau2[i-1]))
  B2=sigma2[i-1]/(sigma2[i-1]+(n2*tau2[i-1]))
  B=tau2[i-1]/(tau2[i-1]+(k*sigma2mu[i-1]))
  thetabar=((n1*theta1[i-1])+(n2*theta2[i-1]))/n
  theta1[i]=rnorm(1,m=B1*mu[i-1]+(1-B1)*xbar1,sd=sqrt(B1*tau2[i-1]))
  theta2[i]=rnorm(1,m=B2*mu[i-1]+(1-B2)*xbar2,sd=sqrt(B2*tau2[i-1]))
  mu[i]=rnorm(1,m=B*mu0+(1-B)*thetabar,sd=sqrt(B*sigma2mu[i-1]))
  ra1=0.5*(sum((x1-theta1[i])^2)+sum((x2-theta2[i])^2))+b1
  ra2=0.5*(((theta1[i]-mu[i])^2)+((theta2[i]-mu[i])^2))+b2
  ra3=0.5*((mu[i]-mu0)^2)+b3
  sigma2[i]=1/rgamma(1,shape=sh1,rate=ra1)
  tau2[i]=1/rgamma(1,shape=sh2,rate=ra2)
  sigma2mu[i]=1/rgamma(1,shape=sh3,rate=ra3)
  par(mfrow=c(2,3))
  hist(mu,xlab=expression(mu),nclass=30,col='gray')
  hist(theta1,xlab=expression(theta[1]),nclass=30,col='gray')
  hist(theta2,xlab=expression(theta[2]),nclass=30,col='gray')
  hist(sqrt(sigma2mu),xlab=expression(sigma[mu]),nclass=30,col='sienna')
  hist(sqrt(tau2),xlab=expression(tau),nclass=30,col='sienna')
  hist(sqrt(sigma2),xlab=expression(sigma),nclass=30,col='sienna')
}

