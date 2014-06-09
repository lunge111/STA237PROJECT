#test_______________________________________________

data <- arfima.sim(10000, model = list(phi =0.81, dfrac = -0.47, theta = 0),
                   sigma2=0.1)






RS_test(data,q=200)
RS_test(data,data_dep=T)
RS_test(data,q=0)

acf(data)$acf

SpecRegre(data,0.5)
re1=reg_sim(model="LMSV",ph=0,d=0.49,sig2=0.05)
re1
re2=reg_sim(model="LMSV",ph=0,d=0.47,sig2=0.11)
re2
setwd("D:/14sp237/project")
Bucharest  <- read.csv("INDICES.csv",head=T)[,2]
Bucharest =diff(log(Bucharest ))
Bucharest1=Bucharest^2 #squared return
Bucharest2=Bucharest #return
Bucharest3=log(abs(Bucharest)) #log absolute return 
Bucharest=log(Bucharest ^2) #log squared return

arfima.est.mc(Bucharest)

par(mfrow=c(2,2))

acf(Bucharest2,lag.max=200,ylim=c(0,0.3),main="acf of BUCH ")
acf(Bucharest1,lag.max=200,ylim=c(0,0.3),main="acf of SQBUCH")
acf(Bucharest,lag.max=200,ylim=c(0,0.3),main="acf of log SQBUCH ")
acf(Bucharest3,lag.max=200,ylim=c(0,0.3),main="acf of log ABSBUCH")

options(digits=3)
reg_sim(model="ARSV",ph=0.95,sig2=0.07)
reg_sim(model="LMSV",ph=0.93,d=0.44,sig2=0.003)
save(Bucharest,file="Bucharest.rda")
acf(Bucharest,lag.max=200,ylim=c(0,0.3),main="acf of log SQBUCH with estimated model acf")
lines(acf)

SpecRegre(Bucharest3,0.45)
SpecRegre(Bucharest3,0.5)
SpecRegre(Bucharest3,0.55)
RS_test(Bucharest3,q=0)
RS_test(Bucharest3,data_dep=T)
RS_test(Bucharest3,q=200)


#estimate____________________________________

data <- arfima.sim(4000, model = list(phi = 0.3, dfrac =0.2, theta = 0),
                   sigma2=1)+log(rchisq(4000,1))

system.time(arfima.est(data))

arfima.est.mc(Bucharest,order=c(1,0))

system.time(sig<-arfima.est.sig(Bucharest))
sig2<-arfima.est.sig(data,order=c(0,0))

result1=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(1000, model = list(phi = 0, dfrac = 0.2, theta = 0.),
                     sigma2=1)+log(rchisq(1000,1))
  result[count,]=arfima.est(data,order=c(1,0))
  if (count==100) break
}


set.seed(1111)
result1=matrix(0,100,1)
count=0
repeat{
  count=count+1
  data <- arfima.sim(1024, model = list(phi = 0.00001, dfrac = -0.4, theta = 0),
                     sigma2=1)+log(rchisq(1024,1))
  result1[count,]=arfima.est(data,order=c(0,0))$par
  if (count==100) break
}


result1.1=matrix(0,100,1)
count=0
repeat{
  count=count+1
  data <- arfima.sim(4096, model = list(phi = 0.00001, dfrac = -0.4, theta = 0),
                     sigma2=1)+log(rchisq(4096,1))
  result1.1[count,]=arfima.est(data,order=c(0,0))$par
  if (count==100) break
}
colMeans(result1.1)
sd(result1.1[,1])


result2=matrix(0,100,1)
count=0
repeat{
  count=count+1
  data <- arfima.sim(1024, model = list(phi = 0.00001, dfrac = -0.2, theta = 0),
                     sigma2=1)+log(rchisq(1024,1))
  result2[count,]=arfima.est(data,order=c(0,0))$par
  if (count==100) break
}

colMeans(result2)
sd(result2[,1])


result2.1=matrix(0,100,1)
count=0
repeat{
  count=count+1
  data <- arfima.sim(4096, model = list(phi = 0.00001, dfrac = -0.2, theta = 0),
                     sigma2=1)+log(rchisq(4096,1))
  result2.1[count,]=arfima.est(data,order=c(0,0))$par
  if (count==100) break
}

colMeans(result2.1)
sd(result2.1[,1])


result3=matrix(0,100,1)
count=0
repeat{
  count=count+1
  data <- arfima.sim(1024, model = list(phi = 0.00001, dfrac = 0, theta = 0),
                     sigma2=1)+log(rchisq(1024,1))
  result3[count,]=arfima.est(data,order=c(0,0))$par
  if (count==100) break
}

colMeans(result3)
sd(result3[,1])


result3.1=matrix(0,100,1)
count=0
repeat{
  count=count+1
  data <- arfima.sim(4096, model = list(phi = 0.00001, dfrac = 0, theta = 0),
                     sigma2=1)+log(rchisq(4096,1))
  result3.1[count,]=arfima.est(data,order=c(0,0))$par
  if (count==100) break
}

colMeans(result3.1)
sd(result3.1[,1])


result4=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(1024, model = list(phi = 0.8, dfrac = -0.2, theta = 0),
                     sigma2=1)+log(rchisq(1024,1))
  result4[count,]=arfima.est(data,order=c(1,0))$par
  if (count==100) break
}

colMeans(result4)
sd(result4[,2])


result4.1=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(4096, model = list(phi = 0.8, dfrac = -0.2, theta = 0),
                     sigma2=1)+log(rchisq(4096,1))
  result4.1[count,]=arfima.est(data,order=c(1,0))$par
  if (count==100) break
}

colMeans(result4.1)
sd(result4.1[,2])


result5=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(1024, model = list(phi = 0.00001, dfrac = 0.4, theta = 0),
                     sigma2=1)+log(rchisq(1024,1))
  result5[count,]=arfima.est(data,order=c(0,0))$par
  if (count==100) break
}

colMeans(result5)

sd(result5[,1])


result5.1=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(4096, model = list(phi = 0.00001, dfrac = 0.4, theta = 0),
                     sigma2=1)+log(rchisq(4096,1))
  result5.1[count,]=arfima.est(data,order=c(0,0))$par
  if (count==100) break
}

colMeans(result5.1)
sd(result5.1[,1])


result6=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(1024, model = list(phi = 0.4, dfrac = 0.2, theta = 0),
                     sigma2=1)+log(rchisq(1024,1))
  result6[count,]=arfima.est(data,order=c(1,0))$par
  if (count==100) break
}

colMeans(result6)
apply(result6,2,sd)


result6.1=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(4096, model = list(phi = 0.4, dfrac = 0.2, theta = 0),
                     sigma2=1)+log(rchisq(4096,1))
  result6.1[count,]=arfima.est(data,order=c(1,0))$par
  if (count==100) break
}

colMeans(result6.1)
apply(result6.1,2,sd)

result7=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(1024, model = list(phi = 0.4, dfrac = 0.4, theta = 0),
                     sigma2=1)+log(rchisq(1024,1))
  result7[count,]=arfima.est(data,order=c(1,0))$par
  if (count==100) break
}

colMeans(result7)
apply(result7,2,sd)

result7.1=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(4096, model = list(phi = 0.4, dfrac = 0.4, theta = 0),
                     sigma2=1)+log(rchisq(4096,1))
  result7.1[count,]=arfima.est(data,order=c(1,0))$par
  if (count==100) break
}

colMeans(result7.1)
apply(result7.1,2,sd)

result8=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(1024, model = list(phi = 0.8, dfrac = 0.2, theta = 0),
                     sigma2=1)+log(rchisq(1024,1))
  result8[count,]=arfima.est(data,order=c(1,0))$par
  if (count==100) break
}

colMeans(result8)
apply(result8,2,sd)

result8.1=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(4096, model = list(phi = 0.8, dfrac = 0.2, theta = 0),
                     sigma2=1)+log(rchisq(4096,1))
  result8.1[count,]=arfima.est(data,order=c(1,0))$par
  if (count==100) break
}

colMeans(result8.1)
apply(result8.1,2,sd)

result9=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(1024, model = list(phi = 0.8, dfrac = 0.4, theta = 0),
                     sigma2=1)+log(rchisq(1024,1))
  result9[count,]=arfima.est(data,order=c(1,0))$par
  if (count==100) break
}

colMeans(result9)
apply(result9,2,sd)

result9.1=matrix(0,100,2)
count=0
repeat{
  count=count+1
  data <- arfima.sim(4096, model = list(phi = 0.8, dfrac = 0.4, theta = 0),
                     sigma2=1)+log(rchisq(4096,1))
  result9.1[count,]=arfima.est(data,order=c(1,0))$par
  if (count==100) break
}

colMeans(result9.1)
apply(result9.1,2,sd)

result10=matrix(0,100,3)
count=0
repeat{
  count=count+1
  data <- arfima.sim(1024, model = list(phi = 0.8, dfrac = 0.4, theta = 0.4),
                     sigma2=1)+log(rchisq(1024,1))
  result10[count,]=arfima.est(data,order=c(1,1))$par
  if (count==100) break
}

colMeans(result10)
apply(result10,2,sd)

result10.1=matrix(0,100,3)
count=0
repeat{
  count=count+1
  data <- arfima.sim(4096, model = list(phi = 0.8, dfrac = 0.4, theta = 0.4),
                     sigma2=1)+log(rchisq(4096,1))
  result10.1[count,]=arfima.est(data,order=c(1,1))$par
  if (count==100) break
}
colMeans(result10.1)
apply(result10.1,2,sd)

par(mfrow=c(1,3))
f=factor(c(rep(1,100),rep(2,100)))
plot(f,c(result1[,1]+0.4,result1.1[,1]+0.4),main="phi=0,d=-0.4",ylab="dhat-d")
f1=factor(c(rep(1,100),rep(2,100)))
plot(f1,c(result3[,1],result3.1[,1]),main="phi=0,d=0",ylab="dhat-d")
f3=factor(c(rep(1,100),rep(2,100)))
plot(f3,c(result5[,1]-0.4,result5.1[,1]-0.4),main="phi=0,d=0.4",ylab="dhat-d")

par(mfrow=c(1,2))
f=factor(c(rep(1,100),rep(2,100)))
plot(f,c(result7[,1]-0.4,result7.1[,1]-0.4),main="phi=0.4,d=0.4",ylab="phihat-phi")
f1=factor(c(rep(1,100),rep(2,100)))
plot(f1,c(result7[,2]-0.4,result7.1[,2]-0.4),main="phi=0.4,d=0.4",ylab="dhat-d")

par(mfrow=c(1,1))
m=matrix(nrow=100,ncol=200)
i=0
repeat{
  i=i+1
  data <- arfima.sim(3629, model = list(phi =0.137, dfrac = 0.423, theta = 0),
                     sigma2=0.780)+log(rchisq(3629,1))*0.927
  m[i,]=acf(data,plot=F,lag.max=200)$acf[1:200]
  if(i==100)break
}
acf=colMeans(m)
plot(acf,type="l")
acf
acf(data,lag.max=200,ylim=c(0,0.3))


acf(data,plot=F)$acf[1:200]
acf(data)$acf