library(arfima)
SpecRegre<-function(data,u){
  n=length(data)
  lowb=floor(n^0.1)
  upb=floor(n^u)
  spec=spec.pgram(data,taper=0,log="no",plot=F)
  y=log(spec$spec)[lowb:upb]  
  x=log(Mod(1-exp(-complex(real=0,imaginary=spec$freq[lowb:upb])))^2)
  m=summary(lm(y~x))
  m[[4]][1,2]
  return(c(d=-m[[4]][2,1],pval=m[[4]][2,4]))
}




RS_test<-function(data,q,data_dep=FALSE){
  require(itsmr)
  if(data_dep==TRUE){
    rho1=acf(data,plot=F)$acf[2]
    q=floor((1.5*length(data))^(1/3)*(2*rho1/(1-rho1^2))^(2/3))
  }
  r=max(cumsum(data)-(1:length(data))*mean(data))-min(cumsum(data)
                                                      -(1:length(data))*mean(data))
  if(q!=0){
    i=c(-rev(1:q),0,1:q)
    gamma=acvf(data,h=q)[abs(i)+1]
    w=-abs(i)/q+1
    s=sqrt(sum(w*gamma))
    Q=r/s
    J=log(Q)/log(length(data))
    V=Q/sqrt(length(data))
    if(V<0.809|V>1.862) {resul=1
    }else resul=0
    
    return(c(J=J,result=resul))
  }else {
    s=sqrt(acvf(data)[1]) 
    Q=r/s
    J=log(Q)/log(length(data))
    V=Q/sqrt(length(data))
    if(V<0.809|V>1.862) {resul=1
    }else resul=0
    
    return(c(J=J,result=resul))
  }
}

reg_sim<-function(model=c("ARSV","LMSV"),ph,d,sig2,
                  replicate=100,samplesize=10000){
  m=matrix(ncol=12,nrow=replicate)
  mod=match.arg(model)
  if(ph==0) ph=1e-4
  for(i in 1:replicate){
    if(mod=="ARSV"){
      dat <- arfima.sim(samplesize, model = list(phi = ph, dfrac = 0, theta = 0),
                        sigma2=sig2)+log(rchisq(samplesize,1))
    }else{
      dat <- arfima.sim(samplesize,  model = list(phi = ph, dfrac = d, theta = 0),
                        sigma2=sig2)+log(rchisq(samplesize,1))
    }
    m[i,]=c(SpecRegre(dat,u=0.45),SpecRegre(dat,u=0.5),
            SpecRegre(dat,u=0.55),RS_test(dat,q=0),
            RS_test(dat,q=0,data_dep=T),RS_test(dat,q=200))
  }
  
  dhat=colMeans(m[,c(1,3,5,7,9,11)])
  sderr=apply(m[,c(1,3,5,7,9,11)],2,sd)
  rej=c(colSums(m[,c(2,4,6)]<0.05)/replicate,colSums(m[,c(8,10,12)])/replicate)
  df=matrix(c(dhat,sderr,rej),nrow=3,byrow=T)
  df=as.data.frame(df)
  names(df)<-c("d,u=0.45","d,u=0.5","d,u=0.55","J_q=0","J_q=qstar","J_q=200")
  rownames(df)<-c("mean","sd","rej_proportion")  
  return(df)
  
}


arfima.est<-function(dat,order=c(1,0),FFT=TRUE,startpoint=NULL,sig1=1,sig2=pi^2/2){
  if(!is.null(startpoint)&length(startpoint)!=sum(order)+1)
    stop("wrong number of start points")
  ar=order[1]
  ma=order[2]
  spec_f<-function(phi,theta,d,x,sig1,sig2){
    if(theta==0) {num=1
    }else num=(1+sum(theta*cos((1:length(theta))*x)))^2+(sum(theta*sin((1:length(theta))*x)))^2
    if(phi==0){
      dem=1
    }else dem=(1-sum(phi*cos((1:length(phi))*x)))^2+sum(phi*sin((1:length(phi))*x))^2
    
    sig1*num/((2-2*cos(x))^d*dem*2*pi)+sig2/2/pi
  }

  
  I<-function(dat,freq){
    1/(2*pi*length(dat))*(sum(dat*cos(freq*(1:length(dat))))^2
                          +sum(dat*sin(freq*(1:length(dat))))^2)
  }
  
  fr<-function(x)
  {
    if(order[1]==0){
      phi=0}else phi=x[1:ar]
    
    if(order[2]==0){
      theta=0}else  theta<-x[(ar+1):(ar+ma)]
    d=x[length(x)]
   
    n=floor(length(dat)/2)
    freq=(1:n)*2*pi/length(dat)
    
    if(FFT){
      In=abs(fft(dat)[2:(n+1)])^2/length(dat)/2/pi 
    }else{
      In=sapply(freq,I,dat=dat)
    }
    f=sapply(freq,spec_f,phi=phi,theta=theta,d=d,sig1=sig1,sig2=sig2)
    sum(log(f)+In/f)
  }  
  if(is.null(startpoint)){
    startpoint=c(rep(0.7,ar+ma),SpecRegre(dat,0.55)[1])
  }
  par=optim(startpoint, fr,hessian=T)$par
  list(par=par,mle=fr(par))
  
}


arfima.est.sig<-function(dat,order=c(1,0),FFT=TRUE,startpoint=NULL){
  if(!is.null(startpoint)&length(startpoint)!=sum(order)+1)
    stop("wrong number of start points")
  ar=order[1]
  ma=order[2]
  spec_f<-function(phi,theta,d,x,sig1,sig2){
    if(theta==0) {num=1
    }else num=(1+sum(theta*cos((1:length(theta))*x)))^2+(sum(theta*sin((1:length(theta))*x)))^2
    if(phi==0){
      dem=1
    }else dem=(1-sum(phi*cos((1:length(phi))*x)))^2+sum(phi*sin((1:length(phi))*x))^2
    
    sig1*num/((2-2*cos(x))^d*dem*2*pi)+sig2/2/pi
  }
  
  
  I<-function(dat,freq){
    1/(2*pi*length(dat))*(sum(dat*cos(freq*(1:length(dat))))^2
                          +sum(dat*sin(freq*(1:length(dat))))^2)
  }
  
  fr<-function(x)
  {
    if(order[1]==0){
      phi=0}else phi=x[1:ar]
    
    if(order[2]==0){
      theta=0}else  theta<-x[(ar+1):(ar+ma)]
    d=x[length(x)]
    sig1=x[sum(order)+1]
    sig2=x[sum(order)+2]
    
    n=floor(length(dat)/2)
    freq=(1:n)*2*pi/length(dat)
    
    if(FFT){
      In=abs(fft(dat)[2:(n+1)])^2/length(dat)/2/pi 
    }else{
      In=sapply(freq,I,dat=dat)
    }
    f=sapply(freq,spec_f,phi=phi,theta=theta,d=d,sig1=sig1,sig2=sig2)
    sum(log(f)+In/f)
  }  
  if(is.null(startpoint)){
    startpoint=c(rep(0.7,ar+ma+2),SpecRegre(dat,0.55)[1])
  }
  par=optim(startpoint, fr,hessian=T)$par
  list(par=par,mle=fr(par))
  
}



arfima.est.mc<-function(dat,order=c(1,0),max.iter=2000,start=c(0.5,5)){
  sig1=start[1]
  sig2=start[2]
  sig1.new=sig1+abs(rnorm(1,sd=0.5))
  sig2.new=sig2+abs(rnorm(1,sd=0.5))
  est=arfima.est(dat,sig1=sig1.new,sig2=sig2.new)
mle=est$mle
par=est$par
mle.old=mle-100
i=0
  while(abs(mle.old-mle)>0.001&i<max.iter){
    
sig1.new=sig1+0.1*abs(runif(1)-0.5)
sig2.new=sig2+abs(runif(1)-0.5)
if(mle>mle.old){
  sig1=sig1.new
  sig2=sig2.new
  mle.old=mle
  par=est$par
}
est=arfima.est(dat,startpoint=par,sig1=sig1.new,sig2=sig2.new)
mle=est$mle
i=i+1
  }
c(sig1sq=sig1,sig2sq=sig2,par)

}







