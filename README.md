STA237PROJECT
=============

R CODE REFERENCE AND MAN  
*****
The pro.R is the functions of our project  
run.R is codes generate figures and plots  
****
SpecRegre() is to do spectral regression test  
Arguments:  
data: the data to be tested  
u: control the upper bound of peridogram ordinates n^u  
Usage:  
install.packages("arfima")  
library(arfima)  
data <- arfima.sim(10000, model = list(phi =0.81, dfrac = -0.47, theta = 0),
                   sigma2=0.1)  
SpecRegre(data,0.5)  

*******
RS_test() is to do R/S test  
Arguments  
data: the data to be tested  
q: to control the lag q if data_dep=FALSE  
data_dep: wether use data dependent formula  
usage:  
RS_test(data,q=200)
RS_test(data,data_dep=T)
RS_test(data,q=0)
*******

reg_sim()is to run simulations of  
usage:
Arguments:  
model: wether a LMSV model or ARSV model  
ph,d,sig2: coefficients in ARFIMA(p,d,q) model  
replicate: number of replication    
samplesize: sample size of each realization  
reg_sim(model="ARSV",ph=0.9,d=0.49,sig2=0.5)  
reg_sim(model="LMSV",ph=0,d=0.47,sig2=0.5)
********

arfima.est():　estimate ARFIMA coefficients with setteled sigma's  
arfima.est.sig():  estimate ARFIMA coefficients and sigma's  
arfima.est.mc():  estimate ARFIMA coefficients and sigma's  using MCMC method  
Arguments:  
order: order of ARFIMA coefficients  
dat: data to be estimated  
FFT: wether use fast fourier transform  
startpoint:  startpoints of optim or MCMC  
sig1,sig2(in arfima.est): prespecified sigma_eta and sigma_epsilon  
max.iter(in arfima.est.mc): max iteration 




