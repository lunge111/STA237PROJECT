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
#install.packages("arfima")  
#library(arfima)  
data <- arfima.sim(10000, model = list(phi =0.81, dfrac = -0.47, theta = 0),
                   sigma2=0.1)  
SpecRegre(data,0.5)  

