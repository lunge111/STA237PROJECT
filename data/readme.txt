The data is processed as follow:
setwd("D:/14sp237/project")
Bucharest  <- read.csv("INDICES.csv",head=T)[,2]
Bucharest =diff(log(Bucharest ))
Bucharest1=Bucharest^2 #squared return
Bucharest2=Bucharest #return
Bucharest3=log(abs(Bucharest)) #log absolute return 
Bucharest=log(Bucharest ^2) #log squared return