# 1. Make resist map using species distribution model

#install package
install.packages("rgdal")                                            
install.packages("ggplot2")
install.packages("dplyr")
install.packages("dismo")

library(dismo)
library(dplyr)
library(ggplot2)
library(rgdal)

#upload data

species_data = read.csv("C:/Users/k/Desktop/frog_n.csv", head=T, sep=",") # occurance species data
head(species_data)


#envrionment data 
species_tm9 <-  stack("C:/Users/k/Desktop/wc2.1_2.5m_tmax_2018-09.tif") 
species_pr9 <-  stack("C:/Users/k/Desktop/wc2.1_2.5m_prec_2018-09.tif") 

#merge RasterLayer object
species_curr=stack(species_tm9,species_pr9)
plot(species_curr$wc2.1_2.5m_tmax_09)

# change object names (wc2.1_2.5m_tmax_09 -->> tmax // wc2.1_2.5m_prec_09 -->> precip)
names(species_curr) <- c('tmax','precip')


#projection 
points(species_data[c("POINT_X","POINT_Y")], pch = "+", cex = 0.5)
projection(species_Curr)



#########################################################################################
# set the coordinate reference system (CRS) (define the projection)#(참고)
projection(x) <- "+proj=utm +zone=48 +datum=WGS84"
#########################################################################################



#data preprocessing
species_locations = select(species_data,POINT_X, POINT_Y) 
head(species_locations) 

##Extract mulit value to points
env_curr = extract(species_curr, species_locations) 
species_env = cbind(species_data, env_curr) 
head(species_env)

## modelling

ggplot(species_env,
       mapping = aes(x=tmax, y=precip, color=present))+
  geom_point() 

logistic_regr_model = glm(present ~ tmax + precip,
                          family = binomial(link = "logit"),
                          data = species_env)

summary(logistic_regr_model) 

## ROC curve ##
presence_data = filter(species_env, present ==1)
absence_data = filter(species_env, present ==0)

evaluation = evaluate(presence_data, absence_data, logistic_regr_model)
plot(evaluation, 'ROC')


## prediction model 
predictions = predict(species_curr,
                      logistic_regr_model,
                      type = "response")

plot(predictions, ext=extent(120,133,33,41))
points(presence_data[c("POINT_X","POINT_Y")], pch = "+", cex = 0.5)

plot(predictions > 0.5, ext = extent(120,133,33,41))

tr = threshold(evaluation, stat = 'prevalence')
plot(predictions > tr, ext = extent(120,133,33,41))
points(presence_data[c("POINT_X","POINT_Y")], pch = "+", cex = 0.5)


# 2. metapopulation model using resist mapping through sdm results


install.packages('rgdal')
install.packages('rgeos')
install.packages('sp')
install.packages('raster')
install.packages('Matrix')
install.packages('gdistance')
install.packages('leastcostpath')
install.packages('quadtree')
library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(Matrix)
library(gdistance)
library(leastcostpath)
library(quadtree)




## uploading resist map 
resist = raster("C:/Users/k/Desktop/R/re_2013.tif")

re = transition(resist,mean,4)
par(mfrow = c(1,2))
plot(resist)
plot(raster(re))

## habitat coordinate
sp = read.csv("C:/Users/k/Desktop/R/gorani_coor.csv", head=T, sep=",") 
sp = data.matrix(sp)
head(sp)
sp1 = sp
sp2 = sp


## Least Cost path 
Cost_D = costDistance(re,sp1,sp2)
Cost_Dmm = nor_minmax(Cost_D)
Cost_Dsd = nor_sd(Cost_D)
head(Cost_D)


## 2. data
install.packages("readxl")
library(readxl)

setwd("C:/Users/k/Desktop/R")
fritty=read_excel("gorani_2012.xlsx") 
head(fritty)
summary(fritty)

# Figure 1 Basic plot of occupied and empty meta-population patches.
attach(fritty)
plot(x.crd, y.crd, asp=1, xlab="x.crd", ylab="y.crd", cex=sqrt(A*5), pch=21,col=p+1, bg=5*p) 

## 3. fitting
## 3.2 Fitting a snapshot in R
#dij
d=dist(cbind(x.crd,y.crd)) 
alpha=1
edis=as.matrix(exp(-alpha*d))
diag(edis)=0
edis=sweep(edis,2,A,"*")
S=rowSums(edis[,p>0])


mod=glm(p~offset(2*log(S))+log(A), family = binomial)
summary(mod)
AIC(mod)

beta=coef(mod)
(xhat=beta[2])

(AO=min(A[p > 0]))

(ey = exp(-beta[1]))
(etilde=AO^xhat)

(ytilde=ey/etilde)

## 3.3 Fitting data from two surveys
P=p+p2
S=rowSums(sweep(edis, 2, P/2, "*"))
mod2=glm(cbind(P, 2-P) ~ offset(2*log(S))+
           log(A), family = binomial)
summary(mod2)
AIC(mod2)
# Figure 2. Fitted incidences.
col=heat.colors(100)[99 * (1-fitted(mod))+1]
plot(x.crd, y.crd, asp=1, xlab="", ylab="", pch=21, col="blue", bg=col, cex=sqrt(A*5))

beta=coef(mod2)
xhat=beta[2]
ey=exp(-beta[1])
etilde=min(A[P>0])^xhat
ytilde=ey/etilde
par=c(xhat, etilde, ytilde)
names(par)=c("x", "e", "y")
par


## 3.4 Separating e and y with two surveys

f=function(x) sum(1/(S*S+x) * (S*S*(1- pmean) + ey * pmean/A^xhat)) - Tpar
(Tpar = sum(p !=p2))

beta=coef(mod2)
xhat=beta[2]
ey=exp(-beta[1])
pmean=P/2
(sol=uniroot(f, c(0, 50)))

etilde=ey/sol$root
par2=c(xhat, etilde, sol$root)
names(par2)=c("x", "e", "y")
rbind(par, par2)

##3.5 Covariates
Class=factor(sample(c("A","B"), length(p), replace = TRUE))
modc=glm(p~ offset(2*log(S)) + log(A) + Class-
           + 1, family=binomial)
coef(modc)

modc=glm(p ~ offset(2* log(S)) + Class/log(A)- 1, family = binomial)
coef(modc)

vec=runif(length(p))
modv = glm(p~offset(2* log(S))+log(A) + vec, family = binomial)

(b= coef(modv))

ey=exp(-(b[1] +b[3] * vec))

modv = glm(p~offset(2*log(S)) + log(A) * vec, family = binomial)

(b = coef(modv))

xhat=b[2] + b[4] * vec
ey=exp(-(b[1]+ b[3] * vec))

modx=glm(p~offset(2*log(S))+log(A) * vec-vec, family = binomial)
coef(modx)

modx=glm(p~offset(2*log(S))+ Class/log(A) - Class, family = binomial)
(b=coef(modx))
AO=min(A[p>0])
(etilde=AO^b[2:3])

(ytilde=exp(-b[1])/etilde)

anova(modc, test = "Chisq")


### 3.6 Estimating a ####

alphascan=function(alpha, d, A, p) {
  edis = as.matrix(exp(-alpha * d))
  diag(edis) = 0
  edis = sweep(edis, 2, A, "*")
  S= rowSums(edis[, p>0])
  mod=glm(p~offset(2*log(S))+log(A), family = binomial)
  deviance(mod)
}
(sol=optimize(alphascan, c(0.15, 5), d=d, p=p, A=A))


##3.7 Confidence intervals of estimates
tmp=summary(mod2)$coefficients
tmp
tmp[1, 1] + c(-2, 2) * tmp[1, 2]
tmp[2, 1] + c(-2, 2) * tmp[2, 2]

# Figure 3. Profile deviance of a. 
nseq=21
alpha=seq(0.1, 5, length=nseq)
prof = numeric(nseq)
for(i in 1:nseq){
  prof[i] = alphascan(alpha[i], d=d, A=A, p=p)} 
plot(alpha, prof, ylab = "Deviance", type = "l", col = "blue", lwd =3)
abline(v=sol$minimum)
abline(v=1, lty=2)
abline(h=sol$objective + qchisq(0.95, 1))

tmp[, 1] + t(outer(c(-2,2), tmp[,2]))

install.packages('Mass')
library(MASS)
confint(mod2)

## 4. Simulation
metastep = function(p, edis, E, y){
  p=p>0
  if (any(p)) {
    S = rowSums(edis[, p, drop = FALSE])
    C = S^2/(S^2 +y)
    cond = ifelse(p, (1-C) *E, C)
    p=ifelse(runif(length(p)) < cond, !p, p)
  }
  as.numeric(p)
}

tmp = p
par

C=S^2/S^2+par[3]

E= pmin(par[2]/A^par[1],1)
tmp=metastep(tmp, edis, E, par[3])
occup=matrix(0, nrow = length(p), ncol = 100+1)
occup[,1] = p
for(t in 1:100){
  occup[, t+1]=metastep(occup[,t], edis, E, par[3])
}
# Figure 4 Simulated population size
plot(colSums(occup), type = "l", col = "blue", lwd =2, xlab = "Time", ylab = "Population Size")
abline(h= mean(colSums(occup)), col="red", lty=2)

ps = colSums(occup) #population size on year

write.csv(ps, 'population_size.gorani_LCPSD.csv')
write.csv(edis,'euclide_d_LCPSD.csv')
write.csv(occup, 'gorani_imf_LCPSD.csv') # 데이터 저장하기
write.csv(E,'extinction_gorani_LCPSD.csv')
write.csv(C,'Colonization_gorani_LCPSD.csv')
write.csv(Cost_D,"lcp_distanceSD.csv")

# Figure 5 Simulated incidences against fitted incidence.
plot(rowMeans(occup), fitted(mod), pch =21, col = "red",
     bg = "yellow", xlab = "Simulated incidence", ylab="Fitted incidence")
abline(0, 1, col ="blue")


## 5. Metacommunity capacity

alpha=1
M=outer(A, A) * as.matrix(exp(-alpha * d))
tmp=eigen(M)

lambda.M = tmp$values[1]
lambda.vec = tmp$vector[, 1]^2


N=length(A)
take=sample(N)
tmp=M[take, take]
cap=numeric(N)
for(i in 1:N){
  cap[i] = eigen(tmp[i:N, i:N])$value[1]
}

#Figure 6 metapopulation capacity as a function of number of randomly selected patches
plot(N:1, cap, xlab = "Number of Patches", ylab = "Metapopulation Capacity",
     type = "b", col = "red", pch = 21, bg ="yellow")


round(lambda.vec, 3)

take=rev(order(lambda.vec))
eigen(M[take[1:5], take[1:5]])$value[1]
eigen(M[take[6:N], take[6:N]])$value[1]

#Figure 7 Simulation results separately for the sites with highest metapopulation capacity("5 Best") and other 45 sites("Rest")
best = matrix(0, nrow=5, ncol=101)
rest = matrix(0, N-5, ncol = 101)
best[, 1] = occup[take[1:5], 101]
rest[, 1] = occup[take[6:N], 101]
i = take[6:N]
for (t in 1:100){
  best[, t+1] = metastep(best[, t], edis[i, i], E[i], par[3])
}
i = take[6:N]
for (t in 1:100){
  rest[, t+1]= metastep(rest[,t], edis[i, i], E[i], par[3])
}
bestline = c(colSums(occup[1:5,]), colSums(best[, -1]))
restline = c(colSums(occup[6:N,]), colSums(rest[,- 1]))
matplot(1:201, cbind(bestline, restline), xlab="Time",
        ylab = "Occupied patches", type = "l", lwd =2, lty =1)
abline(v=101)
legend(150, 0.8 * max(restline), c("5 Best", "Rest"), lty=1, col=1:2, lwd =2)


#######################################################################################################################################################

# E N D #

#######################################################################################################################################################


