#Sample data: Separate file ?site1.csv? 
#This data file contains three columns: plot.id (sample identity), year (time stamp), and dist (Bray-Curtis distance)
 
#Computer code to perform model fitting and comparison. 

###################################################################################
# Bagchi et al. Quantifying long-term trajectories of plant community using animal movement #models: implications for ecological resilience
# R code for fitting the movement models. Updated 2016-07-18
# Author: Sumanta Bagchi & Navinder J Singh.
###################################################################################
setwd("PATH") # set path
list.files() 
#save the above sample data as a file called ?site1.csv? and arrange into three continuous columns
site<-read.csv("/?/site1.csv",h=T) # read the accompanying csv file
#site1.csv is the input datafile with dissimilarity (displacement) values
#it should contain plot-identity, year, and dissimilarity as columns
##################################################################################
# The steps below require the net distance be already calculated as 'dist'
# column in the data frame
##################################################################################
head(site)
site$plot.id<-as.factor(site$plot.id) # Convert plot.id to a factor
# Convert year into time steps as 1,2,3......n.
nYrs <- vector()  
for(n in 1:length(levels(site$plot.id))) {
  dataID <- subset(site, subset = c(plot.id == levels(site$plot.id)[n]),drop=TRUE) 
  nYrsA <- with(dataID,year - min(year))  
  nYrs <- c(nYrs,nYrsA)
  rm(nYrsA,dataID)
}

site$nYrs <- nYrs
###########################################################################
library(nlme)
########################################################
########################################################   

##################  
## GRADUAL LINEAR DYNAMICS  
## Separate regressions are fitted to each plot

site$plot.id<-as.factor(site$plot.id)
plots<-unique(site$plot.id)
res_nomad<-list(data=NA)
aic_nomad<-list(data=NA)
for(k in 1:length(plots)){
  quad<-site[site$plot.id == plots[k],]
  nomad.Mod <- NA
  print(k)
  nomad.Mod <- try(gnls(dist ~ C + M*(nYrs), data = quad, 
                        correlation=corAR1(), 
                        start=list(C=0,M=0)))
  print(k)
  print(class(nomad.Mod))
  if(class(nomad.Mod)[1] == "gnls" ){
    res_nomad[[k]]<-fitted(nomad.Mod)
    aic_nomad[[k]]<-AIC(nomad.Mod)
  }
  else{
    res_nomad[[k]]<- rep(NA,nrow(quad))
    aic_nomad[[k]]<-NA
    
  }
}
options(digits=4)
site$nomad_gnls <- unlist(res_nomad)
head(site) 
# fitted values are appended to the existing dataframe "site"
temp<-data.frame(plots)
colnames(temp)<-"plot.id"
temp
temp1<-data.frame(matrix(unlist(aic_nomad),nrow=k, byrow=T))
colnames(temp1)<-c("aic_linear")
temp1
aic_result<-cbind(aic_result, temp1)

###############
## REVERSIBLE DYNAMICS  ## 
## Separate regressions are fitted to each plot

site$plot.id<-as.factor(site$plot.id)
plots<-unique(site$plot.id)
res_mig<-list(data=NA)
aic_mig<-list(data=NA)
for(k in 1:length(plots)){
  quad<-site[site$plot.id == plots[k],]
  migr.Mod <- NA
  print(k)
  migr.Mod <- try(gnls(dist ~ AsymA/(1+exp((xmidA-nYrs)/scal1)) + (-AsymB /(1 + exp((xmidB-nYrs)/scal2))),
                       data = quad, 
                       correlation=corAR1(),
                       start = list(AsymA = 0.6 , AsymB = 0.6, xmidA = 10,xmidB = 20, scal1 = 1, scal2 = 1),     
                       control=nlmeControl(maxIter=200, pnlsMaxIter=200, niterEM=400, returnObject=TRUE)))
  print(k)
  print(class(migr.Mod))
  if(class(migr.Mod)[1] == "gnls" ){
    res_mig[[k]]<-fitted(migr.Mod)
    aic_mig[[k]]<-AIC(migr.Mod)
  }
  else{
    res_mig[[k]]<- rep(NA,nrow(quad))
    aic_mig[[k]]<-NA
  }
}
site$mig_gnls <- unlist(res_mig)
head(site)
# fitted values are appended to the dataframe 
temp2<-data.frame(matrix(unlist(aic_mig),nrow=k, byrow=T))
colnames(temp2)<-c("aic_migr")
aic_result<-cbind(aic_result, temp2)



## STABLE BEHAVIOUR #
## Separate regressions are fitted to each plot

res_stab<-list(data=NA)
aic_stab<-list(data=NA)
for(k in 1:length(plots)){
  quad<-site[site$plot.id == plots[k],]
  null.mod<-NA  
  asym.HRmod<-NA
  stab.Mod <- NA
  print(k)
  
  ### Step 1 Null Model
  null.mod <- try(gnls(dist ~ A, data = quad, 
                       correlation = corAR1(),
                       start = list(A = mean(quad[,'dist'])), 
                       control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=TRUE))
  ## Step 2 Asymp model
  
  asym.HRmod <- try(gnls(dist ~ Asym*(1-exp(lrc*nYrs)), data = quad, 
                         correlation = corAR1(),
                         start = list(Asym =summary(null.mod)$tTable[1],lrc=-0.059), 
                         control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=TRUE))
  ## Step 3 Full stability model
  
  stab.Mod <- try(gnls(dist ~ (Asym)*(1 - exp(lrc*(nYrs))),data = quad,
                       correlation = corAR1(),
                       start = c(Asym = summary(asym.HRmod)$tTable[1], lrc = summary(asym.HRmod)$tTable[2]),
                       control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=TRUE))
  #print(class(stab.Mod))
  
  if(class(stab.Mod)[1] == "gnls"){
    res_stab[[k]]<-fitted(stab.Mod)
    aic_stab[[k]]<-AIC(stab.Mod)
  }
  else{
    res_stab[[k]]<- rep(NA,nrow(quad))
    aic_stab[[k]]<-NA
  }
}
site$stab_gnls <- unlist(res_stab)
head(site)
# fitted values are appended to the dataframe
temp3<-data.frame(matrix(unlist(aic_stab),nrow=k, byrow=T))
colnames(temp3)<-c("aic_stab")
aic_result<-cbind(aic_result, temp3)
aic_result


## ABRUPT NONLINEAR BEHAVIOUR  
## Separate regressions are fitted to each plot

res_abrup<-list(data=NA)
aic_abrup<-list(data=NA)

for(k in 1:length(plots)){
  quad<-site[site$plot.id == plots[k],]
  null.mod<-NA  
  asym.HRmod<-NA
  disp.Mod <- NA
  print(k)
  
  ## Step 1 NULL model
  null.mod <- try(gnls(dist ~ A, data = quad, 
                       correlation = corAR1(),
                       start = c(A = mean(quad[,'dist'])), 
                       control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=F))
  ## Step 2 Asymp model
  asym.HRmod <- try(gnls(dist ~ Asym*(1-exp(lrc*nYrs)), data = quad, 
                         correlation = corAR1(),
                         start = c(Asym =summary(null.mod)$tTable[1],lrc=-0.059), 
                         control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=F))
  ## Step 3 Disp model 
  disp.Mod <- try(gnls(dist ~ (Asym)/(1 + exp((xmid-nYrs)/scal)),data = quad,
                       correlation = corAR1(), 
                       na.action = na.exclude,
                       start = c(Asym = summary(asym.HRmod)$tTable[1], xmid = 10, scal = 0.5),
                       control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=F))
  print(k)
  print(class(disp.Mod))
  if(class(disp.Mod)[1] == "gnls"){
    res_abrup[[k]]<-fitted(disp.Mod)
    aic_abrup[[k]]<-AIC(disp.Mod)
  }
  
  else{
    res_abrup[[k]]<- rep(NA,nrow(quad))
    aic_abrup[[k]]<-NA
  }
}

site$abrupt_gnls <- unlist(res_abrup)
head(site)
# fitted values are appended to the dataframe
temp4<-data.frame(matrix(unlist(aic_abrup),nrow=k, byrow=T))
colnames(temp4)<-c("aic_abrup")
aic_result<-cbind(aic_result, temp4)
aic_result
#View(aic_result)
#####################################################
### Calculation of CC-score
## create a new object called "siteCC" based on object "site" above
## site is the updated dataframe from the gnls (nlme) fitted values 
## set x = siteCC$dist # this is observed distance
## set y_i = siteCC$fitted # this is model-fitted value
## CC-score will be caclulated from x and y 
## use y1 for linear, y2 for reversible, y3 for stable, y4 for abrupt

dat<-site
head(dat)
siteCC <- array(NA, dim=c(nrow(dat),6))
colnames(siteCC)<-c("plot.id", "x", "y1", "y2", "y3", "y4")
siteCC<-data.frame(unlist(siteCC))
head(siteCC)
dim(siteCC)

names(dat)
siteCC$plot.id<-dat$plot.id
siteCC$x<-dat$dist

#### select y from the columns, one at a time
siteCC$y1<-dat$nomad_gnls 
siteCC$y2<-dat$mig_gnls
siteCC$y3<-dat$stab_gnls
siteCC$y4<-dat$abrupt_gnls

##CC score for gradual linear
siteCC$plot.id<-as.factor(siteCC$plot.id) 
plots<-unique(siteCC$plot.id) 
CC_lin<-list(data=NA)
for(k in 1:length(plots)){
  quad<-siteCC[siteCC$plot.id == plots[k],]
  num<-NA
  denom1<-NA
  denom2<-NA
  denom3<-NA
  print(k) 
  num<-sum((quad$x-quad$y1)^2)
  denom1<-sum((quad$x-mean(quad$y1))^2)
  denom2<-sum((quad$y1-mean(quad$y1))^2)
  denom3<-length(quad$x)*(mean(quad$x)-mean(quad$y1))^2
  CC_lin[[k]]<- 1-((num)/(denom1+denom2+denom3))
  print(CC_lin[[k]])
}
temp5<-data.frame(matrix(unlist(CC_lin),nrow=k, byrow=T))
colnames(temp5)<-c("CC_linear")
good_fit<-cbind(good_fit, temp5)


##CC score for reversible 
siteCC$plot.id<-as.factor(siteCC$plot.id) 
plots<-unique(siteCC$plot.id) 
CC_mig<-list(data=NA)
for(k in 1:length(plots)){
  quad<-siteCC[siteCC$plot.id == plots[k],]
  num<-NA
  denom1<-NA
  denom2<-NA
  denom3<-NA
  print(k) 
  num<-sum((quad$x-quad$y2)^2)
  denom1<-sum((quad$x-mean(quad$y2))^2)
  denom2<-sum((quad$y2-mean(quad$y2))^2)
  denom3<-length(quad$x)*(mean(quad$x)-mean(quad$y2))^2
  CC_mig[[k]]<- 1-((num)/(denom1+denom2+denom3))
  print(CC_mig[[k]])
}
temp6<-data.frame(matrix(unlist(CC_mig),nrow=k, byrow=T))
colnames(temp6)<-c("CC_migr")
good_fit<-cbind(good_fit, temp6)


##CC score for stable curve
siteCC$plot.id<-as.factor(siteCC$plot.id) 
plots<-unique(siteCC$plot.id) 
CC_stab<-list(data=NA)
for(k in 1:length(plots)){
  quad<-siteCC[siteCC$plot.id == plots[k],]
  num<-NA
  denom1<-NA
  denom2<-NA
  denom3<-NA
  print(k) 
  num<-sum((quad$x-quad$y3)^2)
  denom1<-sum((quad$x-mean(quad$y3))^2)
  denom2<-sum((quad$y3-mean(quad$y3))^2)
  denom3<-length(quad$x)*(mean(quad$x)-mean(quad$y3))^2
  CC_stab[[k]]<- 1-((num)/(denom1+denom2+denom3))
  print(CC_stab[[k]])
}
temp7<-data.frame(matrix(unlist(CC_stab),nrow=k, byrow=T))
colnames(temp7)<-c("CC_stab0")
good_fit<-cbind(good_fit, temp7)

##CC score for abrupt 
siteCC$plot.id<-as.factor(siteCC$plot.id) 
plots<-unique(siteCC$plot.id) 
CC_abrup<-list(data=NA)
for(k in 1:length(plots)){
  quad<-siteCC[siteCC$plot.id == plots[k],]
  num<-NA
  denom1<-NA
  denom2<-NA
  denom3<-NA
  print(k) 
  num<-sum((quad$x-quad$y4)^2)
  denom1<-sum((quad$x-mean(quad$y4))^2)
  denom2<-sum((quad$y4-mean(quad$y4))^2)
  denom3<-length(quad$x)*(mean(quad$x)-mean(quad$y4))^2
  CC_abrup[[k]]<- 1-((num)/(denom1+denom2+denom3))
  print(CC_abrup[[k]])
}
temp8<-data.frame(matrix(unlist(CC_abrup),nrow=k, byrow=T))
colnames(temp8)<-c("CC_abrup0")
good_fit<-cbind(good_fit, temp8)
View(good_fit) # This is the end of the analysis 
# Results in object ?good_fit? contain goodness-of-fit and parsimony information

