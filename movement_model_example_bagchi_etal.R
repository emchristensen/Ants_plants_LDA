# Adapted from:
###################################################################################
# Bagchi et al. Quantifying long-term trajectories of plant community using animal movement #models: implications for ecological resilience
# R code for fitting the movement models. Updated 2016-07-18
# Author: Sumanta Bagchi & Navinder J Singh.
###################################################################################

#Sample data: Separate file site1.csv 
#This data file contains three columns: plot.id (sample identity), year (time stamp), and dist (Bray-Curtis distance)
 
#Computer code to perform model fitting and comparison. 

library(nlme)
source('movement_model_functions.R')


# read in sample data:
site<-read.csv("site1.csv",h=T) # read the accompanying csv file
site$plot.id<-as.factor(site$plot.id) # Convert plot.id to a factor

# list of unique "plots" for analysis
plots<-unique(site$plot.id)


# loop through 6 different plots contained in site file

allquads = data.frame()
aictable = data.frame()
for(k in 1:length(plots)){
  quad<-site[site$plot.id == plots[k],]
  quad$nYrs <- quad$year - min(quad$year)
  # gradual linear dynamics
  nomad.Mod <- gradual_linear_dyn(quad)
  quad$nomad_gnls <- unlist(fitted(nomad.Mod))
  # reversible dynamics
  migr.Mod <- reversible_dyn(quad)
  quad$mig_gnls <- unlist(fitted(migr.Mod))
  # stable dynamics
  stab.Mod = stable_dyn(quad)
  quad$stab_gnls <- unlist(fitted(stab.Mod))
  # abrupt shift
  disp.Mod = abrupt_dyn(quad)
  quad$abrupt_gnls <- unlist(fitted(disp.Mod))
  
  # look at data with 4 different model fits
  plot_fits(quad)
  
  # put together results into one big dataframe
  allquads = rbind(allquads,quad)
  
  # obtain aic from each model
  quad_aic = data.frame(plot.id=plots[k],
                        aic_nomad=AIC(nomad.Mod),
                        #aic_migr=AIC(migr.Mod),
                        aic_stab=AIC(stab.Mod),
                        aic_disp=AIC(disp.Mod))
  aictable = rbind(aictable,quad_aic)
}


aictable


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
#good_fit<-cbind(good_fit, temp5)
good_fit<-temp5


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

