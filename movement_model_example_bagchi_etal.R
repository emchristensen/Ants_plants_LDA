# Adapted from:
###################################################################################
# Bagchi et al. Quantifying long-term trajectories of plant community using animal movement #models: implications for ecological resilience
# R code for fitting the movement models. Updated 2016-07-18
# Author: Sumanta Bagchi & Navinder J Singh.
###################################################################################

#Sample data: Separate file site1.csv 
#This data file contains three columns: plot.id (sample identity), year (time stamp), and dist (Bray-Curtis distance)
 
#Computer code to perform model fitting and comparison. 

source('movement_model_functions.R')


# read in sample data:
site<-read.csv("site1.csv",h=T) # read the accompanying csv file
site$plot.id<-as.factor(site$plot.id) # Convert plot.id to a factor

# list of unique "plots" for analysis
plots<-unique(site$plot.id)


# loop through 6 different plots contained in site file

allquads = data.frame()
aictable = data.frame()
cctable = data.frame()
for(k in seq(1:length(plots))){
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
  
  # calculate cc-score from each model
  quad_cc = data.frame(plot.id=plots[k],
                       cc_nomad=calc_cc_score(quad,'nomad_gnls'),
                       #cc_migr=calc_cc_score(quad,'mig_gnls'),
                       cc_stab=calc_cc_score(quad,'stab_gnls'),
                       cc_disp=calc_cc_score(quad,'abrupt_gnls'))
  cctable = rbind(cctable,quad_cc)
}

quad_aic
quad_cc
aictable
cctable

