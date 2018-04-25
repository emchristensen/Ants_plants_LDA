# Using methods from Bagchi et al 2017 -- using movement models to categorize change in species composition over time (plants)
# "When depicted as a time series,
# the distance between each successive sample, relative to
# initial conditions, corresponds to community displacement
# in multivariate space, or the trajectory of change"

# "Gradual linear and abrupt nonlinear dynamics
# would be most common in the annual grasslands that
# are thought to be dominated by nonequilibrium behaviors
# (Bartolome et al. 2009); the warm season community
# was expected to be more changeable than the cool
# season community. Abrupt nonlinear changes related to
# thresholds (or tipping points) should also be represented
# in the desert grasslands that have experienced shrub
# encroachment (Ratajczak et al. 2014) and/or prolonged
# droughts (Bestelmeyer et al. 2011a)."



#source('../portalr/R/PlantData.R')
#source('../portalr/R/plant_data_processing.R')
source('create_data_tables.R')
source('movement_model_functions.R')


# do analysis with rodent data (there's basically no change if you remove rare species)
rodentsummer = rodent_dist_summers(selected_plots=c(2,4,8,11,12,14,17,22),dist_method = 'bray')


rodentsummer$nYrs <- rodentsummer$year - min(rodentsummer$year)
# gradual linear dynamics
nomad.Mod <- gradual_linear_dyn(rodentsummer)
rodentsummer$nomad_gnls <- unlist(fitted(nomad.Mod))
# reversible dynamics
migr.Mod <- reversible_dyn(rodentsummer)
rodentsummer$mig_gnls <- unlist(fitted(migr.Mod))
# stable dynamics
stab.Mod = stable_dyn(rodentsummer)
rodentsummer$stab_gnls <- unlist(fitted(stab.Mod))
# abrupt shift
disp.Mod = abrupt_dyn(rodentsummer)
rodentsummer$abrupt_gnls <- unlist(fitted(disp.Mod))

# look at data with 4 different model fits
plot_fits(rodentsummer)

# obtain aic from each model
rodent_aic = data.frame(plot.id=NA,
                      aic_nomad=AIC(nomad.Mod),
                      #aic_migr=AIC(migr.Mod),
                      aic_stab=AIC(stab.Mod),
                      aic_disp=AIC(disp.Mod))
rodent_aic
# calculate cc-score from each model
rodent_cc = data.frame(plot.id=NA,
                     cc_nomad=calc_cc_score(rodentsummer,'nomad_gnls'),
                     #cc_migr=calc_cc_score(rodentsummer,'mig_gnls'),
                     cc_stab=calc_cc_score(rodentsummer,'stab_gnls'),
                     cc_disp=calc_cc_score(rodentsummer,'abrupt_gnls'))

rodent_cc

# -----------------------------------------------------------------
# portal plant data

# control plots for everything 1977-present
#select_plots = c(11,14)
# control plots w.r.t ants and rodents 1977-2015
select_plots = c(2,11,14,22)
# plots that are rodent controls and not "annual removals"
select_plots = c(8,11,12,14,17)

#summertable = summer_annual_byplot(select_plots)
wintertable = winter_annual_byplot(select_plots,dist_method='bray')
plots = unique(wintertable$plot)

allquads = data.frame()
aictable = data.frame()
cctable = data.frame()
for(k in seq(1:length(plots))){
  quad<-wintertable[wintertable$plot == plots[k],]
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
                        aic_migr=AIC(migr.Mod),
                        aic_stab=AIC(stab.Mod)) #,
                        aic_disp=AIC(disp.Mod))
  aictable = rbind(aictable,quad_aic)
  
  # calculate cc-score from each model
  quad_cc = data.frame(plot.id=plots[k],
                       cc_nomad=calc_cc_score(quad,'nomad_gnls'),
                       cc_migr=calc_cc_score(quad,'mig_gnls'),
                       cc_stab=calc_cc_score(quad,'stab_gnls')) #,
                       cc_disp=calc_cc_score(quad,'abrupt_gnls'))
  cctable = rbind(cctable,quad_cc)
}





#write.csv(wintertable,'winterannuals.csv',row.names=F)

# Used Bray-Curtis (except for one dataset was just presence/absence so they used Sorensen)
# I think they did plot-level aggregate with Portal
# Distance from each sample relative to initial

#summerdist = vegdist(summertable[,-c(1:4)],method='bray')
#summer_ts = summerdist[1:length(summertable$year-1)]
#plot(unique(summertable$year),summer_ts,type='b')


