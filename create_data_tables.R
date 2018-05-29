library(dplyr)
library(portalr)

# ====================================================================================
# Examples: write to csvs

# control plots w.r.t. just rodents: [2,4,8,11,12,14,17,22]
# control plots w.r.t. just ants: [1,2,5,6,7,9,11,14,16,18,22,24]
# control plots w.r.t. ants and rodents: [2,11,14,22]
# control plots w.r.t. plants (and also all three taxa): [11,14]

# rod_table = rodent_table_julys()
# write.csv(rod_table,'Rodent_julys.csv',row.names=F,quote=F)

# rodent_summer_table = rodent_table_summers(selected_plots=c(2,4,8,11,12,14,17,22))
# write.csv(rodent_summer_table,'Rodent_summer_table.csv',row.names=F)

# rodent_yr_table = rodent_table_yearly()
# write.csv(rodent_yr_table,'Rodent_yearly_avg.csv',row.names=F)

# ant_table = ant_colony_presence(selected_plots=c(2,11,14,22))
# write.csv(ant_table,'Ant_colony_presence.csv',row.names=F)

# winterannuals = seasonal_plant_table(selected_plots=c(2,4,8,11,12,14,17,22),
#                                       plant_community='Annuals',
#                                       summer_winter='winter')
# write.csv(winterannuals,'WinterAnnuals.csv',row.names=F)

# summerannuals = seasonal_plant_table(selected_plots=c(2,4,8,11,12,14,17,22),
#                                       plant_community='Annuals',
#                                       summer_winter='summer')
# write.csv(summerannuals,'SummerAnnuals.csv',row.names=F)

# perennials = seasonal_plant_table(selected_plots=c(2,4,8,11,12,14,17,22),
#                                   plant_community='Perennials',
#                                   summer_winter='both')
# write.csv(perennials,'Perennials.csv',row.names=F)

# =====================================================================================
#' @title create table of rodent data: just July censuses
#'
#' @param
#'
#' @return table of rodent species counts by year
#'
rodent_table_julys = function() {
  rodents = read.csv(text=RCurl::getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"),
                     na.strings=c(""), colClasses=c('tag'='character'), stringsAsFactors = FALSE)
  Julys = c(1,13,24,36,47,60,70,80,91,101,113,124,136,149,161,173,185,197,210,222,233,245,256,266,279,290,302,314,325,338,350,363,375,386,395,407,418,427)

  # restrict rodent data to July censuses, control plots, target species, and years in which ant data was collected
  rod = filter(rodents, period %in% Julys,
               plot %in% c(2,11,14,22),
               species %in% c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'),
               year %in% c(1977:1986,1988:1994,1998:2009))
  rod$x = rep(1)

  rod_byplot = aggregate(rod$x,by=list(year = rod$year,
                                       plot = rod$plot,
                                       species = rod$species),FUN=sum)
  rod_byplot$index = paste(rod_byplot$year,rod_byplot$plot,sep='-')

  rod_dat = select(rod_byplot,index,species,x)
  rod_table = make_crosstab(rod_dat,variable_name = 'x')
  rod_table[is.na(rod_table)] = 0

  # remove species that have only one capture ever -- so extremely rare species don't have too much influence on results
  rod_table = rod_table[,!names(rod_table) %in% c('PL','SF')]

  # add row of zeros for plot 22 1985
  plt221985 = c('1985-22',rep(0,14))
  rod_table = rbind(rod_table,plt221985)

  # put rows in order
  rod_table = rod_table[order(rod_table$index),]
  return(rod_table)
}


#' @title create table of rodent count data: averaged over 4 summer months
#'
#' @param selected_plots vector of plots to include in analysis
#'
#' @return table of rodent species counts by year (1978-2015), average for summer months, control plots only
#'
rodent_table_summers = function(selected_plots) {
  rodents = abundance('..',level='Plot',time='date',shape='flat',clean=F,min_plots=23,na_drop=T)
  rodent_control = filter(rodents,plot %in% selected_plots)
  rodent_control$month = format(rodent_control$censusdate,'%m')
  rodent_control$year = format(rodent_control$censusdate,'%Y')
  rodent_control$summer = rep(NA)
  rodent_control$summer[rodent_control$month %in% c('06','07','08','09')]=1

  rodent_summer = filter(rodent_control,summer==1,year>1977,year<2015)
  rodent_summer_tot = aggregate(rodent_summer$abundance,by=list(censusdate=rodent_summer$censusdate,species=rodent_summer$species,
                                                                year=rodent_summer$year),FUN=sum)
  rodent_summer_avg = aggregate(rodent_summer_tot$x,by=list(species=rodent_summer_tot$species,
                                                            year=rodent_summer_tot$year),FUN=mean)

  rodent_summer_table = portalr::make_crosstab(rodent_summer_avg,variable_name='x')

  # remove species that have only one capture ever -- so extremely rare species don't have too much influence on results
  rodent_summer_table = rodent_summer_table[,!names(rodent_summer_table) %in% c('PH','PI','PL','RF','RO','SO')]
  

  # put rows in order
  rodent_summer_table = rodent_summer_table[order(rodent_summer_table$year),]
  #rodent_summer_table = rodent_summer_table[order(rodent_summer_table$year,rodent_summer_table$plot),]

  # put columns in order
  #rodent_summer_table = rodent_summer_table[,c(18,3:17)]
  # round to nearest integer
  rodent_summer_table[,-1] = round(rodent_summer_table[,-1])
  return(rodent_summer_table)
}

#' @title rodent distance: summer only
#'
#' @description extracts rodent data from summer months and long term control plots only,
#' creates data frame of year and community distance compared to first year (baseline)
#'
#' @param selected_plots list of plots to include (long-term rodent controls 1977-2015 are 2,11,14,22)
#' @param dist_method distance method (from vegan) to be used (default Bray-Curtis)
#' @param remove_rare_sp T/F remove species with extremely low capture rates from analysis? -----WIP -------
#'
#' @return data frame of year by distance
#'
rodent_dist_summers = function(selected_plots=c(2,11,14,22),dist_method='bray',remove_rare_sp=F) {
  rodents = abundance('..',level='Plot',time='date',shape='flat',incomplete=T)
  rodent_control = filter(rodents,plot %in% selected_plots)
  rodent_control$month = format(rodent_control$censusdate,'%m')
  rodent_control$year = format(rodent_control$censusdate,'%Y')
  rodent_control$summer = rep(NA)
  rodent_control$summer[rodent_control$month %in% c('06','07','08','09')]=1

  rodent_summer = filter(rodent_control,summer==1, year>1977, year<2015)
  rodent_summer_avg = aggregate(rodent_summer$abundance,by=list(species=rodent_summer$species,
                                                                year=rodent_summer$year),FUN=mean,na.rm=T)

  rodent_summer_table = make_crosstab(rodent_summer_avg,variable_name='x')

  # remove species with very low capture rates? -- this is a WIP!!
  if (remove_rare_sp == T) {
    rodent_summer_table = rodent_summer_table[,!names(rodent_summer_table) %in% c('BA','PH','PI','PL','RF','RO','SO')]
  }

  # put rows in order
  rodent_summer_table = rodent_summer_table[order(rodent_summer_table$year),]

  # calculate distance for each year relative to first (1978)
  rodent_distance = data.frame(yr = c(), d = c())
  init_yr = rodent_summer_table[1,-1]
  for (yr in rodent_summer_table$year) {
    d = vegan::vegdist(rbind(init_yr, rodent_summer_table[rodent_summer_table$year==yr,-1]),method=dist_method)
    rodent_distance = rbind(rodent_distance,c(as.numeric(yr),d))
  }
  names(rodent_distance) = c('year','dist')

  return(rodent_distance)
}


#' @title create table of rodent data: averaged over whole year
#'
#' @param
#'
#' @return table of rodent species counts by year
#'
rodent_table_yearly = function() {
  rodents = abundance('..',level='Plot',time='date',shape='flat',incomplete=T)
  rodent_control = filter(rodents,plot %in% c(2,11,14,22))
  rodent_control$year = format(rodent_control$censusdate,'%Y')
  rodent_control = filter(rodent_control,year %in% c(1977:1986,1988:1994,1998:2009))

  rodent_yr_avg = aggregate(rodent_control$abundance,by=list(plot=rodent_control$plot,
                                                             species=rodent_control$species,
                                                             year=rodent_control$year),FUN=mean,na.rm=T)

  rodent_yr_table = make_crosstab(rodent_yr_avg,variable_name='x')
  rodent_yr_table$index = rep(NA)
  for (n in 1:length(rodent_yr_table$index)) {
    rodent_yr_table$index[n] = paste0(rodent_yr_table$year[n],'-',rodent_yr_table$plot[n])
  }

  # remove species that are too rare
  #rodent_yr_table = rodent_yr_table[,!names(rodent_yr_table) %in% c('PH','PI','PL','RF','RO','SO')]

  # put rows in order
  rodent_yr_table = rodent_yr_table[order(rodent_yr_table$year,rodent_yr_table$plot),]

  # put columns in order
  rodent_yr_table = rodent_yr_table[,c(24,3:23)]
  return(rodent_yr_table)
}


#' @title ant species presence/absence
#'
#' @description get ant species presence/absence from colony data -- plot level
#'
#' @param selected_plots which plots to include in analysis
#'
ant_colony_presence = function(selected_plots) {
  colony_stake = colony_presence_absence(level='Stake',rare_sp=T)

  controls = filter(colony_stake,plot %in% selected_plots)
  control_agg = aggregate(controls$presence,by=list(year = controls$year,
                                                    plot = controls$plot,
                                                    species = controls$species), FUN=sum,na.rm=T)
  control_agg$index = paste(control_agg$year,control_agg$plot,sep='-')
  control_avg = aggregate(control_agg$x,by=list(year=control_agg$year,
                                                species=control_agg$species),FUN=mean)

  dat_table = reshape(control_avg,idvar='year',timevar='species',direction='wide')

  # remove phei yell
  dat_table = dat_table[,names(dat_table) != 'x.phei yell']

  # put rows in order
  dat_table = dat_table[order(dat_table$year),]
  # round to nearest integer
  dat_table[,-1] = round(dat_table[,-1])
  return(dat_table)
}


#' @title get summer annual plant data
#'
#' @param selected_plots plot numbers: 1-24
#'
#' @return table of plant counts by species, at plot level
#'
summer_annual_byplot = function(selected_plots) {
  plant_data = plant_abundance('..',level='Plot',type='Annuals',
                               correct_sp=T,unknowns=F,length='all',
                               shape='flat')
  # remove data before 1983 to avoid having to adjust by quadrat area per plot
  summer_data = filter(plant_data,season=='summer',year>1982)

  # only included species that occurred in >2% of samples (plots)
  nsamples_S = select(summer_data,year,season,plot) %>% unique() %>% nrow()
  sp_table_S = table(summer_data$species)
  splist_S = sp_table_S[sp_table_S > (nsamples_S * .02)] %>% names()
  summerplants = filter(summer_data,species %in% splist_S,plot %in% select_plots)
  summertable = make_crosstab(summerplants,variable_name='abundance')
  summertable[is.na(summertable)] <- 0

  return(summertable)
}

#' @title get winter or summer plant data: crosstab
#'
#' @param selected_plots plot numbers: 1-24
#' @param plant_community which plant community (e.g. annuals or perennials; input 'type' in plant_abundance function)
#' @param summer_winter 'summer' or 'winter' census; or "both" for perennials
#'
#' @return table of plant counts by species
#'
seasonal_plant_table = function(selected_plots,plant_community,summer_winter) {
  plant_data = plant_abundance('..',level='Plot',type=plant_community,
                               correct_sp=T,unknowns=F,length='all',
                               shape='flat')
  # remove data before 1983 to avoid having to adjust by quadrat area per plot; stop at 2015 to avoid plot switch

  season_data = filter(plant_data,year>1982,year<2015)
  if (tolower(summer_winter) %in% c('summer','winter')) {
    season_data = dplyr::filter(season_data,season==summer_winter)
  }

  # find species that occurred in >10% of years
  transients = find_transient_species(season_data,threshhold=.1)

  # filter based on selected_plots, remove transient species
  seasonplants = filter(season_data,plot %in% selected_plots, !(species %in% transients$species))
  # sum by year (effort doesn't vary by year)
  seasontotal = aggregate(seasonplants$abundance,by=list(year=seasonplants$year,season=seasonplants$season,species=seasonplants$species),FUN=sum)
  seasontable = portalr::make_crosstab(seasontotal,variable_name='x')
  seasontable[is.na(seasontable)] <- 0
  return(seasontable)
}


#' @title find transient species
#'
#' @description identifies species found in <=10% of sampling events
#'
#' @param data_table table of data including species and year
#' @param threshhold value for determining transient. Default 10%
#'
find_transient_species = function(data_table,threshhold=0.1) {
  ntimesteps = dplyr::select(data_table,year) %>% unique() %>% nrow()
  sp_time = dplyr::select(data_table,year,species) %>% unique()
  yrs_present = dplyr::count(sp_time,species)
  transient_names = dplyr::filter(yrs_present,n <= threshhold*ntimesteps) %>% select(species)
  return(transient_names)
}

#' @title get winter annual plant data: 'distance' from initial year
#'
#' @param selected_plots plot numbers: 1-24
#' @param dist_method distance method for vegdist function (vegan package)
#'
#' @return table of plant counts by species, at plot level
#'
winter_annual_distance = function(selected_plots, dist_method='bray') {
  plant_data = plant_abundance('..',level='Plot',type='Annuals',
                               correct_sp=T,unknowns=F,length='all',
                               shape='flat')
  # remove data before 1983 to avoid having to adjust by quadrat area per plot
  winter_data = filter(plant_data,season=='winter',year>1982)

  # find species that occurred in >2% of samples (plots)
  nsamples_W = dplyr::select(winter_data,year,season,plot) %>% unique() %>% nrow()
  sp_table_W = table(winter_data$species)
  splist_W = sp_table_W[sp_table_W > (nsamples_W * .02)] %>% names()

  # filter based on species and selected_plots
  winterplants = filter(winter_data,species %in% splist_W,plot %in% select_plots)
  wintertable = make_crosstab(winterplants,variable_name='abundance')
  wintertable[is.na(wintertable)] <- 0

  # calculate distance for each plot
  winter_distance = data.frame()
  for (plt in selected_plots) {
    plt_dat = filter(wintertable, plot==plt)
    init_yr = plt_dat[plt_dat$year==min(plt_dat$year),-c(1:4)]
    D=data.frame()
    for (yr in plt_dat$year) {
      d = vegan::vegdist(rbind(init_yr, plt_dat[plt_dat$year==yr,-c(1:4)]),method='bray')
      D = rbind(D,c(plt,yr,d))
    }
    names(D) = c('plot','year','dist')
    winter_distance = rbind(winter_distance,D)
  }
    return(winter_distance)
}

#' @title get perennial plant data
#'
#' @param selected_plots plot numbers: 1-24
#' @param dist_method distance method for vegdist function (vegan package)
#'
#' @return table of plant counts by species, at plot level
#'
perennial_byplot = function(selected_plots, dist_method='bray') {
  plant_data = plant_abundance('..',level='Plot',type='Non-woody',
                               correct_sp=T,unknowns=F,length='all',
                               shape='flat')
  # remove data before 1983 to avoid having to adjust by quadrat area per plot
  winter_data = filter(plant_data,season=='winter',year>1982)

  # find species that occurred in >2% of samples (plots)
  nsamples_W = dplyr::select(winter_data,year,season,plot) %>% unique() %>% nrow()
  sp_table_W = table(winter_data$species)
  splist_W = sp_table_W[sp_table_W > (nsamples_W * .02)] %>% names()

  # filter based on species and selected_plots
  winterplants = filter(winter_data,species %in% splist_W,plot %in% select_plots)
  wintertable = make_crosstab(winterplants,variable_name='abundance')
  wintertable[is.na(wintertable)] <- 0

  # calculate distance for each plot
  winter_distance = data.frame()
  for (plt in selected_plots) {
    plt_dat = filter(wintertable, plot==plt)
    init_yr = plt_dat[plt_dat$year==min(plt_dat$year),-c(1:4)]
    D=data.frame()
    for (yr in plt_dat$year) {
      d = vegan::vegdist(rbind(init_yr, plt_dat[plt_dat$year==yr,-c(1:4)]),method='bray')
      D = rbind(D,c(plt,yr,d))
    }
    names(D) = c('plot','year','dist')
    winter_distance = rbind(winter_distance,D)
  }
  return(winter_distance)
}



# =====================
# ants: abundance; from bait data, just granivorous species


# ====================================================================================================================

