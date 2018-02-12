library(dplyr)
library(portalr)
library(RCurl)
library(portalr)


# ====================
# rodents
rodents = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"),
                   na.strings=c(""), colClasses=c('tag'='character'), stringsAsFactors = FALSE)
Julys = c(1,13,24,36,47,60,70,80,91,101,113,124,136,149,161,173,185,197,210,222,233,245,256,266,279,290,302,314,325,338,350,363,375,386,395,407,418,427)

# restrict rodent data to July censuses, control plots, target species, and years in which ant data was collected
rod = filter(rodents, period %in% Julys, 
             plot %in% c(2,11,14,22), 
             species %in% c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'),
             year %in% c(1977:1986,1988:1994,1998:2009))
rod$x = rep(1)

rod_byplot = aggregate(rod$x,by=list(year = rod$year, plot = rod$plot, species = rod$species),FUN=sum)
rod_byplot$index = paste(rod_byplot$year,rod_byplot$plot,sep='-')

rod_dat = select(rod_byplot,index,species,x)
rod_table = reshape(rod_dat,idvar='index',timevar='species',direction='wide')
rod_table[is.na(rod_table)] = 0

# remove species that have only one capture ever -- so extremely rare species don't have too much influence on results
rod_table = rod_table[,!names(rod_table) %in% c('x.PL','x.SF')]

# add row of zeros for plot 22 1985

# put rows in order
rod_table = rod_table[order(rod_table$index),]

write.csv(rod_table,'Rodent_julys.csv',row.names=F)

# ===========================
# rodents; average sp comp per plot, averaged over 4 summer months
rodents = abundance('..',level='Plot',time='date',shape='flat',incomplete=T)
rodent_control = filter(rodents,plot %in% c(2,11,14,22))
rodent_control$month = format(rodent_control$censusdate,'%m')
rodent_control$year = format(rodent_control$censusdate,'%Y')
rodent_control$summer = rep(NA)
rodent_control$summer[rodent_control$month %in% c('06','07','08','09')]=1

rodent_summer = filter(rodent_control,summer==1)
rodent_summer_avg = aggregate(rodent_summer$abundance,by=list(plot=rodent_summer$plot,
                                                              species=rodent_summer$species,
                                                              year=rodent_summer$year),FUN=mean,na.rm=T)

rodent_summer_table = make_crosstab(rodent_summer_avg,variable_name='x')
rodent_summer_table$index = rep(NA)
for (n in 1:length(rodent_summer_table$index)) {
  rodent_summer_table$index[n] = paste0(rodent_summer_table$year[n],'-',rodent_summer_table$plot[n])
}

# remove species that have only one capture ever -- so extremely rare species don't have too much influence on results
rodent_summer_table = rodent_summer_table[,!names(rodent_summer_table) %in% c('PH','PI','PL','RF','RO','SO')]

# put rows in order
rodent_summer_table = rodent_summer_table[order(rodent_summer_table$year,rodent_summer_table$plot),]

# put columns in order
rodent_summer_table = rodent_summer_table[,c(18,3:17)]

write.csv(rodent_summer_table,'Rodent_summer_avg.csv',row.names=F)


# ===========================
# rodents; average sp comp per plot, averaged over whole year
rodents = abundance('..',level='Plot',time='date',shape='flat',incomplete=T)
rodent_control = filter(rodents,plot %in% c(2,11,14,22))
rodent_control$year = format(rodent_control$censusdate,'%Y')
rodent_control = filter(rodent_control,year<2010)

rodent_yr_avg = aggregate(rodent_control$abundance,by=list(plot=rodent_control$plot,
                                                              species=rodent_control$species,
                                                              year=rodent_control$year),FUN=mean,na.rm=T)

rodent_yr_table = make_crosstab(rodent_yr_avg,variable_name='x')
rodent_yr_table$index = rep(NA)
for (n in 1:length(rodent_yr_table$index)) {
  rodent_yr_table$index[n] = paste0(rodent_yr_table$year[n],'-',rodent_yr_table$plot[n])
}

# remove species that have only one capture ever -- so extremely rare species don't have too much influence on results
#rodent_yr_table = rodent_yr_table[,!names(rodent_yr_table) %in% c('PH','PI','PL','RF','RO','SO')]

# put rows in order
rodent_yr_table = rodent_yr_table[order(rodent_yr_table$year,rodent_yr_table$plot),]

# put columns in order
rodent_yr_table = rodent_yr_table[,c(24,3:23)]

write.csv(rodent_yr_table,'Rodent_yearly_avg.csv',row.names=F)

# ===================
# ants: stake level presence
colony_stake = colony_presence_absence(level='Stake',rare_sp=T)

# Only 4 plots are definitely controls with respect to ants and rodents, 1977-2009
controls = filter(colony_stake,plot %in% c(2,11,14,22))
control_agg = aggregate(controls$presence,by=list(year = controls$year, plot = controls$plot, species = controls$species), FUN=sum,na.rm=T)
control_agg$index = paste(control_agg$year,control_agg$plot,sep='-')

dat = select(control_agg,index,species,x)
dat_table = reshape(dat,idvar='index',timevar='species',direction='wide')

# remove phei yell
dat_table = dat_table[,names(dat_table) != 'x.phei yell']

# put rows in order
dat_table = dat_table[order(dat_table$index),]

write.csv(dat_table,'Ant_colony_numstakes.csv',row.names=F)


# =====================
# ants: abundance; from bait data, just granivorous species


# ======================
# plants: stake level presence --WIP