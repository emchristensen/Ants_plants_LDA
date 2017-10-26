library(dplyr)
library(portalr)



# ====================
# rodents
rodents = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"),
                   na.strings=c(""), colClasses=c('tag'='character'), stringsAsFactors = FALSE)
Julys = c(1,13,24,36,47,60,70,80,91,101,113,124,136,149,161,173,185,197,210,222,233,245,256,266,279,290,302,314,325,338,350,363,375,386,395,407,418,427)

rod = filter(rodents, period %in% Julys, plot %in% c(2,11,14,22), species %in% c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'))
rod$x = rep(1)

rod_byplot = aggregate(rod$x,by=list(year = rod$year, plot = rod$plot, species = rod$species),FUN=sum)
rod_byplot$index = paste(rod_byplot$year,rod_byplot$plot,sep='-')

rod_dat = select(rod_byplot,index,species,x)
rod_table = reshape(rod_dat,idvar='index',timevar='species',direction='wide')
rod_table[is.na(rod_table)] = 0
write.csv(rod_table,'Rodent_julys.csv',row.names=F)

# ===================
# ants: stake level presence
colony_stake = colony_presence_absence(level='Stake',rare_sp=T)

# Only 4 plots are definitely controls with respect to ants and rodents, 1977-2009
controls = filter(colony_stake,plot %in% c(2,11,14,22))
control_agg = aggregate(controls$presence,by=list(year = controls$year, plot = controls$plot, species = controls$species), FUN=sum,na.rm=T)
control_agg$index = paste(control_agg$year,control_agg$plot,sep='-')

dat = select(control_agg,index,species,x)
dat_table = reshape(dat,idvar='index',timevar='species',direction='wide')
write.csv(dat_table,'Ant_colony_numstakes.csv',row.names=F)
