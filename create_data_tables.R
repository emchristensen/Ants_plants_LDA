library(dplyr)
library(portalr)


# ===================
# stake level presence
colony_stake = colony_presence_absence(level='Stake',rare_sp=T)

# Only 4 plots are definitely controls with respect to ants and rodents, 1977-2009
controls = filter(colony_stake,plot %in% c(2,11,14,22))
control_agg = aggregate(controls$presence,by=list(year = controls$year, plot = controls$plot, species = controls$species), FUN=sum,na.rm=T)
control_agg$index = paste(control_agg$year,control_agg$plot,sep='-')

dat = select(control_agg,index,species,x)
dat_table = reshape(dat,idvar='index',timevar='species',direction='wide')
write.csv(dat_table,'Ant_colony_numstakes.csv',row.names=F)
