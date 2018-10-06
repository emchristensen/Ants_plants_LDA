# checking out the ant bait data for weirdness
# EMC 9/2017

# Bait data:
#    month      - always July
#    year       - 1988-1994; 1998-2009
#    plot       - should be all of them for all years
#    stake
#    species    - all species are found in Portal_ant_species.csv
#    abundance

library(dplyr)
library(RCurl)
library(ggplot2)
library(vegan)

source('get_data_for_LDA.R')


bait = read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Ants/Portal_ant_bait.csv"),stringsAsFactors = F)
colony = read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Ants/Portal_ant_colony.csv"),stringsAsFactors = F)
antsp = read.csv(text = getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Ants/Portal_ant_species.csv"),stringsAsFactors = F)

# control plots with regard to ants or rodents for the period 1988-2009
controls = c(2, 11, 14, 22)

# =====================================================
# just looking at data: what do we have, is there data missing etc

# plots and stakes that should be censused
plotstake = expand.grid(plot=seq(24),stake=unique(bait$stake))

# just look at what the data look like year by year -- are any stakes missing?
dat = filter(bait,year==1989)
setdiff(unique(dat$plot),seq(24))
unique(bait$abundance)
setdiff(plotstake,select(dat,plot,stake) %>% unique())

# I think the data 1988-1994 and 1998-2009 are ok. Some data sheets exist for <1988 and 1995-1997, but quality is questionable

# ========================================================
# total by species and year
byspecies = aggregate(bait$abundance,by=list(species=bait$species,year = bait$year),FUN=sum)
byspecies
ggplot(byspecies,aes(x=year,y=x,col=species)) +
  geom_point() +
  geom_line()


# ==========================================================
# make a matrix of species by year, just control plots

bait_c = filter(bait,plot %in% controls)

# this is a table of proportion of stakes where each species was found (out of 72)
bait_table = table(bait_c$year,bait_c$species)

# this table is presence/absence of species by year
bait_pres_abs <- bait_table
bait_pres_abs[bait_pres_abs>0] = 1

# table of total numbers of each species per year
bait_table = ant_bait_table(controls)

# =========================================================
# distance metrics
d = as.matrix(vegdist(bait_table),method='bray')
plot(row.names(bait_table)[-1],d[row(d) == col(d) + 1],type='l',xlab='',ylab='distance')



# ======================================================
# table of what species were even counted in each census
yearly = aggregate(colony$colonies,by=list(year=colony$year,species=colony$species),FUN=sum)
yearly_table = reshape(yearly,idvar='year',timevar='species',direction='wide')
yearly_table = yearly_table[order(yearly_table$year),]
write.csv(yearly_table,'Portal_ant_species_counted_bait.csv',row.names = F)

presenceonly=filter(colony,flag %in% c(7,2)) %>% select(year,species,flag) %>% unique()
