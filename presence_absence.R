# attempt to get at presence/absence of species
# EMC 10/19/2017

library(dplyr)
library(RCurl)

source('ant_data_summaries.R')

colony_presence = colony_presence_absence(level='Site',rare_sp=T)

test = reshape(colony_presence,idvar='year',timevar='species',direction='wide')

# ==============================
# how does this compare to presence/absence of bait data
baitpresence = bait_presence_absence(level='Site')

A = filter(colony_presence,year==1990)
B = filter(baitpresence,year==1990)

setdiff(A,B)
setdiff(B,A)
