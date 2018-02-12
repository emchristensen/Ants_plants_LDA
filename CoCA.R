#######################################################################
### Co-Correspondence Analysis: Ant, Plant, & Rodents (2 at a time) ###
########## Finds patterns between two species matrices or uses ########
########## one species matrix to attempt to predict the other #########
###################### Joan Meiners, November 2017 ####################

library(cocorresp)

## load rodent data -- abundances in each of 4 control plots in July census 1977-2009
rodents = read.csv("Rodent_julys.csv", header = TRUE, row.names = 1)
dim(rodents)
View(rodents)

## load ant data -- # stakes presence/absence in annual (July) census of 4 control plots 1977-2009
ants = read.csv("Ant_colony_numstakes.csv", header = TRUE, row.names = 1)
dim(ants)
View(ants)

## load plant data (when it's ready)


## remove row from ant data because no matching row in rodent data
ants = ants[row.names(ants) != '1985-22',]


## SYMMETRIC COCA
## First try in rodent-ant order
# log transform rodent data
rodents_log = log(rodents + 1)
## fit the model to look at rodent-ant influence
rod.ant.sym = coca(rodents_log ~ ., data = ants, method = "symmetric")
rod.ant.sym
summary(rod.ant.sym)
plot(rod.ant.sym)
title("Log(Rodents):Ants July Symmetric CoCA")

## Then compare in ant-rodent order
# log transform rodent data
ants_log = log(ants + 1)
## fit the model to look at rodent-ant influence
ant.rod.sym = coca(ants_log ~ ., data = rodents, method = "symmetric")
ant.rod.sym
summary(ant.rod.sym)
plot(ant.rod.sym)
title("Log(Ants):Rodents July Symmetric CoCA")

## plot rodent and ant influences on each other side by side (still unclear exactly how to interpret this)
quartz(width = 10, height = 6)
layout(matrix(1:2, ncol = 2))
plot(rod.ant.sym, which = "response", main = "Rodents", display = "species")
plot(rod.ant.sym, which = "predictor", main = "Ants", display = "species")
layout(1)

# same as above showing sites instead of species
quartz(width = 10, height = 6)
layout(matrix(1:2, ncol = 2))
plot(rod.ant.sym, which = "response", main = "Rodents", display = "sites", type = "text")
plot(rod.ant.sym, which = "predictor", main = "Ants", display = "sites", type = "text")
layout(1)


## Plot interpretation notes:
# the closer things are to origin, the less distinct they are (probably), so the second plot above tells us less than the first



## PREDICTIVE COCA using SIMPLS and formula interface
ra.pred = coca(rodents ~ ., data = ants)
summary(ra.pred)
ra.pred
