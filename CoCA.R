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

ants = ants[row.names(ants) != '1985-22',]

# plants = 

## SYMMETRIC COCA
# log transform rodent data
rodents_log = log(rodents + 1)  #might consider swapping ant and rodent data positions here, not sure of the influence of either
ants_log = log(ants + 1)

## fit the model
ra.sym = coca(ants_log ~ ., data = rodents, method = "symmetric")
ra.sym
summary(ra.sym)
plot(ra.sym)


## PREDICTIVE COCA using SIMPLS and formula interface
ra.pred = coca(rodents ~ ., data = ants)
summary(ra.pred)
ra.pred
