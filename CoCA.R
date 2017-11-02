############################################################################################################
################# Co-Correspondence Analysis of Ant, Plant, and Rodent data (two at a time) ################
## Finds patterns between two species matrices or uses one species matrix to attempt to predict the other ##
#################################### Joan Meiners, November 2017 ###########################################

library(cocorresp)

## load data 
rodents = read.csv("Rodent_julys.csv", header = TRUE, row.names = 1)
dim(rodents)
View(rodents)

ants = read.csv("Ant_colony_numstakes.csv", header = TRUE, row.names = 1)
dim(ants)
View(ants)

# plants = 

## SYMMETRIC COCA
# log transform rodent data
rodents = log(rodents + 1)  #might consider swapping ant and rodent data positions here, not sure of the influence of either

## fit the model
ra.sym = coca(rodents ~ ., data = ants, method = "symmetric")
ra.sym
summary(ra.sym)
plot(ra.sym)


## PREDICTIVE COCA using SIMPLS and formula interface
ra.pred = coca(rodents ~ ., data = ants)