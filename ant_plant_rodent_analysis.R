
devtools::install_github("weecology/portalr")
devtools::install_github("weecology/LDATS")
library(LDATS)
source('create_data_tables.R')
source('figure_scripts.R')

###################################
# TO DO:
#   1. How to deal with rodent data. Annual mean? Summer mean? Biennial? Keep it monthly?
#   2. should we use the same plots for all taxa, or increase the power of each time series by using all cotrols for that specific taxon?






###############################################
# Data

winterannuals = seasonal_plant_table(selected_plots=c(11,14),
                                     plant_community='Annuals',
                                     summer_winter='winter')
winterannuals = winterannuals[,-2]

summerannuals = seasonal_plant_table(selected_plots=c(11,14),
                                     plant_community='Annuals',
                                     summer_winter='summer')
summerannuals = summerannuals[,-2]

ant_table = ant_colony_presence(selected_plots=c(11,14))

rodent_summer_table = read.csv('Rodent_summer_table.csv')


################################################
# Look at data

plot_pop_ts(winterannuals)
plot_pop_ts(summerannuals)
plot_pop_ts(ant_table)
plot_pop_ts(rodent_summer_table)


################################################
# # LDA -- WIP
# 
# control_time_LDA = LDATS::LDA(data = select(rodent_summer_table, -year), ntopics =  c(2, 3, 4, 5),
#                               nseeds = 4, ncores = 4)
# control_time_LDA_use = LDATS:::LDA_select(lda_models = control_time_LDA, LDA_eval = quote(AIC), correction = TRUE,
#                                           LDA_selector = quote(min))
# plot.LDA(control_time_LDA_use)


################################################
# PCA
# try this https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/

pcaant = PCA(ant_table[,-1])
head(pcaant$x)
head(pcaant$rotation)
plot(pcaant)
plot(pcaant$sdev)
