# LDA and changepoint analysis pipeline on rodent data -- uses VEM method
#
#  1. prepare data
#      - specific script for each data set
#  2a. run LDA with different number of topics and use AIC to select best model
#  2b. run LDA with best number of topics (determined by 2a)
#  3. run changepoint model
#  4. produce figures


library(topicmodels)

library(multipanelfigure)


source('get_data_for_LDA.r')

source('../Extreme-events-LDA/AIC_model_selection.R')
source('../Extreme-events-LDA/LDA_figure_scripts.R')
source('changepointmodel.r')
source('LDA-distance.R')

# ===================================================================
# prepare rodent data
dat = ant_colony_table(1977,2009,c(1,2,5,6,7,9,11,14,15,16,18,21,22))
dat_openings = ant_opening_table(1977,2009,c(1,2,5,6,7,9,11,14,15,16,18,21,22))
dat_abund = ant_bait_table(1977,2009,c(1,2,5,6,7,9,11,14,15,16,18,21,22))

# ===================================
# model parameters:
topic_min = 2
topic_max = 9
nspp=length(dat)
# ==================================================================
# select number of topics

# Fit a bunch of LDA models with different seeds
# Only use every other seed because consecutive seeds give identical results (!?)
seeds = 2*seq(200)

# repeat aic calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic_abund = repeat_VEM(dat_abund,
                         seeds,
                         topic_min=2,
                         topic_max=10)

# plot histogram of how many seeds chose how many topics
hist(best_ntopic[,1],breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),xlab='best # of topics', main='')

# ============================================================
ntopics=2
ldamodel = LDA(dat,ntopics, control = list(seed = 1),method='VEM')
beta1 = community_composition(ldamodel)
plot_community_composition(beta1)
plot_component_communities(ldamodel,ntopics,as.Date(row.names(dat),format='%Y'))
