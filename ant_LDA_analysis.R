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


#source('get_data_for_LDA.r')

source('../Extreme-events-LDA/AIC_model_selection.R')
source('LDA_figure_scripts_antsplants.R')
source('../Extreme-events-LDA/changepointmodel.r')
source('../Extreme-events-LDA/LDA-distance.R')

# ===================================================================
# prepare ant data
dat = read.csv('Ant_colony_numstakes.csv',as.is=T)
dat$date = substr(dat$index,1,4) %>% as.numeric()
dat = dat %>% select(-index)
yrdat = aggregate(.~date,data = dat,sum)
antdat = yrdat %>% select(-date)

# ===================================
# model parameters:
topic_min = 2
topic_max = 9
nspp=length(antdat)
# ==================================================================
# select number of topics

# Fit a bunch of LDA models with different seeds
# Only use every other seed because consecutive seeds give identical results (!?)
seeds = 2*seq(50)

# repeat aic calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic_abund = repeat_VEM(antdat,
                         seeds,
                         topic_min=2,
                         topic_max=10)

# plot histogram of how many seeds chose how many topics
hist(best_ntopic_abund$k,breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10,10.5),xlab='best # of topics', main='')

# 5 topics seems to be best, but 6, 7, 8 are also well-represented

# ============================================================
ntopics=5
ldamodel = LDA(antdat,ntopics, control = list(seed = 1),method='VEM')
beta1 = community_composition(ldamodel)
plot_community_composition(beta1)
plot_component_communities(ldamodel,ntopics,as.Date(as.character(yrdat$date),format='%Y'))
