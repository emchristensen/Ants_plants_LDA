# Main function for doing LDA on different data sets
# using LDATS package

devtools::install_github("weecology/LDATS")

library(LDATS)
source('LDA_figure_scripts_antsplants.R')


# =========================================================================================
# example from repo
data(rodents)
lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
ts_data <- data.frame("time" = data.frame(rodents)[ , "newmoon"])

r_LDA <- LDATS::parLDA(data = lda_data, ntopics = 2:5, nseeds = 2, ncores = 4)
r_LDATS <- LDATS::LDA_TS(lda_data, ts_data, formula = c(~1, ~time),
                         ntopics = 2:5, nseeds = 2, ncores = 4, nit = 100)

# =======================================================================================
# 

# load data
summerann = read.csv('SummerAnnuals.csv')
winterann = read.csv('WinterAnnuals.csv')
wann_lda = dplyr::select(winterann,-c(year))
sann_lda = dplyr::select(summerann,-c(year))

# setup parameters
nseeds = 2
ncores = 4
data = sann_lda
ts_data = summerann$year

# run a bunch of LDA models with 2-5 topics
r_LDA = LDATS::parLDA(data=data, ntopics = 2:4, nseeds = nseeds, ncores = ncores)

# choose optimal number of topics
selected_lda = LDATS::LDA_select(lda_models = r_LDA, LDA_eval = quote(AIC), correction = T, LDA_selector = quote(min))

# plot this
ntopics = selected_lda@k
composition = community_composition(selected_lda)
composition
plot_community_composition_gg(composition)
plot_component_communities(selected_lda,ntopics = ntopics,xticks = ts_data)


# changepoint analysis

