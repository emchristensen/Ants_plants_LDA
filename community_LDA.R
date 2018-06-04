# Main function for doing LDA on different data sets
# using LDATS package

#devtools::install_github("weecology/LDATS")

library(LDATS)
library(dplyr)
source('LDA_figure_scripts_antsplants.R')


# =========================================================================================
# LDATS example from repo
data(rodents)
lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
ts_data <- data.frame("time" = data.frame(rodents)[ , "newmoon"])

r_LDA <- LDATS::parLDA(data = lda_data, ntopics = 2:5, nseeds = 2, ncores = 4)
r_LDATS <- LDATS::LDA_TS(lda_data, ts_data, formula = c(~1, ~time),
                         ntopics = 2:5, nseeds = 2, ncores = 4, nit = 100)

# =======================================================================================
# load data

# rodent data
data(rodents)
lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
ts_data <- data.frame("time" = data.frame(rodents)[ , "newmoon"])

# dates to go with count data
moondat = read.csv(text=RCurl::getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"),stringsAsFactors = F)
moondat$date = as.Date(moondat$censusdate)

period_dates = filter(moondat,period %in% ts_data$time) %>% select(period,date)
dates = period_dates$date

# plant data
summerann = read.csv('SummerAnnuals.csv')
winterann = read.csv('WinterAnnuals.csv')
wann_lda = dplyr::select(winterann,-c(year))
sann_lda = dplyr::select(summerann,-c(year))

# ======================================================================================
# setup parameters
nseeds = 2
ncores = 4
data = lda_data
ts_data = ts_data

# run a bunch of LDA models with 2-5 topics
r_LDA = LDATS::parLDA(data=data, ntopics = 5, nseeds = nseeds, ncores = ncores)

# choose best lda model (number of topics)
selected_lda = LDATS::LDA_select(lda_models = r_LDA, LDA_eval = quote(AIC), correction = T, LDA_selector = quote(min))

# plot this
ntopics = selected_lda@k
composition = community_composition(selected_lda)
composition
plot_community_composition(composition)
plot_component_communities(selected_lda,ntopics = ntopics,xticks = ts_data$time)

# ======================================================================================
# changepoint analysis (from Extreme-events-LDA)

source('changepoint_model.R')
#source('../Extreme-events-LDA/changepointmodel.r')

# set up parameters for model
year_continuous = 1970 + as.integer(julian(dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)

pop_matrix = as.matrix(lda_data)

# timeseries matrix
ts_matrix = selected_lda@gamma
ts_matrix = pop_matrix

# run models with 1, 2, 3, 4, 5 changepoints
cp_results_rodent = changepoint_model(ts_matrix, x, 1, weights = rep(1,length(year_continuous)))
cp_results_rodent2 = changepoint_model(selected_lda, x, 2, weights = rep(1,length(year_continuous)))
cp_results_rodent3 = changepoint_model(selected_lda, x, 3, weights = rep(1,length(year_continuous)))
cp_results_rodent4 = changepoint_model(selected_lda, x, 4, weights = rep(1,length(year_continuous)))
cp_results_rodent5 = changepoint_model(selected_lda, x, 5, weights = rep(1,length(year_continuous)))

# some quick histograms of changepoint model results
hist(year_continuous[cp_results_rodent$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')
annual_hist(cp_results_rodent4,year_continuous)

# turn changepoint results into data frame
df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,])) %>% reshape::melt()
df_4$value = year_continuous[df_4$value]

# find 95% confidence intervals on each changepoint:
quantile(df_4[df_4$variable=='V1','value'],probs=c(.025,.975)) %>% lubridate::date_decimal() %>% format('%d-%m-%Y')
quantile(df_4[df_4$variable=='V2','value'],probs=c(.025,.975)) %>% lubridate::date_decimal() %>% format('%d-%m-%Y')
quantile(df_4[df_4$variable=='V3','value'],probs=c(.025,.975)) %>% lubridate::date_decimal() %>% format('%d-%m-%Y')
quantile(df_4[df_4$variable=='V4','value'],probs=c(.025,.975)) %>% lubridate::date_decimal() %>% format('%d-%m-%Y')

# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))

