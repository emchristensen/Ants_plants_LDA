# Main function for doing LDA on different data sets
# using LDATS package

devtools::install_github("weecology/LDATS")

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
# LDA model -- setup parameters
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
# changepoint analysis (modified from Extreme-events-LDA)

source('changepoint_model.R')


# set up parameters for model -- rodents
year_continuous = 1970 + as.integer(julian(dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)
pop_matrix = as.matrix(lda_data)
# timeseries matrix
#ts_matrix = selected_lda@gamma
ts_matrix = pop_matrix
weights = rep(1,length(year_continuous))


# set up parameters for model -- winter annuals
year_continuous = winterann$year
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)
ts_matrix = winterann[,!names(winterann)=='year'] %>% as.matrix()
weights = rep(1,length(x$year_continuous))



# run models with 1:6 changepoints
cp_results_rodent = list()
cp_results_winterann = list()
i = 1
for (npts in 1:6) {
  cp_results = changepoint_model(ts_matrix, x, npts, weights = weights)
  cp_results_winterann[[i]] = cp_results
  i = i + 1
}



# ============================================================
# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
nvars = dim(ts_matrix)[2]

for (n in 1:length(cp_results_winterann)) {
  npoints = dim(cp_results_winterann[[n]]$saved)[1]
  mean_dev = mean(cp_results_winterann[[n]]$saved_lls * -2) + 2*(3*(nvars-1)*(npoints+1)+(npoints))
  print(c(npoints,mean_dev))
}

# best model is one with lowest mean deviation
cp_rodent = cp_results_rodent[[5]]
cp_winterann = cp_results_winterann[[4]]
# ========================================================================================
# some quick histograms of changepoint model results
hist(year_continuous[cp_rodent$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')
annual_hist(cp_rodent,year_continuous)
annual_hist(cp_winterann,year_continuous)

# turn changepoint results into data frame
df = as.data.frame(t(cp_rodent$saved[,1,])) %>% reshape2::melt()
df$value = year_continuous[df$value]

# find 95% confidence intervals on each changepoint:
quantile(df[df$variable=='V1','value'],probs=c(.025,.975)) %>% lubridate::date_decimal() %>% format('%d-%m-%Y')
quantile(df[df$variable=='V2','value'],probs=c(.025,.975)) %>% lubridate::date_decimal() %>% format('%d-%m-%Y')
quantile(df[df$variable=='V3','value'],probs=c(.025,.975)) %>% lubridate::date_decimal() %>% format('%d-%m-%Y')
quantile(df[df$variable=='V4','value'],probs=c(.025,.975)) %>% lubridate::date_decimal() %>% format('%d-%m-%Y')

# time series plots
cpts = find_changepoint_location(cp_rodent)
cpt_plot = get_ll_non_memoized_plot(ts_matrix,x,cpts,weights=weights)

