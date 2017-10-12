library(RCurl)
library(dplyr)
library(reshape2)

#' Create ant species abundance table -- number of colonies --- WIP!!!!!
#'
#' Processes ant census data so it can be used for LDA analysis
#' This script parallels rodent_data_for_LDA.R in Extreme-events-LDA project
#'
#' @param yr_first first year of data desired
#' @param yr_last last year of data desired
#' @param selected_plots vector of plot numbers to be included
#' 
#' @return Table of colony counts species per period
#' 
#' @example
#' a_table = ant_colony_table(1977,2009,c(1,2,5,6,7,9,11,14,15,16,18,21,22))
#'
ant_colony_table = function(yr_first,yr_last,selected_plots=1:24) {
  
  # retrieve current version of ant data
  ants = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Ants/Portal_ant_colony.csv"),
                  na.strings=c("",-99), stringsAsFactors = FALSE)
  
  # extract desired data by year and plot; remove unkn species and 'sole xylo' which was not censused in every year; remove 1981 because rare species 
# not censused
  A = filter(ants, Year>=yr_first, Year<=yr_last,
             Plot %in% selected_plots) %>% select(Year,Plot,Species,Colonies)
  
  # sum number of colonies per year
  A_colonies = aggregate(A$Colonies,by=list(Year = A$Year, Species = A$Species),FUN=sum)
  A_table = reshape(A_colonies, idvar = 'Year', timevar = 'Species', direction = 'wide')
  
  # change column names
  #names(A_table)[-1] %>% strsplit('\\.') %>% unlist()
  
  #order rows by year
  A_table = A_table[order(A_table$Year),]
  #fill in NA with 0
  A_table[is.na(A_table)] = 0
  row.names(A_table) = A_table$Year
  
  return(A_table[,-1])
}




#' Create ant species abundance table -- number of openings -- WIP!!!
#'
#' Processes ant census data so it can be used for LDA analysis
#' This script parallels rodent_data_for_LDA.R in Extreme-events-LDA project
#'
#' @param yr_first first year of data desired
#' @param yr_last last year of data desired
#' @param selected_plots vector of plot numbers to be included
#' 
#' @return Table of colony opening counts per species per period
#' 
#' @example
#' a_table = ant_opening_table(1977,2009,c(1,2,5,6,7,9,11,14,15,16,18,21,22))
#'
ant_opening_table = function(yr_first,yr_last,selected_plots=1:24) {
  
  # retrieve current version of ant data
  ants = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Ants/Portal_ant_colony.csv"),
                  na.strings=c("",-99), stringsAsFactors = FALSE)
  
  # extract desired data by year and plot; remove unkn species and 'sole xylo' which was not censused in every year; remove 1981 because rare species 
  # not censused
  A = filter(ants, Year>=yr_first, Year<=yr_last,
             Plot %in% selected_plots,
             !(Species %in% c('unkn','sole xylo')),
             Year != 1981) %>% select(Year,Plot,Species,Openings)
  
  # sum number of colonies per year
  A_colonies = aggregate(A$Openings,by=list(Year = A$Year, Species = A$Species),FUN=sum)
  A_table = reshape(A_colonies, idvar = 'Year', timevar = 'Species', direction = 'wide')
  
  # change column names
  #names(A_table)[-1] %>% strsplit('\\.') %>% unlist()
  
  #order rows by year
  A_table = A_table[order(A_table$Year),]
  #fill in NA with 0
  A_table[is.na(A_table)] = 0
  row.names(A_table) = A_table$Year
  
  return(A_table[,-1])
}




#' Create ant species abundance table -- bait data
#'
#' Processes ant census data so it can be used for LDA analysis
#' This script parallels rodent_data_for_LDA.R in Extreme-events-LDA project
#'
#' @param selected_plots vector of plot numbers to be included
#' 
#' @return Table of counts per species per year
#' 
#' @example
#' bait_table = ant_bait_table(selected_plots = c(2,11,14,22))
#' write.csv(bait_table,'ant_bait_table.csv',row.names=T)
#'
ant_bait_table = function(selected_plots=1:24) {
  
  # retrieve current version of ant data
  ants = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Ants/Portal_ant_bait.csv"),
                  stringsAsFactors = FALSE)
  
  # extract desired data by year and plot; remove unkn species 
  A = filter(ants, plot %in% selected_plots) %>% select(year,plot,species,abundance)
  
  # sum abundance per year
  A_abund = aggregate(A$abundance,by=list(year = A$year, species = A$species),FUN=sum)
  A_table = reshape(A_abund, idvar = 'year', timevar = 'species', direction = 'wide')
  
  # change column names
  #names(A_table)[-1] %>% strsplit('\\.') %>% unlist()
  
  #order rows by year
  A_table = A_table[order(A_table$year),]
  #fill in NA with 0
  A_table[is.na(A_table)] = 0
  row.names(A_table) = A_table$year
  
  return(A_table[,-1])
}

