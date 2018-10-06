library(ggplot2)

#' @title plot population time series
#'
#' @param dattable table of species counts by time. first column is time column
#'
plot_pop_ts = function(dattable) {
  dat = tidyr::gather(dattable,species,n,2:dim(dattable)[2])
  ggplot(dat,aes(x=year,y=n,colour=species)) +
    geom_line()
}