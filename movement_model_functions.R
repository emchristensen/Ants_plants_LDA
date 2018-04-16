library(nlme)

#' @title Gradual Linear Dynamics
#' 
#' @param quad data frame of a single timeseries of a quadrat/plot: contains columns "dist" and "nYrs"
#'    
gradual_linear_dyn = function(quad) {
  nomad.Mod <- NA
  nomad.Mod <- try(gnls(dist ~ C + M*(nYrs), data = quad,
                        correlation = corAR1(),
                        start = list(C=0,M=0)))
  return(nomad.Mod)
}

#' @title Reversible Dynamics
#'
#' @param quad data frame of a single timeseries of a quadrat/plot: contains columns "dist" and "nYrs"
#'
reversible_dyn = function(quad) {
  migr.Mod <- NA
  migr.Mod <- try(gnls(dist ~ AsymA/(1+exp((xmidA-nYrs)/scal1)) + (-AsymB /(1 + exp((xmidB-nYrs)/scal2))),
                       data = quad, 
                       correlation=corAR1(),
                       start = list(AsymA = 0.6 , AsymB = 0.6, xmidA = 10,xmidB = 20, scal1 = 1, scal2 = 1),     
                       control=nlmeControl(maxIter=200, pnlsMaxIter=200, niterEM=400, returnObject=TRUE)))
  return(migr.Mod)
}

#' @title Stable Dynamics
#' 
#' @param quad data frame of a single timeseries of a quadrat/plot: contains columns "dist" and "nYrs"
#' 
stable_dyn = function(quad) {
  null.mod<-NA  
  asym.HRmod<-NA
  stab.Mod <- NA
  ### Step 1 Null Model
  null.mod <- try(gnls(dist ~ A, data = quad, 
                       correlation = corAR1(),
                       start = list(A = mean(quad[,'dist'])), 
                       control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=F))
  ## Step 2 Asymp model
  asym.HRmod <- try(gnls(dist ~ Asym*(1-exp(lrc*nYrs)), data = quad, 
                         correlation = corAR1(),
                         start = list(Asym =summary(null.mod)$tTable[1],lrc=-0.059), 
                         control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=F))
  ## Step 3 Full stability model
  stab.Mod <- try(gnls(dist ~ (Asym)*(1 - exp(lrc*(nYrs))),data = quad,
                       correlation = corAR1(),
                       start = c(Asym = summary(asym.HRmod)$tTable[1], lrc = summary(asym.HRmod)$tTable[2]),
                       control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=F))
  return(stab.Mod)
}

#' @title Abrupt dynamics
#' 
#' @param quad data frame of a single timeseries of a quadrat/plot: contains columns "dist" and "nYrs"
#' 
abrupt_dyn= function(quad) {
  null.mod<-NA  
  asym.HRmod<-NA
  disp.Mod <- NA
  
  ## Step 1 NULL model
  null.mod <- try(gnls(dist ~ A, data = quad, 
                       correlation = corAR1(),
                       start = c(A = mean(quad[,'dist'])), 
                       control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=F))
  ## Step 2 Asymp model
  asym.HRmod <- try(gnls(dist ~ Asym*(1-exp(lrc*nYrs)), data = quad, 
                         correlation = corAR1(),
                         start = c(Asym =summary(null.mod)$tTable[1],lrc=-0.059), 
                         control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=F))
  ## Step 3 Disp model 
  disp.Mod <- try(gnls(dist ~ (Asym)/(1 + exp((xmid-nYrs)/scal)),data = quad,
                       correlation = corAR1(), 
                       na.action = na.exclude,
                       start = c(Asym = summary(asym.HRmod)$tTable[1], xmid = 10, scal = 0.5),
                       control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=F))
  return(disp.Mod)
}

#' @title plot_fits
#' @description plot data for single plot as well as the 4 different fits
#' 
#' 
plot_fits = function(df) {
  plot(df$year,df$dist,type='b',ylim=c(0,1))
  lines(df$year,df$nomad_gnls,col='red')
  lines(df$year,df$mig_gnls,col='blue')
  lines(df$year,df$stab_gnls,col='forestgreen')
  lines(df$year,df$abrupt_gnls,col='purple')
  legend(2000,.4, legend=c('data','nomad/linear','migration','stable','abrupt change'),
         col=c('black','red','blue','forestgreen','purple'), lty=1, cex=0.8)
  
}
