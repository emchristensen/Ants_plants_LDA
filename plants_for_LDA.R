# 

library(dplyr)
#reshape

source('../Extreme-Events-LDA/AIC_model_selection.R')
source('LDA_figure_scripts_antsplants.R')


dat = read.csv("C:/Users/EC/Desktop/git/PortalData/Plants/Portal_plant_quadrats.csv",stringsAsFactors = F)
plant_list = read.csv('C:/Users/EC/Desktop/git/PortalData/Plants/Portal_plant_species.csv',stringsAsFactors = F)

# filter out unknown species
perennial = filter(plant_list,duration=='Perennial',commonname!='Unknown')
annual = filter(plant_list,duration=='Annual',commonname!='Unknown')

# only control plot data, separate by winter/summer/perennial community
controls = c(2,4,8,11,12,14,17,22)
#summer = filter(dat,season=='summer', plot %in% controls, year > 1982, species %in% annual$speciescode)
#winter = filter(dat,season=='winter', plot %in% controls, year > 1982, species %in% annual$speciescode)
perenn = filter(dat, plot %in% controls, year > 1988, species %in% perennial$speciescode)
# I can use earlier data if I restrict the quadrats: 1981-1982 only did 8 quads per plot
summer = filter(dat,season=='summer', plot %in% controls, quadrat %in% c(13,15,31,37,51,57,73,75), species %in% annual$speciescode)
winter = filter(dat,season=='winter', plot %in% controls, quadrat %in% c(11,17,33,35,53,55,71,77), species %in% annual$speciescode)

# aggregate to get site-wide total
summer_byplot = aggregate(summer$abundance,by=list(year = summer$year, species = summer$species),FUN=sum)
winter_byplot = aggregate(winter$abundance,by=list(year = winter$year, species = winter$species),FUN=sum)
perenn_byplot = aggregate(perenn$abundance,by=list(year = perenn$year, species = perenn$species),FUN=sum)

# convert to table
summer_table = reshape(summer_byplot,timevar='species',idvar='year',direction='wide')
summer_table[is.na(summer_table)] = 0
row.names(summer_table) = summer_table$year
summer_table = summer_table[order(summer_table$year),-1]

winter_table = reshape(winter_byplot,timevar='species',idvar='year',direction='wide')
winter_table[is.na(winter_table)] = 0
row.names(winter_table) = winter_table$year
winter_table = winter_table[order(winter_table$year),-1]

perenn_table = reshape(perenn_byplot,timevar='species',idvar='year',direction='wide')
perenn_table[is.na(perenn_table)] = 0
row.names(perenn_table) = perenn_table$year
perenn_table = perenn_table[order(perenn_table$year),-1]


# ===========================================================================================
# LDA perennials

seeds = 2*seq(20)

# repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(perenn_table,
                         seeds,
                         topic_min=2,
                         topic_max=24)

SEED = 2
ntopics = 4
ldamodel = LDA(perenn_table,ntopics, control = list(seed = SEED),method='VEM')


#get parameter estimates
z=posterior(ldamodel)
commun.plot=z$topics
commun.spp=z$term

# look at what's in the groups
structure(round(exp(ldamodel@beta), 3), dimnames = list(NULL, ldamodel@terms))

dates = c(1988:2009,2011:2017)

plot_component_communities(ldamodel,ntopics,dates)

# =============================
# LDA summer annuals

seeds = 2*seq(20)

# repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(summer_table,
                         seeds,
                         topic_min=18,
                         topic_max=26)

SEED = 4
ntopics = 5
ldamodel = LDA(summer_table,ntopics, control = list(seed = SEED),method='VEM')

#get parameter estimates
z=posterior(ldamodel)
commun.plot=z$topics
commun.spp=z$term

# look at what's in the groups
structure(round(exp(ldamodel@beta), 3), dimnames = list(NULL, ldamodel@terms))

dates = c(1981:2002,2004:2008,2011,2014:2016)

plot_component_communities(ldamodel,ntopics,dates)

# =============================
# LDA winter annuals

seeds = 2

# repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(winter_table,
                         seeds,
                         topic_min=10,
                         topic_max=20)

SEED = 2
ntopics = 5
ldamodel = LDA(winter_table,ntopics, control = list(seed = SEED),method='VEM')


#get parameter estimates
z=posterior(ldamodel)
commun.plot=z$topics
commun.spp=z$term

# look at what's in the groups
structure(round(exp(ldamodel@beta), 3), dimnames = list(NULL, ldamodel@terms))

dates = c(1983:1995,1997:1999,2001:2009,2012:2017)

plot_component_communities(ldamodel,ntopics,dates)

