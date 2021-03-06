# Ants_plants_LDA
Code and data for project comparing change over time in three taxa recorded at Portal: rodents, ants, and plants.

## Data files
  * __Rodent_summer_table.csv__ table of species counts by year, can be used to run LDA. 1978-2015. Control plots combined: 2,4,8,11,12,14,17,22. Average of summer censuses for each year (June-Sept), rounded to nearest integer.
  * __Ant_colony_presence.csv__ table of species by year, entries are avg # of stakes per plot (0-49) where species was present in that year. 1977-2009. Control plots: 2,11,14,22
  * __WinterAnnuals.csv__ table of species by year, winter annual plant community. 1983-2014. Sum of abundances on all quadrats on all control plots per year (plots 2,4,8,11,12,14,17,22). Rare species removed (present in 10% or fewer of the censuses)
  * __SummerAnnuals.csv__ table of species by year, summer annual plant community. 1983-2014. Same processing as winter annuals.
