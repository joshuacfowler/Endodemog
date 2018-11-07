## Authors: Josh and Tom
## Purpose: Create a script that imports Endodemog data, perform all raw data manipulation,
## and create an .RData object that can be loaded for analysis
## Last update: 23 Oct 2018
######################################################

library(tidyverse)
library(reshape2)
library(lubridate)
library(readxl)

# read in data from POAL
Poal_data <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POALnew complete with 2016 data.xlsx", sheet = "POAL")
Poal_data_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POALnew complete with 2016 data.xlsx", sheet = "POAL (NEW) recruits")
Poal_data_old <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POALold complete with 2016 data.xlsx", sheet = "POAL (OLD)")
Poal_data_old_r <- read_xlsx("/Users/joshuacfowler/Dropbox/EndodemogData/rawdatafilesbyspecies/POALold complete with 2016 data.xlsx", sheet = "POAL (OLD) recruits")

# A list of the columns that are redundant for now, or will be used in our seed estimates
pmain <- Poal_data %>%  
  select(-survive4may08,-`Surv4/11`, -Aphids1, -aphids2, -aphid3, -contains("SB"),-notes2,
         -notes2__1, -notes3, -notes4, -notes5,-`notes5/4/2008`,-notes6, -notes7, 
         -notes8, -notes9, -othernotes3, -data, -`Planting Notes`, -TAG, -Endocheck, 
         -EndoDateCheck, -EndoDateCheck_Day, -EndoDateCheck_Month, -EndoDateCheck_Year,
         -`Coll Date`, -TotTillers11sep08, -tilleradjust, -VisualEST1, -endoyr, 
         -seed2surv, -seed3surv, -seed4surv, -actseed2, -contains("Est"), -contains("Hbv"), 
         -contains("Lvs"), -contains("Infl"), -contains("CB"), -contains("Prop"))

# Combining measurements across years for the “New” data ------------------


## Combining measurements across years for Surv, Growth, and Flowering using melt
## Recoding those measurements for the year they are taken

psurv <- Poal_data %>%
  rename("Birth Year" = "Planted Date") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"),
       measure.var = c("survive1", "survive2", "Survive3", "Survive4", 
                       "Survive5", "survive6", "survive7", "survive8", 
                       "survive9"),
       value.name = "surv") 
psurv$year<- ifelse(psurv$variable == "survive1", 2008, ifelse(psurv$variable  == "survive2", 2009, ifelse(psurv$variable  == "Survive3", 2010, ifelse(psurv$variable  == "Survive4", 2011, ifelse(psurv$variable  == "Survive5", 2012, ifelse(psurv$variable  == "survive6", 2013,ifelse(psurv$variable == "survive7", 2014,ifelse(psurv$variable == "survive8", 2015,ifelse(psurv$variable  == "survive9", 2016, NA)))))))))
# View(psurv)

pgrow <- Poal_data %>% 
  rename("Birth Year" = "Planted Date") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("Tottillers1", "TotTillers2", "TotTillers3",
                       "TotTillers4", "TotTillers5", "TotTillers6", 
                       "TotTillers7", "TotTillers8", "TotTillers9"), 
       value.name = "size") 
pgrow$year<- ifelse(pgrow$variable == "Tottillers1", 2008, ifelse(pgrow$variable  == "TotTillers2", 2009, ifelse(pgrow$variable  == "TotTillers3", 2010, ifelse(pgrow$variable  == "TotTillers4", 2011, ifelse(pgrow$variable  == "TotTillers5", 2012, ifelse(pgrow$variable  == "TotTillers6", 2013, ifelse(pgrow$variable == "TotTillers7", 2014, ifelse(pgrow$variable == "TotTillers8", 2015, ifelse(pgrow$variable  == "TotTillers9", 2016, NA)))))))))
# View(pgrow)

pflw <- Poal_data %>% 
  rename("Birth Year" = "Planted Date") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", "Birth Year", 
                  "TRT", "Plant"), 
       measure.var = c("Flwtillers1", "FlwTillers2", "FlwTillers3", 
                       "FlwTillers4", "FlwTillers5", "FlwTillers6", 
                       "FlwTillers7", "FlwTillers8", "FlwTillers9"), 
       value.name = "flw") 
pflw$year<- ifelse(pflw$variable == "Flwtillers1", 2008, ifelse(pflw$variable  == "FlwTillers2", 2009, ifelse(pflw$variable  == "FlwTillers3", 2010, ifelse(pflw$variable  == "FlwTillers4", 2011, ifelse(pflw$variable  == "FlwTillers5", 2012, ifelse(pflw$variable  == "FlwTillers6", 2013,ifelse(pflw$variable == "FlwTillers7", 2014,ifelse(pflw$variable == "FlwTillers8", 2015,ifelse(pflw$variable  == "FlwTillers9", 2016, NA)))))))))
# View(pflw)

pmerge_sg <- merge(psurv, pgrow, by = c( "plot", "pos", "tag", "Endo", 
                                              "Loc'n", "Birth Year", "TRT",
                                              "Plant", "year"))
# View(pmerge_sg)

pmerge_sgf <- merge(pmerge_sg, pflw, by = c( "plot", "pos", "tag", "Endo", 
                                             "Loc'n", "Birth Year", "TRT",
                                             "Plant", "year"))
# View(pmerge_sgf)

# getting a dataframe with t and t_1
pmerge_t1 <-pmerge_sgf %>%
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(pmerge_t1)

pmerge_t <-pmerge_sgf %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, surv_t = surv, size_t = size, flw_t = flw) 
# View(pmerge_t)
## merge and set origin, coded as 0 for original plants and 1 for recruits
pmerge <- pmerge_t1 %>% 
  full_join(pmerge_t, by = c("plot", "pos", "tag", "Endo", 
                             "Loc'n", "Birth Year", "TRT",
                             "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
# View(pmerge)


# Combining measurements across years for the "New" recruits data ---------------


## Combining measurements across years for the recruits data
## recoding for the year of measurement
## merging these measurements into one dataframe
rsurv <- Poal_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(survive10 = NA) %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"), 
       measure.var = c("survive10", "Survive11", "Survive12", "Survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
rsurv$year<- ifelse(rsurv$variable == "survive10", 2010, ifelse(rsurv$variable == "Survive11", 2011, ifelse(rsurv$variable  == "Survive12", 2012, ifelse(rsurv$variable  == "Survive13", 2013, ifelse(rsurv$variable  == "Survive14", 2014, ifelse(rsurv$variable  == "Survive15", 2015, ifelse(rsurv$variable  == "Survive16", 2016, NA)))))))
# View(rsurv)

rgrow <- Poal_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TOTtiller10","TOTtiller11", "TOTtiller12", 
                       "TOTtiller13", "TOTtiller14", "TOTtiller15", 
                       "TOTtiller16"),
       value.name = "size") 
rgrow$year<- ifelse(rgrow$variable == "TOTtiller10", 2010, ifelse(rgrow$variable == "TOTtiller11", 2011, ifelse(rgrow$variable  == "TOTtiller12", 2012, ifelse(rgrow$variable  == "TOTtiller13", 2013, ifelse(rgrow$variable  == "TOTtiller14", 2014, ifelse(rgrow$variable  == "TOTtiller15", 2015, ifelse(rgrow$variable  == "TOTtiller16", 2016, NA)))))))
# View(rgrow)

rflw <- Poal_data_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
rflw$year<- ifelse(rflw$variable == "FLWtiller10", 2010, ifelse(rflw$variable == "FLWtiller11", 2011, ifelse(rflw$variable  == "FLWtiller12", 2012, ifelse(rflw$variable  == "FLWtiller13", 2013, ifelse(rflw$variable  == "FLWtiller14", 2014, ifelse(rflw$variable  == "FLWtiller15", 2015, ifelse(rflw$variable  == "FLWtiller16", 2016, NA)))))))
# View(rflw)

rmerge_sg <- merge(rsurv, rgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(rmerge_sg)

rmerge_sgf <- merge(rmerge_sg, rflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(rmerge_sgf)

## getting a dataframe with time t and t_1
rmerge_t1 <-rmerge_sgf %>%
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(rmerge_t1)

rmerge_t <-rmerge_sgf %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, surv_t = surv, size_t = size, flw_t = flw) 
# View(rmerge_t)

rmerge <- rmerge_t1 %>% 
  full_join(rmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
# View(rmerge)




# Combining measurements across years for the “Old” data ------------------


## Combining data across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe
poldsurv <- Poal_data_old %>%
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>%
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", 
                  "Birth Year", "TRT", "Plant"),
       measure.var = c("survive1", "Survive3", "survive4", 
                       "survive5", "survive6", "survive7", 
                       "survive8", "survive9", "survive10"),
       value.name = "surv") 
poldsurv$year<- ifelse(poldsurv$variable == "survive1", 2008, ifelse(poldsurv$variable  == "Survive3", 2009, ifelse(poldsurv$variable  == "survive4", 2010, ifelse(poldsurv$variable  == "survive5", 2011, ifelse(poldsurv$variable  == "Survive6", 2012, ifelse(poldsurv$variable  == "survive7", 2013,ifelse(poldsurv$variable == "survive8", 2014,ifelse(poldsurv$variable == "survive9", 2015,ifelse(poldsurv$variable  == "survive10", 2016, NA)))))))))
# View(poldsurv)

poldgrow <- Poal_data_old %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n", 
                  "Birth Year", "TRT", "Plant"), 
       measure.var = c("TotTillers1", "TotTillers3","TotTillers4", 
                       "TotTillers5", "TotTillers6","TotTillers7", 
                       "TotTillers8", "TotTillers9", "TotTillers10"), 
       value.name = "size") 
poldgrow$year<- ifelse(poldgrow$variable == "TotTillers1", 2008, ifelse(poldgrow$variable  == "TotTillers3", 2009, ifelse(poldgrow$variable  == "TotTillers4", 2010, ifelse(poldgrow$variable  == "TotTillers5", 2011, ifelse(poldgrow$variable  == "TotTillers6", 2012, ifelse(poldgrow$variable  == "TotTillers7", 2013, ifelse(poldgrow$variable == "TotTillers8", 2014, ifelse(poldgrow$variable == "TotTillers9", 2015, ifelse(poldgrow$variable  == "TotTillers10", 2016, NA)))))))))
# View(poldgrow)

poldflw <- Poal_data_old %>% 
  rename("Birth Year" = "Date", "plot" = "PLOT", "pos" = "POS", "tag" = "TAG") %>% 
  melt(id.var = c("plot","pos", "tag", "Endo", "Loc'n",
                  "Birth Year", "TRT", "Plant"), 
       measure.var = c("Flwtillers1", "FlwTillers3", "FlwTillers4", 
                       "FlwTillers5", "FlwTillers6", "FlwTillers7", 
                       "FlwTillers8", "FlwTillers9", "FlwTillers10"), 
       value.name = "flw") 
poldflw$year<- ifelse(poldflw$variable == "Flwtillers1", 2008, ifelse(poldflw$variable  == "FlwTillers3", 2009, ifelse(poldflw$variable  == "FlwTillers4", 2010, ifelse(poldflw$variable  == "FlwTillers5", 2011, ifelse(poldflw$variable  == "FlwTillers6", 2012, ifelse(poldflw$variable  == "FlwTillers7", 2013,ifelse(poldflw$variable == "FlwTillers8", 2014,ifelse(poldflw$variable == "FlwTillers9", 2015,ifelse(poldflw$variable  == "FlwTillers10", 2016, NA)))))))))
# View(poldflw)


poldmerge_sg <- merge(poldsurv, poldgrow, by = c( "plot","pos", "tag", "Endo", 
                                         "Loc'n", "Birth Year", "TRT",
                                         "Plant", "year"))
# View(poldmerge_sg)

poldmerge_sgf <- merge(poldmerge_sg, poldflw, by = c( "plot","pos", "tag", "Endo", 
                                             "Loc'n", "Birth Year", "TRT",
                                             "Plant", "year"))
# View(poldmerge_sgf)

# getting a dataframe with t and t_1
poldmerge_t1 <-poldmerge_sgf %>%
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(pmerge_t1)

poldmerge_t <-poldmerge_sgf %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, surv_t = surv, size_t = size, flw_t = flw) 
# View(pmerge_t)

poldmerge <- poldmerge_t1 %>% 
  full_join(poldmerge_t, by = c("plot","pos", "tag", "Endo", 
                             "Loc'n", "Birth Year", "TRT",
                             "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 0) %>% 
  mutate(`Birth Year` = year(`Birth Year`))
# View(poldmerge)


# Combining measurements across years for the “Old” recruits data --------

## Combining data for recruits across years from the "Old" excel sheet
## recoding the values for year
## Merging these into one dataframe
roldsurv <- Poal_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  mutate(survive09 = "NA") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("survive09", "Survive10", "Survive11", "Survive12", "Survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
roldsurv$year<- ifelse(roldsurv$variable == "survive09", 2009, ifelse(roldsurv$variable == "Survive10", 2010, ifelse(roldsurv$variable == "Survive11", 2011, ifelse(roldsurv$variable  == "Survive12", 2012, ifelse(roldsurv$variable  == "Survive13", 2013, ifelse(roldsurv$variable  == "Survive14", 2014, ifelse(roldsurv$variable  == "Survive15", 2015, ifelse(roldsurv$variable  == "Survive16", 2016, NA))))))))
# View(roldsurv)


roldgrow <- Poal_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("TOTtiller09", "TOTtiller10","TOTtiller11", "TOTtiller12", 
                       "TOTtiller13", "TOTtiller14", "TOTtiller15", 
                       "TOTtiller16"),
       value.name = "size") 
roldgrow$year<- ifelse(roldgrow$variable == "TOTtiller09", 2009, ifelse(roldgrow$variable == "TOTtiller10", 2010, ifelse(roldgrow$variable == "TOTtiller11", 2011, ifelse(roldgrow$variable  == "TOTtiller12", 2012, ifelse(roldgrow$variable  == "TOTtiller13", 2013, ifelse(roldgrow$variable  == "TOTtiller14", 2014, ifelse(roldgrow$variable  == "TOTtiller15", 2015, ifelse(roldgrow$variable  == "TOTtiller16", 2016, NA))))))))
# View(roldgrow)

roldflw <- Poal_data_old_r %>%
  rename("Birth Year" = "Date", "tag" = "Tag", "plot" = "Plot", "pos" = "RecruitNo") %>% 
  melt(id.var = c("plot", "pos", "tag", "Endo", "Birth Year"),
       measure.var = c("FLWTiller09", "FLWtiller10","FLWtiller11", "FLWtiller12", 
                       "FLWtiller13", "FLWtiller14", "FLWtiller15",
                       "FLWtiller16"),
       value.name = "flw") 
roldflw$year<- ifelse(roldflw$variable == "FLWTiller09", 2009, ifelse(roldflw$variable == "FLWtiller10", 2010, ifelse(roldflw$variable == "FLWtiller11", 2011, ifelse(roldflw$variable  == "FLWtiller12", 2012, ifelse(roldflw$variable  == "FLWtiller13", 2013, ifelse(roldflw$variable  == "FLWtiller14", 2014, ifelse(roldflw$variable  == "FLWtiller15", 2015, ifelse(roldflw$variable  == "FLWtiller16", 2016, NA))))))))
# View(roldflw)

roldmerge_sg <- merge(roldsurv, roldgrow, by = c( "plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(roldmerge_sg)

roldmerge_sgf <- merge(roldmerge_sg, roldflw, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year"))
# View(roldmerge_sgf)

## getting a dataframe with time t and t_1
roldmerge_t1 <-roldmerge_sgf %>%
  rename(year_t1 = year, surv_t1 = surv, size_t1 = size, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1)
# View(roldmerge_t1)

roldmerge_t <-roldmerge_sgf %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, surv_t = surv, size_t = size, flw_t = flw) 
# View(roldmerge_t)

roldmerge <- roldmerge_t1 %>% 
  full_join(roldmerge_t, by = c("plot", "pos", "tag", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = 1) %>% 
  mutate(`Loc'n` = NA) %>% 
  mutate(TRT = NA) %>% 
  mutate(Plant = NA)
# View(roldmerge)










# Combining the old and new and original and recruit dataframes ---------
pmerge <- pmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                     "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                     "year_t", "surv_t", "size_t", "flw_t")]
poldmerge <- poldmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                         "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                         "year_t", "surv_t", "size_t", "flw_t")]
rmerge <- rmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                   "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                   "year_t", "surv_t", "size_t", "flw_t")]
roldmerge <- roldmerge[c("plot", "pos", "tag", "Endo", "origin", "Loc'n", "Birth Year",
                         "TRT", "Plant", "year_t1", "surv_t1", "size_t1", "flw_t1",
                         "year_t", "surv_t", "size_t", "flw_t")]

Poal <- pmerge %>% 
  rbind(poldmerge) %>% 
  rbind(roldmerge) %>% 
  rbind(rmerge) %>% 
  mutate(species = "POAL")
View(Poal)
# Pulling out the seed production estimates -------------------------------

# Pulling out the seed production estimates
pseed <- Poal_data %>% 
  select()

# What we want our final data frame to look like
LTREB_endodemog <- data.frame(colnames("quad", "species", "origin", "plot", "pos", "id", "surv_t1", "size_t1", "seed_t1", "flower_t1", "size_t", "endo", "birth", "year_t", "year_t1"))

dim(LTREB_endodemog)

