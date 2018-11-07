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

pmerge <- pmerge_t1 %>% 
  full_join(pmerge_t, by = c("plot", "pos", "tag", "Endo", 
                             "Loc'n", "Birth Year", "TRT",
                             "Plant", "year_t"), all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = "O")
View(pmerge)

## Recruits data
rsurv <- Poal_data_r %>%
  rename("Birth Year" = "Date") %>% 
  melt(id.var = c("Tag", "Plot", "Endo", "Birth Year"),
       measure.var = c("Survive11", "Survive12", "Survive13", "Survive14", 
                       "Survive15", "Survive16"),
       value.name = "surv") 
rsurv$year<- ifelse(rsurv$variable == "Survive11", 2011, ifelse(rsurv$variable  == "Survive12", 2012, ifelse(rsurv$variable  == "Survive13", 2013, ifelse(rsurv$variable  == "Survive14", 2014, ifelse(rsurv$variable  == "Survive15", 2015, ifelse(rsurv$variable  == "Survive16", 2016, NA))))))
# View(rsurv)

rgrow <- Poal_data_r %>%
  rename("Birth Year" = "Date") %>% 
  melt(id.var = c("Tag", "Plot", "Endo", "Birth Year"),
       measure.var = c("TOTtiller10","TOTtiller11", "TOTtiller12", "TOTtiller13", "TOTtiller14", 
                       "TOTtiller15", "TOTtiller16"),
       value.name = "size") 
rgrow$year<- ifelse(rgrow$variable == "TOTtiller10", 2010, ifelse(rgrow$variable == "TOTtiller11", 2011, ifelse(rgrow$variable  == "TOTtiller12", 2012, ifelse(rgrow$variable  == "TOTtiller13", 2013, ifelse(rgrow$variable  == "TOTtiller14", 2014, ifelse(rgrow$variable  == "TOTtiller15", 2015, ifelse(rgrow$variable  == "TOTtiller16", 2016, NA)))))))
# View(rgrow)

rflw <- Poal_data_r %>%
  rename("Birth Year" = "Date") %>% 
  melt(id.var = c("Tag", "Plot", "Endo", "Birth Year"),
       measure.var = c("FLWtiller10","FLWtiller11", "FLWtiller12", "FLWtiller13", "FLWtiller14", 
                       "FLWtiller15", "FLWtiller16"),
       value.name = "flw") 
rflw$year<- ifelse(rflw$variable == "FLWtiller10", 2010, ifelse(rflw$variable == "FLWtiller11", 2011, ifelse(rflw$variable  == "FLWtiller12", 2012, ifelse(rflw$variable  == "FLWtiller13", 2013, ifelse(rflw$variable  == "FLWtiller14", 2014, ifelse(rflw$variable  == "FLWtiller15", 2015, ifelse(rflw$variable  == "FLWtiller16", 2016, NA)))))))
# View(rflw)

rmerge_sg <- merge(rsurv, rgrow, by = c( "Tag", "Plot", "Endo", "Birth Year", "year"))
# View(rmerge_sg)

rmerge_sgf <- merge(rmerge_sg, rflw, by = c("Tag", "Plot", "Endo", "Birth Year", "year"))
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
  full_join(rmerge_t, by = c("Tag", "Plot", "Endo", "Birth Year", "year_t"),
            all.x = all, all.y = all) %>% 
  select(-contains("variable")) %>% 
  mutate(origin = "R")
View(rmerge)














# Pulling out the seed production estimates
pseed <- Poal_data %>% 
  select()

# What we want our final data frame to look like
LTREB_endodemog <- data.frame(colnames("quad", "species", "origin", "plot", "pos", "id", "surv_t1", "size_t1", "seed_t1", "flower_t1", "size_t", "endo", "birth", "year_t", "year_t1"))

dim(LTREB_endodemog)

