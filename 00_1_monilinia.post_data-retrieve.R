## ******************************************************************* ##
## ***** DATA LOAD / PRE-PROCESSING AND STATISTICAL ANALYSIS  ******** ##
## ******************************************************************* ##


## libraries
{
library(plyr); library(dplyr)  
library(tidyr); library(ggplot2); library(broom)
library(Rmisc)
library(graphics)
library(svglite)
library(stats); library(lme4)
library(yaml)
library(deSolve)
library(FME)
library(minpack.lm)
library(GA)
library(directlabels)
library(ggforce)
library(popbio)
library(gridExtra)
library(MASS); library(wesanderson)
library(survival); library(flexsurv); library(survminer); library(mpr); library(survAUC); library(Hmisc) #library(tidysurv)
library(ranger); library(GGally)
library(treemap)
library(pec)
library(data.table)
library(splitstackshape)
#library(simPH)
library(rsample)
library(ROSE)
library(caret)
library(parallel)
library(paletteer); library(wesanderson)
library(splitstackshape)
library(multcompView)
}

rm(list = ls())
gc()
cat('\014')

source("user-defined_functions/f_get_dd.R")  # calculate degree-days 
source("user-defined_functions/f_get_cdd.R") # calculate cumulative degree-days


## ************************************************************************
# 1) Let's consider the first 2 days at 2°C ? (storage.start = 0, otherwise = 2)
storage.start = 0

# 2) I do not consider the last days of experiments in 2019 ? storage.length2019 = 14, otherwise = 21)  
storage.length2019 = 21


# 3) In ** monilinia2018 ** only Syst = S2 (ECO) is considered


# 4) caliber: strong or mild regroup?
strong.cal = T
## *************************************************************************



# full bloom date             : 17-03-2018
# 1st harvest date            : 12-07-2018
# 2nd harvest date            : 19-07-2018
# 3rd harvest date            : 27-07-2018
# thinning                    : 22-05-2018
# irrigation differentiation  : 18-06-2018

date_bloom2018_julian = as.numeric( format(as.Date("2018-03-17"),"%j") ) # julian date of full blooming day
dafb_thinning2018     = as.numeric( format(as.Date("2018-05-22"),"%j") ) - date_bloom2018_julian # thinning day (in DAFB)
dafb_irrigation2018   = as.numeric( format(as.Date("2018-06-18"),"%j") ) - date_bloom2018_julian # start of differentiation in irrigation regimes (in DAFB)
dafb_1harvest2018     = as.numeric( format(as.Date("2018-07-12"),"%j") ) - date_bloom2018_julian
dafb_2harvest2018     = as.numeric( format(as.Date("2018-07-19"),"%j") ) - date_bloom2018_julian
dafb_3harvest2018     = as.numeric( format(as.Date("2018-07-27"),"%j") ) - date_bloom2018_julian

# full bloom date             : 10-03-2019
# 1st harvest date            : 22-07-2019
# 2nd harvest date            : 29-07-2019

date_bloom2019 = as.Date("19-03-10")
date_bloom2019_julian = as.numeric( format(as.Date("2019-03-10"),"%j") )
dafb_1harvest2019     = as.numeric( format(as.Date("2019-07-22"),"%j") ) - date_bloom2019_julian
dafb_2harvest2019     = as.numeric( format(as.Date("2019-07-29"),"%j") ) - date_bloom2019_julian




## *************************************************** 
##        PRE-HARVEST ENVIRONMENTAL CONDITIONS
## ***************************************************
{
  
  environconditions_pre_2018 = read.table("data/Avignon2018_daily_INRA_STATION_84007004.csv", dec=".", sep=";", header=TRUE, skip = 13)%>%
    mutate( date_j = as.numeric( format(as.Date(paste0(year,"-",month,"-",day)),"%j") ), 
            date_dafb = date_j - date_bloom2018_julian,
            date = format(as.Date(paste0( day,'/',month,'/',substr(year,3,4) ), format = "%d/%m/%y" ), "%d/%m/%Y"), date=as.character(date)   ) %>%
    filter(date_dafb>=0)
  
  dd18 = vector()
  for (i in 1:dim(environconditions_pre_2018)[1]) { dd18[i] =  f_get_dd(Tmax = environconditions_pre_2018$t_max[i],  Tmin=environconditions_pre_2018$t_min[i]) }
  
  environconditions_pre_2018 = cbind(environconditions_pre_2018, deg_day = dd18)%>%
    mutate( cum_deg_day = cumsum(dd18) )

  
  # create different meteorological index (with different lag times) to test for the statistical models
  environconditions_pre_2018 = environconditions_pre_2018 %>%
    mutate( 
      
      l.wd2=lag(wetness_duration_h, 2),
      l.wd3=lag(wetness_duration_h, 3),
      l.wd4=lag(wetness_duration_h, 4),
      l.wd5=lag(wetness_duration_h, 5),
      l.wd6=lag(wetness_duration_h, 6),
      l.wd7=lag(wetness_duration_h, 7),
      
      wd1 = lag(wetness_duration_h, 1),
      wd2 = (wd1 + l.wd2)/2,
      wd3 = (wd1 + l.wd2 + l.wd3)/3,
      wd4 = (wd1 + l.wd2 + l.wd3 + l.wd4)/4,
      wd5 = (wd1 + l.wd2 + l.wd3 + l.wd4 + l.wd5)/5,
      wd6 = (wd1 + l.wd2 + l.wd3 + l.wd4 + l.wd5 + l.wd6)/6,
      wd7 = (wd1 + l.wd2 + l.wd3 + l.wd4 + l.wd5 + l.wd6 + l.wd7 )/7,
      
      rh1 = (lag(RH, 1)),
      rh2 = (rh1 + lag(RH, 2))/2,
      rh3 = (rh2 + lag(RH, 3))/2,
      rh4 = (rh3 + lag(RH, 4))/2,
      rh5 = (rh4 + lag(RH, 5))/2,
      rh6 = (rh5 + lag(RH, 6))/2,
      rh7 = (rh6 + lag(RH, 7))/2,   
      
      rain1 = lag(rain_mm, 1),
      # rain2 = rain1 + lag(rain_mm, 2),
      # rain3 = rain2 + lag(rain_mm, 3),
      # rain4 = rain3 + lag(rain_mm, 4),
      # rain5 = rain4 + lag(rain_mm, 5),
      # rain6 = rain5 + lag(rain_mm, 6),
      # rain7 = rain6 + lag(rain_mm, 7),
      rain2 = (rain1 + lag(rain_mm, 2))/2,
      rain3 = (rain2 + lag(rain_mm, 3))/2,
      rain4 = (rain3 + lag(rain_mm, 4))/2,
      rain5 = (rain4 + lag(rain_mm, 5))/2,
      rain6 = (rain5 + lag(rain_mm, 6))/2,
      rain7 = (rain6 + lag(rain_mm, 7))/2,      
      
      t1 = (lag(t_mean, 1)),
      t2 = (t1 + lag(t_mean, 2))/2,
      t3 = (t2 + lag(t_mean, 3))/2,
      t4 = (t3 + lag(t_mean, 4))/2,
      t5 = (t4 + lag(t_mean, 5))/2,
      t6 = (t5 + lag(t_mean, 6))/2,
      t7 = (t6 + lag(t_mean, 7))/2        )
    
  
  temp.max.list = shift(environconditions_pre_2018$t_max, n=1:7, type="lag")
  temp.min.list = shift(environconditions_pre_2018$t_min, n=1:7, type="lag")
    
  
  environconditions_pre_2018$tmin1 =  apply( t(matrix(unlist(temp.min.list[1]), ncol = dim(environconditions_pre_2018)[1],  byrow = TRUE)), 
                                            1, min, na.rm=T) 
  environconditions_pre_2018$tmin2 =  apply( t(matrix(unlist(temp.min.list[1:2]), ncol = dim(environconditions_pre_2018)[1],  byrow = TRUE)), 
                                             1, min, na.rm=T) 
  environconditions_pre_2018$tmin3 =  apply( t(matrix(unlist(temp.min.list[1:3]), ncol = dim(environconditions_pre_2018)[1],  byrow = TRUE)), 
                                             1, min, na.rm=T) 
  environconditions_pre_2018$tmin4 =  apply( t(matrix(unlist(temp.min.list[1:4]), ncol = dim(environconditions_pre_2018)[1],  byrow = TRUE)), 
                                             1, min, na.rm=T) 
  environconditions_pre_2018$tmin5 =  apply( t(matrix(unlist(temp.min.list[1:5]), ncol = dim(environconditions_pre_2018)[1],  byrow = TRUE)), 
                                             1, min, na.rm=T) 
  environconditions_pre_2018$tmin6 =  apply( t(matrix(unlist(temp.min.list[1:6]), ncol = dim(environconditions_pre_2018)[1],  byrow = TRUE)), 
                                             1, min, na.rm=T) 
  environconditions_pre_2018$tmin7 =  apply( t(matrix(unlist(temp.min.list), ncol = dim(environconditions_pre_2018)[1],  byrow = TRUE)), 
                                             1, min, na.rm=T) 
  
  environconditions_pre_2018$tmax1 =  apply( t(matrix(unlist(temp.max.list[1]), ncol = dim(environconditions_pre_2018)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2018$tmax2 =  apply( t(matrix(unlist(temp.max.list[1:2]), ncol = dim(environconditions_pre_2018)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2018$tmax3 =  apply( t(matrix(unlist(temp.max.list[1:3]), ncol = dim(environconditions_pre_2018)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2018$tmax4 =  apply( t(matrix(unlist(temp.max.list[1:4]), ncol = dim(environconditions_pre_2018)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2018$tmax5 =  apply( t(matrix(unlist(temp.max.list[1:5]), ncol = dim(environconditions_pre_2018)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2018$tmax6 =  apply( t(matrix(unlist(temp.max.list[1:6]), ncol = dim(environconditions_pre_2018)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2018$tmax7 =  apply( t(matrix(unlist(temp.max.list), ncol = dim(environconditions_pre_2018)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  
  
  
  
  ## ****** 2019
  environconditions_pre_2019 = read.table("data/Avignon2019_daily_INRA_STATION_84007004.csv", dec=".", sep=";", header=TRUE, skip = 13)%>%
    mutate( date_j = as.numeric( format(as.Date(paste0(year,"-",month,"-",day)),"%j") ),
            date_dafb = date_j - date_bloom2019_julian,
            date = format(as.Date(paste0( day,'/',month,'/',substr(year,3,4) ), format = "%d/%m/%y" ), "%d/%m/%Y"), date=as.character(date)     )%>%
    filter(date_dafb>=0)
  
  dd19 = vector()
  for (i in 1:dim(environconditions_pre_2019)[1]) { dd19[i] =  f_get_dd(Tmax = environconditions_pre_2019$t_max[i],  Tmin=environconditions_pre_2019$t_min[i]) }
  
  environconditions_pre_2019 = cbind(environconditions_pre_2019, deg_day = dd19)%>%
    mutate( cum_deg_day = cumsum(dd19) )
  
  # create different index to test for the statistical models
  environconditions_pre_2019 = environconditions_pre_2019 %>%
    mutate( 

      l.wd2=lag(wetness_duration_h, 2),
      l.wd3=lag(wetness_duration_h, 3),
      l.wd4=lag(wetness_duration_h, 4),
      l.wd5=lag(wetness_duration_h, 5),
      l.wd6=lag(wetness_duration_h, 6),
      l.wd7=lag(wetness_duration_h, 7),
      
      wd1 = lag(wetness_duration_h, 1),
      wd2 = (wd1 + l.wd2)/2,
      wd3 = (wd1 + l.wd2 + l.wd3)/3,
      wd4 = (wd1 + l.wd2 + l.wd3 + l.wd4)/4,
      wd5 = (wd1 + l.wd2 + l.wd3 + l.wd4 + l.wd5)/5,
      wd6 = (wd1 + l.wd2 + l.wd3 + l.wd4 + l.wd5 + l.wd6)/6,
      wd7 = (wd1 + l.wd2 + l.wd3 + l.wd4 + l.wd5 + l.wd6 + l.wd7 )/7,
      
      rh1 = (lag(RH, 1)),
      rh2 = (rh1 + lag(RH, 2))/2,
      rh3 = (rh2 + lag(RH, 3))/2,
      rh4 = (rh3 + lag(RH, 4))/2,
      rh5 = (rh4 + lag(RH, 5))/2,
      rh6 = (rh5 + lag(RH, 6))/2,
      rh7 = (rh6 + lag(RH, 7))/2,        
            
      rain1 = lag(rain_mm, 1),
      # rain2 = rain1 + lag(rain_mm, 2),
      # rain3 = rain2 + lag(rain_mm, 3),
      # rain4 = rain3 + lag(rain_mm, 4),
      # rain5 = rain4 + lag(rain_mm, 5),
      # rain6 = rain5 + lag(rain_mm, 6),
      # rain7 = rain6 + lag(rain_mm, 7),
      rain2 = (rain1 + lag(rain_mm, 2))/2,
      rain3 = (rain2 + lag(rain_mm, 3))/2,
      rain4 = (rain3 + lag(rain_mm, 4))/2,
      rain5 = (rain4 + lag(rain_mm, 5))/2,
      rain6 = (rain5 + lag(rain_mm, 6))/2,
      rain7 = (rain6 + lag(rain_mm, 7))/2,      
      
      
      t1 = (lag(t_mean, 1)),
      t2 = (t1 + lag(t_mean, 2))/2,
      t3 = (t2 + lag(t_mean, 3))/2,
      t4 = (t3 + lag(t_mean, 4))/2,
      t5 = (t4 + lag(t_mean, 5))/2,
      t6 = (t5 + lag(t_mean, 6))/2,
      t7 = (t6 + lag(t_mean, 7))/2
            
            )
  
  temp.max.list = shift(environconditions_pre_2019$t_max, n=1:7, type="lag")
  temp.min.list = shift(environconditions_pre_2019$t_min, n=1:7, type="lag")
  
  environconditions_pre_2019$tmin1 =   apply( t(matrix(unlist(temp.min.list[1]), ncol = dim(environconditions_pre_2019)[1],  byrow = TRUE)), 
                                             1, min, na.rm=T) 
  environconditions_pre_2019$tmin2 =  apply( t(matrix(unlist(temp.min.list[1:2]), ncol = dim(environconditions_pre_2019)[1],  byrow = TRUE)), 
                                            1, min, na.rm=T) 
  environconditions_pre_2019$tmin3 =  apply( t(matrix(unlist(temp.min.list[1:3]), ncol = dim(environconditions_pre_2019)[1],  byrow = TRUE)), 
                                            1, min, na.rm=T) 
  environconditions_pre_2019$tmin4 =  apply( t(matrix(unlist(temp.min.list[1:4]), ncol = dim(environconditions_pre_2019)[1],  byrow = TRUE)), 
                                            1, min, na.rm=T) 
  environconditions_pre_2019$tmin5 =  apply( t(matrix(unlist(temp.min.list[1:5]), ncol = dim(environconditions_pre_2019)[1],  byrow = TRUE)), 
                                            1, min, na.rm=T) 
  environconditions_pre_2019$tmin6 =  apply( t(matrix(unlist(temp.min.list[1:6]), ncol = dim(environconditions_pre_2019)[1],  byrow = TRUE)), 
                                            1, min, na.rm=T) 
  environconditions_pre_2019$tmin7 =  apply( t(matrix(unlist(temp.min.list), ncol = dim(environconditions_pre_2019)[1],  byrow = TRUE)), 
                                             1, min, na.rm=T) 
  
  environconditions_pre_2019$tmax1 =  apply( t(matrix(unlist(temp.max.list[1]), ncol = dim(environconditions_pre_2019)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2019$tmax2 =  apply( t(matrix(unlist(temp.max.list[1:2]), ncol = dim(environconditions_pre_2019)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2019$tmax3 =  apply( t(matrix(unlist(temp.max.list[1:3]), ncol = dim(environconditions_pre_2019)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2019$tmax4 =  apply( t(matrix(unlist(temp.max.list[1:4]), ncol = dim(environconditions_pre_2019)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2019$tmax5 =  apply( t(matrix(unlist(temp.max.list[1:5]), ncol = dim(environconditions_pre_2019)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2019$tmax6 =  apply( t(matrix(unlist(temp.max.list[1:6]), ncol = dim(environconditions_pre_2019)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  environconditions_pre_2019$tmax7 =  apply( t(matrix(unlist(temp.max.list), ncol = dim(environconditions_pre_2019)[1], byrow = TRUE)),
                                             1, max, na.rm=T) 
  
  rm('temp.max.list','temp.min.list')
  
  
}


## *************************************************** 
##        STORAGE CONDITIONS
## ***************************************************
{
  
  ## --------------- 2018 --------------- ##   
  environconditions_post_2018 <- read.table("data/environconditions_post_2018.csv", dec=".", sep=";", header=TRUE) %>%
    mutate( harvest_date = as.character(as.Date(harvest_date, "%Y-%m-%d" )) )
  
  
  
   ## --------------- 2019 --------------- ##

  capteurs <- read.table("data/capteurs_total_postharvest2019.csv", dec=".", sep=";", header=TRUE) %>%
    mutate(  date_time_gmt = as.POSIXct(as.character(date_time_gmt_plus2), format= "%d/%m/%y %H:%M:%S", tz="GMT+2"),
             date= as.Date(date), 
             sensor_HOBO = as.factor(sensor_HOBO), 
             T_factor=as.factor(T_factor) ) 
    # filter( (date >= date_1harvest19 & date <= j10_harvest1) )
  
  
  ## ------------------- DAILY MEAN VALUES ------------------- ##
  {
    capteurs_daily = capteurs %>%
      group_by( date, T_factor ) %>% dplyr::summarize( T_C_mean = mean(T_C, na.rm=T), T_C_sd = sd(T_C, na.rm=T), Ha_mean = mean(Ha, na.rm=T), Ha_sd = sd(Ha, na.rm=T) )%>%
      ungroup()%>% droplevels()
    
    capteurs_hour = capteurs %>%
      group_by( date, time, T_factor ) %>% dplyr::summarize( T_C = mean(T_C, na.rm=T), Ha = mean(Ha, na.rm=T) )%>%
      ungroup()%>% droplevels()
    
    
    }
  
  # start and end date of the storage period (2019 PH)
  storage_start_date = as.Date("19-07-22")
  storage_differentiation_date = as.Date("19-07-24")
  storage_end_date = as.Date("19-08-01")
  
  # set specific conditions
  clim_2C_beforedifferentiation = capteurs_daily %>% filter( date < storage_differentiation_date  & T_factor=="2")
  clim_15_afterdifferentiation = capteurs_daily %>% filter( (date <= storage_end_date & date >= storage_differentiation_date ) & T_factor=="15") 
  clim_25_afterdifferentiation = capteurs_daily %>% filter( (date <= storage_end_date & date >= storage_differentiation_date ) & T_factor=="25") 
  
  
  ## --------storage conditions (3 temperature levels) --------------------- ##
  clim_02C = capteurs_daily %>% filter( date <= storage_end_date & T_factor=="2") %>% mutate(temp = "02", T_factor=NULL, T_cumul = cumsum(T_C_mean)) 
  clim_15C = rbind(clim_2C_beforedifferentiation, clim_15_afterdifferentiation )  %>% mutate(temp = "15", T_factor=NULL, T_cumul = cumsum(T_C_mean))
  clim_25C = rbind(clim_2C_beforedifferentiation, clim_25_afterdifferentiation )  %>% mutate(temp = "25", T_factor=NULL, T_cumul = cumsum(T_C_mean))
  
  ## daily mean storage conditions for 2019 postharvest biochemical experimentations
  environconditions_post_2019 = rbind(clim_02C, clim_15C, clim_25C) %>% mutate(date=as.character(date))
  ## ----------------------------------------------------------------------- ##
  
  rm("clim_02C", "clim_15C", "clim_25C", "clim_15_afterdifferentiation", "clim_25_afterdifferentiation", "clim_2C_beforedifferentiation",
     "capteurs", "capteurs_hour")
  
  
}



## *************************************************** 
##        PREHARVEST ENVIRONMENTAL CONDITIONS - plots
## ***************************************************
{

# p.temp = ggplot()+
#     geom_line(environconditions_pre_2018, mapping = aes(date_dafb, t_mean), color='red')+
#     geom_line(environconditions_pre_2018, mapping = aes(date_dafb, t_max), color='red', alpha=0.3)+
#     geom_line(environconditions_pre_2018, mapping = aes(date_dafb, t_min), color='red', alpha=0.3)+
#     geom_line(environconditions_pre_2019, mapping = aes(date_dafb, t_mean), color='blue')+
#     geom_line(environconditions_pre_2019, mapping = aes(date_dafb, t_max), color='blue', alpha=0.3)+
#     geom_line(environconditions_pre_2019, mapping = aes(date_dafb, t_min), color='blue', alpha=0.3)+
#     geom_vline(xintercept = dafb_1harvest2018, color='red')+
#     geom_vline(xintercept = dafb_2harvest2018, color='red')+
#     geom_vline(xintercept = dafb_3harvest2018, color='red')+
#     geom_vline(xintercept = dafb_1harvest2019, color='blue')+
#     geom_vline(xintercept = dafb_2harvest2019, color='blue')+
#     theme_pubr()
# p.temp
# 
# p.rh = ggplot()+
#   ylim(c(30,100))+
#   geom_line(environconditions_pre_2018, mapping = aes(date_dafb, RH), color='red')+
#   geom_line(environconditions_pre_2019, mapping = aes(date_dafb, RH), color='blue')+
#   geom_vline(xintercept = dafb_1harvest2018, color='red')+
#   geom_vline(xintercept = dafb_2harvest2018, color='red')+
#   geom_vline(xintercept = dafb_3harvest2018, color='red')+
#   geom_vline(xintercept = dafb_1harvest2019, color='blue')+
#   geom_vline(xintercept = dafb_2harvest2019, color='blue')+
#   theme_pubr()
# p.rh
# 
# p.rain = ggplot()+
#   geom_line(environconditions_pre_2018, mapping = aes(date_dafb, rain_mm), color='red')+
#   geom_line(environconditions_pre_2019, mapping = aes(date_dafb, rain_mm), color='blue')+
#   geom_vline(xintercept = dafb_1harvest2018, color='red')+
#   geom_vline(xintercept = dafb_2harvest2018, color='red')+
#   geom_vline(xintercept = dafb_3harvest2018, color='red')+
#   geom_vline(xintercept = dafb_1harvest2019, color='blue')+
#   geom_vline(xintercept = dafb_2harvest2019, color='blue')+
#   theme_pubr()
# p.rain
# 
# 
# p.wd = ggplot()+
#   geom_line(environconditions_pre_2018, mapping = aes(date_dafb, wetness_duration_h), color='red')+
#   geom_line(environconditions_pre_2019, mapping = aes(date_dafb, wetness_duration_h), color='blue')+
#   geom_vline(xintercept = dafb_1harvest2018, color='red')+
#   geom_vline(xintercept = dafb_2harvest2018, color='red')+
#   geom_vline(xintercept = dafb_3harvest2018, color='red')+
#   geom_vline(xintercept = dafb_1harvest2019, color='blue')+
#   geom_vline(xintercept = dafb_2harvest2019, color='blue')+
#   theme_pubr()
# p.wd
# 
# ggplot()+
#     geom_point(environconditions_pre_2018, mapping = aes(RH, wetness_duration_h), color='red')+
#     geom_point(environconditions_pre_2019, mapping = aes(RH, wetness_duration_h), color='blue')+
#     theme_pubr()
# 
# 
# p.cdd = ggplot()+
#   geom_line(environconditions_pre_2018, mapping = aes(date_dafb, cum_deg_day), color='red')+
#   geom_line(environconditions_pre_2019, mapping = aes(date_dafb, cum_deg_day), color='blue')+
#   geom_vline(xintercept = dafb_1harvest2018, color='red')+
#   geom_vline(xintercept = dafb_2harvest2018, color='red')+
#   geom_vline(xintercept = dafb_3harvest2018, color='red')+
#   geom_vline(xintercept = dafb_1harvest2019, color='blue')+
#   geom_vline(xintercept = dafb_2harvest2019, color='blue')+
#   theme_pubr()
# p.cdd

}



## **************************************************************** 
##    DATAFRAMES HARVEST - Ecopêche 18-19 Avignon (Daniel Plenet)
## ***************************************************************
{
## ********** 2018 ##

harvest_caliber18 <- read.table("./data/Harvest_caliber_2018.csv", dec=".", sep=";", header=TRUE) %>%
  mutate( Harvest=as.factor(Harvest), Tree=as.factor(Tree), 
          perc_rot=round( (N_rot_cal/N_tot_cal)*100, 1 ),
          Syst=as.factor(substr(Mod,1,2)), Exp=as.factor(substr(Mod,3,4)), 
          Irr= as.factor(ifelse( Exp=="M1" | Exp=="M2", '100%',
                                 ifelse(Exp=="M3" | Exp=='M4', '50%', NA) ) ),
          FL= as.factor(ifelse( Exp=="M1" | Exp=="M3", '400',
                                ifelse(Exp=="M2" | Exp=='M4', '200', NA) ) )
  ) %>%
  group_by(Harvest, Syst, Mod, Tree) %>%
  mutate( perc_cal = round( N_tot_cal/sum(N_tot_cal)*100, 1 ), 
          perc_rot_rel = round((perc_rot/100)*(perc_cal/100),2)  # (percentual of rotten fruit for caliber i) X (number of fruits from caliber i over harvested total)
  )%>%
  ungroup()

levels(harvest_caliber18$Cal) <- c("C","B","A","2A","3A")
harvest_caliber18$perc_rot_rel[is.nan(harvest_caliber18$perc_rot_rel)]=NA
harvest_caliber18$perc_rot[is.nan(harvest_caliber18$perc_rot)]=NA

harvest_rot18 = harvest_caliber18 %>% group_by(Syst, Harvest) %>% 
  dplyr::summarize( n_tot = sum(N_sain_cal,na.rm = T), n_rot = sum(N_rot_cal , na.rm = T)  )  %>% mutate( perc_rot = n_rot/n_tot*100  ) %>% droplevels()


## ********** 2019 ##
harvest_caliber19 <- read.table("./data/Harvest_caliber_2019.csv", dec=".", sep=";", header=TRUE) 

harvest_rot19 = harvest_caliber19 %>% group_by(Syst, Harvest) %>%
  dplyr::summarize( n_tot = sum(N_sain_cal,na.rm = T), n_rot = sum(N_rot_cal , na.rm = T)  )  %>% mutate( perc_rot = n_rot/n_tot*100  ) %>% droplevels()
}



## *************************************************** 
##     LOAD BROWN ROT DATA
## ***************************************************  
{

 ## ************************************************************** ##
 ##                             2018
 ## ************************************************************** ##
 {
   
  monilinia18 = read.table("data/nectarlove_monilinia_post_2018.csv", dec=".", sep=";", header=TRUE) %>%
    mutate( 
      box = as.factor(box), n_fruit=as.factor(n_fruit), dah = as.factor(dah), harvest_n=as.factor(harvest_n),
           
      harvest_date = as.character(harvest_date),
           
      id.fruit = paste0(harvest_n, ".", mod, ".", box, ".", n_fruit),
           
      Syst = as.factor(substr(mod,1,2)),
      Exp = as.factor(substr(mod,3,4)),
      Irr= as.factor(ifelse( Exp=="M1" | Exp=="M2", '100%',
                             ifelse(Exp=="M3" | Exp=='M4', '50%', NA) ) ),
      FL= as.factor(ifelse( Exp=="M1" | Exp=="M3", '400',
                            ifelse(Exp=="M2" | Exp=='M4', '200', NA) ) )
      
      )

  # length( unique(monilinia18[monilinia18$harvest_n=="1" & monilinia18$Syst=="S2",]$id.fruit ) )
  # length( unique(monilinia18[monilinia18$harvest_n=="2" & monilinia18$Syst=="S2",]$id.fruit ) )
  # length( unique(monilinia18[monilinia18$harvest_n=="3" & monilinia18$Syst=="S2",]$id.fruit ) )
    
  monilinia18 <- monilinia18[monilinia18$infection!="Absent ",]
  monilinia18 <- monilinia18[monilinia18$infection!="Absent",]
  
  # length( unique(monilinia18[monilinia18$harvest_n=="1" & monilinia18$Syst=="S2",]$id.fruit ) )
  # length( unique(monilinia18[monilinia18$harvest_n=="2" & monilinia18$Syst=="S2",]$id.fruit ) )
  # length( unique(monilinia18[monilinia18$harvest_n=="3" & monilinia18$Syst=="S2",]$id.fruit ) )
  

  ## ***************************** 
  ## TABLE FRUIT SIZE   
  poidsmon = read.table("./data/PH_PoidsMon_2018.csv", dec=".", sep=";", header=TRUE) %>%
      mutate(
        Box=as.factor(Box), Tree=as.factor(Tree), Fruit=as.factor(Fruit), Harvest=as.factor(Harvest),
        # Weight=as.numeric(as.character(Weight)),
        id.fruit = paste0(Harvest, ".", Mod, ".", Box, ".", Fruit),
        # Syst=substr(Mod,1,2),
        H_Date=NULL)
  
    
  colnames(poidsmon)[colnames(poidsmon)=="X0" ] <- 0
  colnames(poidsmon)[colnames(poidsmon)=="X4" ] <- 4
  colnames(poidsmon)[colnames(poidsmon)=="X6" ] <- 6
  colnames(poidsmon)[colnames(poidsmon)=="X12"] <- 12
  
    
    
    ## *****************************
    ## TABLE MERGE BROWN ROT - FRUIT SIZE
    monilinia18 <- merge(monilinia18, poidsmon, by.x ="id.fruit", by.y ="id.fruit"  ) %>%
      mutate(
        weight_g = as.numeric(Weight), Weight = NULL, H_Date=NULL,
        
        dah.int = as.numeric(as.character(dah)),
        temp = 20,
        
        caliber = ifelse( weight_g >= 65 &  weight_g < 85,   "D",
                          ifelse( weight_g >= 85 &  weight_g < 105,  "C",  
                                  ifelse( weight_g >= 105 &  weight_g < 135, "B",         
                                          ifelse( weight_g >= 135 &  weight_g < 180, "A",
                                                  ifelse( weight_g >= 180 &  weight_g < 220, "2A",
                                                          ifelse( weight_g >= 220 &  weight_g < 300, "3A",
                                                                  ifelse( weight_g >= 300 ,                  "4A",NA  ))))))),
        caliber = as.factor(caliber),
        # caliber = car::recode(caliber, "c('B','C','D') = 'D+C+B' ; c('2A','3A','4A') = '2A+3A+4A'  "), 
        caliber = car::recode(caliber, "c('C','D') = 'D+C' ; c('3A','4A') = '3A+4A'  ")

        
      ) %>%
      droplevels()
    
    monilinia18 = monilinia18 %>% filter( !is.na(caliber)  ) %>% droplevels()
    # monilinia18$caliber = factor(monilinia18$caliber, levels=c("D+C+B","A", "2A+3A+4A") )
    monilinia18$caliber = factor(monilinia18$caliber, levels=c("D+C","B","A", "2A", "3A+4A") ) 
    
    rm('poidsmon')
    
    
    monilinia18 = 
      merge(monilinia18, environconditions_pre_2018, by.x = "harvest_date", by.y="date") %>% droplevels() %>%
      merge( harvest_rot18, by.x = c("harvest_n","Syst"), by.y = c("Harvest","Syst") ) %>% droplevels() %>%
      # filter( dah.int >= storage.start ) %>% mutate( dah.int = dah.int - storage.start  )
      filter( dah.int >= storage.start & Syst == "S2") %>% 
      mutate( dah.int = dah.int - storage.start  )
    

    # length( unique(monilinia18[monilinia18$harvest_n=="1" & monilinia18$Syst=="S2",]$id.fruit ) )
    # length( unique(monilinia18[monilinia18$harvest_n=="2" & monilinia18$Syst=="S2",]$id.fruit ) )
    # length( unique(monilinia18[monilinia18$harvest_n=="3" & monilinia18$Syst=="S2",]$id.fruit ) )
    
    
  ## *****************************
  ##  TABLE SURVIVAL ANALYSIS

  un_id = unique(monilinia18$id.fruit)
  
  mon18 = monilinia18 %>% dplyr::select(c('id.fruit','harvest_n', 'temp', 'weight_g', 'year',
                                          'cum_deg_day', 'RH', "t_mean", 'date_dafb', 'harvest_date', 'perc_rot', 'Irr', 'FL',
                                          
                                          'rain_mm', 'rain1', 'rain2', 'rain3', 'rain4', 'rain5', 'rain6', 'rain7', 
                                          't1', 't2','t3','t4', 't5', 't6', 't7', 
                                          'rh1', 'rh2', 'rh3', 'rh4', 'rh5', 'rh6', 'rh7', 
                                          'tmin1', 'tmin2', 'tmin3', 'tmin4', 'tmin5', 'tmin6', 'tmin7',
                                          'tmax1', 'tmax2', 'tmax3', 'tmax4', 'tmax5', 'tmax6', 'tmax7',
                                          'wetness_duration_h', 'wd1', 'wd2', 'wd3', 'wd4', 'wd5', 'wd6', 'wd7')) %>% distinct()
  time_status = c()
  
  for ( i in 1:length(unique(monilinia18$id.fruit)) ){
      
      m = monilinia18[ monilinia18$id.fruit == un_id[i], ]   
      m = m[order(m$dah.int),]
      mm = min(which(m$infection==1))
      
      if (mm<Inf){
        tim = m[mm,'dah.int']
        status = 1
      } else {
        tim = max(m$dah.int)
        status = 0
      }
      time_status = rbind(time_status, c('id.fruit' = un_id[i],
                                         'time' = tim,
                                         'status' = status))
    }
  
  mon18 = merge(mon18, time_status, by='id.fruit')
  mon18$time = as.numeric(mon18$time)
  mon18$status = as.numeric(mon18$status)
  mon18$id.fruit = as.factor(mon18$id.fruit)
  
  mon.survival18 = mon18
  
  length( unique(mon.survival18[mon.survival18$harvest_n=="1",]$id.fruit ) )
  length( unique(mon.survival18[mon.survival18$harvest_n=="2",]$id.fruit ) )
  length( unique(mon.survival18[mon.survival18$harvest_n=="3",]$id.fruit ) )
  
  length( unique(mon.survival18$id.fruit ) )
  
  
  rm("mon18", "m")
  
  
  ## *****************************
  ## PH model calibration DB
  
  monilinia_stat.caliber18 = monilinia18 %>%

    group_by(harvest_n, dah, caliber, temp, year, harvest_date) %>%

    dplyr::summarize( mean_w = mean( weight_g, na.rm = T ),
                      # incidence.test = sum( infection=="1" | infection=="2", na.rm = T ),
                      n_tot  =  length(infection),
                      incidence = sum( infection=="1" | infection=="2", na.rm = T ) /( length(infection)),
                      n_I = sum( infection=="1" | infection=="2", na.rm = T ),

                      wd4   = mean(wd4,na.rm=T),
                      temp4 = mean(t4) ,
                      rot = mean(perc_rot),
                      cum_deg_day = mean(cum_deg_day, na.rm=T)

    ) %>%

    ungroup() %>%
    mutate (  dah_int = as.numeric(as.character(dah)),
              temp = as.factor(temp) ) %>% droplevels()
    
  

    

  # ggplot(monilinia_stat.caliber18)+
  #   scale_colour_brewer(palette = "Set1")+
  #   # scale_colour_brewer(palette = wes_palette(name="Zissou1", type = "continuous"))+
  #   geom_line(mapping = aes(dah_int, n_tot - n_I), size=1.5) +
  #   geom_line(mapping = aes(dah_int, n_tot), size = 1.2)+
  #   ylab('Fruit number')+
  #   xlab('Days after harvest')+
  #   facet_grid(as.factor(cum_deg_day)~caliber)+
  #   theme_bw(base_size = 18)+
  #   theme(legend.position="bottom")+
  #   guides(shape = guide_legend(title = "Storage temperature"),
  #          color = guide_legend(title = "Storage temperature") )
  
 }
  
  
  
  ## ************************************************************** ##
  ##                             2019
  ## ************************************************************** ##
 {
   
  monilinia19 = read.table("./data/nectarlove_monilinia_post_2019.csv", dec=".", sep=";", header=TRUE) %>%
  mutate(temp=as.factor(temp), irr=as.factor(irr), fruit_n=as.factor(fruit_n), 
         dah.int = dah, dah = as.factor(dah),
         harvest_n = as.factor(harvest_n), 
         weight_g = as.numeric(as.character(weight)), weight = NULL, 
         mod = as.factor(paste0(temp,"_",irr)),
         
         # fruit unique id
         id.fruit = paste0(harvest_n, ".", temp, ".", irr, ".", box, ".", fruit_n),
         
         caliber = ifelse( weight_g >= 65 &  weight_g < 85,   "D",
                           ifelse( weight_g >= 85 &  weight_g < 105,  "C",  
                                   ifelse( weight_g >= 105 &  weight_g < 135, "B",         
                                           ifelse( weight_g >= 135 &  weight_g < 180, "A",
                                                   ifelse( weight_g >= 180 &  weight_g < 220, "2A",
                                                           ifelse( weight_g >= 220 &  weight_g < 300, "3A",
                                                                   ifelse( weight_g >= 300 ,                  "4A",NA  ))))))),
         caliber = as.factor(caliber),
         # caliber = car::recode(caliber, "c('B','C','D') = 'D+C+B' ; c('2A','3A','4A') = '2A+3A+4A'  "), 
         caliber = car::recode(caliber, "c('C','D') = 'D+C' ; c('3A','4A') = '3A+4A'  ") 
         
         )

 monilinia19 = monilinia19 %>% filter( !is.na(caliber)  ) %>% droplevels()
 # monilinia19$caliber  = factor(monilinia19$caliber, levels=c("D+C+B","A", "2A+3A+4A") )
 monilinia19$caliber  = factor(monilinia19$caliber, levels=c("D+C","B","A", "2A", "3A+4A") ) 
 
 
 length( unique(monilinia19[monilinia19$harvest_n=="1",]$id.fruit ) )
 length( unique(monilinia19[monilinia19$harvest_n=="2",]$id.fruit ) )
 
 length( unique(monilinia19$id.fruit ) ) 
 
  monilinia19 = merge(monilinia19, environconditions_pre_2019, by.x = "harvest_date", by.y="date") %>%
    merge( harvest_rot19, by.x = "harvest_n", by.y = "Harvest" ) %>%
    filter( dah.int >= storage.start & dah.int <= storage.length2019 ) %>% droplevels() %>%
    mutate( dah.int = dah.int - storage.start )
  

  ## *****************************
  ##  dataframe for survival analysis
   
   un_id = unique(monilinia19$id.fruit)
   
   mon19 = monilinia19 %>% dplyr::select(c('id.fruit', 'harvest_n', 'temp','weight_g', 'year',
                                           'cum_deg_day', 'RH', "t_mean", 'date_dafb', 'harvest_date', 'perc_rot', 'irr',
                                           
                                           'rain_mm', 'rain1', 'rain2', 'rain3', 'rain4', 'rain5', 'rain6', 'rain7', 
                                           't_mean', 't1', 't2','t3','t4', 't5', 't6', 't7', 
                                           'rh1', 'rh2', 'rh3', 'rh4', 'rh5', 'rh6', 'rh7',
                                           'tmin1', 'tmin2', 'tmin3', 'tmin4', 'tmin5', 'tmin6', 'tmin7',
                                           'tmax1', 'tmax2', 'tmax3', 'tmax4', 'tmax5', 'tmax6', 'tmax7',
                                           'wetness_duration_h', 'wd1', 'wd2', 'wd3', 'wd4', 'wd5', 'wd6', 'wd7')) %>% distinct() 
   time_status = c()
   
   for ( i in 1:length(unique(monilinia19$id.fruit)) ){
    
    m = monilinia19[ monilinia19$id.fruit == un_id[i], ]   
    m = m[order(m$dah.int),]
    mm = min(which(m$infection==1))
    
    if (mm<Inf){
      tim = m[mm,'dah.int']
      status = 1
    } else {
      tim = max(m$dah.int)
      status = 0
    }
    time_status = rbind(time_status, c('id.fruit' = un_id[i],
                                       'time' = tim,
                                       'status' = status))
  }
   
   mon19 = merge(mon19, time_status, by='id.fruit')
   mon19$time =   as.numeric(mon19$time)
    mon19$status = as.numeric(mon19$status)
   mon19$id.fruit = as.factor(mon19$id.fruit)
   mon.survival19 = mon19%>%
     mutate( temp = as.numeric(as.character(substr(temp,2,3))) )
   
   rm("mon19","m")
  
   
   length( unique(mon.survival19[mon.survival19$harvest_n=="1",]$id.fruit ) )
   length( unique(mon.survival19[mon.survival19$harvest_n=="2",]$id.fruit ) )
   
   ## *****************************
   ## PH model calibration DB
   
   monilinia_stat.caliber19 = monilinia19 %>%
     group_by(harvest_n, dah, caliber, temp, year, harvest_date) %>%
     dplyr::summarize( 
       
       mean_w = mean( weight_g, na.rm = T ), 
       # incidence.test = sum( infection=="1" | infection=="2", na.rm = T ),
       n_tot  =  length(infection),
       incidence = sum( infection=="1" | infection=="2", na.rm = T ) /( length(infection)),
       n_I = sum( infection=="1" | infection=="2", na.rm = T ),
       
       wd4   = mean(wd4,na.rm=T),
       temp4 = mean(t4) ,
       rot = mean(perc_rot),
       cum_deg_day = mean(cum_deg_day, na.rm=T)
       
     ) %>%
     ungroup() %>%
     mutate (  dah_int = as.numeric(as.character(dah)),
               temp    = as.factor(substr(temp,2,3))        ) %>% droplevels()
   

   # ggplot(monilinia_stat.caliber19)+
   #   scale_colour_brewer(palette = "Set1")+
   #   # scale_colour_brewer(palette = wes_palette(name="Zissou1", type = "continuous"))+
   #   geom_line(mapping = aes(dah_int, n_tot - n_I, colour=temp), size=1.5) +
   #   geom_line(mapping = aes(dah_int, n_tot, colour=temp), size = 1.2)+
   #   ylab('Fruit number')+
   #   xlab('Days after harvest')+
   #   facet_grid(as.factor(cum_deg_day)~caliber)+
   #   theme_bw(base_size = 18)+
   #   theme(legend.position="bottom")+
   #   guides(shape = guide_legend(title = "Storage temperature"),
   #          color = guide_legend(title = "Storage temperature") )
   
   
   

  
  

  
  ## *****************************
  ## PH model calibration DB
  
  # monilinia_stat.caliber = rbind(monilinia_stat.caliber19, monilinia_stat.caliber18) %>%
  #   mutate (  dah_int = as.numeric(as.character(dah)) )
  # # mutate(temp  = as.numeric(as.character(temp)) )
 }
  
  
  
  }


## ************************************************************** ##
##     MONILINIA SURVIVAL DATAFRAME   //// 2018 + 2019 ////
## ************************************************************** ##
{
  
  enviropost19 = environconditions_post_2019 %>%
    group_by(temp)%>%dplyr::summarize(Ha = mean(Ha_mean)/100)%>%ungroup() %>%
    mutate(temp = as.numeric(as.character(temp)))
  
  mon.survival19 = merge(mon.survival19, enviropost19)  
  
  mon.survival18$Ha=NA
  
  
# mon.survival = rbind(mon.survival18, mon.survival19) %>% 
mon.survival = rbind(mon.survival18%>%dplyr::select(-Irr,-FL), 
                     mon.survival19%>%dplyr::select(-irr) ) %>% 
  mutate( 
    caliber = ifelse( weight_g >= 65 &  weight_g < 85,   "D",
                      ifelse( weight_g >= 85 &  weight_g < 105,  "C",  
                              ifelse( weight_g >= 105 &  weight_g < 135, "B",         
                                      ifelse( weight_g >= 135 &  weight_g < 180, "A",
                                              ifelse( weight_g >= 180 &  weight_g < 220, "2A",
                                                      ifelse( weight_g >= 220 &  weight_g < 300, "3A",
                                                              ifelse( weight_g >= 300 ,                  "4A",NA  ))))))),
    caliber = as.factor(caliber),
    # caliber = car::recode(caliber, "c('B','C','D') = 'D+C+B' ; c('2A','3A','4A') = '2A+3A+4A'  "), 
    caliber = car::recode(caliber, "c('C','D') = 'D+C' ; c('3A','4A') = '3A+4A'  ") , 
    
    ## interaction term 'mean temperature before harvest' x 'cumulative wetness duration' on the # days before harvest
    wd_x_t = wetness_duration_h*t_mean,
    wd_x_t2 = wd2*t2,
    wd_x_t3 = wd3*t3,
    wd_x_t4 = wd4*t4,
    wd_x_t7 = wd7*t7,
    
    temp.factor = as.factor(temp),
    cum_deg_day.factor = as.factor(cum_deg_day)
    
  ) %>% droplevels() 
# filter( temp != "2" ) %>% droplevels()

# mon.survival$caliber = factor(mon.survival$caliber, levels=c("D+C+B","A", "2A+3A+4A") )
mon.survival$caliber = factor(mon.survival$caliber, levels=c("D+C","B","A", "2A", "3A+4A") ) 


mon.survival$rep = mon.survival$time-1

mon.survival$wd7=round(mon.survival$wd7,2)


mon.survival$tempK = mon.survival$temp + 273.15
mon.survival$tempK_1 = (1/mon.survival$tempK)
mon.survival$tempK_1_vh = (1/293.15) - (1/mon.survival$tempK)


mon.survival$weight_g.scale = scale(mon.survival$weight_g)
mon.survival$rh7.scale      = scale(mon.survival$rh7)
mon.survival$t7.scale       = scale(mon.survival$t7)
mon.survival$tmin7.scale    = scale(mon.survival$tmin7)
mon.survival$tmax7.scale    = scale(mon.survival$tmax7)
mon.survival$wd7.scale      = scale(mon.survival$wd7)


wg.mu.scale = attributes(mon.survival$weight_g.scale)[[2]]
wg.si.scale = attributes(mon.survival$weight_g.scale)[[3]]

mon.survival$weight_g.scale = round(as.numeric( mon.survival$weight_g.scale ),2)


mon.survival$cum_deg_day.scale = scale(mon.survival$cum_deg_day)
GDD.mu.scale = attributes(mon.survival$cum_deg_day.scale)[[2]]
GDD.si.scale = attributes(mon.survival$cum_deg_day.scale)[[3]]

mon.survival$cum_deg_day.scale = round(as.numeric( mon.survival$cum_deg_day.scale ),2)
}

# ## verify unique values ##
# unique(mon.survival$cum_deg_day.factor)
# unique(mon.survival$wd1)
# unique(mon.survival$wd2)
# unique(mon.survival$wd3)
# unique(mon.survival$wd4)
# unique(mon.survival$wd7)
# unique(mon.survival$perc_rot)


mon.survival18$weight_g.scale =    (mon.survival18$weight_g-wg.mu.scale) / wg.si.scale
mon.survival18$cum_deg_day.scale = (mon.survival18$cum_deg_day-GDD.mu.scale) / GDD.si.scale

mon.survival19$weight_g.scale =    (mon.survival19$weight_g-wg.mu.scale) / wg.si.scale
mon.survival19$cum_deg_day.scale = (mon.survival19$cum_deg_day-GDD.mu.scale) / GDD.si.scale
 
rm("enviropost19")



## *****************************************
##   Nectarine weight - ANOVA and Tukey test
## *****************************************

generate_label_df <- function(TUKEY, variable){

  # TUKEY = TukeyWeight18
  # variable = "id.hort"

  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])

  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}


mon18Weight = monilinia18 %>% group_by(harvest_n, box, tree, Irr, FL, n_fruit) %>%
  dplyr::summarise( weight_g = mean(weight_g, na.rm=T) ) %>% mutate( id.hort = as.factor(paste0(harvest_n,"_", Irr,"_", FL) ) ) %>% ungroup()


summarySE(mon18Weight, measurevar = "weight_g", groupvars = "id.hort")

aovWeight18 = aov(weight_g ~ id.hort, data = mon18Weight)
TukeyWeight18 = TukeyHSD(aovWeight18)
labelTukeyWeight18 <- generate_label_df(TukeyWeight18 , "id.hort"); labelTukeyWeight18




mon19Weight = monilinia19 %>% group_by(harvest_n, box, fruit_n, irr, temp) %>%
  dplyr::summarise( weight_g = mean(weight_g, na.rm=T) ) %>% mutate( id.hort = as.factor(paste0(harvest_n,"_", irr,"_", temp) ) ) %>% ungroup()

summarySE(mon19Weight, measurevar = "weight_g", groupvars = "id.hort")

aovWeight19 = aov(weight_g ~ id.hort, data = mon19Weight)
TukeyWeight19 = TukeyHSD(aovWeight19)
labelTukeyWeight19 <- generate_label_df(TukeyWeight19 , "id.hort"); labelTukeyWeight19


## *****************************************************************************************************************
##   Cumulative prevalence at the end of storage - Kruskal-Wallis and Conover-Iman test (ad-hoc multiple comparison)
## *****************************************************************************************************************

library(FSA)

mon18.prevalence = monilinia18 %>% 
  filter(Syst == "S2" & harvest_date != "3") %>%
  group_by(harvest_n, Irr, FL, box, dah, Exp) %>%
  dplyr::summarize( incidence = round(sum( infection=="1" )*100 /( length(infection)),2) ) %>% ungroup() %>%
  mutate( exp_harvest = paste0(Exp,"_",harvest_n) )

mon19.prevalence = monilinia19 %>% 
  group_by(harvest_n, irr, temp, box, dah, mod) %>%
  dplyr::summarize( incidence = sum(  infection=="1" | infection=="2", na.rm = T )*100 /( length(infection)), 
                    weight_mean = mean(weight_g, na.rm=T), weight_sd = sd(weight_g, na.rm=T) )  # cumulated incidence 



## 2018 - 1 harvest
mon18.prevalence.1h = mon18.prevalence %>% filter( dah=="14"  & harvest_n=="1")
kw_18.1h = kruskal.test( incidence ~ Exp, data = mon18.prevalence.1h)
conover.test_18.1h = conover.test::conover.test(mon18.prevalence.1h$incidence, g=mon18.prevalence.1h$Exp, method='bh', kw=F, label = T)
rcompanion::cldList(P.adjusted ~ comparison,
                    data      = conover.test_18.1h,
                    threshold = 0.05)

## 2018 - 2 harvest
mon18.prevalence.2h = mon18.prevalence %>% filter( dah=="14"  & harvest_n=="2")
kw_18.2h = kruskal.test( incidence ~ Exp, data = mon18.prevalence.2h)
conover.test_18.2h = conover.test::conover.test(mon18.prevalence.2h$incidence, g=mon18.prevalence.2h$Exp, method='bh', kw=F, label = T)
rcompanion::cldList(P.adjusted ~ comparison,
                    data      = conover.test_18.2h,
                    threshold = 0.05)

## 2019 - 1 harvest
mon19.prevalence.1h = mon19.prevalence %>% filter( dah=="17"  & harvest_n=="1")
kw_19.1h = kruskal.test( incidence ~ mod, data = mon19.prevalence.1h)
conover.test_19.1h = conover.test::conover.test(mon19.prevalence.1h$incidence, g=mon19.prevalence.1h$mod, method='bh', kw=F, label = T)
rcompanion::cldList(P.adjusted ~ comparison,
                    data      = conover.test_19.1h,
                    threshold = 0.05)

## 2019 - 2 harvest
mon19.prevalence.2h = mon19.prevalence %>% filter( dah=="21"  & harvest_n=="2")
kw_19.2h = kruskal.test( incidence ~ mod, data = mon19.prevalence.2h)
conover.test_19.2h = conover.test::conover.test(mon19.prevalence.2h$incidence, g=mon19.prevalence.2h$mod, method='bh', kw=F, label = T)
rcompanion::cldList(P.adjusted ~ comparison,
                    data      = conover.test_19.2h,
                    threshold = 0.05)




