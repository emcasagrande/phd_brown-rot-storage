## libraries
{
  library(broom); library(tidyr); library(dplyr); library(ggplot2)
  library(Rmisc); library(stats)
  library(graphics)
  library(svglite)
  library(stats ); library(lme4 )
  library(yaml )
  library(deSolve)
  library(FME)
  library(minpack.lm)
  library(GA)
  library(plotly)
  library(directlabels)
  library(ggforce)
  library(popbio)
  library(gridExtra); library(ggpubr); library(plot3D); library(scatterplot3d); library(car); library(svglite)
  library(MASS); library(ggplot2); library(nlme)
}

rm(list = ls())
gc()
cat('\014')

date_1harvest19 = as.Date("19/07/22")
j10_harvest1 = as.Date("19/08/01")
date_2harvest19 = as.Date("19/07/29")
j10_harvest2 = as.Date("19/08/08")

capteurs <- read.table("capteurs_total_postharvest2019.csv", dec=".", sep=";", header=TRUE) %>%
  mutate(  date_time_gmt = as.POSIXct(as.character(date_time_gmt_plus2), format= "%d/%m/%y %H:%M:%S", tz="GMT+2"),
           date=as.Date(date), 
           sensor_HOBO = as.factor(sensor_HOBO), 
           T_factor=as.factor(T_factor) ) %>%
  filter( (date >= date_1harvest19 & date <= j10_harvest1) |
            (date >= date_2harvest19 & date <= j10_harvest2))


summary(capteurs)

unique(capteurs$date)
unique(capteurs$sensor_HOBO)

date_j0_h1_cond = as.POSIXct("2019-07-22 10:00:00")
date_j4_h1_cond = as.POSIXct("2019-07-26 10:00:00")
date_j7_h1_cond = as.POSIXct("2019-07-29 10:00:00")
date_j10_h1_cond = as.POSIXct("2019-08-01 10:00:00")

date_j0_h2_cond = as.POSIXct("2019-07-29 10:00:00")
date_j4_h2_cond = as.POSIXct("2019-08-02 10:00:00")
date_j7_h2_cond = as.POSIXct("2019-08-05 10:00:00")
date_j10_h2_cond = as.POSIXct("2019-08-08 10:00:00")

date_measure_cond_h1 = c(date_j0_h1_cond, date_j4_h1_cond, date_j7_h1_cond, date_j10_h1_cond) 
date_measure_cond_h2 = c(date_j0_h2_cond, date_j4_h2_cond, date_j7_h2_cond, date_j10_h2_cond) 



# search for specific position in the T-Ha time series
time_start = as.POSIXct("2019-07-26 10:45:00", tz="GMT+2")
time_end   = as.POSIXct("2019-07-26 18:10:00", tz="GMT+2")

# time_start = as.POSIXct("2019-08-03 00:00:00", tz="GMT+2")
# time_end   = as.POSIXct("2019-08-03 09:00:00", tz="GMT+2")

  
capteurs_specific = capteurs %>%
  mutate(T_factor = as.factor(T_factor)) %>%
  filter( date_time_gmt >= time_start & date_time_gmt <= time_end & T_factor=="15")

pl_hum =
  ggplot(capteurs_specific, aes(date_time_gmt, Ha, colour=sensor_HOBO) )+
  geom_line()+
  geom_point()+
  geom_smooth(se=F)+
  # facet_grid(~T_factor)+
  theme_bw()+
  ggtitle("Relative humidity [%]")
pl_tem =
  ggplot(capteurs_specific, aes(date_time_gmt, T_C, colour=sensor_HOBO) )+
  geom_line()+
  geom_point()+
  geom_smooth(se=F)+
  # facet_grid(~T_factor)+
  theme_bw()+
  ggtitle("Temperature [°C]")
ggarrange(pl_hum, pl_tem, nrow=2, common.legend = T, legend = "bottom")


  ggplot(capteurs_specific, aes(x = date_time_gmt, colour=sensor_HOBO) )+
  geom_line( aes(y = T_C, colour="Temperature") )+
  geom_point( aes(y = T_C, colour="Temperature") )+
  geom_line( aes(y = Ha/10, colour="RH") )+
  geom_point( aes(y = Ha/10, colour="RH") )+
  scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Relative humidity [%]"))+
  scale_colour_manual(values = c("blue", "red"))+
  facet_grid(sensor_HOBO~.)
  


plot(capteurs_specific$date_time_gmt,capteurs_specific$T_C,type="l",col="red")
par(new=T)
plot(capteurs_specific$date_time_gmt,capteurs_specific$Ha,ax=F,type="l",col="blue")
axis(4)


## ------------------- PLOTS ------------------- ##
{
# temperature
ggplot(capteurs, aes(date_time_gmt, T_C, colour=sensor_HOBO) )+
  geom_line()+
  geom_vline(xintercept = date_measure_cond_h1,  linetype="dashed", color = "red")+
  geom_vline(xintercept = date_measure_cond_h2,  linetype="dashed", color = "blue")+
  geom_smooth(se=F)+
  facet_grid(T_factor~.)+
  theme_bw()

# relative humidity
ggplot(capteurs, aes(date_time_gmt, Ha, colour=sensor_HOBO) )+
  geom_line()+
  geom_vline(xintercept = date_measure_cond_h1,  linetype="dashed", color = "red")+
  geom_vline(xintercept = date_measure_cond_h2,  linetype="dashed", color = "blue")+
  geom_smooth(se=F)+
  facet_grid(T_factor~.)+
  theme_bw()+
  ggtitle("Relative humidity [%]")


# relative humidity vs temperature
ggplot(capteurs, aes(T_C, Ha, colour=sensor_HOBO) )+
  geom_line()+
  facet_grid(T_factor~.)+
  theme_bw()+
  ggtitle("Relative humidity vs temperature")
}


## ------------------- DAILY MEAN VALUES ------------------- ##
{
  capteurs_daily = capteurs %>%
    group_by( date, T_factor ) %>% dplyr::summarize( T_C_mean = mean(T_C, na.rm=T), T_C_sd = sd(T_C, na.rm=T), Ha_mean = mean(Ha, na.rm=T), Ha_sd = sd(Ha, na.rm=T) )%>%
    ungroup()%>% droplevels()
  
  capteurs_hour = capteurs %>%
    group_by( date, time, T_factor ) %>% dplyr::summarize( T_C = mean(T_C, na.rm=T), Ha = mean(Ha, na.rm=T) )%>%
    ungroup()%>% droplevels()
  
  a=ggplot(data=capteurs_daily, aes(x=date, y=T_C_mean, colour=T_factor))+
    # geom_vline(xintercept = date_measure_cond_h1,  linetype="dashed", color = "red")+
    # geom_vline(xintercept = date_measure_cond_h2,  linetype="dashed", color = "blue")+
    geom_line()+
    geom_ribbon( aes(ymin= T_C_mean-T_C_sd , ymax=T_C_mean+T_C_sd , colour=T_factor), alpha=0.1 )+
    geom_point()+
    ylab("Temperature [°C]")+
    theme_bw(base_size = 12)+
    labs( colour = "Temperature level" )
  
  b=ggplot(capteurs_daily, aes(date, Ha_mean, colour=T_factor))+
    # geom_vline(xintercept = date_measure_cond_h1,  linetype="dashed", color = "red")+
    # geom_vline(xintercept = date_measure_cond_h2,  linetype="dashed", color = "blue")+
    geom_line()+
    geom_ribbon( aes(ymin= Ha_mean-Ha_sd , ymax=Ha_mean+Ha_sd, colour=T_factor ), alpha=0.1 )+
    geom_point()+
    ylab("Relative humidity [%]")+
    theme_bw(base_size = 12)+
    labs( colour = "Temperature level" )
  
  ggarrange(a,b, nrow=2, common.legend = T)

  
}
