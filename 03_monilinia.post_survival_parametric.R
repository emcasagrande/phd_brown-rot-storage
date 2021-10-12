## ******************************************************************* ##
## ******************* MONILINIA SURVIVAL ANALYSIS ******************* ##
## ******************************************************************* ##

# in this version of the script I made the hp. 
# *** System REF (S1) 2018 was not taken in consideration (fongicide)

## **************
## load observed data + plots (check filter informations before run analysis) 
source("00_monilinia.post_data-retrieve.R")
## **************

# plot monilinia together with parametric
mon18.db.plot = monilinia18 %>%
  group_by(harvest_n, harvest_date, dah, Irr, FL, Syst, Exp, dah.int) %>%
  dplyr::summarize( incidence = round(sum( infection=="1" )*100 /( length(infection)),2), 
                    weight_mean = round(mean(weight_g, na.rm=T),2),  weight_sd = sd(weight_g, na.rm=T)) %>% # cumulated incidence 
  ungroup() %>%
  group_by(harvest_n, harvest_date, Irr, FL, Syst, Exp) %>%
  mutate(max.dah = max(dah.int) ) %>% ungroup() %>%
  filter(Syst == "S2")

mon19.db.plot = monilinia19 %>%
  group_by(harvest_n, dah, temp, irr, mod, dah.int) %>%
  dplyr::summarize( incidence = sum(  infection=="1" | infection=="2", na.rm = T )*100 /( length(infection)), 
                    weight_mean = mean(weight_g, na.rm=T), weight_sd = sd(weight_g, na.rm=T) ) %>% # cumulated incidence 
  ungroup()%>%
  group_by(harvest_n, irr, temp, mod) %>%
  mutate(max.dah = max(dah.int),
         temp.int = as.numeric(as.character(substr(temp,2,3))), temp=as.factor(temp.int)  
         ) %>% ungroup() 




### **** https://link.springer.com/article/10.1007/s11157-017-9443-0

# CTMI model - Lobry 1991

# Tmin = -2
# Tmax = 38
# Topt = 23

Tmin = -9.6 # calibrated
Tmax = 35.4
Topt = 25.4

mon.survival$temp.bell = ( (mon.survival$temp - Tmax) * ((mon.survival$temp-Tmin)^2) ) /
                         ( (Topt-Tmin) * ( ((Topt-Tmin)*(mon.survival$temp-Topt)) - ((Topt-Tmax)*(Topt+Tmin-2*mon.survival$temp)) ) ) # CTMI

mon.survival$temp.bell = ((Tmax-mon.survival$temp)/(Tmax - Topt))*(((mon.survival$temp-Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt))) # Yan
mon.survival$temp.bell=round(mon.survival$temp.bell,2)

unique(mon.survival$temp.bell)



mon.survival$rain7.bool = as.factor(1*(mon.survival$rain7>0))


## ******************************************************************* ##
##                PARAMETRIC REGRESSION 18-19 
## ******************************************************************* ##
fit_gomp.mpr.step = mpr( Surv(time, status) ~ list(~ 1,
                                                   ~ 1),
                         data = mon.survival,
                         # data = mon.survival[mon.survival$temp!=20,],
                         family = "Gompertz"); summary(fit_gomp.mpr.step)


step.mpt_gomp = stepmpr(fit_gomp.mpr.step, 
                        # scope = list(lower = ~ 1, upper = ~ wd7 + weight_g.scale + temp.bell + rain7 + perc_rot),
                        scope = list(lower = ~ 1, upper = ~ wd7 + weight_g.scale + temp + rain7.bool + perc_rot  + tmin7 + tmax7 + t7 + cum_deg_day.scale ),
                        comp = 1,
                        # comp =1:fit_gomp.mpr.step$ncomp,
                        direction = "both",
                        joint = T,
                        aic=F
                        ); summary(step.mpt_gomp)

count_combinations = sum(choose(n=9,k=1:9))


# # **************** best model (temp) with covariate-dependent scale and location
# fit_gomp.bestmpr = flexsurvreg( Surv(time, status) ~ temp + weight_g.scale + wd7 + perc_rot + t7,
#                                 data = mon.survival,
#                                 anc = list(shape = ~ t7 ),
#                                 dist = "gompertz"); fit_gomp.bestmpr


# # best model (temp) with covariate-dependent scale (WEIBULL AFT)
fit_weib.mpr = mpr( Surv(time, status) ~ list(~ temp + perc_rot + weight_g.scale + wd7 ,
                                              ~ 1),
                    data = mon.survival,
                    family = "WeibullAFT"); summary(fit_weib.mpr)
fit_weib.bestmpr = flexsurvreg( Surv(time, status) ~ temp + perc_rot + weight_g.scale + wd7 + rain7 ,
                                data = mon.survival,
                                # anc = list(shape = ~ 1 ),
                                dist = "weibull"); fit_weib.bestmpr



# # best model (temp) with covariate-dependent scale (GOMPERTZ)
fit_gomp.mpr = mpr( Surv(time, status) ~ list(~ temp + weight_g.scale + wd7 + perc_rot  ,
                                              ~ 1),
                    data = mon.survival,
                    family = "Gompertz"); summary(fit_gomp.mpr)
fit_gomp.bestmpr = flexsurvreg( Surv(time, status) ~ temp + weight_g.scale + wd7 + perc_rot ,
                                data = mon.survival,
                                # anc = list(shape = ~ 1 ),
                                dist = "gompertz"); fit_gomp.bestmpr
fit_gomp.bestmpr = flexsurvreg( Surv(time, status) ~ temp + weight_g + wd7 + perc_rot   ,
                                data = mon.survival,
                                #anc = list(shape = ~ 1 ),
                                dist = "gompertz"); fit_gomp.bestmpr

results_fitgomp = as.data.frame(fit_gomp.bestmpr$res[3:6,])
results_fitgomp$est * fit_gomp.bestmpr$datameans
tidy(fit_gomp.bestmpr, transform = "coefs.exp")


# summary.mpr(fit_gomp.mpr, newdata = new.data, type = "percentile")

# fit_gomp.bestmpr$coefficients["temp"]     * mean(mon.survival$temp) 
# fit_gomp.bestmpr$coefficients["weight_g"] * mean(mon.survival$weight_g)
# fit_gomp.bestmpr$coefficients["wd7"]      * mean(mon.survival$wd7)
# fit_gomp.bestmpr$coefficients["perc_rot"] * mean(mon.survival$perc_rot)
# c(fit_gomp.bestmpr$coefficients["temp"]     * min(mon.survival$temp),     fit_gomp.bestmpr$coefficients["temp"]     * max(mon.survival$temp))
# c(fit_gomp.bestmpr$coefficients["weight_g"] * min(mon.survival$weight_g), fit_gomp.bestmpr$coefficients["weight_g"] * max(mon.survival$weight_g))
# c(fit_gomp.bestmpr$coefficients["wd7"]      * min(mon.survival$wd7),      fit_gomp.bestmpr$coefficients["wd7"]      * max(mon.survival$wd7))
# c(fit_gomp.bestmpr$coefficients["perc_rot"] * min(mon.survival$perc_rot), fit_gomp.bestmpr$coefficients["perc_rot"] * max(mon.survival$perc_rot))

### ************************************************************************* ##  
###                         TEST PARAMETERS IMPORTANCE
### ************************************************************************* ##
{

  ## WEIGHT
  new.data =  mon.survival %>% ## plot solo valori medi weight
    mutate( temp=mean(temp),  wd7=mean(wd7), perc_rot=mean(perc_rot) ) %>% distinct( weight_g.scale , .keep_all = T )

  length(unique(mon.survival$weight_g.scale))

  surv_gomp = summary(model.surv.g, newdata = new.data, type = "median", tidy = TRUE )

  ggplot(surv_gomp, aes(weight_g.scale, est) )+
    geom_point()


  ## TEMPERATURE
  new.data =  mon.survival %>% ## plot solo valori medi weight
    mutate( weight_g.scale=mean(weight_g.scale),  wd7=mean(wd7), perc_rot=mean(perc_rot) ) %>% distinct(  temp, .keep_all = T )

  length(unique(mon.survival$temp))

  surv_gomp = summary(model.surv.g, newdata = new.data, type = "median", tidy = TRUE )

  ggplot(surv_gomp, aes(temp, est) )+
    geom_point()


  ## PERC ROT
  new.data =  mon.survival %>% ## plot solo valori medi weight
    mutate( weight_g.scale=mean(weight_g.scale),  wd7=mean(wd7), temp=mean(temp) ) %>% distinct(  perc_rot, .keep_all = T )

  length(unique(mon.survival$perc_rot))

  surv_gomp = summary(model.surv.g, newdata = new.data, type = "median", tidy = TRUE )

  ggplot(surv_gomp, aes(perc_rot, est) )+
    geom_point()



}


### ************************************************************************* ##  
###                                 MODEL PLOT
### ************************************************************************* ##

model.surv.w = fit_weib.bestmpr #flexsurvreg (weibull)
model.surv.g = fit_gomp.bestmpr #flexsurvreg


## *********** 2018 *************
{
new.data18 = mon.survival %>% ## plot solo valori medi weight
  group_by( temp, cum_deg_day.scale, wd7, perc_rot, rain7, t7, year, Ha) %>%
  dplyr::summarize( weight_g = mean(weight_g, na.rm=T),
                    weight_g.scale = mean(weight_g.scale, na.rm=T) ) %>% ungroup() %>%
  filter(year == 2018)

  
  surv_gomp18 = summary(model.surv.g, newdata = new.data18, type = "survival", tidy = TRUE, t = c(0:max(mon.survival18$time)) ) 
  # surv_weib18 = summary(model.surv.w, newdata = new.data18, type = "survival", tidy = TRUE, t = c(0:max(mon.survival18$time)) ) 
  
  # surv_gomp18 = summary(model.surv.g, newdata = new.data18, type = "survival", tidy = TRUE, t = c(0:max(mon.survival18$time)) ) %>%
  # mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
  # surv_weib18 = summary(model.surv.w, newdata = new.data18, type = "survival", tidy = TRUE, t = c(0:max(mon.survival18$time)) ) %>%
  # mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )


km_fit18 = survfit(Surv(time, status) ~ wd7 , data = mon.survival18, conf.int=0.95 )
names(km_fit18$strata) <- gsub("wd7=", "",  c("0.56","3.09", "3.26") )

# col=rep("black", length(unique(mon.survival$wd7)))
km.plot18 = ggsurvplot(km_fit18, xlim =c(-0.3,(max(mon.survival18$time)+0.3)), conf.int = F, legend.title="Harvest date (dafb)",
                       # palette = c("#009E73", "#E69F00", "#0072B2",
                       #             "#009E73", "#E69F00", "#0072B2"),
                       palette = c("#009E73", "#E69F00", "#0072B2"),
                       linetype = "twodash",
                       # fun = "event",   
                       legend="none")

# km.plot18$plot$data$surv = 1- km.plot18$plot$data$surv
# km.plot18$plot$data$surv = km.plot18$plot$data$surv*100

km.plot18$plot +
  # ggplot() +
  ylab("Survival probability (%)")+
  xlab("Storage time (d)")+
  # geom_line(surv_gomp18, mapping = aes(x = time, y = 1-est, colour=as.factor(wd7)), size = 0.9, )+
  geom_line(surv_gomp18, mapping = aes(x = time, y = est, colour=as.factor(wd7)), size = 0.9, )+
  # geom_hline(yintercept = 0.5, colour="black", linetype="dashed")+
  # scale_y_continuous( limits = c(0,1.02), expand = c(0,0), breaks = seq(0,1,0.2) )+
      scale_y_continuous( limits = c(-0.01,1.02), expand = c(0,0), breaks = seq(0,1,0.25) )+
  scale_x_continuous( expand = c(0,0), breaks = seq(0,14,2) )+
  theme_bw(base_size = 16) + 
  theme(
    axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),
  # legend.position = "none",
  legend.position = "bottom",
  # Hide panel borders and remove grid lines
  panel.grid.major = element_blank())

}

## *********** 2018 with irrigation *************
{
  new.data18 = mon.survival18 %>% ## plot solo valori medi weight
    group_by( temp, cum_deg_day.scale, wd7, perc_rot, rain7, t7, year, Ha, Irr) %>%
    dplyr::summarize( weight_g = mean(weight_g, na.rm=T),
                      weight_g.scale = mean(weight_g.scale, na.rm=T) ) %>% ungroup() %>%
    filter(year == 2018)
  
  
  surv_gomp18 = summary(model.surv.g, newdata = new.data18, type = "survival", tidy = TRUE, t = c(0:max(mon.survival18$time)) ) 
  surv_gomp18.test = merge(surv_gomp18, new.data18 ) 
  
  
  km_fit18 = survfit(Surv(time, status) ~ wd7 + Irr , data = mon.survival18, conf.int=0.95 )
  
  # col=rep("black", length(unique(mon.survival$wd7)))
  km.plot18 = ggsurvplot(km_fit18, xlim =c(-0.3,(max(mon.survival18$time)+0.3)), conf.int = F, legend.title="Harvest date (dafb)",
                         # palette = c("#009E73", "#E69F00", "#0072B2",
                         #             "#009E73", "#E69F00", "#0072B2"),
                         palette = c("#009E73", "#E69F00", "#0072B2", 
                           "#009E73", "#009E73", "#E69F00","#E69F00", "#0072B2", "#0072B2"),
                         linetype = "twodash",
                         # fun = "event",   
                         legend="none")
  
  # km.plot18$plot$data$surv = 1- km.plot18$plot$data$surv
  # km.plot18$plot$data$surv = km.plot18$plot$data$surv*100
  
  km.plot18$plot + facet_grid(.~Irr)+
    # ggplot() +
    ylab("Survival probability (%)")+
    xlab("Storage time (d)")+
    # geom_line(surv_gomp18, mapping = aes(x = time, y = 1-est, colour=as.factor(wd7)), size = 0.9, )+
    geom_line(surv_gomp18.test, mapping = aes(x = time, y = est, colour=as.factor(wd7)), size = 0.9, )+
    # geom_hline(yintercept = 0.5, colour="black", linetype="dashed")+
    # scale_y_continuous( limits = c(0,1.02), expand = c(0,0), breaks = seq(0,1,0.2) )+
    scale_y_continuous( limits = c(-0.01,1.02), expand = c(0,0), breaks = seq(0,1,0.25) )+
    scale_x_continuous( expand = c(0,0), breaks = seq(0,14,2) )+
    theme_bw(base_size = 16) + 
    theme(
      axis.text.x=element_text(colour="black"),
      axis.text.y=element_text(colour="black"),
      # legend.position = "none",
      legend.position = "bottom",
      strip.background = element_blank(), strip.text.x = element_blank(), 
      # Hide panel borders and remove grid lines
      panel.grid.major = element_blank(),  panel.grid.minor = element_blank())
  
}

## *********** 2019 *************
{
  
  new.data19 = mon.survival19 %>% ## plot solo valori medi weight
    group_by( temp, cum_deg_day.scale, wd7, perc_rot, rain7, t7, year, Ha) %>%
    # group_by( temp, cum_deg_day.scale, wd7, perc_rot, rain7, t7, year, irr) %>%
    dplyr::summarize( weight_g = mean(weight_g, na.rm=T),
                      weight_g.scale = mean(weight_g.scale, na.rm=T) ) %>% ungroup() %>%
    filter(year == 2019)
  
  surv_gomp19 = summary(model.surv.g, newdata = new.data19, type = "survival", tidy = TRUE, t = c(0:max(mon.survival19$time))) 
  surv_gomp19.test = merge(surv_gomp19, new.data19 ) 

  
  # surv_gomp19.all = summary(model.surv.g, newdata = mon.survival19, type = "survival", tidy = TRUE, t = c(0:max(mon.survival19$time)))

  # surv_gomp19.all90 = surv_gomp19.all %>% filter( quantile(weight_g.scale, 0.9) < weight_g.scale)
  # surv_gomp19.all50 = surv_gomp19.all %>% filter( quantile(weight_g.scale, 0.5) < weight_g.scale)
  # surv_gomp19.all10 = surv_gomp19.all %>% filter( quantile(weight_g.scale, 0.1) > weight_g.scale)
  
  
  
  
  
  mon.survival19$wd7 = round(mon.survival19$wd7,2)
  
  
  km_fit19 = survfit(Surv(time, status) ~ wd7 + temp , data = mon.survival19)
  km.plot19 = ggsurvplot(km_fit19,xlim =c(-0.3,(max(mon.survival19$time)+0.3)),  conf.int = F, legend.title="Storage temperature (°C)",
                         palette = c("blue3","darkorange","firebrick3",
                                     "blue3","darkorange","firebrick3",
                                     "blue3","darkorange","firebrick3"),
                         # 
                         # palette = c("#5C5A5A","grey80","#000000","#5C5A5A","grey80","#000000",
                         #             "#5C5A5A","grey80","#000000","#5C5A5A","grey80","#000000"),
                         
                         linetype = "twodash", 
                         # fun="event",
                         legend="none"
                         )
  
  labels <- data.frame(wd7=unique(km.plot19$plot$data$wd7), label = c("a", "b"))
  km.plot19$plot + facet_grid(.~wd7)+
    # ggplot() + facet_grid(wd7~.)+
    ylab("Survival probability (%)")+
    xlab("Storage time (d)")+
    
    # geom_ribbon(surv_gomp19, mapping = aes(x = time, ymin = lcl, ymax= ucl, fill=as.factor(temp)), inherit.aes = FALSE, alpha=0.3 )+
    geom_line(surv_gomp19, mapping = aes(x = time, y = est, colour=as.factor(temp)), size = 0.9)+
   
    # stat_summary(surv_gomp19.all, mapping = aes(x = time, y = est, colour=as.factor(temp)), fun.data = "mean_sd", geom = "line", size = 0.9)+
    # stat_summary(surv_gomp19.all10, mapping = aes(x = time, y = est, colour=as.factor(temp)), fun.data = "mean_sd", geom = "line", size = 0.3)+
    # stat_summary(surv_gomp19.all90, mapping = aes(x = time, y = est, colour=as.factor(temp)), fun.data = "mean_sd", geom = "line", size = 0.3)+
    geom_text(x = 2, y = 0.1, aes(label = label), data = labels, fontface = "bold", size = 7)+
    # scale_colour_manual(values = c("grey80","#5C5A5A","#000000" ), name = "Storage temperature (°C)")+
    # scale_linetype_discrete(name = "Storage temperature (°C)")+
    # geom_line(surv_weib18, mapping = aes(x = time, y = est, colour=wd7), size = 0.6, colour="red")+
    # geom_hline(yintercept = 0.5, colour="black", linetype="dashed")+
    scale_y_continuous( limits = c(-0.01,1.02), expand = c(0,0), breaks = seq(0,1,0.25) )+
    scale_x_continuous( expand = c(0,0),
                        breaks = seq(0,21,3) )+
    theme_bw(base_size = 16) + 
    theme(
      axis.text.x=element_text(colour="black"),
      axis.text.y=element_text(colour="black"),
      # legend.position = "none",
      legend.position = "bottom",
      # Hide panel borders and remove grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),  
      strip.background = element_blank(), strip.text.x = element_blank())+
    guides(fill = guide_colourbar( ticks=element_line(color='black'), border=element_line(color='black'), ))
  
}

## *********** 2019 with irrigation *************
{
  
  new.data19 = mon.survival19 %>% ## plot solo valori medi weight
    group_by( temp, cum_deg_day.scale, wd7, perc_rot, rain7, t7, year, irr) %>%
    dplyr::summarize( weight_g = mean(weight_g, na.rm=T),
                      weight_g.scale = mean(weight_g.scale, na.rm=T) ) %>% ungroup() %>%
    filter(year == 2019)
  
  surv_gomp19 = summary(model.surv.g, newdata = new.data19, type = "survival", tidy = TRUE, t = c(0:max(mon.survival19$time))) 
  surv_gomp19.test = merge(surv_gomp19, new.data19 ) 
  
    # surv_gomp19.all = summary(model.surv.g, newdata = mon.survival19, type = "survival", tidy = TRUE, t = c(0:max(mon.survival19$time)))
  # surv_gomp19.all90 = surv_gomp19.all %>% filter( quantile(weight_g.scale, 0.9) < weight_g.scale)
  # surv_gomp19.all50 = surv_gomp19.all %>% filter( quantile(weight_g.scale, 0.5) < weight_g.scale)
  # surv_gomp19.all10 = surv_gomp19.all %>% filter( quantile(weight_g.scale, 0.1) > weight_g.scale)
  
  
  mon.survival19$wd7 = round(mon.survival19$wd7,2)
  
  
  km_fit19 = survfit(Surv(time, status) ~ wd7 + temp + irr , data = mon.survival19)
  km.plot19 = ggsurvplot(km_fit19,xlim =c(-0.3,(max(mon.survival19$time)+0.3)),  conf.int = F, legend.title="Storage temperature (°C)",
                         palette = c("blue3", "darkorange", "firebrick3",
                           "blue3","blue3", "darkorange", "darkorange", "firebrick3", "firebrick3",
                                     "blue3","blue3", "darkorange", "darkorange", "firebrick3", "firebrick3"),
                         # 
                         # palette = c("#5C5A5A","grey80","#000000","#5C5A5A","grey80","#000000",
                         #             "#5C5A5A","grey80","#000000","#5C5A5A","grey80","#000000"),
                         
                         linetype = "twodash", 
                         # fun="event",
                         # legend="none"
  )
  
  labels <- data.frame(wd7=unique(km.plot19$plot$data$wd7), label = c("a", "b"))
  km.plot19$plot + facet_grid(irr~wd7)+
    # ggplot() + facet_grid(wd7~.)+
    ylab("Survival probability (%)")+
    xlab("Storage time (d)")+
    
    # geom_ribbon(surv_gomp19, mapping = aes(x = time, ymin = lcl, ymax= ucl, fill=as.factor(temp)), inherit.aes = FALSE, alpha=0.3 )+
    # geom_line(surv_gomp19, mapping = aes(x = time, y = est, colour=as.factor(temp)), size = 0.9)+
    
    stat_summary(surv_gomp19, mapping = aes(x = time, y = est, colour=as.factor(temp)), fun.data = "mean_sd", geom = "line", size = 0.9)+
    # stat_summary(surv_gomp19.all10, mapping = aes(x = time, y = est, colour=as.factor(temp)), fun.data = "mean_sd", geom = "line", size = 0.3)+
    # stat_summary(surv_gomp19.all90, mapping = aes(x = time, y = est, colour=as.factor(temp)), fun.data = "mean_sd", geom = "line", size = 0.3)+
    geom_text(x = 2, y = 0.1, aes(label = label), data = labels, fontface = "bold", size = 7)+
    # scale_colour_manual(values = c("grey80","#5C5A5A","#000000" ), name = "Storage temperature (°C)")+
    # scale_linetype_discrete(name = "Storage temperature (°C)")+
    # geom_line(surv_weib18, mapping = aes(x = time, y = est, colour=wd7), size = 0.6, colour="red")+
    # geom_hline(yintercept = 0.5, colour="black", linetype="dashed")+
    scale_y_continuous( limits = c(-0.01,1.02), expand = c(0,0), breaks = seq(0,1,0.25) )+
    scale_x_continuous( expand = c(0,0),
                        breaks = seq(0,21,3) )+
    theme_bw(base_size = 16) + 
    theme(
      axis.text.x=element_text(colour="black"),
      axis.text.y=element_text(colour="black"),
      # legend.position = "none",
      legend.position = "bottom",
      # Hide panel borders and remove grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),  
      strip.background = element_blank(), strip.text.x = element_blank())+
    guides(fill = guide_colourbar( ticks=element_line(color='black'), border=element_line(color='black'), ))
  
}


## *********** test18 ************* 
## survival analysis per experimental treatment
{

  new.data18 = mon.survival18 %>% ## plot solo valori medi weight
    group_by( temp, cum_deg_day.scale, wd7, perc_rot, rain7, t7, Irr, FL, harvest_n) %>%
    dplyr::summarize( weight_g.scale = mean(weight_g.scale, na.rm=T) ) %>% ungroup() 
  
  model.surv.g = fit_gomp.bestmpr #flexsurvreg (weibull)
  # model.surv.w = fit_weib.bestmpr #flexsurvreg
  
  surv_gomp18 = summary(model.surv.g, newdata = new.data18, type = "survival", tidy = TRUE, t = c(0:max(mon.survival18$time))) %>%
    mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
  surv_gomp18.test = merge(surv_gomp18, new.data18 ) 
  
  # surv_weib18 = summary(model.surv.w, newdata = new.data18, type = "survival", tidy = TRUE, t = c(0:max(mon.survival18$time))) %>%
  #   mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
  # surv_weib18.test = merge(surv_weib18, new.data18 )
  
  
  tidy(model.surv.g, pars = "coefs")


  mon.survival18$wd7 = round(mon.survival18$wd7,2) 
  
  # km_fit18 = survfit(Surv(time, status) ~ Irr + FL + harvest_n, data = mon.survival18)
  # km.plot18 = ggsurvplot(km_fit18, xlim =c(0,15), legend="none", fun="event", ylim=c(0,1), size=0.7, linetype = "dashed")
  #                        # palette = c("#5C5A5A","grey80","#000000",
  #                        #             "#5C5A5A","grey80","#000000",
  #                        #             "#5C5A5A","grey80","#000000"),
  #                        # linetype = c(11,6,7,11,6,7,11,6,7)  )
  # 
  # 
  # km.plot18$plot + facet_grid(FL~Irr)+
  #   # ggplot() + facet_grid(wd7~.)+
  #   ylab("Cumulative incidence (%)")+
  #   geom_line(surv_gomp18.test, mapping = aes(x = time, y = 1-est, colour=harvest_n), size = 0.8)+
  #   # scale_colour_manual(values = c("grey80","#5C5A5A","#000000" ), name = "Storage temperature (°C)")+
  #   # scale_linetype_discrete(name = "Storage temperature (°C)")+
  #   # geom_line(surv_weib18, mapping = aes(x = time, y = est, colour=wd7), size = 0.6, colour="red")+
  #   # geom_hline(yintercept = 0.5, colour="black", linetype="dashed")+
  #   theme_pubr() + theme(legend.position = "none")

  
  ggplot()+
    geom_point(mon18.db.plot, mapping = aes(dah.int, incidence, colour=harvest_n, shape=harvest_n), size=1.8 )+
    geom_line(mon18.db.plot, mapping = aes(dah.int, incidence, colour=harvest_n), linetype="dashed", size=0.4, alpha=0.8 )+
    geom_line( surv_gomp18.test, mapping = aes(x = time, y = (1-est)*100, colour=harvest_n), size = 0.8)+
    geom_ribbon( surv_gomp18.test, mapping = aes(x = time, ymin = (1-lcl)*100, ymax = (1-ucl)*100, fill=harvest_n), alpha = 0.3)+
    
    ylab("Cumulative incidence (%)")+
    # theme(legend.position = "none")+
    facet_grid(FL~Irr)+
    theme_pubr() 
    
  ggplot()+
    stat_summary(mon18.db.plot, mapping = aes(dah.int, incidence, colour=harvest_n), fun.data ="mean_cl_boot", geom="point", size=2.5 )+
    stat_summary(mon18.db.plot, mapping = aes(dah.int, incidence, colour=harvest_n), fun.data ="mean_cl_boot", geom="line", linetype="dashed", size=0.4, alpha=0.8 )+

    stat_summary( surv_gomp18.test, mapping = aes(x = time, y = (1-est)*100, colour=harvest_n), fun.data ="mean_cl_boot", geom="line")+
    stat_summary( surv_gomp18.test, mapping = aes(x = time, y = (1-lcl)*100, colour=harvest_n), fun.data ="mean_cl_boot", geom="line", alpha=0.3)+
    stat_summary( surv_gomp18.test, mapping = aes(x = time, y = (1-ucl)*100, colour=harvest_n), fun.data ="mean_cl_boot", geom="line", alpha=0.3)+
    
    scale_colour_manual(values=c("royalblue","orangered3", "red4"))+
    # scale_color_brewer(name="Harvest date", palette = "YlOrRd")+
    
    ylab("Cumulative incidence (%)")+
    # theme(legend.position = "none")+
    facet_grid(.~Irr)+
    theme_pubr() 
  
  
 
}

## *********** test19 *************
## survival analysis per experimental treatment
{
  new.data19 = mon.survival19 %>% ## plot solo valori medi weight
    group_by( temp, cum_deg_day.scale, wd7, perc_rot, rain7, t7, irr, harvest_n) %>%
    dplyr::summarize( weight_g.scale = mean(weight_g.scale, na.rm=T) ) %>% ungroup() 
  
  model.surv.g = fit_gomp.bestmpr #flexsurvreg (weibull)
  # model.surv.w = fit_weib.bestmpr #flexsurvreg
  
  surv_gomp19 = summary(model.surv.g, newdata = new.data19, type = "survival", tidy = TRUE, t = c(0:max(mon.survival19$time))) %>%
    mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
  surv_gomp19.test = merge(surv_gomp19, new.data19 )
  
  # surv_weib19 = summary(model.surv.w, newdata = new.data19, type = "survival", tidy = TRUE, t = c(0:max(mon.survival19$time))) %>%
  #   mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
  # surv_weib19.test = merge(surv_weib19, new.data19 ) 
    
  
  
  mon.survival19$wd7 = round(mon.survival19$wd7,2) 
  
  
  # km_fit19 = survfit(Surv(time, status) ~ irr + harvest_n + temp, data = mon.survival19)
  # km.plot19 = ggsurvplot(km_fit19, xlim =c(0,21), legend="none", fun="event", ylim=c(0,1), size=0.7, linetype = "dashed")
  # # palette = c("#5C5A5A","grey80","#000000",
  # #             "#5C5A5A","grey80","#000000",
  # #             "#5C5A5A","grey80","#000000"),
  # # linetype = c(11,6,7,11,6,7,11,6,7)  )
  # 
  # km.plot19$plot + facet_grid(irr~harvest_n)+
  #   # ggplot() + facet_grid(wd7~.)+
  #   ylab("Survival probability")+
  #   geom_line(surv_gomp19.test, mapping = aes(x = time, y = 1-est, colour=as.factor(temp)), size = 0.8)+
  #   # scale_colour_manual(values = c("grey80","#5C5A5A","#000000" ), name = "Storage temperature (°C)")+
  #   # scale_linetype_discrete(name = "Storage temperature (°C)")+
  #   # geom_line(surv_weib18, mapping = aes(x = time, y = est, colour=wd7), size = 0.6, colour="red")+
  #   # geom_hline(yintercept = 0.5, colour="black", linetype="dashed")+
  #   theme_pubr() + theme(legend.position = "none")
  
  
  ggplot()+
    geom_point(mon19.db.plot, mapping = aes(dah.int, incidence, colour=temp, shape=temp), size=1.8 )+
    geom_line(mon19.db.plot, mapping = aes(dah.int, incidence, colour=temp), linetype="dashed", size=0.2, alpha=0.8 )+
    
    geom_line(surv_gomp19.test, mapping = aes(x = time, y = (1-est)*100, colour=as.factor(temp)), size = 0.8)+
    # geom_ribbon( surv_gomp19.test, mapping = aes(x = time, ymin = (1-lcl)*100, ymax = (1-ucl)*100, fill=as.factor(temp)), alpha = 0.3)+
    
    ylab("Cumulative incidence (%)")+
    # theme(legend.position = "none")+
    facet_grid(irr~harvest_n)+
    theme_pubr() 
  
  
  
  mon19.db.plot = mon19.db.plot %>% filter(temp!="2")
  # surv_gomp19.test = surv_gomp19.test %>% filter(temp!="2")
  
  
  ggplot()+
    stat_summary(mon19.db.plot, mapping = aes(dah.int, incidence, colour=harvest_n), fun.data ="mean_cl_boot" )+
    geom_jitter(mon19.db.plot, mapping = aes(dah.int, incidence, colour=harvest_n))+
    stat_summary(mon19.db.plot, mapping = aes(dah.int, incidence, colour=harvest_n), fun.data ="mean_cl_boot", geom="line", linetype="dashed", size=0.4, alpha=0.8 )+
    
    stat_summary( surv_gomp19.test, mapping = aes(x = time, y = (1-est)*100, colour=harvest_n), fun.data ="mean_cl_boot", geom="line")+
    stat_summary( surv_gomp19.test, mapping = aes(x = time, y = (1-lcl)*100, colour=harvest_n), fun.data ="mean_cl_boot", geom="line", alpha=0.3)+
    stat_summary( surv_gomp19.test, mapping = aes(x = time, y = (1-ucl)*100, colour=harvest_n), fun.data ="mean_cl_boot", geom="line", alpha=0.3)+
    
    scale_colour_manual(values=c("royalblue","orangered3"))+
    # scale_color_brewer(name="Harvest date", palette = "YlOrRd")+
    
    ylab("Cumulative incidence (%)")+
    # theme(legend.position = "none")+
    facet_grid(temp~irr)+
    theme_pubr() 
  
  
  
}



## *********** test18 ************* 
## cumulative incidence per experimental treatment
{
  new.data18 = mon.survival18 %>% dplyr::select(harvest_n, Irr, FL, temp, wd7, perc_rot, cum_deg_day.scale, weight_g.scale, id.fruit, rain7)
  
  length(unique(new.data18$weight_g.scale))
  n_occur <- data.frame(table(new.data18$weight_g.scale))
  n_occur[n_occur$Freq > 1,]
  
  model.surv.g = fit_gomp.bestmpr #flexsurvreg (gompertz)
  # model.surv.g = fit_weib.bestmpr #flexsurvreg (weibull)
  
  surv_gomp18 = summary(model.surv.g, newdata = new.data18, type = "mean", tidy = TRUE) %>%
    mutate( weight_g.h = (weight_g.scale*wg.si.scale) + wg.mu.scale )
  
  # cumulated incidence
  surv_gomp18.cuminc = base::merge(surv_gomp18, new.data18, all.y = FALSE) %>% mutate(rep=15)
  
  rep.vec = c()
  for(i in 1:length(surv_gomp18.cuminc$est) ){
    rep.vec =  c(rep.vec,0:14)  }
  
  surv_gomp18.cuminc = expandRows(surv_gomp18.cuminc, "rep")
  


    surv_gomp18.cuminc$time   = rep.vec
    surv_gomp18.cuminc$infection   = ifelse(surv_gomp18.cuminc$time >= surv_gomp18.cuminc$est, "1", "0")
    
    surv_gomp18.cuminc = surv_gomp18.cuminc %>%
      group_by(harvest_n, time, Irr, FL) %>%
      dplyr::summarize( incidence = round(sum( infection=="1" )*100 /( length(infection)),2), 
                        weight_mean = round(mean(weight_g.h, na.rm=T),2),  weight_sd = sd(weight_g.h, na.rm=T)) %>% # cumulated incidence 
      ungroup()   
    
  ggplot()+
    geom_line( surv_gomp18.cuminc, mapping=aes(time,    incidence, colour=harvest_n) )+
    # geom_smooth( surv_gomp18.cuminc, mapping=aes(time,    incidence, colour=harvest_n), se=F )+
    # geom_point(surv_gomp18.cuminc, mapping=aes(time,    incidence, colour=harvest_n) )+
    # geom_line( mon18.db.plot,      mapping=aes(dah.int, incidence, colour=harvest_n ), linetype=2 )+
    # geom_point(mon18.db.plot,      mapping=aes(dah.int, incidence, colour=harvest_n ) )+
    # geom_line(surv_gomp18.cuminc, mapping=aes(time,incidence, colour=harvest_n) )+
    facet_grid(FL~Irr)+
    theme_bw()
}



## *********** test19 ************* 
## cumulative incidence per experimental treatment
{
  new.data19 = mon.survival19 
  
  length(unique(new.data19$weight_g.scale))
  n_occur <- data.frame(table(new.data19$weight_g.scale))
  n_occur[n_occur$Freq > 1,]
  
  model.surv.g = fit_gomp.bestmpr #flexsurvreg (gompertz)
  # model.surv.g = fit_weib.bestmpr #flexsurvreg (weibull)
  
  surv_gomp19 = summary(model.surv.g, newdata = new.data19, type = "mean", tidy = TRUE) %>%
    mutate( weight_g.h = (weight_g.scale*wg.si.scale) + wg.mu.scale )
  
  # cumulated incidence
  surv_gomp19.cuminc = base::merge(surv_gomp19, new.data19, all.y = FALSE) %>% mutate(rep=21)
  
  rep.vec = c()
  for(i in 1:length(surv_gomp19.cuminc$est) ){
    rep.vec =  c(rep.vec,0:20)  }
  
  surv_gomp19.cuminc = expandRows(surv_gomp19.cuminc, "rep")
  
  
  
  surv_gomp19.cuminc$time   = rep.vec
  surv_gomp19.cuminc$infection   = ifelse(surv_gomp19.cuminc$time >= surv_gomp19.cuminc$est, "1", "0")
  
  surv_gomp19.cuminc = surv_gomp19.cuminc %>%
    group_by(harvest_n, time, irr, temp) %>%
    dplyr::summarize( incidence = round(sum( infection=="1" )*100 /( length(infection)),2), 
                      weight_mean = round(mean(weight_g.h, na.rm=T),2),  weight_sd = sd(weight_g.h, na.rm=T)) %>% # cumulated incidence 
    ungroup()   
  
  ggplot()+
    geom_line( surv_gomp19.cuminc, mapping=aes(time,    incidence, colour=as.factor(temp)) )+
    geom_line( mon19.db.plot,      mapping=aes(dah.int, incidence, colour=as.factor(temp) ), linetype=2 )+
    # geom_point(mon18.db.plot,      mapping=aes(dah.int, incidence, colour=harvest_n ) )+
    # geom_line(surv_gomp18.cuminc, mapping=aes(time,incidence, colour=harvest_n) )+
    facet_grid(irr~harvest_n)+
    theme_bw()
}


## **********************************************************************************************

new.data.article = mon.survival %>% ## plot solo valori medi weight
  group_by( cum_deg_day.scale, wd7, perc_rot, rain7, t7, cum_deg_day, year, harvest_n) %>%
  dplyr::summarize( weight_g.scale = mean(weight_g.scale, na.rm=T), 
                    weight_sd = sd(weight_g, na.rm=T),
                    weight_g = mean(weight_g, na.rm=T) ) %>% ungroup()



new.data = mon.survival %>% ## plot solo valori medi weight
  group_by( temp, cum_deg_day.scale, wd7, perc_rot, rain7, t7, cum_deg_day, year, harvest_n) %>%
  dplyr::summarize( weight_g.scale = mean(weight_g.scale, na.rm=T), weight_g = mean(weight_g, na.rm=T) ) %>% ungroup()
# new.data = mon.survival ## plot tutti valori weight


haz_gomp3  = summary(model.surv.g, newdata = new.data, type = "hazard",   tidy = TRUE, t = c(0:21)) %>%
  mutate( weight_g.h = (weight_g.scale*wg.si.scale) + wg.mu.scale )
surv_gomp3 = summary(model.surv.g, newdata = new.data, type = "survival", tidy = TRUE, t = c(0:21)) %>%
  mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
haz_weib3  = summary(model.surv.w, newdata = new.data, type = "hazard",   tidy = TRUE, t = c(0:21)) %>%
  mutate( weight_g.h = (weight_g.scale*wg.si.scale) + wg.mu.scale )
surv_weib3 = summary(model.surv.w, newdata = new.data, type = "survival", tidy = TRUE, t = c(0:21)) %>%
  mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )

summary(model.surv.g, newdata = new.data, type = "median", tidy = TRUE, t = c(0:21)) %>%
  mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )



km_fit = survfit(Surv(time, status) ~ temp + wd7 , data = mon.survival)
col=rep("black", length(unique(mon.survival$temp))*length(unique(mon.survival$wd7)))
km.plot = ggsurvplot(km_fit, xlim =c(0,22), palette = col, legend="none"  )
# km.plot = ggsurvplot(km_fit) # statitsticam. !=

km.plot$plot + facet_grid(temp~wd7, labeller = label_both)+
# ggplot() + facet_grid(temp~wd7, labeller = label_both)+
  ylab("Survival probability")+
  geom_line(surv_gomp3, mapping = aes(x = time, y = est), size = 0.6, colour="blue")+
  geom_line(surv_weib3, mapping = aes(x = time, y = est), size = 0.6, colour="red")+
  geom_hline(yintercept = 0.5, colour="black", linetype="dashed")+
  geom_label( surv_gomp3, mapping = aes(x = 5.5, y = 0.25, label=paste0("mean.weight = ",weight_g.h)), size=2.5)+
  geom_label( surv_gomp3, mapping = aes(x = 5.5, y = 0.10, label=paste0("perc_rot= ",round(perc_rot,2))), size=2.5)+
  theme_bw()




### test perc_rot
cum_deg_day.scale = mean(mon.survival$cum_deg_day.scale)
weight_g.scale  = mean(mon.survival$weight_g.scale);
temp = unique(mon.survival$temp)
wd7 = mean(mon.survival$wd7)
t7=25
perc_rot = c(unique(mon.survival$perc_rot),10,20)
rain7 =  mean(mon.survival$rain7)
  
new.data.list = list(temp = temp, weight_g.scale =weight_g.scale, cum_deg_day.scale=cum_deg_day.scale, 
                     wd7=wd7, t7=t7, perc_rot=perc_rot, rain7=rain7)
new.data = expand.grid(new.data.list)
new.data$temp.bell   = ( (new.data$temp - Tmax) * ((new.data$temp-Tmin)^2) ) /
  ( (Topt-Tmin) * ( ((Topt-Tmin)*(new.data$temp-Topt)) - ((Topt-Tmax)*(Topt+Tmin-2*new.data$temp)) ) ); 

model.surv = fit_gomp.bestmpr #flexsurvreg

haz_gomp3  = summary(model.surv, newdata = new.data, type = "hazard",   tidy = TRUE, t = c(0:21))%>%
  # mutate( GDD.h = round((cum_deg_day.scale*GDD.si.scale) + GDD.mu.scale,2),
          # weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
  mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
surv_gomp3 = summary(model.surv, newdata = new.data, type = "survival", tidy = TRUE, t = c(0:21))%>%
  # mutate( GDD.h = round((cum_deg_day.scale*GDD.si.scale) + GDD.mu.scale,2),
  mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )


ggplot() + facet_grid(perc_rot~., labeller = label_both) +
  geom_line(haz_gomp3, mapping = aes(x = time, y = est, group=temp, colour=temp), size = 0.6)+
  # scale_color_paletteer_c("harrypotter::ronweasley", direction = -1, name = "mean daily wetness duration (h)")+
  scale_color_paletteer_c("harrypotter::ronweasley", direction = -1, name = "storage temperature (°C)")+
  theme_bw() 


ggplot() +  facet_grid(.~round(perc_rot,2), labeller = label_both) + ylab("infection probability (%)")+
    geom_line(surv_gomp3, mapping = aes(x = time, y = est, group=temp, colour=temp), size = 1)+
  # scale_color_paletteer_c("harrypotter::ronweasley", direction = -1, name = "mean daily wetness duration (h)")+
  scale_color_paletteer_c("pals::coolwarm", direction = 1, name = "storage temperature (°C)")+
  theme_bw(base_size = 16) + theme(legend.position = "bottom")
  



### **************
###   test weight
### **************

# cum_deg_day.scale = mean(mon.survival$cum_deg_day.scale)
cum_deg_day.scale = max(mon.survival$cum_deg_day.scale)


weight_g.test  = c(100,150,200,250);
weight_g.scale =(weight_g.test - wg.mu.scale)/wg.si.scale 

# temp = unique(mon.survival$temp)
temp = c(2,10,15,20,25)

# wd7 = mean(mon.survival$wd7)
wd7 = max(mon.survival$wd7)


t7=25

# perc_rot = mean(mon.survival$perc_rot)
perc_rot = max(mon.survival$perc_rot)


## *****
new.data.list = list(temp = temp, weight_g.scale =weight_g.scale, cum_deg_day.scale=cum_deg_day.scale, 
                     wd7=wd7, t7=t7, perc_rot=perc_rot)
new.data = expand.grid(new.data.list)
new.data$temp.bell   = ( (new.data$temp - Tmax) * ((new.data$temp-Tmin)^2) ) /
  ( (Topt-Tmin) * ( ((Topt-Tmin)*(new.data$temp-Topt)) - ((Topt-Tmax)*(Topt+Tmin-2*new.data$temp)) ) ); 

model.surv = fit_gomp.bestmpr #flexsurvreg

haz_gomp3  = summary(model.surv, newdata = new.data, type = "hazard",   tidy = TRUE, t = c(0:21))%>%
  mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
surv_gomp3 = summary(model.surv, newdata = new.data, type = "survival", tidy = TRUE, t = c(0:21))%>%
  mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
mean_gomp3 = summary(model.surv, newdata = new.data, type = "mean", tidy = TRUE)%>%
  mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )


ggplot() + facet_grid(.~temp, labeller = label_both) + ylab("infection probability (%)")+
  geom_line(haz_gomp3, mapping = aes(x = time, y = est, group=weight_g.h, colour=weight_g.h), size = 1)+
  # scale_color_paletteer_c("harrypotter::ronweasley", direction = -1, name = "mean daily wetness duration (h)")+
  scale_color_paletteer_c("pals::coolwarm", direction = 1, name = "fruit weight (g)")+
  # theme(legend.position = "bottom")+
  theme_bw(base_size = 16) 


ggplot() + facet_grid(.~temp, labeller = label_both) + ylab("infection probability (%)")+
  geom_line(surv_gomp3, mapping = aes(x = time, y = 1-est, group=weight_g.h, colour=weight_g.h), size = 1)+
  # scale_color_paletteer_c("harrypotter::ronweasley", direction = -1, name = "mean daily wetness duration (h)")+
  # scale_color_paletteer_c("pals::coolwarm", direction = 1, name = "fruit weight (g)")+
  geom_hline(yintercept=0.5, linetype=2)+
  theme_bw(base_size = 16) +   theme(legend.position = "bottom")
  # theme_bw(base_size = 16)


# ### test
# cum_deg_day.scale = mean(mon.survival$cum_deg_day.scale)
# weight_g.scale  = mean(mon.survival$weight_g.scale);
# temp = unique(mon.survival$temp)
# wd7 = mean(mon.survival$wd7)
# t7=c(20,25,30)
# perc_rot = c(0,2,10,20)
# 
# new.data.list = list(temp = temp, weight_g.scale =weight_g.scale, cum_deg_day.scale=cum_deg_day.scale, 
#                      wd7=wd7, t7=t7, perc_rot=perc_rot)
# new.data = expand.grid(new.data.list)
# new.data$temp.bell   = ( (new.data$temp - Tmax) * ((new.data$temp-Tmin)^2) ) /
#   ( (Topt-Tmin) * ( ((Topt-Tmin)*(new.data$temp-Topt)) - ((Topt-Tmax)*(Topt+Tmin-2*new.data$temp)) ) ); 
# 
# model.surv = fit_gomp.bestmpr #flexsurvreg
# 
# haz_gomp3  = summary(model.surv, newdata = new.data, type = "hazard",   tidy = TRUE, t = c(0:21))%>%
#   # mutate( GDD.h = round((cum_deg_day.scale*GDD.si.scale) + GDD.mu.scale,2),
#   # weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
#   mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
# surv_gomp3 = summary(model.surv, newdata = new.data, type = "survival", tidy = TRUE, t = c(0:21))%>%
#   # mutate( GDD.h = round((cum_deg_day.scale*GDD.si.scale) + GDD.mu.scale,2),
#   mutate( weight_g.h = round((weight_g.scale*wg.si.scale) + wg.mu.scale,2) )
# 
# 
# ggplot() + facet_grid(perc_rot~temp, labeller = label_both) +
#   geom_line(haz_gomp3, mapping = aes(x = time, y = est, group=t7, colour=t7), size = 0.6)+
#   # scale_color_paletteer_c("harrypotter::ronweasley", direction = -1, name = "mean daily wetness duration (h)")+
#   scale_color_paletteer_c("harrypotter::ronweasley", direction = -1, name = "storage temperature (°C)")+
#   theme_bw() 
# 
# 
# ggplot() +  facet_grid(perc_rot~temp, labeller = label_both) + ylab("infection probability (%)")+
#   geom_line(surv_gomp3, mapping = aes(x = time, y = 1-est, group=t7, colour=t7), size = 1)+
#   # scale_color_paletteer_c("harrypotter::ronweasley", direction = -1, name = "mean daily wetness duration (h)")+
#   scale_color_paletteer_c("pals::coolwarm", direction = 1, name = "temperature (°C)")+
#   theme_bw(base_size = 16) + theme(legend.position = "bottom")

