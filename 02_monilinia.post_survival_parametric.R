## ******************************************************************* ##
## ******************* MONILINIA SURVIVAL ANALYSIS ******************* ##
## ******************************************************************* ##

# in this version of the script I made the hp. 
# *** System REF (S1) 2018 was not taken in consideration (fongicide)

## **************
## load observed data + plots (check filter informations before run analysis) 
source("00_1_monilinia.post_data-retrieve.R")
## **************

## ******************************************************************* ##
##    PARAMETRIC REGRESSION 18-19 - best model from cross-validation 
## ******************************************************************* ##

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



### ************************************************************************* ##  
###                                 MODEL PLOT
### ************************************************************************* ##

model.surv.g = fit_gomp.bestmpr #flexsurvreg

## *********** 2018 *************
{
    new.data18 = mon.survival %>% 
      group_by( temp, cum_deg_day.scale, wd7, perc_rot, rain7, t7, year, Ha) %>%
      dplyr::summarize( weight_g = mean(weight_g, na.rm=T),
                        weight_g.scale = mean(weight_g.scale, na.rm=T) ) %>% ungroup() %>%
      filter(year == 2018)
    
      
    surv_gomp18 = summary(model.surv.g, newdata = new.data18, type = "survival", tidy = TRUE, t = c(0:max(mon.survival18$time)) ) 
      
    
    km_fit18 = survfit(Surv(time, status) ~ wd7 , data = mon.survival18, conf.int=0.95 )
    names(km_fit18$strata) <- gsub("wd7=", "",  c("0.56","3.09", "3.26") )
    
    km.plot18 = ggsurvplot(km_fit18, xlim =c(-0.3,(max(mon.survival18$time)+0.3)), conf.int = F, legend.title="Harvest date (dafb)",
                           palette = c("#009E73", "#E69F00", "#0072B2"),
                           linetype = "twodash",
                           legend="none")
    
    km.plot18$plot +
      ylab("Survival probability (%)")+
      xlab("Storage time (d)")+
      geom_line(surv_gomp18, mapping = aes(x = time, y = est, colour=as.factor(wd7)), size = 0.9, )+
      scale_y_continuous( limits = c(-0.01,1.02), expand = c(0,0), breaks = seq(0,1,0.25) )+
      scale_x_continuous( expand = c(0,0), breaks = seq(0,14,2) )+
      theme_bw(base_size = 16) + 
      theme(
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        legend.position = "bottom",
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
  km.plot18 = ggsurvplot(km_fit18, xlim =c(-0.3,(max(mon.survival18$time)+0.3)), conf.int = F, legend.title="Harvest date (dafb)",
                         palette = c("#009E73", "#E69F00", "#0072B2", 
                           "#009E73", "#009E73", "#E69F00","#E69F00", "#0072B2", "#0072B2"),
                         linetype = "twodash",
                         legend="none")
  
  
  km.plot18$plot + facet_grid(.~Irr)+
    ylab("Survival probability (%)")+
    xlab("Storage time (d)")+
    geom_line(surv_gomp18.test, mapping = aes(x = time, y = est, colour=as.factor(wd7)), size = 0.9, )+
    scale_y_continuous( limits = c(-0.01,1.02), expand = c(0,0), breaks = seq(0,1,0.25) )+
    scale_x_continuous( expand = c(0,0), breaks = seq(0,14,2) )+
    theme_bw(base_size = 16) + 
    theme(
      axis.text.x=element_text(colour="black"),
      axis.text.y=element_text(colour="black"),
      legend.position = "bottom",
      strip.background = element_blank(), strip.text.x = element_blank(), 
      panel.grid.major = element_blank(),  panel.grid.minor = element_blank())
  
}

## *********** 2019 *************
{
  
  new.data19 = mon.survival19 %>% 
    group_by( temp, cum_deg_day.scale, wd7, perc_rot, rain7, t7, year, Ha) %>%
    dplyr::summarize( weight_g = mean(weight_g, na.rm=T),
                      weight_g.scale = mean(weight_g.scale, na.rm=T) ) %>% ungroup() %>%
    filter(year == 2019)
  
  surv_gomp19 = summary(model.surv.g, newdata = new.data19, type = "survival", tidy = TRUE, t = c(0:max(mon.survival19$time))) 
  surv_gomp19.test = merge(surv_gomp19, new.data19 ) 
  
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
    ylab("Survival probability (%)")+
    xlab("Storage time (d)")+
    geom_line(surv_gomp19, mapping = aes(x = time, y = est, colour=as.factor(temp)), size = 0.9)+
    geom_text(x = 2, y = 0.1, aes(label = label), data = labels, fontface = "bold", size = 7)+
    scale_y_continuous( limits = c(-0.01,1.02), expand = c(0,0), breaks = seq(0,1,0.25) )+
    scale_x_continuous( expand = c(0,0),
                        breaks = seq(0,21,3) )+
    theme_bw(base_size = 16) + 
    theme(
      axis.text.x=element_text(colour="black"),
      axis.text.y=element_text(colour="black"),
      legend.position = "bottom",
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

