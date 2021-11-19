## ******************************************************************* ##
## ******************* BROWN ROT PROGRESS PLOT *********************** ##
## ******************************************************************* ##

## **************
## load observed data + plots (check filter informations before run analysis) 
source("00_1_monilinia.post_data-retrieve.R")
## **************

# for plot 
max.time= 22 # max storage time (days)  
max.incidence = 100  # %


## ***** 2018
{
mon18.db.plot = monilinia18 %>%
  group_by(harvest_n, harvest_date, dah, Irr, FL, Syst, Exp, dah.int) %>%
  dplyr::summarize( incidence = round(sum( infection=="1" )*100 /( length(infection)),2), 
                    weight_mean = round(mean(weight_g, na.rm=T),2),  weight_sd = sd(weight_g, na.rm=T)) %>% # cumulated incidence 
  ungroup() %>%
  group_by(harvest_n, harvest_date, Irr, FL, Syst, Exp) %>%
  mutate(max.dah = max(dah.int) ) %>% ungroup() %>%
  filter(Syst == "S2")


mon18.db.box = monilinia18 %>%
  group_by(harvest_n, harvest_date, dah, Irr, FL, Syst, Exp, dah.int, Box) %>%
  dplyr::summarize( incidence = round(sum( infection=="1" )*100 /( length(infection)),2), 
                    weight_mean = round(mean(weight_g, na.rm=T),2),  weight_sd = sd(weight_g, na.rm=T)) %>% # cumulated incidence 
  ungroup() %>%
  group_by(harvest_n, harvest_date, Irr, FL, Syst, Exp, Box) %>%
  mutate(max.dah = max(dah.int), 
         dafbHarvest = as.numeric( format(as.Date(harvest_date,tryFormats = "%d/%m/%Y"),"%j") ) - date_bloom2018_julian) %>% ungroup() %>%
  filter(Syst == "S2")

mon18.db.box %>% group_by(harvest_n, harvest_date, dah, Irr, FL, Syst, Exp, dah.int) %>%
  dplyr::summarize( i_sd = sd(incidence, na.rm=T) ) %>%
  filter(dah.int==14 | dah.int==12)


## ************************* cumulative incidence ~ irrigation + harvest date 
ggplot(mon18.db.plot) +
  stat_summary( aes(dah.int, incidence, colour=harvest_n ), fun.data = "mean_cl_boot", geom = "line")+
  stat_summary( aes(dah.int, incidence, colour=harvest_n), fun.data = "mean_cl_boot", position=position_dodge(width=0.5))+
  facet_grid(.~Irr, labeller = label_both)+
  scale_color_grey(name="Harvest date", start=0.6, end=0.2)+
  scale_linetype_discrete(name="Harvest date")+
  xlab("Time (d)")+
  ylab("Brown rot incidence (%)")+
  xlim(c(0,max.time)) + ylim(c(0,max.incidence)) +
  theme_bw(base_size = 14)+  theme(legend.position="bottom")


## ************************* cumulative incidence ~ irrigation + fruit load + harvest date 
ggplot(mon18.db.box) +
  stat_summary( aes(dah.int, incidence, colour=FL, linetype=FL ), fun.data = "mean_sd", geom = "line")+
  stat_summary( aes(dah.int, incidence, colour=FL), fun.data = "mean_sd", position=position_dodge(width=0.5))+
  facet_grid(harvest_n~Irr, labeller = label_both)+
  scale_color_grey(name="Fruit Load", start=0.6, end=0.2)+
  scale_linetype_discrete(name="Fruit Load")+
  xlab("Time (d)")+
  ylab("Brown rot incidence (%)")+
  xlim(c(0,max.time)) + ylim(c(0,max.incidence)) +
  theme_bw(base_size = 14)+  theme(legend.position="bottom")



labels <- data.frame(Irr = c("100%", "50%"), label = c("a", "b"))
ggplot(mon18.db.box) +
  stat_summary( aes(dah.int, incidence, colour=as.factor(dafbHarvest) ), fun.data = "mean_cl_boot", geom = "line", linetype='dotted')+
  stat_summary(  aes(dah.int, incidence, colour=as.factor(dafbHarvest) ), fun.data = "mean_cl_boot")+
  facet_grid(.~Irr)+
  # geom_text(x = min(mon18.db.box$dah.int)+0.4, y = max(mon18.db.box$incidence), aes(label = label), data = labels, fontface = "bold", size=5)+
  scale_linetype_discrete(name="Harvest date")+
  scale_color_manual(values=c("#009E73", "#E69F00", "#0072B2"), name="Harvest date (dafb)")+
  xlab("Storage time (d)")+
  ylab("Brown rot cumulative incidence (%)")+
  scale_x_continuous( limits = c(-0.5,15), expand = c(0,0),
                      breaks = seq(0,14,2) )+
  theme_bw(base_size = 14)+  theme(legend.position="bottom",
                                   axis.text.x=element_text(colour="black"),
                                   axis.text.y=element_text(colour="black"),
                                   panel.grid.major = element_blank(),
                                   strip.background = element_blank(),
                                         strip.text.x = element_blank())
  

}

## ***** 2019
{
mon19.db.plot = monilinia19 %>%
  group_by(harvest_n, dah, temp, irr, mod, dah.int) %>%
  dplyr::summarize( incidence = sum(  infection=="1" | infection=="2", na.rm = T )*100 /( length(infection)),
                    weight_mean = mean(weight_g, na.rm=T), weight_sd = sd(weight_g, na.rm=T) ) %>% # cumulated incidence 
  ungroup()%>%
  group_by(harvest_n, irr, temp, mod) %>%
  mutate(max.dah = max(dah.int),
         temp.int = as.numeric(as.character(substr(temp,2,3))) ) %>% ungroup() 

mon19.db.box = monilinia19 %>%
  group_by(harvest_n, dah, temp, irr, mod, dah.int, box, date_dafb) %>%
  dplyr::summarize( incidence = sum(  infection=="1" | infection=="2", na.rm = T )*100 /( length(infection)), 
                    weight_mean = mean(weight_g, na.rm=T), weight_sd = sd(weight_g, na.rm=T) ) %>% # cumulated incidence 
  ungroup()%>%
  group_by(harvest_n, irr, temp, mod, box) %>%
  mutate(max.dah = max(dah.int),
         temp.int = as.numeric(as.character(substr(temp,2,3))), 
         dateHarvest=date_dafb, date_dafb=NULL) %>% ungroup() 
mon19.db.box %>% group_by(harvest_n, dah, temp, irr, mod, dah.int) %>%
dplyr::summarize( i_sd = sd(incidence, na.rm=T) ) %>%
  filter(dah.int==17 | dah.int==21)


ggplot(mon19.db.plot) +
  geom_line(  aes(dah.int, incidence, colour=as.factor(temp.int)) )+
  geom_point( aes(dah.int, incidence, colour=as.factor(temp.int), shape=as.factor(temp.int)), size=1.8 )+
  facet_grid(harvest_n~irr, labeller = label_both)+
  xlab("Time (d)")+
  ylab("Brown rot incidence (%)")+
  scale_color_grey(name="Storage temperature (°C)", start=0.75, end=0.3)+
  scale_shape_discrete(name="Storage temperature (°C)")+
  # scale_color_manual(values=c("#003399","#FFCC00","#FF6600"), name="Temperature (°C)")+
  xlim(c(0,max.time)) + ylim(c(0,max.incidence)) +
  theme_bw(base_size = 14)+  theme(legend.position="bottom")

ggplot(mon19.db.plot) + 
  geom_line(  aes(dah.int, incidence, colour=as.factor(temp.int),linetype=harvest_n ), size=0.8 )+
  geom_point( aes(dah.int, incidence, colour=as.factor(temp.int), shape=harvest_n), size=3 )+
  facet_grid(.~irr, labeller = label_both)+  
  xlab("Time (d)")+
  ylab("Brown rot incidence (%)")+
  scale_color_grey(name="Storage temperature (°C)", start=0.75, end=0.3)+
  scale_shape_discrete(name="harvest_n")+
  # scale_color_manual(values=c("#003399","#FFCC00","#FF6600"), name="Temperature (°C)")+
  xlim(c(0,max.time)) + ylim(c(0,max.incidence)) +
  theme_bw(base_size = 14)+  theme(legend.position="bottom")



labels <- data.frame(dateHarvest=c(134,134, 141,141), irr=c("100%","50%","100%","50%"), label = c("a", "b", "c", "d"))
ggplot(mon19.db.box) + 
  stat_summary(  aes(dah.int, incidence, colour=as.factor(temp.int) ), fun.data = "mean_cl_boot", geom='line', linetype='dotted' )+
  stat_summary(  aes(dah.int, incidence, colour=as.factor(temp.int) ), fun.data = "mean_cl_boot", position=position_jitter(width=0.3,height=0) )+
  # geom_jitter(  aes(dah.int, incidence, colour=as.factor(temp.int) ) )+
  facet_grid(dateHarvest~irr, labeller = label_both)+
  geom_text(x = min(mon19.db.box$dah.int)+1, y = max(mon19.db.box$incidence)-5, aes(label = label), data = labels, fontface = "bold", size = 5)+
  xlab("Storage time (d)")+
  ylab("Brown rot cumulative incidence (%)")+
  scale_color_manual(values=c("blue3","darkorange","firebrick3"), name="Temperature (°C)")+
  scale_x_continuous( limits = c(-0.5,22), expand = c(0,0),
                      breaks = seq(0,21,3) )+
  theme_bw(base_size = 14)+  theme(legend.position="bottom",
                                   axis.text.x=element_text(colour="black"),
                                   axis.text.y=element_text(colour="black"),
                                   panel.grid.major = element_blank(),
                                   strip.background = element_blank(),
                                         strip.text.x = element_blank(),  strip.text.y = element_blank())
 

## ************************* harvest 
ggplot(mon19.db.plot) + 
  stat_summary(aes(dah.int, incidence, colour=as.factor(temp.int)), fun.data= mean_se, geom = "line" )+
  # stat_summary(aes(dah.int, incidence, colour=as.factor(temp.int)), fun.data = mean_se)+  
  geom_point(aes(dah.int, incidence, colour=as.factor(temp.int)))+  
  
  facet_grid(.~harvest_n, labeller = label_both)+  
  xlab("Time (d)")+
  ylab("Brown rot incidence (%)")+
  scale_color_grey(name="Storage temperature (°C)", start=0.75, end=0.3)+
  scale_shape_discrete(name="harvest_n")+
  # scale_color_manual(values=c("#003399","#FFCC00","#FF6600"), name="Temperature (°C)")+
  xlim(c(0,max.time)) + ylim(c(0,max.incidence)) +
  theme_bw(base_size = 14)+  theme(legend.position="bottom")

}
