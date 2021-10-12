## ******************************************************************* ##
## ******************* MONILINIA SURVIVAL ANALYSIS ******************* ##
## ******************************************************************* ##

# in this version of the script I made the hp. 
# *** Modalité REF 2018 was not taken in consideratation (fungicide)

## **************
## load observed data + plots (check filter informations before run analysis) 
source("00_monilinia.post_data-retrieve.R")
## **************

mon.survival$rain7.bool = as.factor(1*(mon.survival$rain7>0))


# ### **** https://link.springer.com/article/10.1007/s11157-017-9443-0
# 
# # CTMI model - Lobry 1991
# 
# # Tmin = -2
# # Tmax = 38
# # Topt = 23
# 
# Tmin = -9.6 # calibrated
# Tmax = 35.4
# Topt = 25.4
# 
# mon.survival$temp.bell = ( (mon.survival$temp - Tmax) * ((mon.survival$temp-Tmin)^2) ) /
#                          ( (Topt-Tmin) * ( ((Topt-Tmin)*(mon.survival$temp-Topt)) - ((Topt-Tmax)*(Topt+Tmin-2*mon.survival$temp)) ) ) # CTMI
# 
# mon.survival$temp.bell = ((Tmax-mon.survival$temp)/(Tmax - Topt))*(((mon.survival$temp-Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt))) # Yan
# mon.survival$temp.bell=round(mon.survival$temp.bell,2)
# 
# unique(mon.survival$temp.bell)



## ********************* models to test
cov.string = c("wd7", "weight_g.scale", "temp", "rain7.bool", "perc_rot", "cum_deg_day.scale", "t7", "tmin7", "tmax7")
n_cov = length(cov.string)

id_cov = unlist(  lapply(1:n_cov, function(i) combn(1:n_cov, i, simplify=FALSE) ), recursive=FALSE )

formula_cov = sapply(id_cov, function(i)   paste("Surv(time, status) ~",paste(cov.string[i],collapse="+")) )
formula_cov.mpr = sapply(id_cov, function(i)   paste("Surv(time, status) ~ list(~",paste(cov.string[i],collapse="+"),",~1)") )
## ****************



## **************** cross-validation folds
 k.cv = 6
 
 # cv.folds = createFolds(as.factor(mon.survival$temp), k = k.cv) # create index for crossvalidation
 cv.folds = createFolds(as.factor(mon.survival$status), k = k.cv) 
 
 
 df.train = lapply(X = 1:k.cv, FUN = function(X) mon.survival[-cv.folds[[X]], ])
 df.test  = lapply(X = 1:k.cv, FUN = function(X) mon.survival[ cv.folds[[X]], ])
## **************** ##
 
 
 # **************** parallel calibration function
 funct_cv.gompertz = function(x, k.cv){  
   
   # x = 1 ## test
   
   ## ******************************************************************* ##
   ##                             AIC - BIC
   ## ******************************************************************* ##
   mod.mpr = mpr( as.formula(formula_cov.mpr[x]),
                       data = mon.survival,
                       family = "Gompertz" ,
                  hessian=F); 
   
   
   
   npar  = mod.mpr$model$npar
   model = mod.mpr$model$family
   
   aic     = summary(mod.mpr)$model["aic"]
   bic     = summary(mod.mpr)$model["bic"]
   loglike = summary(mod.mpr)$model["loglike"]
   
   ## ******************************************************************* ##
   ##                Cross-Validation
   ## ******************************************************************* ##
   

   # ConcordanceUno = NULL
   Dxy         = NULL
   Concordance = NULL
   
   for (i in 1:k.cv) {
     
     # i=1 ## test
     
     train.fit = flexsurvreg( as.formula(formula_cov[x])  ,
     # train.fit = flexsurvreg( Surv(time, status) ~ temp + weight_g.scale + wd7 + perc_rot   ,
                              data = df.train[[i]],
                              # anc = list(shape = ~ 1 ),
                              dist = "gompertz") 
     

     lpnew.calc = summary(train.fit, newdata = df.test[[i]], type = "median", tidy = TRUE) %>% mutate(epx.lp = exp(est))
     lpnew = lpnew.calc$est
     
     Surv.rsp.new = Surv( df.test[[i]]$time,  df.test[[i]]$status)
     
     # ConcordanceUno[i] = UnoC(Surv.rsp, Surv.rsp.new, lpnew); print(ConcordanceUno[i]) 
     rcorr = rcorr.cens(lpnew, Surv.rsp.new)
     Dxy[i]          = rcorr["Dxy"];     # print(Dxy[i])
     Concordance[i]  = rcorr["C Index"]; # print(Concordance[i])
     
    }
   
   db.res = data.frame(model = model,
                       mod.formula = formula_cov[x],
                       n.par = npar,
                       aic = aic, bic = bic, loglike = loglike,
                       D.xy_somers = mean(Dxy, na.rm=T),
                       D.xy_sd =       sd(Dxy, na.rm=T),
                       C.index = mean(Concordance, na.rm=T),
                       k.crossval = k.cv ) 
   
   return(db.res)
   
 }
 funct_cv.weibull  = function(x, k.cv){  
   
   x = 2 ## test
   
   ## ******************************************************************* ##
   ##                             AIC - BIC
   ## ******************************************************************* ##
   mod.mpr = mpr( as.formula(formula_cov.mpr[x]),
                  data = mon.survival,
                  family = "WeibullAFT"); 
   
   npar  = mod.mpr$model$npar
   model = mod.mpr$model$family
   
   n_coef = length(coef(mod.mpr)$beta) + 1
   
   aic = summary(mod.mpr)$model["aic"]
   bic = summary(mod.mpr)$model["bic"]
   loglike = summary(mod.mpr)$model["loglike"]
   
   ## ******************************************************************* ##
   ##                Cross-Validation
   ## ******************************************************************* ##
   
   
   # ConcordanceUno = NULL
   Dxy         = NULL
   Concordance = NULL
   
   for (i in 1:k.cv) {
     
     # i=1 ## test
     
     train.fit = flexsurvreg( as.formula(formula_cov[x])  ,
                              # train.fit = flexsurvreg( Surv(time, status) ~ temp + weight_g.scale + wd7 + perc_rot   ,
                              data = df.train[[i]],
                              # anc = list(shape = ~ 1 ),
                              dist = "weibull") 
     
     
     lpnew.calc = summary(train.fit, newdata = df.test[[i]], type = "median", tidy = TRUE) %>% mutate(epx.lp = exp(est))
     lpnew = lpnew.calc$est
     
     # Surv.rsp     = Surv(df.train[[i]]$time, df.train[[i]]$status)
     Surv.rsp.new = Surv( df.test[[i]]$time,  df.test[[i]]$status)
     
     # ConcordanceUno[i] = UnoC(Surv.rsp, Surv.rsp.new, lpnew); print(ConcordanceUno[i])
     rcorr = rcorr.cens(lpnew, Surv.rsp.new)
     Dxy[i]          = rcorr["Dxy"];     # print(Dxy[i])
     Concordance[i]  = rcorr["C Index"]; # print(Concordance[i])
     
   }
   
   db.res = data.frame(model = model,
                       mod.formula = formula_cov[x],
                       n.par = npar,
                       aic = aic, bic = bic, loglike = loglike,
                       D.xy_somers = mean(Dxy, na.rm=T),
                       D.xy_sd =       sd(Dxy, na.rm=T),
                       C.index = mean(Concordance, na.rm=T),
                       k.crossval = k.cv ) 
   
   return(db.res)
   
 }
 ## **************** ## 
 
 
 # res.cv = mclapply(X = 1:2, funct_cv.gompertz, mc.cores = 4, k.cv = k.cv )
 # res.cv.list.gomp = mclapply(X = 1:length(formula_cov), funct_cv.gompertz, mc.cores = 4, k.cv = k.cv )
 # res.cv.list.weib = mclapply(X = 1:length(formula_cov), funct_cv.weibull,  mc.cores = 4, k.cv = k.cv )
 # 
 # results.crossvalidation_gompertz.db = do.call(rbind.data.frame, res.cv.list.gomp)
 # results.crossvalidation_weibull.db  = do.call(rbind.data.frame,  res.cv.list.weib)
 # results.crossvalidation.db = rbind(results.crossvalidation_gompertz.db, results.crossvalidation_weibull.db )
 # save(results.crossvalidation.db, file = paste0("results.crossvalidation.db_", Sys.Date(),"_k",k.cv,"_median.RData") )

 load("results.crossvalidation.db_2020-11-24_k6_median.RData")

 results.crossvalidation.db$AICsurv=results.crossvalidation.db$aic + ((2*((results.crossvalidation.db$n.par-2)+2)*((results.crossvalidation.db$n.par-2)+3))/(1560-(results.crossvalidation.db$n.par-2)-3))
 
 utopia = c(bic = min(results.crossvalidation.db$bic) , D.xy_somers = min(results.crossvalidation.db$D.xy_somers) ) 
 
 # euclidean distance between the Pareto solutions and the utopian point
 euclidean = sqrt( (utopia[1] - results.crossvalidation.db$bic)^2 + (utopia[2] - results.crossvalidation.db$D.xy_somers)^2 )
 
 results.crossvalidation.db.utopia = results.crossvalidation.db %>%
   mutate( obj_euclidean_dist = euclidean   )
 
 par_opt_utopia = results.crossvalidation.db.utopia[ with(results.crossvalidation.db.utopia, which(obj_euclidean_dist==min(obj_euclidean_dist) ) ),]
 
 library(rPref)
 a = psel(results.crossvalidation.db, low(bic) * high(D.xy_somers) )  
 
 ggplot(results.crossvalidation.db) +
   geom_point(aes(bic, D.xy_somers)) +
   xlim(c(3200,4000))+
   theme_bw()
 
 # results.crossvalidation.db.order = results.crossvalidation.db[order(results.crossvalidation.db$loglike),]%>%
   # mutate( id = as.numeric(rownames(results.crossvalidation.db.order)) )
 
 # ggplot(results.crossvalidation.db.order) +
   # geom_point( aes(reorder(id, loglike), -loglike) ) +
   # theme_bw()
 
 
 ### ---------------------- ##
 ### plot training vs. predictive error
 ### ---------------------- ##
 
 # attach(results.crossvalidation.db.order)
 # 
 # ## add extra space to right margin of plot within frame
 # par(mar=c(5, 4, 4, 6) + 0.1)
 # 
 # ## Plot first set of data and draw its axis
 # plot(reorder(id, loglike), -loglike, pch=16, axes=FALSE, ylim=c(1500,2200), xlab="", ylab="", 
 #      type="b",col="black", main=" Performance")
 # axis(2, ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
 # mtext("- Loglike",side=2,line=2.5)
 # box()
 # 
 # ## Allow a second plot on the same graph
 # par(new=TRUE)
 # # dev.off()
 # 
 # 
 # ## Plot the second plot and put axis scale on right
 # plot(reorder(id, loglike), D.xy_somers, pch=99,  xlab="", ylab="", ylim=c(0.4,0.7), 
 #      axes=FALSE, type="p", col="yellow")
 # 
 # ## a little farther out (line=4) to make room for labels
 # mtext("D.xy Somers", side=4, col="red",line=4) 
 # axis(4, ylim=c(1500,2200), col="red",col.axis="red",las=1)
 # 
 # ## Draw the time axis
 # axis(1,pretty(range(id),10))
 # mtext("id model",side=1,col="black",line=2.5)  
 # 
 # ## Add Legend
 # legend("topleft",legend=c("- Loglike","id model"),
 #        text.col=c("black","red"),pch=c(16,15),col=c("black","red"))
 # 
 # detach(results.crossvalidation.db.order)
