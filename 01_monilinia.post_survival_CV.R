## ******************************************************************* ##
## ******************* MONILINIA SURVIVAL ANALYSIS ******************* ##
## ******************************************************************* ##

## **************
## load observed data 
source("00_1_monilinia.post_data-retrieve.R")
## **************


## ********************* models to test
cov.string = c("wd7", "weight_g.scale", "temp", "rain7", "perc_rot", "cum_deg_day.scale", "t7")
n_cov = length(cov.string)

id_cov = unlist(  lapply(1:n_cov, function(i) combn(1:n_cov, i, simplify=FALSE) ), recursive=FALSE )

formula_cov = sapply(id_cov, function(i)   paste("Surv(time, status) ~",paste(cov.string[i],collapse="+")) )
formula_cov.mpr = sapply(id_cov, function(i)   paste("Surv(time, status) ~ list(~",paste(cov.string[i],collapse="+"),",~1)") )
## ****************


## **************** CREATE CROSS-VALIDATION FOLDS
 k.cv = 6
 
 cv.folds = createFolds(as.factor(mon.survival$status), k = k.cv) # create index for crossvalidation
 
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
   Dxy         = NULL # Somer's rank correlation for censored response variable
   Concordance = NULL # c index  
   
   for (i in 1:k.cv) {
     
     # i=1 ## test
     
     train.fit = flexsurvreg( as.formula(formula_cov[x])  ,
                              data = df.train[[i]],
                              dist = "gompertz") 
     
     lpnew.calc = summary(train.fit, newdata = df.test[[i]], type = "median", tidy = TRUE) %>% mutate(epx.lp = exp(est))
     lpnew = lpnew.calc$est
     
     Surv.rsp.new = Surv( df.test[[i]]$time,  df.test[[i]]$status)
     
     # ConcordanceUno[i] = UnoC(Surv.rsp, Surv.rsp.new, lpnew); print(ConcordanceUno[i]) 
     rcorr           = rcorr.cens(lpnew, Surv.rsp.new)
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
 
 
 res.cv = mclapply(X = 1:2, funct_cv.gompertz, mc.cores = 4, k.cv = k.cv )
 results.crossvalidation.db = mclapply(X = 1:length(formula_cov), funct_cv.gompertz, mc.cores = 4, k.cv = k.cv )

 save(results.crossvalidation.db, file = paste0("results.crossvalidation.db_", Sys.Date(),"_k",k.cv,"_median.RData") )
 # load("results.crossvalidation.db_............_median.RData")
 
 
 results.crossvalidation.db = results.crossvalidation_gompertz.db
 
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
