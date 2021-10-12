# Degree-days calculation 
## reference: "Zalom, F.G, P.B. Goodell, L.T. Wilson, W.W. Barnett, and W.J. Bentley, 1983.
##            Degree-days: The calculation and use of heat units in pest management"
#############################################################

# 1) methode 'single sine': it estimates the degrees cumulated in a certain period
#--------------------------
f_get_cdd <- function(Tu=35,Tl=7,Tmax,Tmin)
{
  degD <- vector()
  
  for (i in 1:(length(Tmax)-1) )
  {
    
	meanT <- (Tmax[i]+Tmin[i])/2  # mean temp DAILY VALUES
	alpha <- (Tmax[i]-Tmin[i])/2  # deviation from mean day temp
	
	if ( is.na(Tmin[i]) | is.na(Tmax[i]) ){ #0 (NA)

	  degD[i] <- NA
	  ###############
	  }else

	if ( Tmin[i]>=Tu ){ #I (max accumulate is 35-7=28) 
	  
	   degD[i] <- Tu-Tl 
	###############
	  }else            
	
	if ( Tmax[i]<=Tl ){ #II (min simulated)
	  
	  degD[i] <- 0 
	###############  
	  }else                         
	
	if ( Tmin[i] >=Tl & Tmax[i]<=Tu ){ #III
	  
	  degD[i] <- meanT - Tl 
	###############  
	  }else     
		
	if ( Tmin[i] < Tl & Tmax[i] <= Tu ) {  #IV  
		
		teta1 <- asin((Tl - meanT)/alpha)
		
		degD[i] <- (((meanT-Tl)*((pi/2)-teta1)) + (alpha*cos(teta1)))/pi
	###############	
	  }else
		
	if ( Tmin[i]>=Tl & Tmax[i]>Tu ) {  #V
	
		teta2 <- asin((Tu - meanT)/alpha)
		
		degD[i] <- (((meanT-Tl) * ((pi/2)+teta2)) + ((Tu-Tl)*((pi/2)-teta2)) -	(alpha*cos(teta2)))/pi
	###############	
		}else
	
	if ( Tmin[i]<Tl & Tmax[i]>Tu ) {  #VI

		teta1 <- asin((Tl - meanT)/alpha)
		teta2 <- asin((Tu - meanT)/alpha)
		
		degD[i] <- (((meanT-Tl)*(teta2-teta1)) +(alpha*(cos(teta1)-cos(teta2))) +((Tu-Tl)*((pi/2)-teta2)))/pi
	}
	
  }
	
	cumulateddegD = sum(degD, na.rm = T)
	
  return(cumulateddegD)
	
}
