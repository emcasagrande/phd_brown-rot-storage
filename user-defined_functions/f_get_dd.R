# Degree-days calculation 
## reference: "Zalom, F.G, P.B. Goodell, L.T. Wilson, W.W. Barnett, and W.J. Bentley, 1983.
##            Degree-days: The calculation and use of heat units in pest management"
#############################################################

# 1) methode 'single sine': it estimates the degrees cumulated in a given day 
#--------------------------
f_get_dd <- function( Tu=35, Tl=7, Tmax, Tmin)
{
  
  # Tmax = T_max; Tmin = T_min
  
	meanT <- (Tmax+Tmin)/2  # mean day temp
	alpha <- (Tmax-Tmin)/2  # deviation from mean day temp
	
	if ( is.na(Tmin)==T | is.na(Tmax)==T  ){ #0 NA 
	  
	  0
	  ###############
	}else 
	
	if (Tmin>=Tu){ #I (max accumulate is 35-7=28) 
	  
	    Tu-Tl 
	###############
	  }else            
	
	if (Tmax<=Tl){ #II (min simulated)
	  
	    0 
	###############  
	  }else                         
	
	if (Tmin >=Tl & Tmax<=Tu){ #III
	  
	   meanT - Tl 
	###############  
	  }else     
		
	if (Tmin < Tl & Tmax <=Tu) {  #IV  
		
		teta1 <- asin((Tl - meanT)/alpha)
		
		(((meanT-Tl)*((pi/2)-teta1)) + (alpha*cos(teta1)))/pi
	###############	
	  }else
		
	if (Tmin>=Tl & Tmax>Tu) {  #V
	
		teta2 <- asin((Tu - meanT)/alpha)
		
		(((meanT-Tl) * ((pi/2)+teta2)) + ((Tu-Tl)*((pi/2)-teta2)) -	(alpha*cos(teta2)))/pi
	###############	
		}else
	
	if (Tmin<Tl & Tmax>Tu) {  #VI

		teta1 <- asin((Tl - meanT)/alpha)
		teta2 <- asin((Tu - meanT)/alpha)
		(((meanT-Tl)*(teta2-teta1)) +(alpha*(cos(teta1)-cos(teta2))) +((Tu-Tl)*((pi/2)-teta2)))/pi
	}

	
	}
