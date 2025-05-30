library(dplyr)
library(gllvm)
data("microbialdata")
Ysoil <- microbialdata$Y
Xenv <- microbialdata$Xenv
dim(Ysoil)
names(Xenv)
unique(Xenv$Site)
unique(Xenv$Region)
unique(Xenv$Soiltype)
View(Xenv)



# Using dplyr's filter function
subset_data <- Xenv %>% filter(Region == "Aus", Site == "A", Soiltype == "B")
print(subset_data)

#"SOM"       "pH"       "Phosp"  fixed effects

# "Region"   "Site=sample=independen"     "Soiltype=group"
#   rr(0 + taxon | Site,2) counts with site are correlated
#   us(1 + Soiltype|taxon) 

meanY <- apply(Ysoil,2, mean)
varY <- apply(Ysoil,2, var)
plot(log(meanY),varY, log = "y", main = "Species mean-variance relationship")

View(Xenv)
#'  For easiness to compare we just choose a simple model. 
#'  glmmTMB can take on the more complex model but the other models might 
#'  not be able to take them on. Thefore we choose simple models. 
#' We consider 
us(1 + group|taxon) +  
  rr(0 + taxon | subject,2)