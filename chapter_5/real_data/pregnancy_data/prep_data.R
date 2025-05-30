library(NBZIMM)

data(Romero)
names(Romero)

otu = Romero$OTU; dim(otu)
sam = Romero$SampleData; dim(sam)
colnames(sam)
names(sam)
dim(otu)
###########################################################
path   =   paste0(getwd(),"/real_data/")
source(paste0(path,"/fun.R"))
###########################################################
data("microbialdata")
data <- (microbialdata$Y); met_dd <- microbialdata$Xenv
dd  <- data[, apply(data, 2, function(col) any(col == 0))]
###########################################################
meta_dd    =   met_dd %>% 
  select(Site, Soiltype) %>% 
  rownames_to_column() %>% 
  setNames(c("sample","site","group")) 
#' group variable is the soiltype
#' subjects are the soil samples
##############################################################
dd_long    =   df_long(dd, otu_names = "OTU", subject_name = "sample")
dd_long    =   left_join(dd_long, meta_dd, by ="sample")  
###########################################################
pp         =    deseqfun(dd, meta_dd, alpha_level=0.1,
                         design   = ~group + site,
                         ref_name = "B")

countdata     =    pp$data$countdata
meta_dd       =    pp$data$meta_data
normalize_fac =    sizeFactors(pp$object)
meta_dd$normalizer  =   normalize_fac
###########################################################
normalizer =   data.frame(normalizer = normalize_fac) %>% 
  rownames_to_column("sample")

df    =   left_join(dd_long, normalizer, by ="sample")  
################################################################
gprior  <- data.frame(prior = "gamma(2, 2.5)",
                      class = "theta_sd",
                      coef = "")
 
