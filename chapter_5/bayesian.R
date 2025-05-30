library(lattice)
library(nlme)
library(BayesRegMicrobiome)
library(tibble)
library(DESeq2)
library(Matrix)
library(huge)

source("func2.R")
source("initial_param0.R")

path = paste0(getwd(),"/",nsubj,"_",ntaxa,"/")
path

# install.packages("BayesRegMicrobiome_1.0.tar.gz", repos = NULL, type="source")
####################################################################
data	  =   readRDS(paste0(path,"otu_meta_list_withzi_taxa.rds"))
################################################################
# cc	=   commandArgs(trailingOnly  = TRUE)
# i	=   as.integer(cc[1])
i =1 
dd	=   data[[i]]
#######################################################
countdata  =  as.matrix(dd$countdata)
met_dd     =   dd$met_data

time    =  1:nsubj 
#rep(1,nsubj); 
ncovariates =  1
hyper <- fn.hyper(countdata, nsubj, time, max(time), ncovariates, ntaxa) 

met_dd = data.frame(index = met_dd$subject,time =  time,
                    group = as.numeric(factor(met_dd$group)) )

Miss_ind = data.frame(index = 1:nsubj, time = 1:nsubj, group = 0)
#######################################################
Dat =  list(Y= countdata, Tot_N = nsubj, Z_1= met_dd, 
            P = ncovariates, J=ntaxa, n = nsubj)


###  MCMC parameters
NN_Burn <- 10  ## No of iterations for Burn-in period
NN_sam <- 100  ## Posterior sample size

#debug(fn.run.MCMC.1)
MCMC_sam <- Bayes_Reg_Microbime(hyper, Dat, NN_Burn, NN_sam)

ini_sam <- fn.initialize(hyper, Dat$Y, Dat$Z_1, Dat$Tot_N, 
                         Dat$P, Dat$J)

names(Dat)
SIM_dat$P
names(SIM_dat)
SIM_dat$Z_1[,2] = 1
hyper <- fn.hyper(SIM_dat$Y, SIM_dat$Tot_N, SIM_dat$Z_1[,2], 
                  max(SIM_dat$Z_1[,2]), SIM_dat$P, SIM_dat$J)

names(hyper)



#200_600
#57 us
#rr 12,35,240,316,369