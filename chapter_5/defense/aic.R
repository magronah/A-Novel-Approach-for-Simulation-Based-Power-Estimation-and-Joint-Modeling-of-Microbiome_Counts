path  =  paste0("defense/",nsubj,"_",ntaxa,"/")
delta  =  c()
for(i in 1:10){
  rr     =  readRDS(paste0(path,"/rr_mod",i,".rds"))
  us     =  readRDS(paste0(path,"/us_mod",i,".rds"))
  delta[i] =  AIC(rr) - AIC(us)
}
plot(delta)
