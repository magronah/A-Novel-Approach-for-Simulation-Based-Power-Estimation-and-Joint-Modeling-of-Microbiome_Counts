pkgs <- scan("reproducible/single_time/pkgs.txt", what = character(1),
             comment.char = "#")
lapply(pkgs, library, character.only = TRUE)

#setwd(here())
theme_set(theme_bw())


# num_cores = 5
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
