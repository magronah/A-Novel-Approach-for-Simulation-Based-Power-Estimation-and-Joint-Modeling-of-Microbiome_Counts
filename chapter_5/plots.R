
dd_sub_list  =  plt  = list()
numbers <- 1:50
group_size <- 10
number_list <- split(numbers, ceiling(seq_along(numbers) / group_size))
position_dodge_width <- 0.5

for(i in 1:length(number_list)){
  
  elem     =  number_list[[i]]  
  dd_sub   =  dd_full[dd_full$param_name %in% paste0("taxon",elem),]
  plt[[i]] =  ggplot(dd_sub, aes(x = param_name, y = true_param, group = model_type, color = model_type)) +
    geom_pointrange(aes(ymin = lwr, ymax = upr),
                    position = position_dodge(width = position_dodge_width)) +
    geom_point(position = position_dodge(width = position_dodge_width), color="black",
               size = 3) 
}

(plt[[1]]|plt[[2]])/(plt[[3]]|plt[[4]]) +  plot_layout(guides = "collect")


  
