setwd("/project/6006158/agronahm/Michael-n-Ben-Repo/")
#source("reproducible/simp_single/Load_Packages.R")
library(tidyverse)

dd  =  readRDS(paste0(getwd(), "/reproducible/simp_single/data/wald_conf.rds"))
true_param  =  readRDS(paste0(getwd(), "/reproducible/simp_single/data/true_param.rds"))

dd  =   dd %>% arrange(sim)

df =  left_join(true_param,dd, by= "param_name")

cov_dd <- df %>%
  mutate(coverage = ifelse(true_param >= lwr & true_param <= upr, 1, 0)) %>%
  group_by(param_name) %>%
  summarise(cov = mean(coverage), true_param=true_param[1])


ggplot(cov_dd, aes(x=cov)) +
  geom_histogram()


dd_cov =  cov_dd %>%
         mutate(param_name = factor(param_name))  %>%
          arrange(true_param)

ggplot(dd_cov, aes(x=true_param,y=cov))+
  geom_point() +
  theme_bw()

ggplot(dd_cov, aes(x = param_name, y = cov)) +
  geom_point() +
  geom_line(aes(group = 1)) +  # Connect the points with lines
  geom_text(data = filter(dd_cov, cov > 90), aes(label = param_name), vjust = -0.5) +  # Add labels for taxa with cov > 90
  theme_bw() +
  labs(title = "Coverage by Parameter Name",
       x = "Parameter Name",
       y = "Coverage") +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(dd_cov, aes(x = reorder(param_name, cov), y = cov)) +
  geom_point() +
  # geom_text(aes(label = ifelse(cov > 90, as.character(param_name), "")), vjust = -0.5) +
 # coord_flip() +  
  theme_bw() +
  labs(title = "Coverage by Parameter Name",
       x = "Parameter Name",
       y = "Coverage") +
  theme(plot.title = element_text(hjust = 0.5))

median(dd_cov$cov)
####
##do power calculation and then compare the average 
##What are the ranges of mean abundance for those taxa with high coverage compared 
## with those with low coverage?

##SIMULATION PHASE
## Next steps
##' MSE, Power, Bias 
##' when will the reduced ranked model work best in com


library(tidyverse)

dd  =  readRDS(paste0(getwd(), "/reproducible/simp_single/data/wald_conf.rds"))
true_param  =  readRDS(paste0(getwd(), "/reproducible/simp_single/data/true_param.rds"))

dd <- dd %>% arrange(sim)
df <- left_join(true_param, dd, by = "param_name")

# Calculate coverage
cov_dd <- df %>%
  mutate(coverage = ifelse(true_param >= lwr & true_param <= upr, 1, 0)) %>%
  group_by(param_name) %>%
  summarise(cov = mean(coverage), true_param = true_param[1])

# Join cov_dd and df
plot_data <- left_join(cov_dd, df, by = c("param_name","true_param"))

# Filter cov < 0.1
plot_data <- plot_data %>% filter(cov < 0.1)



# Group by param_name and add ranking
plot_data <- plot_data %>%
  ungroup() %>%
  mutate(param_name = reorder(factor(param_name), true_param)) %>%
  group_by(param_name) %>%
  mutate(r = rank(est_param)) 

str(plot_data)
# Create the caterpillar plot
ggplot(plot_data, aes(x = r, y = est_param)) +
  geom_linerange(aes(ymin = lwr, ymax = upr)) +
  geom_hline(aes(yintercept = true_param), color="red") +
  geom_hline(yintercept = 0,lty=2) +
  geom_line() +
  facet_wrap(~param_name) +
  theme_minimal() 
######################################################


