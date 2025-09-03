
########## PLOT POPULATION PARAMETERS ##########

library(tidyverse)

rick <- readRDS(file.path("RDS_Files", "rick.rds"))
spp_list <- readRDS(file.path("RDS_Files", "spp_lists"))

# get into arrays

lams_df <- as.data.frame(rick)[,grepl("lam", names(as.data.frame(rick)))]
lams_array <- array(as.numeric(unlist(lams_df)), dim = c(nrow(lams_df), 13))

alphas_df <- as.data.frame(rick)[,grepl("alpha", names(as.data.frame(rick)))]
alphas_array <- array(as.numeric(unlist(alphas_df)), dim = c(nrow(lams_df), 13, 13))

# get medians 
apply(lams_array, 2, median)
apply(alphas_array, c(2,3), median)

ggplot(data = cbind.data.frame(med = apply(lams_array, 2, median), 
                               up = apply(lams_array, 2, quantile, 0.975), 
                               low = apply(lams_array, 2, quantile, 0.025), 
                               sp = spp_list[[1]][,2]), 
       aes(x = sp, y = med, ymin = low, ymax = up)) + 
  geom_pointrange() + ylab("λ")

# some of these look real small...

ggplot(data = cbind.data.frame(med = c(apply(alphas_array, c(2,3), median)), 
                               up = c(apply(alphas_array, c(2,3), quantile, 0.975)), 
                               low = c(apply(alphas_array, c(2,3), quantile, 0.975)), 
                               sp = rep(spp_list[[1]][,2], 13), 
                               comp = rep(spp_list[[2]][,2], each = 13)), 
       aes(x = paste(sp, comp, sep = "-"), y = med, ymin = low, ymax = up)) + 
  geom_pointrange() + ylab("α") + 
  geom_hline(yintercept = 0, linetype = "dashed")




