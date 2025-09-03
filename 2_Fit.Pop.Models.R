
########## FITTING POPULATION MODELS ##########

set.seed(6)

library(tidyverse)
library(rstan)
library(parallel)
library(tidybayes)

source("1_Data.Prep.R")

# Prep data for Stan ------------------------------------------------------

dat2 <- dat %>%
  filter(competitor != "MGLO")

# get focal species into numbers
foc_sp_list <- cbind.data.frame(num = 1:length(unique(dat2$focal_species)), 
                                sp = sort(unique(dat2$focal_species)))
comp_sp_list <- cbind.data.frame(num = 1:length(unique(dat2$competitor)), 
                                 sp = sort(unique(dat2$competitor)))
#idk y pn is in there
foc_sp_list <- foc_sp_list[which(foc_sp_list$sp != "PN"),]
comp_sp_list <- comp_sp_list[which(comp_sp_list$sp != "PN"),]

saveRDS(list(foc_sp_list, comp_sp_list), file.path("RDS_Files", "spp_lists.RDS"))

SLA_vec <- traits %>% filter(sp != "MGLO") %>% 
  select(SLA) %>% unlist() %>% unname() %>% scale() %>% as.numeric()

stan_dat <- list(N = nrow(dat2), 
                 FS = nrow(foc_sp_list),
                 CS = nrow(foc_sp_list),
                 sp = foc_sp_list$num[match(dat2$focal_species, foc_sp_list$sp)], 
                 comp = comp_sp_list$num[match(dat2$competitor, comp_sp_list$sp)], 
                 count = dat2$N_neighbors, 
                 fec = round(dat2$fitness), 
                 SLA = outer(SLA_vec, SLA_vec, FUN = function(x,y){abs(x-y)}))


# Fit model ---------------------------------------------------------------

rick <- stan(file="Ricker.stan", data=stan_dat, 
             cores = detectCores(), chains = 4, 
             pars = c("lambda", "alpha", "log_lik"))

saveRDS(rick, file.path("RDS_Files", "rick.RDS"))


