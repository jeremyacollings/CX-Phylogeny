
########## FITTING PATH ANALYSES ##########

set.seed(6)

library(tidyverse)
library(tidybayes)
library(brms)
library(GGally)

source("1_Data.Prep.R")

dat2 <- dat %>%
  filter(competitor != "MGLO")


# Calculate Coexistence Metrics -------------------------------------------

mod <- readRDS(file.path("RDS_Files","rick.RDS"))
spp <- readRDS(file.path("RDS_Files","spp_lists.RDS"))

# extract paramter values

lam.df <-  mod %>% 
  spread_draws(lambda[s])

alpha.df <- mod %>%
  spread_draws(alpha[s,c])

N.species <- n_distinct(lam.df$s)
N.draws <- n_distinct(lam.df$.draw)

lambdas <- lam.df %>%
  ungroup() %>%
  select(lambda) %>% unlist() %>%
  matrix(ncol = N.species, nrow = N.draws)

alphas <- alpha.df %>%
  ungroup() %>%
  select(alpha) %>% unlist() %>%
  array(dim = c(N.draws, N.species, N.species))

# calculate metrics

ND <- CR1 <- CR2 <- CR <- DR1 <- DR2 <- DR <- FI <- array(NA, dim = dim(alphas))
CX <- winners <- array(NA, dim = dim(alphas))
for(i in 1:N.species){
  for(j in 1:N.species){
    CR1[,i,j] <- sqrt((alphas[,j,i]*alphas[,j,j])/(alphas[,i,i]*alphas[,i,j]))
    CR2[,i,j] <- sqrt((alphas[,i,i]*alphas[,i,j])/(alphas[,j,i]*alphas[,j,j]))
    DR1[,i,j] <- lambdas[,i]/lambdas[,j]
    DR2[,i,j] <- lambdas[,j]/lambdas[,i]
    FI[,i,j] <- ifelse(DR1[,i,j]*CR1[,i,j] > DR2[,i,j]*CR2[,i,j], 
                         DR1[,i,j]*CR1[,i,j], DR2[,i,j]*CR2[,i,j])
    CR[,i,j] <- ifelse(DR1[,i,j]*CR1[,i,j] > DR2[,i,j]*CR2[,i,j],
                         CR1[,i,j], CR2[,i,j])
    DR[,i,j] <- ifelse(DR1[,i,j]*CR1[,i,j] > DR2[,i,j]*CR2[,i,j],
                         DR1[,i,j], DR2[,i,j])
    ND[,i,j] <- 1 - sqrt((alphas[,i,j]*alphas[,j,i])/
                             (alphas[,i,i]*alphas[,j,j]))
    
    # need to write code to say 1. is ND = 1 & FI = 0... then NA
    # 2. do ND and FI suggest CX... then 1
    # 3. do ND and FI suggest EX... then 2
    # 4. do ND and FI suggest PE... then 3
    
    # then a seperate array of the winners... 
    # if first array != 2... then NA
    # if first array == 2...
    # who won? 1, 2, or 3? 
    
    CX[,i,j] <- ifelse(ND[,i,j] == 1 & FI[,i,j] == 0, NA, 
                         ifelse(1 - ND[,i,j] < FI[,i,j] &
                                  FI[,i,j] < 1 / (1 - ND[,i,j]), 1, # coexistence
                                ifelse(1 - ND[,i,j] < FI[,i,j] &
                                         FI[,i,j] > 1 / (1 - ND[,i,j]), 2, # exclusions
                                       ifelse(1 - ND[,i,j] > FI[,i,j] &
                                                FI[,i,j] > 1 / (1 - ND[,i,j]), 3, NA)))) # priority effects
    
    winners[,i,j] <- ifelse(CX[,i,j] %in% c(NA, 1, 3), NA, 
                              ifelse(CX[,i,j] == 2 &
                                       DR1[,i,j]*CR1[,i,j] > DR2[,i,j]*CR2[,i,j], i, 
                                     ifelse(CX[,i,j] == 2 &
                                              DR1[,i,j]*CR1[,i,j] < DR2[,i,j]*CR2[,i,j], j, NA)))
  }
}

# compile into single dataframe
traits_red <- traits[which(traits$sp != "MGLO"),]

get_distance <- function(x){
  mat <- outer(x, x, function(x, y) abs(x - y)) # get matrix of distances
  mat[lower.tri(mat)] # extract lower triangle of matrix
}

# through trait distances into df
traits_dist <- cbind.data.frame(SLA = get_distance(traits_red$SLA), 
                                RMF = get_distance(traits_red$RMF), 
                                height = get_distance(traits_red$height), 
                                LA = get_distance(traits_red$LA), 
                                LDMC = get_distance(traits_red$LDMC))

# add in medians and lower & upper quantiles of coexistence metrics
traits_dist$ND_med <- apply(ND, 2:3, median)[lower.tri(apply(ND, 2:3, median))]
traits_dist$ND_low <- apply(ND, 2:3, function(x) quantile(x, .025))[lower.tri(apply(ND, 2:3, median))]
traits_dist$ND_up <- apply(ND, 2:3, function(x) quantile(x, .975))[lower.tri(apply(ND, 2:3, median))]
traits_dist$CR_med <- apply(CR, 2:3, median)[lower.tri(apply(CR, 2:3, median))]
traits_dist$CR_low <- apply(CR, 2:3, function(x) quantile(x, .025))[lower.tri(apply(CR, 2:3, median))]
traits_dist$CR_up <- apply(CR, 2:3, function(x) quantile(x, .975))[lower.tri(apply(CR, 2:3, median))]
traits_dist$DR_med <- apply(DR, 2:3, median)[lower.tri(apply(DR, 2:3, median))]
traits_dist$DR_low <- apply(DR, 2:3, function(x) quantile(x, .025))[lower.tri(apply(DR, 2:3, median))]
traits_dist$DR_up <- apply(DR, 2:3, function(x) quantile(x, .975))[lower.tri(apply(DR, 2:3, median))]
traits_dist$FI_med <- apply(FI, 2:3, median)[lower.tri(apply(FI, 2:3, median))]
traits_dist$FI_low <- apply(FI, 2:3, function(x) quantile(x, .025))[lower.tri(apply(FI, 2:3, median))]
traits_dist$FI_up <- apply(FI, 2:3, function(x) quantile(x, .975))[lower.tri(apply(FI, 2:3, median))]
traits_dist$CX <- apply(CX, 2:3, function(x) sum(x == 1)/length(x))[lower.tri(apply(CX, 2:3, median))]

traits_dist$pair <- c(outer(spp[[1]][,2], spp[[2]][,2], function(x, y) paste(x, y, sep = "-"))[lower.tri(outer(spp[[1]][,2], spp[[2]][,2], function(x, y) paste(x, y, sep = "-")))])
traits_dist$dist <- dists$dist[match(traits_dist$pair, paste(dists$sp2, dists$sp1, sep = "-"))]

# haphazardly scaling everything for the sake of throwing into path analysis..
# probably don't wanna scale some of these things for viz or other sorts of models
traits_dist[,-19] <- lapply(traits_dist[,-19], function(x) c(scale(x)))

# Fit Path Analyses -------------------------------------------------------

# quick look at correlations
ggpairs(traits_dist[,c(1:6, 9, 12, 15, 18, 20)])

## Median Estimates -------------------------------------------------------

full_mod <- brm(bf(SLA ~ dist) + 
                  bf(height ~ dist) + 
                  bf(RMF ~ dist) + 
                  bf(LA ~ dist) + 
                  bf(LDMC ~ dist) + 
                  bf(ND_med ~ SLA + height + RMF + LA + LDMC) + 
                  bf(DR_med ~ SLA + height + RMF + LA + LDMC) + 
                  bf(CR_med ~ SLA + height + RMF + LA + LDMC) + 
                  set_rescor(FALSE), data = traits_dist,
                cores = 4)

## Lower Estimates --------------------------------------------------------

full_mod_low <- brm(bf(SLA ~ dist) + 
                  bf(height ~ dist) + 
                  bf(RMF ~ dist) + 
                  bf(LA ~ dist) + 
                  bf(LDMC ~ dist) + 
                  bf(ND_low ~ SLA + height + RMF + LA + LDMC) + 
                  bf(DR_low ~ SLA + height + RMF + LA + LDMC) + 
                  bf(CR_low ~ SLA + height + RMF + LA + LDMC) + 
                  set_rescor(FALSE), data = traits_dist,
                cores = 4)

## Upper Estimates --------------------------------------------------------

full_mod_up <- brm(bf(SLA ~ dist) + 
                      bf(height ~ dist) + 
                      bf(RMF ~ dist) + 
                      bf(LA ~ dist) + 
                      bf(LDMC ~ dist) + 
                      bf(ND_up ~ SLA + height + RMF + LA + LDMC) + 
                      bf(DR_up ~ SLA + height + RMF + LA + LDMC) + 
                      bf(CR_up ~ SLA + height + RMF + LA + LDMC) + 
                      set_rescor(FALSE), data = traits_dist,
                    cores = 4)
