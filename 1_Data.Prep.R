
########## PREPPING THE DATA ##########

library(tidyverse)

dat <- read.csv(file.path("Data","raw_data.csv"), header=TRUE,na.strings=c("","NA"))

trait_dat0 <- read.csv(file.path("Data", "trait_data_EC.csv"))

trait_dat <- trait_dat0 %>%
  filter(Region == "South") %>%
  select(unique_ID, gensp, focal_species, SLA, RMF, plant_height = plant_height_mm, 
         LA, LDMC)

traits <- group_by(trait_dat, focal_species)%>%
  summarise(SLA=mean(SLA, na.rm = TRUE),
            height=mean(plant_height, na.rm = TRUE),
            RMF = mean(RMF, na.rm = TRUE), 
            LA = mean(LA, na.rm = TRUE), 
            LDMC = mean(LDMC, na.rm = TRUE))%>%
  rename(sp = focal_species)

species.list <- unique(dat$focal_species)
for(s in species.list){
  # s<-species.list[1]
  temp0 <- dat %>% filter(focal_species==s, category %in% c("Intraspecific","Interspecific"))
  alone_data <- dat %>% filter(focal_species==s, category=="Alone")
  comp.list <- unique(temp0$competitor)
  for(s2 in comp.list) {
    temp1 <- mutate(alone_data, competitor = s2)
    dat <- rbind(dat, temp1)
  }}

dat <- dat %>% drop_na(competitor) %>%
  droplevels()

dat <- dat %>% mutate(fitness = ifelse(focal_species=="CG", Dry_weight, N_flowers))

dat <- dat %>% drop_na(fitness) %>% drop_na(N_neighbors)
dat <- dat %>% filter(focal_species != "MGLO" & 
                        focal_species != "PN" &
                        competitor != "PN")

dists <- read.csv(file.path("Data", "distances.csv"), header=TRUE,na.strings=c("","NA"))
dists <- dists[!duplicated(dists),]
