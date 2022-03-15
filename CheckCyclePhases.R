# Determine Phases - Keith (1990) delineation of phases as cited in Oli et al (2020)
# Jenilee Gobin
# 14 Mar 2022


# 1 - CHECK/LOAD PACKAGES ----- 

# function to look for & install multiple packages
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#  vector of needed packages
analysis.packages <- c("tidyverse", "viridis", "data.table", "cowplot", "ggplot2")

# apply function to packages
check.packages(analysis.packages)

# 2 - LOAD DATA -----

# get density data

kluane_density_dat <- fread(file = 'Kluane_density_data.csv')

kluane_density_dat <- kluane_density_dat %>% 
  mutate(year = c(1987:1996), 
         winter = year +1) %>% 
  relocate(year, .after = winter) 

# get difference between years (e.g., 1988 would be difference between 87-88 and 88-89) and express this as a proportion relative to the first year (i.e., e.g., 88)
hare_phases_1 <- kluane_density_dat[1:nrow(kluane_density_dat)-1,] %>% 
  select(c(winter, hare_density)) %>% 
  cbind(kluane_density_dat$hare_density[2:nrow(kluane_density_dat)], diff(kluane_density_dat$hare_density)) %>% 
  rename(hare_density_t1 = hare_density, hare_density_t2 = V2) %>% 
  mutate(prop = abs(hare_density_t2/hare_density_t1), 
         phase = case_when(prop < 0.44 ~ 'decline', 
                           prop > 1.89 ~ 'increase', 
                           prop > 0.43 & prop < 1.90 ~ 'peak/low')
  )
  
  
  
 





