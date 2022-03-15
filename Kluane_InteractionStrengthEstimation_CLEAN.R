# Estimating food web interaction strengths using type II and best-fit functiona response relationships
# Jenilee Gobin
# 18 November 2020

# 1. LOAD LIBRARIES -----

check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#  vector of needed packages
analysis.packages <- c("tidyverse", "viridis", "data.table")

# apply function to packages
check.packages(analysis.packages)

# 2. LOAD DATA -----

# get density and parameter data

kluane_density_dat <- fread(file = 'Kluane_density_data.csv')

kluane_density_dat <- kluane_density_dat %>% 
  mutate (year = c(1987:1996)) %>% 
  relocate(year, .after = winter)

kluane_parameters <- fread(file = 'Kluane_FR_parameters.csv')

kluane_parameters_wide <- kluane_parameters %>% 
  dplyr::select (-source) %>% 
  pivot_wider(names_from = parameter, values_from = value)

# 3. ESTIMATE INTERACTION STRENGTHS -----

# start with type II functional response, which makes the form of the equation the same but employs different parameters for each predator and prey

calc_type_II_IS <- function(a, h, N){
  
  IS <- a*N/(1 + a*h*N)
  
}


type_II_parameters <- kluane_parameters_wide %>% 
  filter (FR_relationship == 'type_II')


hare_lynx_param <-type_II_parameters %>% 
  subset (prey =='hare' & predator == 'lynx')

rsq_lynx_param <-type_II_parameters %>% 
  subset (prey =='rsq' & predator == 'lynx')

hare_coyote_param <-type_II_parameters %>% 
  subset (prey =='hare' & predator == 'coyote')

rsq_coyote_param <-type_II_parameters %>% 
  subset (prey =='rsq' & predator == 'coyote')

hare_owl_param <-type_II_parameters %>% 
  subset (prey =='hare' & predator == 'owl')

get_IS_dat <- function(prey_name, pred_name, params, density_dat){
  
  IS_dat <- data.frame(year = kluane_density_dat$year,
                       prey = prey_name, 
                       predator = pred_name, 
                       IS = calc_type_II_IS(params$a, params$h, density_dat)
  )
  
  return (IS_dat)
  
}

hare_lynx_IS <- get_IS_dat('hare', 'lynx', hare_lynx_param, kluane_density_dat$hare_density)
rsq_lynx_IS <- get_IS_dat('rsq', 'lynx', rsq_lynx_param, kluane_density_dat$rsq_density)
hare_coyote_IS <- get_IS_dat('hare', 'coyote', hare_coyote_param, kluane_density_dat$hare_density)
rsq_coyote_IS <- get_IS_dat('rsq', 'coyote', rsq_coyote_param, kluane_density_dat$rsq_density)
hare_owl_IS <- get_IS_dat('hare', 'owl', hare_owl_param, kluane_density_dat$hare_density)

all_type_II_IS <- rbind(hare_lynx_IS, rsq_lynx_IS, hare_coyote_IS, rsq_coyote_IS, hare_owl_IS) %>% 
  rename(resource = prey, consumer = predator, interaction_strength = IS) 

all_type_II_IS$resource<-recode(all_type_II_IS$resource, 'rsq' = 'red_squirrel')

# adapt to best fit
# for hare_owl, this won't change so will use existing type_II output
# for hare_lynx, this will change to a type III (sigmoidal) prey-dependent FR

calc_type_III_IS <- function(a, h, N){
  
  IS <- a*N^2/(1 + a*h*N^2)
  
}

# rsq_coyote and hare_coyote are predator-dependent linear

calc_RD_linear_IS <- function(a, N, P, theta){
  
  IS <- a*N*P^-theta
  
}

# rsq_lynx is inverse sigmoidal

calc_inv_typeIII_IS <- function(b, c, N){
  
  IS <- c/(1+b*N^2)
  
}

best_fit_parameters <- kluane_parameters_wide %>% 
  filter (FR_relationship == 'best_fit')

hare_lynx_param_bestfit <- best_fit_parameters %>% 
  subset (prey == 'hare' & predator == 'lynx')

rsq_lynx_param_bestfit <-best_fit_parameters %>% 
  subset (prey =='rsq' & predator == 'lynx')

hare_coyote_param_bestfit <-best_fit_parameters %>% 
  subset (prey =='hare' & predator == 'coyote')

rsq_coyote_param_bestfit <-best_fit_parameters %>% 
  subset (prey =='rsq' & predator == 'coyote')

get_IS_dat_TypeIII <- function(prey_name, pred_name, params, density_dat){
  
  IS_dat_TypeIII <- data.frame(year = kluane_density_dat$year,
                       prey = prey_name, 
                       predator = pred_name, 
                       IS = calc_type_III_IS(params$a, params$h, density_dat)
  )
  
  return (IS_dat_TypeIII)
  
}

get_IS_dat_RD <- function(prey_name, pred_name, params, prey_density, pred_density){
  
  IS_dat_RD <- data.frame(year = kluane_density_dat$year,
                       prey = prey_name, 
                       predator = pred_name, 
                       IS = calc_RD_linear_IS(params$a, prey_density, pred_density, 1)
  )
  
  return (IS_dat_RD)
  
}

get_IS_dat_invTypeIII <- function(prey_name, pred_name, params, prey_density){
  
  IS_dat_invTypeIII <- data.frame(year = kluane_density_dat$year,
                          prey = prey_name, 
                          predator = pred_name, 
                          IS = calc_inv_typeIII_IS(params$b, params$c, prey_density)
  )
  
  return (IS_dat_invTypeIII)
  
}


hare_lynx_IS_bestfit <- get_IS_dat_TypeIII('hare', 'lynx', hare_lynx_param_bestfit, kluane_density_dat$hare_density)
rsq_lynx_IS_bestfit <- get_IS_dat_invTypeIII('rsq', 'lynx', rsq_lynx_param_bestfit, kluane_density_dat$hare_density)
hare_coyote_IS_bestfit <- get_IS_dat_RD('hare', 'coyote', hare_coyote_param_bestfit, kluane_density_dat$hare_density, kluane_density_dat$coyote_density)
rsq_coyote_IS_bestfit <- get_IS_dat_RD('rsq', 'coyote', rsq_coyote_param_bestfit, kluane_density_dat$rsq_density, kluane_density_dat$coyote_density)
hare_owl_IS_bestfit <- get_IS_dat('hare', 'owl', hare_owl_param, kluane_density_dat$hare_density)

all_bestfit_IS <- rbind(hare_lynx_IS_bestfit, rsq_lynx_IS_bestfit, hare_coyote_IS_bestfit, rsq_coyote_IS_bestfit, hare_owl_IS_bestfit) %>% 
  rename(resource = prey, consumer = predator, interaction_strength = IS) 

all_bestfit_IS$resource<-recode(all_bestfit_IS$resource, 'rsq' = 'red_squirrel')


# 4. OUTPUT TO SEPARATE FILES -----
# will output to separate files from each year for now that can can be copy/paste into appropriate folders

for (i in 1987:1996){
  
  filename1 <- paste('type_II', i, sep = '_')
  filename2 <- paste(filename1, '.csv', sep="")
  
  annual_dat <- subset(all_type_II_IS, year == i) %>% 
    dplyr::select (-year) %>% 
    write_csv(filename2)
  
}

for (i in 1987:1996){
  
  filename1 <- paste('bestfit', i, sep = '_')
  filename2 <- paste(filename1, '.csv', sep="")
  
  annual_dat <- subset(all_bestfit_IS, year == i) %>% 
    dplyr::select (-year) %>% 
    write_csv(filename2)
  
}

# 5. SENSITIVITY ANALYSIS -----
# estimate IS strengths separately for when hare density estimates are low and high

# hares
# type II low
hare_lynx_IS_TypeII_hare_low <- get_IS_dat('hare', 'lynx', hare_lynx_param, kluane_density_dat$hare_CI_low)
hare_coyote_IS_TypeII_hare_low <- get_IS_dat('hare', 'coyote', hare_coyote_param, kluane_density_dat$hare_CI_low)
hare_owl_IS_TypeII_hare_low <- get_IS_dat('hare', 'owl', hare_owl_param, kluane_density_dat$hare_CI_low)

all_type_II_IS_hare_low <- rbind(hare_lynx_IS_TypeII_hare_low, rsq_lynx_IS, hare_coyote_IS_TypeII_hare_low, rsq_coyote_IS, hare_owl_IS_TypeII_hare_low) %>% 
  rename(resource = prey, consumer = predator, interaction_strength = IS) 

all_type_II_IS_hare_low$resource<-recode(all_type_II_IS_hare_low$resource, 'rsq' = 'red_squirrel')

# hares
# type II high
hare_lynx_IS_TypeII_hare_high <- get_IS_dat('hare', 'lynx', hare_lynx_param, kluane_density_dat$hare_CI_high)
hare_coyote_IS_TypeII_hare_high <- get_IS_dat('hare', 'coyote', hare_coyote_param, kluane_density_dat$hare_CI_high)
hare_owl_IS_TypeII_hare_high <- get_IS_dat('hare', 'owl', hare_owl_param, kluane_density_dat$hare_CI_high)

all_type_II_IS_hare_high <- rbind(hare_lynx_IS_TypeII_hare_high, rsq_lynx_IS, hare_coyote_IS_TypeII_hare_high, rsq_coyote_IS, hare_owl_IS_TypeII_hare_high) %>% 
  rename(resource = prey, consumer = predator, interaction_strength = IS) 

all_type_II_IS_hare_high$resource<-recode(all_type_II_IS_hare_high$resource, 'rsq' = 'red_squirrel')

# next to do best fit... 

# hares
# best fit low
hare_lynx_IS_bestfit_hare_low <- get_IS_dat_TypeIII('hare', 'lynx', hare_lynx_param_bestfit, kluane_density_dat$hare_CI_low)
rsq_lynx_IS_bestfit_hare_low <- get_IS_dat_invTypeIII('rsq', 'lynx', rsq_lynx_param_bestfit, kluane_density_dat$hare_CI_low)
hare_coyote_IS_bestfit_hare_low <- get_IS_dat_RD('hare', 'coyote', hare_coyote_param_bestfit, kluane_density_dat$hare_CI_low, kluane_density_dat$coyote_density)
hare_owl_IS_bestfit_hare_low <- get_IS_dat('hare', 'owl', hare_owl_param, kluane_density_dat$hare_CI_low)

all_bestfit_IS_hare_low <- rbind(hare_lynx_IS_bestfit_hare_low, rsq_lynx_IS_bestfit, hare_coyote_IS_bestfit_hare_low, rsq_coyote_IS_bestfit, hare_owl_IS_bestfit_hare_low) %>% 
  rename(resource = prey, consumer = predator, interaction_strength = IS) 

all_bestfit_IS_hare_low$resource<-recode(all_bestfit_IS_hare_low$resource, 'rsq' = 'red_squirrel')

# hares
# best fit high
hare_lynx_IS_bestfit_hare_high <- get_IS_dat_TypeIII('hare', 'lynx', hare_lynx_param_bestfit, kluane_density_dat$hare_CI_high)
rsq_lynx_IS_bestfit_hare_high <- get_IS_dat_invTypeIII('rsq', 'lynx', rsq_lynx_param_bestfit, kluane_density_dat$hare_CI_high)
hare_coyote_IS_bestfit_hare_high <- get_IS_dat_RD('hare', 'coyote', hare_coyote_param_bestfit, kluane_density_dat$hare_CI_high, kluane_density_dat$coyote_density)
hare_owl_IS_bestfit_hare_high <- get_IS_dat('hare', 'owl', hare_owl_param, kluane_density_dat$hare_CI_high)

all_bestfit_IS_hare_high <- rbind(hare_lynx_IS_bestfit_hare_high, rsq_lynx_IS_bestfit, hare_coyote_IS_bestfit_hare_high, rsq_coyote_IS_bestfit, hare_owl_IS_bestfit_hare_high) %>% 
  rename(resource = prey, consumer = predator, interaction_strength = IS) 

all_bestfit_IS_hare_high$resource<-recode(all_bestfit_IS_hare_high$resource, 'rsq' = 'red_squirrel')

# output data for sensitivity analysis to separate files

# hares low
for (i in 1987:1996){
  
  filename3 <- paste('type_II_hare_low', i, sep = '_')
  filename4 <- paste(filename3, '.csv', sep="")
  
  annual_dat <- subset(all_type_II_IS_hare_low, year == i) %>% 
    dplyr::select (-year) %>% 
    write_csv(filename4)
  
}

for (i in 1987:1996){
  
  filename3 <- paste('bestfit_hare_low', i, sep = '_')
  filename4 <- paste(filename3, '.csv', sep="")
  
  annual_dat <- subset(all_bestfit_IS_hare_low, year == i) %>% 
    dplyr::select (-year) %>% 
    write_csv(filename4)
  
}

# hares high
for (i in 1987:1996){
  
  filename3 <- paste('type_II_hare_high', i, sep = '_')
  filename4 <- paste(filename3, '.csv', sep="")
  
  annual_dat <- subset(all_type_II_IS_hare_high, year == i) %>% 
    dplyr::select (-year) %>% 
    write_csv(filename4)
  
}

for (i in 1987:1996){
  
  filename3 <- paste('bestfit_hare_high', i, sep = '_')
  filename4 <- paste(filename3, '.csv', sep="")
  
  annual_dat <- subset(all_bestfit_IS_hare_high, year == i) %>% 
    dplyr::select (-year) %>% 
    write_csv(filename4)
  
}


# 6. PLOT INTERACTION STRENGTHS (FIG. A1) -----

# only need IS (per capita kill rate) for best fit


str(all_bestfit_IS)

# info for cycle phases
cycle <- tibble(phase = c('increase', 'peak', 'decline', 'low', 'increase'),
                start = c(1988, 1989, 1991, 1993, 1994), 
                end = c(1989, 1991, 1993, 1994, 1997)) %>% 
  mutate(phase = factor(phase, levels = c('increase', 'peak', 'decline', 'low')), 
         center = (start+end)/2)

grey_palette <- gray.colors(n=4, start = 0.75, end = 0.95, alpha = 0.5, rev = T)

v_palette <- viridis(5, 
                     alpha = 1, begin = 0, end = 1, direction = 1, 
                     option = "D") # generate colour palette; could use this to manually set colours for each species/node

# text to annotate phases of the cycle
ann_text <- data.frame(winter = cycle$center, 
                       interaction_strength = 0.02, 
                       label = c('increase', 'peak', 'decline', 'low', 'increase'))

(all_bestfit_IS_plot <- all_bestfit_IS %>% 
    mutate(winter = year +1, 
           pair = paste(consumer, resource, sep = '-'),
           pair = factor(pair, levels = c('lynx-hare', 'lynx-red_squirrel', 'coyote-hare', 'coyote-red_squirrel', 'owl-hare'))) %>% 
  ggplot()
  + geom_rect(aes(xmin = start, xmax = end, fill = as.factor(phase),
                  ymin = 0, ymax = Inf), lty = 0, data = cycle)
  + geom_point(aes(x = winter, y = interaction_strength, col = pair, pch = pair), cex = 4, stroke = 2)
  + geom_path(aes(x = winter, y = interaction_strength, col = pair, lty = pair), lwd = 1.5)
  + scale_color_manual(values=c(v_palette[c(1,1,3,3,4)]))
  + scale_shape_manual(values = c(16,1,15,0,17))
  + scale_linetype_manual(values = c(1, 2, 1, 2, 1))
  + theme_classic()
  + scale_x_continuous(name = 'winter', limits = c(1988,1997), breaks = c(1988:1997)) 
  + scale_y_continuous(name = 'interaction strength (per capita kill rate)')
  + guides(fill = 'none', color = guide_legend(ncol = 3))
  + theme(legend.position = 'top', 
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.key.width = unit(3, 'cm'), 
          legend.key.height = unit(1, 'cm'),
          legend.background = element_blank(), 
          axis.text = element_text(size = 12, face = 'bold'), 
          axis.title = element_text(size = 15, face = 'bold'))
  + scale_fill_manual(values = c(grey_palette[2], grey_palette[1], grey_palette[3], grey_palette[4], grey_palette[2]))
  + geom_text(aes(x = winter, y = interaction_strength), data = ann_text, label = ann_text$label, col = 'black', fontface = 'bold', size = 5)
)

ggsave('A1_interactionstrengths.tiff', all_bestfit_IS_plot, dpi = 300, width = 11, height = 8)




