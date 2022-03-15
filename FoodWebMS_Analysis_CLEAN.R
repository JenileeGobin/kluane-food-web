# Investigating the impact of functional response relationships on food web node and network properties. 
# Jenilee Gobin
# 14 Aug 2020

# To run this code, ensure the community data are stored in a folder named "communities" within the working/project directory

# 1 - CHECK/LOAD PACKAGES ----- 

# function to look for & install multiple packages
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#  vector of needed packages
analysis.packages <- c("cheddar", "tidyverse", "viridis", "data.table", "cowplot", "ggplot2", "png", "grid", "cvequality")

# apply function to packages
check.packages(analysis.packages)

# 2 - LOAD COLLECTION OF COMMUNITIES -----

# load collection from working directory; food webs for each year will be a separate community within the collection
f_path <- getwd() # get working directory path
kluane_collection <- LoadCollection(paste(f_path))

# check the data
length(kluane_collection) # verify number of communities in the collection
names(kluane_collection) # veify the names of the communities in the collection
sapply(kluane_collection, 'NumberOfTrophicLinks') # number of trophic links in each community
CollectionCPS(kluane_collection) # CPS = collection community properties
CollectionNPS(kluane_collection) # collection node properties
CollectionTLPS(kluane_collection) # collection trophic link properties

# compute some properties for communities in the collection
CollectionCPS(kluane_collection, c(
  'NumberOfNodes',
  'NumberOfTrophicLinks',
  'DirectedConnectance',
  'NvMSlope')) 
# clearly 'DirectedConnectance' is not weighted but based solely on the number of links relative to the total number of possible links

# 3. GET NODE & NETWORK QUANTITATIVE DESCRIPTORS -----
# node and network quantitative descriptors are generated for weighted networks

# how to reference a single community in the collection
kluane_collection[["kluane_87_88"]] 

# make vectors for info that needs to be added to node quantitative descriptor output for the collection 
community <- rep(names(kluane_collection), each = 5)
node <- CollectionNPS(kluane_collection)$node
FR <- rep(CollectionCPS(kluane_collection)$FR, each = 5)
CI <- rep(CollectionCPS(kluane_collection)$CI_dat, each = 5)
CI_species <- rep(CollectionCPS(kluane_collection)$CI_species, each = 5)
CI_type <- rep(CollectionCPS(kluane_collection)$CI_type, each = 5)

# get node quantitative descriptors
NQD_all <- kluane_collection %>% 
  lapply(NodeQuantitativeDescriptors, weight = 'interaction_strength') %>% 
  do.call(rbind,.) %>% # remove highest list level (i.e., community names)
  as.data.frame() %>%  # make into dataframe
  cbind(community, node, FR, CI, CI_species, CI_type) %>% # add community names, nodes, FR, and CI info
  separate(community, into = c(NA, 'fall', 'winter'), sep = '_', remove = F) %>% # make continuous variables for years
  mutate(fall = as.numeric(fall), winter = as.numeric(winter))

# extract predation matrices to incorporate mass

PM_all <- kluane_collection %>% 
  lapply(PredationMatrix, weight = 'interaction_strength') %>% 
  do.call(rbind,.) %>% # remove highest list level (i.e., community names)
  as.data.frame() %>% 
  cbind(community, node, FR, CI, CI_species, CI_type) %>% # add community names, nodes, FR, and CI info
  separate(community, into = c(NA, 'fall', 'winter'), sep = '_', remove = F) %>% # make continuous variables for years
  mutate(fall = as.numeric(fall), winter = as.numeric(winter))

# make vectors for info that needs to be added to network quantitative descriptor output for the collection 
community <- rep(names(kluane_collection), each = 19)
FR <- rep(CollectionCPS(kluane_collection)$FR, each = 19)
CI <- rep(CollectionCPS(kluane_collection)$CI_dat, each = 19)
CI_species <- rep(CollectionCPS(kluane_collection)$CI_species, each = 19)
CI_type <- rep(CollectionCPS(kluane_collection)$CI_type, each = 19)

# get network quantitative descriptors
QD_all <- kluane_collection %>% 
  lapply(QuantitativeDescriptors, weight = 'interaction_strength') %>% 
  do.call(rbind,.) %>% # remove highest list level (i.e., community names)
  as.data.frame() %>% 
  mutate(metric = gsub('.[0-9]+', '', row.names(.))) %>% 
  cbind(community, FR, CI, CI_species, CI_type) %>% 
  separate(community, into = c(NA, 'fall', 'winter'), sep = '_', remove = F) %>% # make continuous variables for years
  mutate(fall = as.numeric(fall), winter = as.numeric(winter))

# reshape data
NQD_all_long <- NQD_all %>% 
  mutate(FR = factor(FR, levels = c('hyperbolic', 'best_fit')), 
         node = factor(node, levels = c('hare', 'red_squirrel', 'lynx', 'coyote', 'owl')), 
         panel = factor(if_else(node == 'hare' | node == 'red_squirrel', 'prey', 'consumer'), levels = c('prey', 'consumer'))) %>% 
  pivot_longer(cols = c(bIn, bOut), names_to = 'biomass_flow', values_to = 'biomass')

NQD_all_wide <- NQD_all_long %>% 
  dplyr::select(winter, panel, node, FR, CI_species, CI_type, biomass, biomass_flow, g, v) %>% 
  pivot_wider(values_from = c(biomass, g, v), names_from = c(CI_species, CI_type)) 

QD_all_wide <- QD_all %>% 
  mutate(FR = factor(FR, levels = c('hyperbolic', 'best_fit'))) %>% 
  filter(metric == 'Connectance'|metric == 'Generality'| metric =='Vulnerability') %>% 
  dplyr::select(-c('community', 'CI')) %>% 
  pivot_wider(values_from = c(Qualitative, Unweighted, Weighted), names_from = c(metric, CI_species, CI_type)) # can redo this to leave out metric...

# incorporate predator biomass values into bIN and bOUT

mass_values <- CollectionNPS(kluane_collection[1]) %>%  # get mass values from collection NPS
  select(-community)

true_biomass <- PM_all %>% 
  merge(mass_values) %>% 
  mutate(true_lynx = lynx*M, 
         true_coyote = coyote*M, 
         true_owl = owl*M) %>% # M values are matched to prey in each row of the predation matrix, therefore multiplying these and then summing (next mutate) yields predation rates that account for differences in the biomass of prey
  dplyr::select(-c(fall, community)) %>% # remove extra columns
  group_by(winter, CI_species, CI_type, FR) %>% 
  mutate(true_bIN_lynx = sum(true_lynx), 
         true_bIN_coyote = sum(true_coyote), 
         true_bIN_owl = sum(true_owl))
  
true_biomass_long <- true_biomass %>% 
  pivot_longer(cols = c('true_bIN_lynx', 'true_bIN_coyote', 'true_bIN_owl'), names_to = 'consumer', values_to = 'true_bIN') %>% 
  mutate(FR = factor(FR, levels = c('hyperbolic', 'best_fit')), 
         CI_species = factor(CI_species, levels = c('hare', 'rsq')),
         CI_type = factor(CI_type, levels = c('high', 'low'))
         )

true_biomass_wide <- true_biomass_long %>% 
  dplyr::select(FR, CI_species, CI_type, winter, consumer, true_bIN) %>% 
  separate(consumer, into = c(NA, NA, 'consumer'), sep = '_', remove = T) %>% 
  pivot_wider(values_from = true_bIN, names_from = c(CI_species, CI_type), values_fn = list(true_bIN = mean)) %>% 
  rename(true_bIN = NA_NA) %>% 
  mutate(consumer = factor(consumer,levels = c('lynx', 'coyote', 'owl')))

# incorporate biomass for prey (which is simply the product of the interaction strength and the mass of each prey type)

NQD_all_wide_prey_biomass <- NQD_all_long %>% 
  merge(CollectionNPS(kluane_collection)[,2:3]) %>% 
  filter(panel == 'prey', biomass_flow == 'bOut') %>% 
  mutate(true_biomass = biomass*M) %>% 
  dplyr::select(winter, panel, node, FR, CI_species, CI_type, true_biomass, g, v) %>% 
  pivot_wider(values_from = c(true_biomass, g, v), names_from = c(CI_species, CI_type), values_fn = list(true_biomass = mean, g = mean, v = mean))

# 4. CALCULATE DIFFERENCES BETWEEN BEST FIT AND HYPERBOLIC FR -----

# for cumulative interaction strengths

str(NQD_all_wide)

biomass_diff_df <- NQD_all_wide %>% 
  dplyr::select(c('winter', 'node', 'FR', 'biomass_NA_NA')) %>% 
  filter(biomass_NA_NA >0) %>% 
  pivot_wider(names_from = FR, values_from = biomass_NA_NA) %>% 
  mutate(abs_diff = hyperbolic-best_fit, 
         rel_diff = abs_diff/best_fit*100) %>% 
  arrange(winter)

head(biomass_diff_df)
write_csv(biomass_diff_df, 'cumulativeinteractionstrengths.csv')

# for biomass

# predator biomass
str(true_biomass_wide)

truebIN_diff_df <- true_biomass_wide %>%
  dplyr::select(c('FR', 'winter', 'consumer', 'true_bIN')) %>% 
  pivot_wider(names_from = FR, values_from = true_bIN) %>% 
  mutate(abs_diff = hyperbolic-best_fit, 
         rel_diff = abs_diff/best_fit*100) %>% 
  arrange(winter)

head(truebIN_diff_df)
write_csv(truebIN_diff_df, 'predator_bIN.csv')

# prey biomass
str(NQD_all_wide_prey_biomass)

truebOUT_diff_df <- NQD_all_wide_prey_biomass %>% 
  dplyr::select(c('winter', 'node', 'FR', 'true_biomass_NA_NA')) %>% 
  pivot_wider(names_from = FR, values_from = true_biomass_NA_NA) %>% 
  mutate(abs_diff = hyperbolic-best_fit, 
         rel_diff = abs_diff/best_fit*100) %>% 
  arrange(winter)

head(truebOUT_diff_df)
write_csv(truebOUT_diff_df, 'prey_bOUT.csv')

# estimate coefficients of variation for cumulative interaction strengths and biomass converted values 

# cv for cumulative interaction strengths (i.e., kill rates)
cv_node <- NQD_all_long %>% 
  filter(CI == 'no') %>% 
  dplyr::select(node, FR, biomass, v, g) %>% 
  rename(killrate = biomass) %>% 
  group_by(node, FR) %>% 
  summarize(across(everything(), cv, .names = "{.col}")) %>% 
  pivot_longer(cols = c(killrate, v, g)) %>% 
  pivot_wider(names_from = FR, values_from = value) %>% 
  rename(metric = name)


# cv for biomass converted values
str(NQD_all_wide_prey_biomass)
str(true_biomass_wide)

true_biomass_pred <- true_biomass_wide %>% 
  select(consumer, FR, true_bIN, winter) %>% 
  rename(true_biomass = true_bIN, node = consumer)

true_biomass_cv <- NQD_all_wide_prey_biomass %>% 
  filter(node == 'hare' | node == 'red_squirrel') %>% 
  select(node, FR, true_biomass_NA_NA, winter) %>% 
  rename(true_biomass = true_biomass_NA_NA) %>% 
  rbind(true_biomass_pred) %>% 
  group_by(node, FR) %>% 
  summarise(biomass_cv = cv(true_biomass)) %>% 
  pivot_wider(names_from = FR, values_from = biomass_cv) %>%
  mutate(metric = 'biomass') %>% 
  relocate(metric, .after = node)

cv_node_all <- rbind(cv_node, true_biomass_cv) %>% 
  relocate(node, .after = metric) %>% 
  filter(hyperbolic >0 & best_fit >0) %>% 
  mutate(metric = as.character(metric), 
         metric = replace(metric, metric == 'v', 'vulnerability'), 
         metric = replace(metric, metric == 'g', 'generality'), 
         metric = factor(metric, levels = c('killrate', 'biomass', 'vulnerability', 'generality'))) %>% 
  arrange(metric, node)

write.csv(cv_node_all, 'CV_node.csv')

# if wanted to test for significant differences in CV

cv_equality_test <- NQD_all_long %>% 
  filter(CI == 'no') %>% 
  dplyr::select(node, FR, biomass, v, g) %>% 
  rename(killrate = biomass) 

with(cv_equality_test, asymptotic_test(g, FR))

# estimate differences for node-level generality and vulnerability

# generality
generality_diff_df <- NQD_all_wide %>% 
  dplyr::select(c('winter', 'node', 'FR', 'g_NA_NA')) %>% 
  filter(g_NA_NA >0) %>% 
  pivot_wider(names_from = FR, values_from = g_NA_NA, values_fn = list(g_NA_NA = unique)) %>% 
  mutate(abs_diff = hyperbolic-best_fit, 
         rel_diff = abs_diff/best_fit*100) %>% 
  arrange(winter)

head(generality_diff_df)
write_csv(generality_diff_df, 'nodelevelgenerality.csv')

# vulnerability
vulnerability_diff_df <- NQD_all_wide %>% 
  dplyr::select(c('winter', 'node', 'FR', 'v_NA_NA')) %>% 
  filter(v_NA_NA >0) %>% 
  pivot_wider(names_from = FR, values_from = v_NA_NA, values_fn = list(v_NA_NA = unique)) %>% 
  mutate(abs_diff = hyperbolic-best_fit, 
         rel_diff = abs_diff/best_fit*100) %>% 
  arrange(winter)

head(vulnerability_diff_df)
write_csv(vulnerability_diff_df, 'nodelevelvulnerability.csv')

# check owl generality
owl_NQD <- NQD_all_long %>% 
  filter(node == 'owl' & CI == 'no') %>% 
  arrange(winter) %>% 
  mutate(conv_factor = g.prime/g) %>% 
  relocate(conv_factor, .after = g)

# estimate differences for network-level metrics 

QD_best_fit <- QD_all_wide %>% 
  dplyr::select(FR, winter, Weighted_Connectance_NA_NA, Weighted_Generality_NA_NA, Weighted_Vulnerability_NA_NA) %>% 
  filter(FR == 'best_fit')

QD_type_II <- QD_all_wide %>% 
  dplyr::select(FR, winter, Weighted_Connectance_NA_NA, Weighted_Generality_NA_NA, Weighted_Vulnerability_NA_NA) %>% 
  filter(FR == 'hyperbolic')

QD_abs_diff_dat <- (QD_type_II[,-(1:2)] - QD_best_fit[,-(1:2)]) %>% 
  cbind(.,QD_best_fit[,2]) %>% 
  relocate(winter, .before = Weighted_Connectance_NA_NA)

QD_rel_diff_dat <- (QD_abs_diff_dat[,-1]/QD_best_fit[,-(1:2)])*100

names (QD_abs_diff_dat) <- c('winter', 'connectance_abs', 'generality_abs', 'vulnerability_abs')
names (QD_rel_diff_dat) <- c('connectance_rel', 'generality_rel', 'vulnerability_rel')

QD_diff_dat <- cbind(QD_abs_diff_dat, QD_rel_diff_dat)

QD_values <- QD_best_fit %>% 
  rename(Weighted_Connectance_bestfit = Weighted_Connectance_NA_NA, 
         Weighted_Generality_bestfit = Weighted_Generality_NA_NA, 
         Weighted_Vulnerability_bestfit = Weighted_Vulnerability_NA_NA) %>% 
  cbind(QD_type_II[,-c(1:2)]) %>% 
  rename(Weighted_Connectance_typeII = Weighted_Connectance_NA_NA, 
         Weighted_Generality_typeII = Weighted_Generality_NA_NA, 
         Weighted_Vulnerability_typeII = Weighted_Vulnerability_NA_NA)

write_csv(QD_values, 'networkquantitativedescriptors.csv') # output values

write_csv(QD_diff_dat, 'networkquantitativedescriptors_differences.csv') # output differences

# cv for network level metrics
cv_network_all <- QD_combine_figs %>% 
  dplyr::select(metric, FR, winter, metric_type, value) %>% 
  mutate(metric_type = factor(metric_type, levels = c('Weighted', 'Unweighted'))) %>% 
  group_by(metric, metric_type, FR) %>% 
  summarize(cv = cv(value)) %>% 
  filter(metric_type != 'Qualitative') %>% 
  pivot_wider(names_from = FR, values_from = cv)%>% 
  ungroup() %>% 
  mutate(metric = as.character(metric), 
         metric_type = as.character(metric_type), 
         across(where(is.character), tolower)) 

write.csv(cv_network_all, 'CV_network.csv')

# 5. LOAD AND RESHAPE DENSITY DATA 
kluane_density_dat <- fread(file = 'Kluane_density_data.csv') %>% 
  mutate (winter = c(88:97), fall = c(87:96)) %>% 
  relocate(fall, .after = winter) %>% 
  pivot_longer(cols = c(lynx_density, 
                        coyote_density, 
                        hare_density, 
                        rsq_density, 
                        owl_density), 
               values_to = 'density', 
               names_to = 'species') %>% 
  separate(species, into = c('node', NA), sep = '_', remove = F) %>% 
  mutate(panel = factor(if_else(node == 'hare' | node == 'rsq', 'prey', 'consumer'), levels = c('prey', 'consumer')),
         rsq_CI_high = if_else(node == 'rsq', density+rsq_SE, NA_real_),
         rsq_CI_low = if_else(node == 'rsq', density-rsq_SE, NA_real_), 
         hare_CI_high = if_else(node == 'hare', hare_CI_high, NA_integer_), 
         hare_CI_low = if_else(node == 'hare', hare_CI_low, NA_integer_), 
         node = ifelse(node == 'rsq', 'red_squirrel', node),
         node = factor(node, levels = c('hare', 'red_squirrel', 'lynx', 'coyote', 'owl'))
  )

# 6. GENERATE PLOTS ----- 

# function to add images to plot

annotation_custom2 <- 
  function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){ 
    layer(data = data, stat = StatIdentity, position = PositionIdentity, 
          geom = ggplot2:::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob, 
                                            xmin = xmin, xmax = xmax, 
                                            ymin = ymin, ymax = ymax)
    )
  }

# get images for plots - all images taken from phylopic with no copyright (public domain dedication 1.0 license)

hare_image <- readPNG('hare.png')
red_squirrel_image <- readPNG('rsq.png')
lynx_image <- readPNG('lynx.png')
coyote_image <- readPNG('coyote.png')
owl_image <- readPNG('owl.png') 

# info for cycle phases
cycle <- tibble(phase = c('increase', 'peak', 'decline', 'low', 'increase'),
                start = c(87, 89, 91, 93, 94), 
                end = c(89, 91, 93, 94, 97)) %>% 
  mutate(phase = factor(phase, levels = c('increase', 'peak', 'decline', 'low')), 
         center = (start+end)/2)

# create grey palette for phases of the cycle
grey_palette <- gray.colors(n=4, start = 0.75, end = 0.95, alpha = 0.5, rev = T)

# create colour palette for nodes
v_palette <- viridis(length(unique(kluane_density_dat$node)), 
                     alpha = 1, begin = 0, end = 1, direction = 1, 
                     option = "D") # generate colour palette; could use this to manually set colours for each species/node

# text to annotate phases of the cycle
ann_text <- data.frame(winter = cycle$center, 
                       density = 5, 
                       node = 'owl', 
                       label = c('increase', 'peak', 'decline', 'low', 'increase')) %>% 
  mutate(node = factor(node, levels = c('hare', 'red_squirrel', 'lynx', 'coyote', 'owl')))

# annotation text version 2
ann_text2 <- data.frame(winter = cycle$center, 
                        density = 0.01, 
                        node = 'owl', 
                        label = c('increase', 'peak', 'decline', 'low', 'increase')) %>% 
  mutate(node = factor(node, levels = c('hare', 'red_squirrel', 'lynx', 'coyote', 'owl')))

#  annotation text version 3
ann_text3 <- ann_text2 %>% 
  mutate(consumer = 'owl') %>% 
  mutate(consumer = factor(consumer, levels = c('hare', 'red_squirrel', 'lynx', 'coyote', 'owl')))

# annotation text version 4
ann_text4 <- ann_text %>% 
  mutate (density = 0.025, metric = 'Connectance') %>% 
  dplyr::select(-node) 

# Figure 3 - Abundance over time -----

# specify positions for images

a_hare <- annotation_custom2(rasterGrob(hare_image, interpolate = T), 
                             xmin = 95, xmax = 97, ymin = 12500, ymax = Inf, 
                             data = kluane_density_dat[3,])
a_red_squirrel <- annotation_custom2(rasterGrob(red_squirrel_image, interpolate = T), 
                                     xmin = 95, xmax = 97, ymin = 7500, ymax = 15000, 
                                     data = kluane_density_dat[4,])
a_lynx <- annotation_custom2(rasterGrob(lynx_image, interpolate = T), 
                             xmin = 94.5, xmax = 97, ymin = 7, ymax = Inf, 
                             data = kluane_density_dat[1,])
a_coyote <- annotation_custom2(rasterGrob(coyote_image, interpolate = T), 
                               xmin = 94.5, xmax = 97, ymin = 4, ymax = Inf, 
                               data = kluane_density_dat[2,])
a_owl <- annotation_custom2(rasterGrob(owl_image, interpolate = T), 
                            xmin = 95, xmax = 97, ymin = 30, ymax = Inf, 
                            data = kluane_density_dat[5,])

# make plot
(abund_over_time <- ggplot(kluane_density_dat, aes(col = node))
  + scale_fill_manual(values = c(grey_palette[2], grey_palette[1], grey_palette[3], grey_palette[4], grey_palette[2]))
  + geom_rect(aes(xmin = start, xmax = end, fill = as.factor(phase),
                  ymin = 0, ymax = Inf, col = NA), lty = 0, data = cycle)
  + geom_point(aes(x = winter, y = density, pch = panel), cex = 3, stroke = 2)
  + geom_line(aes(x = winter, y = density), lwd = 1.5)
  + geom_errorbar(aes(x = winter, y = density, ymin = hare_CI_low, ymax = hare_CI_high), width = 0.2)
  + geom_errorbar(aes(x = winter, y = density, ymin = rsq_CI_low, ymax = rsq_CI_high), width = 0.2)
  + facet_grid(node~., scales = 'free')
  + scale_x_continuous(name = expression(bold(winter)), breaks = c(87:97))
  + scale_y_continuous(name = expression(bold("density"~"(individuals/100"~km^2*")")))
  + theme_classic()
  + theme(legend.position = 'none', 
          axis.text = element_text(size = 12, face = 'bold'), 
          axis.title = element_text(size = 18, face = 'bold'),
          strip.text = element_blank(), 
          strip.background = element_blank())
  + scale_color_viridis_d(option = 'viridis')
  + scale_shape_manual(values = c(1,2))
  + geom_text(aes(x = winter, y = density), data = ann_text, label = ann_text$label, col = 'black', fontface = 'bold', size = 5)
  + a_hare
  + a_red_squirrel
  + a_lynx
  + a_coyote
  + a_owl
)


# export image
ggsave('3_AbundanceOverTime.tiff', abund_over_time, dpi = 300, width = 8, height = 11)

# Figure 4 - Cumulative per capita daily kill rates -----

# specify positions for images

a_hare <- annotation_custom2(rasterGrob(hare_image, interpolate = T), 
                             xmin = 94.5, xmax = 95.5, ymin = 2, ymax = Inf, 
                             data = kluane_density_dat[3,])
a_red_squirrel <- annotation_custom2(rasterGrob(red_squirrel_image, interpolate = T), 
                                     xmin = 95.5, xmax = 96.5, ymin = 2.5, ymax = 4, 
                                     data = kluane_density_dat[4,])
a_lynx <- annotation_custom2(rasterGrob(lynx_image, interpolate = T), 
                             xmin = 95, xmax = 97, ymin = 2.5, ymax = 5.5, 
                             data = kluane_density_dat[1,])
a_coyote <- annotation_custom2(rasterGrob(coyote_image, interpolate = T), 
                               xmin = 94, xmax = 96, ymin = 1.75, ymax = 2.75, 
                               data = kluane_density_dat[2,])
a_owl <- annotation_custom2(rasterGrob(owl_image, interpolate = T), 
                            xmin = 95.25, xmax = 96.75, ymin = 0.03, ymax = 0.09, 
                            data = kluane_density_dat[5,])

(per_capita_kill_rate_images <- ggplot(subset(NQD_all_wide, biomass_NA_NA != 0), aes(col = node))
  + geom_rect(aes(xmin = start, xmax = end, fill = as.factor(phase),
                  ymin = 0, ymax = Inf, col = NA), lty = 0, data = cycle)
  + geom_text(aes(x = winter, y = density), data = ann_text2, label = ann_text2$label, col = 'black', fontface = 'bold', size = 5)
  + geom_point(aes(x = winter, y = biomass_NA_NA, shape = FR), cex = 3, stroke = 2)
  + geom_errorbar(aes(x = winter, y = biomass_NA_NA, ymin = biomass_hare_low, ymax = biomass_hare_high), width = 0.5)
  + facet_grid(node~., scales = 'free')
  + theme_classic()
  + scale_y_continuous(name = "")
  + scale_x_continuous(breaks = c(87:97))
  + scale_color_viridis_d()
  + scale_shape_manual(values = c(1,2))
  + scale_fill_manual(values = c(grey_palette[2], grey_palette[1], grey_palette[3], grey_palette[4], grey_palette[2]))
  + guides(colour = 'none', fill = 'none')
  + theme(legend.position = 'top', 
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1.5, 'cm'),
          axis.text = element_text(size = 12, face = 'bold'), 
          axis.title = element_text(size = 15, face = 'bold'),
          strip.text = element_blank(), 
          strip.background = element_blank())
  + a_hare
  + a_red_squirrel
  + a_lynx
  + a_coyote
  + a_owl
)

per_capita_kill_text <- ggdraw() + draw_text('cumulative per capita daily kill rate', fontface = 'bold', size = 17, angle = 90)

per_capita_kill_rate_final <- plot_grid(per_capita_kill_text, per_capita_kill_rate_images, nrow = 1, ncol = 2, rel_widths = c(1, 30))

ggsave('4_CulumlativeKillRates.tiff', per_capita_kill_rate_final, dpi = 300, width = 8, height = 11)

# Figure 5 - Incorporating biomass into kill rates -----

# specify positions for images

# for prey

a_hare <- annotation_custom2(rasterGrob(hare_image, interpolate = T), 
                             xmin = 94.5, xmax = 95.5, ymin = 3.5, ymax = Inf, 
                             data = kluane_density_dat[3,])

a_red_squirrel <- annotation_custom2(rasterGrob(red_squirrel_image, interpolate = T), 
                                     xmin = 95.5, xmax = 96.5, ymin = 0.55, ymax = 0.9, 
                                     data = kluane_density_dat[4,])
# for predators

kluane_density_dat_pred <-rbind(kluane_density_dat[1,], kluane_density_dat[2,], kluane_density_dat[5,]) %>% 
  rename(consumer = node)

a_lynx <- annotation_custom2(rasterGrob(lynx_image, interpolate = T), 
                             xmin = 95.5, xmax = 97, ymin = 0, ymax = 1.25, 
                             data = kluane_density_dat_pred[1,])
a_coyote <- annotation_custom2(rasterGrob(coyote_image, interpolate = T), 
                               xmin = 94, xmax = 96, ymin = 2.25, ymax = 3.75, 
                               data = kluane_density_dat_pred[2,])
a_owl <- annotation_custom2(rasterGrob(owl_image, interpolate = T), 
                            xmin = 95.25, xmax = 96.75, ymin = 0.04, ymax = 0.13, 
                            data = kluane_density_dat_pred[3,])

# plot for prey

(bOut_prey_images <- ggplot(NQD_all_wide_prey_biomass, aes(col = node))
 + geom_rect(aes(xmin = start, xmax = end, fill = as.factor(phase),
                 ymin = 0, ymax = Inf, col = NA),lty = 0, data = cycle)
 #+ geom_text(aes(x = winter, y = density), data = ann_text, label = ann_text$label, col = 'black', fontface = 'bold', size = 5)
 + geom_point(aes(x = winter, y = true_biomass_NA_NA, shape = FR), cex = 3, stroke = 2)
 + geom_errorbar(aes(x = winter, y = true_biomass_NA_NA, ymin = true_biomass_hare_low, ymax = true_biomass_hare_high), width = 0.5)
 + facet_grid(node~., scales = 'free')
 + theme_classic()
 + scale_y_continuous(name = "")
 + scale_x_continuous(breaks = c(87:97))
 + scale_color_manual(values = v_palette[c(2,5)])
 + scale_shape_manual(values = c(1,2))
 + scale_fill_manual(values = c(grey_palette[2], grey_palette[1], grey_palette[3], grey_palette[4], grey_palette[2]))
 + guides(colour = 'none', fill = 'none')
 + theme(legend.position = 'top', 
         legend.title = element_blank(),
         legend.text = element_text(size = 12),
         legend.key.size = unit(1.5, 'cm'),
         axis.text = element_text(size = 12, face = 'bold'), 
         axis.title = element_text(size = 15, face = 'bold'),
         strip.text = element_blank(), 
         strip.background = element_blank(), 
         axis.title.x = element_blank(), 
         axis.text.x = element_blank(), 
         axis.ticks.x = element_blank())
  + a_hare
  + a_red_squirrel
)

# plot for predators

(bIn_pred_plot_images <- ggplot(true_biomass_wide, aes(col = consumer))
  + geom_rect(aes(xmin = start, xmax = end, fill = as.factor(phase),
                  ymin = 0, ymax = Inf, col = NA), lty = 0, data = cycle)
  + geom_text(aes(x = winter, y = density), data = ann_text3, label = ann_text3$label, col = 'black', fontface = 'bold', size = 5)
  + geom_point(aes(x = winter, y = true_bIN, shape = FR), cex = 3, stroke = 2)
  + geom_errorbar(aes(x = winter, y = true_bIN, ymin = hare_low, ymax = hare_high), width = 0.5)
  + facet_grid(consumer~., scales = 'free')
  + theme_classic()
  + scale_y_continuous(name = "")
  + scale_x_continuous(breaks = c(87:97))
  + scale_color_manual(values=c(v_palette[c(1,3,4)]))
  + scale_shape_manual(values = c(1,2))
  + scale_fill_manual(values = c(grey_palette[2], grey_palette[1], grey_palette[3], grey_palette[4], grey_palette[2]))
  + guides(colour = 'none', fill = 'none')
  + theme(legend.position = 'none',
          # legend.title = element_blank(),
          # legend.text = element_text(size = 12),
          # legend.key.size = unit(1.5, 'cm'),
          axis.text = element_text(size = 12, face = 'bold'),
          axis.title = element_text(size = 15, face = 'bold'),
          strip.text = element_blank(),
          strip.background = element_blank())
  + a_lynx
  + a_coyote
  + a_owl
)

# combine prey and predator plots

biomass_text <- ggdraw() + draw_text('biomass (kg)', fontface = 'bold', size = 17, angle = 90)

prey_pred_biomass_converted <- plot_grid(biomass_text, plot_grid(bOut_prey_images, bIn_pred_plot_images, nrow = 2, ncol = 1, rel_heights = c(2.25, 3)), nrow = 1, ncol = 2, rel_widths = c(0.075,1))

ggsave('5_bInOut_preypred_biomassconversion.tiff', prey_pred_biomass_converted, dpi = 300, width = 8, height = 11)

# Figure 6 - node level vulnerability and generality -----

# specify positions for images

a_hare <- annotation_custom2(rasterGrob(hare_image, interpolate = T), 
                             xmin = 93.5, xmax = 94.5, ymin = 3, ymax = Inf, 
                             data = kluane_density_dat[3,])
a_red_squirrel <- annotation_custom2(rasterGrob(red_squirrel_image, interpolate = T), 
                                     xmin = 96, xmax = 97, ymin = 1.5, ymax = 2.5, 
                                     data = kluane_density_dat[4,])
a_lynx <- annotation_custom2(rasterGrob(lynx_image, interpolate = T), 
                             xmin = 95.5, xmax = 97, ymin = 0, ymax = 2, 
                             data = kluane_density_dat[1,])
a_coyote <- annotation_custom2(rasterGrob(coyote_image, interpolate = T), 
                               xmin = 95.5, xmax = 97, ymin = 0.95, ymax = 2.25, 
                               data = kluane_density_dat[2,])
a_owl <- annotation_custom2(rasterGrob(owl_image, interpolate = T), 
                            xmin = 95, xmax = 97, ymin = 0.15, ymax = 0.25, 
                            data = kluane_density_dat[5,])

# plot generality of predators and vulnerability of prey change over time (at the node/species level)

(generalityvulnerability_plot_images<- ggplot()
  + geom_rect(aes(xmin = start, xmax = end, fill = as.factor(phase),
                  ymin = 0, ymax = Inf), lty = 0, data = cycle)
  + geom_point(data = subset(NQD_all_wide, g_NA_NA!=0), aes(x = winter, y = g_NA_NA, col = node, shape = FR), cex = 3, stroke = 2)
  + geom_errorbar(data = subset(NQD_all_wide, g_NA_NA!=0), aes(x = winter, y = g_NA_NA, ymin = g_hare_low, ymax = g_hare_high, col = node), width = 0.5)
  + geom_point(data = subset(NQD_all_wide, v_NA_NA!=0), aes(x = winter, y = v_NA_NA, col = node, shape = FR), cex = 3, stroke = 2)
  + geom_errorbar(data = subset(NQD_all_wide, v_NA_NA!=0), aes(x = winter, y = v_NA_NA, ymin = v_hare_low, ymax = v_hare_high, col = node), width = 0.5)
  + facet_grid(node~., scales = 'free')
  + theme_classic()
  + scale_y_continuous(name = "")
  + scale_x_continuous(breaks = c(87:97))
  + scale_color_viridis_d()
  + scale_shape_manual(values = c(1,2))
  + geom_text(aes(x = winter, y = density), data = ann_text2, label = ann_text2$label, col = 'black', fontface = 'bold', size = 5)
  + scale_fill_manual(values = c(grey_palette[2], grey_palette[1], grey_palette[3], grey_palette[4], grey_palette[2]))
  + guides(colour = 'none', fill = 'none')
  + theme(legend.position = 'top', 
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1.5, 'cm'),
          axis.text = element_text(size = 12, face = 'bold'), 
          axis.title = element_text(size = 15, face = 'bold'),
          strip.text = element_blank(), 
          strip.background = element_blank())
  + a_hare
  + a_red_squirrel
  + a_lynx
  + a_coyote
  + a_owl
)

generality_text <- ggdraw() + draw_text('generality', fontface = 'bold', size = 17, angle = 90)
vulnerability_text <- ggdraw() + draw_text('vulnerability', fontface = 'bold', size = 17, angle = 90)

generality_vulnerability_plot_final <- plot_grid(plot_grid(NULL, vulnerability_text, generality_text, nrow = 3, ncol = 1, rel_heights = c(0.3,2,3)), generalityvulnerability_plot_images, nrow = 1, ncol = 2, rel_widths = c(1, 30))

ggsave('6_GeneralityVulnerability_NOpanels.tiff', generality_vulnerability_plot_final, dpi = 300, width = 8, height = 11)

# Network quantitative descriptors

# reshape data to combine all network-level metrics into a single paneled figure

QD_combine_figs <- QD_all %>% 
  mutate(FR = factor(FR, levels = c('hyperbolic', 'best_fit'))) %>% 
  filter(metric == 'Connectance'|metric == 'Generality'| metric =='Vulnerability') %>% 
  dplyr::select(-c('community', 'CI')) %>% 
  pivot_longer(cols = c(Qualitative:Weighted), names_to = 'metric_type', values_to = 'value') %>% 
  pivot_wider(values_from = value, names_from = c(CI_species, CI_type)) %>% 
  rename(value = NA_NA) %>% 
  mutate (metric = factor(metric, levels = c('Vulnerability', 'Generality', 'Connectance'))) 

# plot weighted metrics

(weighted_network_metrics <- QD_combine_figs %>% 
    filter(metric_type == 'Weighted') %>% 
    ggplot()
  + geom_rect(aes(xmin = start, xmax = end, fill = as.factor(phase),
                  ymin = 0, ymax = Inf), lty = 0, data = cycle)
  + geom_point(aes(x = winter, y = value, pch = FR), cex = 3, stroke = 2)
  + geom_path(aes(x = winter, y = value, lty = FR))
  + geom_errorbar(aes(x = winter, y = value, ymin = hare_low, ymax = hare_high), width = 0.5)
  + facet_grid(as.factor(metric)~., scales = 'free')
  + theme_classic()
  + scale_shape_manual(values = c(1,2))
  + scale_x_continuous(breaks = c(87:97))
  + geom_text(aes(x = winter, y = density), data = ann_text4, label = ann_text4$label, col = 'black', fontface = 'bold', size = 5)
  + scale_fill_manual(values = c(grey_palette[2], grey_palette[1], grey_palette[3], grey_palette[4], grey_palette[2]))
  + guides(colour = 'none', fill = 'none')
  + theme(legend.position = 'top', 
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1.5, 'cm'),
          axis.text = element_text(size = 12, face = 'bold'), 
          axis.title = element_text(size = 15, face = 'bold'), 
          axis.title.y = element_blank())
)

ggsave('7_WeightedNQD.tiff', weighted_network_metrics, dpi = 300, width = 6, height = 8)

# plot unweighted metrics 

(unweighted_network_metrics <- QD_combine_figs %>% 
    filter(metric_type == 'Unweighted') %>% 
    ggplot()
  + geom_rect(aes(xmin = start, xmax = end, fill = as.factor(phase),
                  ymin = 0, ymax = Inf), lty = 0, data = cycle)
  + geom_point(aes(x = winter, y = value, pch = FR), cex = 3, stroke = 2)
  + geom_path(aes(x = winter, y = value, lty = FR))
  + geom_errorbar(aes(x = winter, y = value, ymin = hare_low, ymax = hare_high), width = 0.5)
  + facet_grid(as.factor(metric)~., scales = 'free')
  + theme_classic()
  + scale_shape_manual(values = c(1,2))
  + scale_x_continuous(breaks = c(87:97))
  + geom_text(aes(x = winter, y = density), data = ann_text4, label = ann_text4$label, col = 'black', fontface = 'bold', size = 5)
  + scale_fill_manual(values = c(grey_palette[2], grey_palette[1], grey_palette[3], grey_palette[4], grey_palette[2]))
  + guides(colour = 'none', fill = 'none')
  + theme(legend.position = 'top', 
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1.5, 'cm'),
          axis.text = element_text(size = 12, face = 'bold'), 
          axis.title = element_text(size = 15, face = 'bold'), 
          axis.title.y = element_blank())
)

ggsave('7b_UnweightedNQD.tiff', unweighted_network_metrics, dpi = 300, width = 6, height = 8)

# ----- SESSION INFO -----

sessionInfo()





























