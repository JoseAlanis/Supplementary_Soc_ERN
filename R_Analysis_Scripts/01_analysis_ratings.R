# --- author: Jose C. Garcia Alanis
# --- encoding: utf-8
# --- r version: 3.4.4 (2018-03-15) -- "Someone to Lean On"
# --- content: analysis of motivational ratings
# --- version: Fri Aug 31 11:56:50 2018

# --- 1) Set paths  -----------------------------------------------------------
project_path <- '/Users/Josealanis/Documents/Experiments/soc_ftask/'
upload_path <- '/Users/Josealanis/Documents/GitHub/Supplementary_Soc_ERN/'
setwd(project_path)


# --- 2) Get helper functions and packages ------------------------------------
source(paste(upload_path, 'r_functions/getPacks.R', sep = ''))

# Install and load multiple r packages necessary for analysis
pkgs <- c('dplyr', 'reshape2', 'tidyr', 
          'car',
          'ggplot2', 'viridis', 'cowplot')
getPacks(pkgs)
rm(pkgs)


# --- 3) Read-in raw data -----------------------------------------------------
# Get all behavioral
ratings <- read.table('./revision/rev_data/motiv_ratings.txt', 
                        header = T)


# --- 4) Analyse ratings data ------------------------------------------------
# # (UNCOMMENT to filter out outliers)
# ratings <- filter(ratings, start_motivation > 0)

# Set up model
contrasts(ratings$group) <- contr.sum(2); contrasts(ratings$group)
mod_ratings <- lm(data = ratings, start_motivation ~ group)
# Anova table
car::Anova(mod_ratings, type = 3)

# To long format for plot
ratings_long <- tidyr::gather(ratings,
                              stage, rating, start_motivation:end_motivation)
# # separate colums and reorder factors for plot
ratings_long <- ratings_long %>% separate(stage, 
                                          c('stage', 'item'), '_')
# Reorder and rename levels of factor
ratings_long$stage <- as.factor(ratings_long$stage)
ratings_long$stage <- factor(ratings_long$stage, 
                             levels(ratings_long$stage)[c(2,1)])
ratings_long$stage <- recode_factor(ratings_long$stage, 
                                    `1` = 'Start', `2` = 'End')
# Create variable that tracks the change in motivation for every person
# (i.e. drop vs. increase)
ratings_long$up_down <- ratings_long$diff_motivation >= 0


# --- 6) Plot changes in motivation -------------------------------------------
# Test wheathermotivation ratintig are significantly diferent,
# depending on context group (competiton vs. cooperation) and 
# stage of experiment (start vs. end).

# # (UNCOMMENT to filter out outliers)
# ratings <- filter(ratings, start_motivation > 0)

# Set up model
contrasts(ratings_long$group) <- contr.sum(2); contrasts(ratings_long$group)
contrasts(ratings_long$stage) <- contr.sum(2); contrasts(ratings_long$stage)
mod_ratings <- lm(data = ratings_long, rating ~ stage * group)

# Anova table
car::Anova(mod_ratings, type = 3)
# Model diagnostics
plot(mod_ratings)


# Plot rating values
ratings_plot <- ggplot(ratings_long, 
                       aes(x = stage, y = rating, fill = group, color = group)) + 
  
  coord_cartesian(ylim = c(-.5, 10)) +
  
  stat_summary(aes(y = rating),
               fun.y=mean,
               position = position_dodge(.7),
               geom="point",
               shape=95,
               size = 20) +
  
  geom_dotplot(binaxis='y', 
               stackdir='center', 
               binwidth = .2, 
               position = position_dodge(.7), 
               color = 'black') +
  
  scale_fill_viridis(option = 'D', begin = .1, end = .5, discrete = T ) +
  scale_color_viridis(option = 'D', begin = .1, end = .5, discrete = T ) +
  
  labs(x ='Stage of experiment', y = 'Mean motivation rating') +
  
  geom_segment(aes(x = 'start', y = -Inf, xend = 'end', yend = -Inf), 
               color = 'black', size = 1) +
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 10), 
               color = 'black', size = 1) +
  
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_text(color = 'black', face = 'bold', size = 13),
        axis.title.x = element_text(face = 'bold', color = 'black', size = 14, 
                                    margin = margin(t = 20)),
        axis.title.y = element_text(face = 'bold', color = 'black', size = 15, 
                                    margin = margin(r = 15)),
        legend.text = element_text(color = 'black', size = 12),
        legend.key.size = unit(.3, 'cm'),
        legend.title = element_blank(),
        legend.position = 'bottom'); ratings_plot

# Save plot to .pdf
cowplot::save_plot('./revision/figures/ratings.pdf', 
                   ratings_plot, 
                   base_height = 5, base_width = 4.5)