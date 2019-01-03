# --- author: jose c. garcia alanis
# --- encoding: utf-8
# --- r version: 3.5.1 (2018-07-02) -- "Feather Spray"
# --- script version: Dez 2018
# --- content: create erp figures

# --- 1) set paths and get workflow functions ----------------------------------
# path to project
setwd('/Volumes/TOSHIBA/manuscrips_and_data/soc_ern/')

# workflow functions
source('./r_functions/getPacks.R')
source('./r_functions/stdResid.R')
source('./r_functions/dataSummary.R')
# source('./r_functions/topoplot.R') # from package eegUtils


# Get helper function
source('./R_Functions/getPacks.R')
source('~/Documents/r_functions/topoplot.R')

# load multiple packages necessary for analysis
getPacks(c('dplyr', 'plyr', 'reshape2', 
           'ggplot2', 'viridis', 'cowplot', 'RColorBrewer'))

# Get electrode locations
chanlocs <- read.table('./meta_dat/chanlocs_ftask.txt',
                       header = T, 
                       sep = ',')
# Rename column
names(chanlocs)[1] <- 'Electrode'


# --- 2) Read in the data  -----------------------------------------------------
load('./data_for_r/Incomp_ERPs.RData')


# --- 3) COMPUTE ∆ERN (incorrect-correct) --------------------------------------
# ----- From long to wide
ERP_wide <- dcast(ERP, Subject + Group + Electrode + Time ~ Reaction,
                  value.var = c('Amplitude'))

# --- Compute ERN
ERP_wide <- ERP_wide %>% mutate(ERN = Incorrect - Correct)



# --- 4) Create data sets for plots --------------------------------------------

#  ****** This section will take a while
#  ****** go get a coffee

#  ****** AVERAGE for electrode and GFP plot
Ave_Elect <- plyr::ddply(ERP_wide, 
                       c('Electrode', 'Time'), 
                       
                       dplyr::summarise,
                       
                       N = sum(!is.na(ERN)),
                       M_Amp = mean(ERN, na.rm = T),
                       sd   = sd(ERN, na.rm = T),
                       se   = sd / sqrt(N),
                       # Compute t-staistic for confidence interval 
                       # (i.e., quantile t, with N-1 degrees of freedom)
                       # and use it to compute individual ci 
                       # (i.e., multiplier for se):
                       ci   = se * qt(.95/2 + .5, N-1) )



# ******* AVERAGE for ∆ERN plot
Ave_ERN <- plyr::ddply(filter(ERP_wide,
                              Electrode == 'FCz' | Electrode == 'Cz'), 
                       c('Group', 'Time'), 
                       
                       dplyr::summarise,
                       
                       N = sum(!is.na(ERN)),
                       M_Amp = mean(ERN, na.rm = T),
                       sd   = sd(ERN, na.rm = T),
                       se   = sd / sqrt(N),
                       # Compute t-staistic for confidence interval 
                       # (i.e., quantile t, with N-1 degrees of freedom)
                       # and use it to compute individual ci 
                       # (i.e., multiplier for se):
                       ci   = se * qt(.95/2 + .5, N-1) )


# ******* AVERAGE for ∆ERN ERP plot
Subj_ERN <- plyr::ddply(filter(ERP_wide, 
                               Electrode == 'FCz' | Electrode == 'Cz'), 
                       c('Group', 'Subject', 'Time'), 
                       
                       dplyr::summarise,
                       
                       N = sum(!is.na(ERN)),
                       M_Amp = mean(ERN, na.rm = T),
                       sd   = sd(ERN, na.rm = T),
                       se   = sd / sqrt(N),
                       # Compute t-staistic for confidence interval 
                       # (i.e., quantile t, with N-1 degrees of freedom)
                       # and use it to compute individual ci 
                       # (i.e., multiplier for se):
                       ci   = se * qt(.95/2 + .5, N-1) )



# ******* AVERAGE ∆ERN topoplot
Topo_ERN <- plyr::ddply(ERP_wide, 
                       c('Group', 'Time', 'Electrode'), 
                       
                       dplyr::summarise,
                       
                       N = sum(!is.na(ERN)),
                       M_Amp = mean(ERN, na.rm = T),
                       sd   = sd(ERN, na.rm = T),
                       se   = sd / sqrt(N),
                       # Compute t-staistic for confidence interval 
                       # (i.e., quantile t, with N-1 degrees of freedom)
                       # and use it to compute individual ci 
                       # (i.e., multiplier for se):
                       ci   = se * qt(.95/2 + .5, N-1) )



# ******* AVERAGE for ERN and CRN plot
Ave_ERP <- plyr::ddply(ERP, 
                       c('Group', 'Reaction', 'Electrode', 'Time'), 
                       
                       dplyr::summarise,
                       
                       N = sum(!is.na(Amplitude)),
                       M_Amp = mean(Amplitude, na.rm = T),
                       sd   = sd(Amplitude, na.rm = T),
                       se   = sd / sqrt(N),
                       # Compute t-staistic for confidence interval 
                       # (i.e., quantile t, with N-1 degrees of freedom)
                       # and use it to compute individual ci 
                       # (i.e., multiplier for se):
                       ci   = se * qt(.95/2 + .5, N-1) )



# --- 5) Plot ∆ERN by electrode ------------------------------------------------
# ----- Data to wide
Ave_Elect_wide <- dcast(Ave_Elect, Time ~ Electrode,
                  value.var = c('M_Amp'))
# ----- Compute GFP
Ave_Elect_wide$GFP <- apply(Ave_Elect_wide[, 2:65], 1, sd)
Ave_Elect_long <- melt(Ave_Elect_wide, id.vars=c("Time"))
names(Ave_Elect_long)[2:3] <- c('Electrode', 'Amplitude')

# ----- Create ∆ERN wave by electrode plot
Ave_ro <- Ave_Elect_long %>% filter(Time >= 0 & Time <= 100) %>% 
  group_by(Electrode) %>% 
  dplyr::summarise(Mean_Amp = mean(Amplitude)) %>% 
  arrange(Mean_Amp)

Ave_Elect_long$Electrode <- factor(Ave_Elect_long$Electrode, 
                                   levels = as.character(Ave_ro$Electrode))

elec_p <- ggplot(filter(Ave_Elect_long, 
                        Electrode %in% c('Fz', 'FCz', 'Cz','CPz', 'Pz', 'GFP'),
                        Time >= -500), 
                 
                 aes(Time, Amplitude, color = Electrode, size = Electrode)) +
  
  annotate("rect", xmin = 0, xmax = 100, 
           ymin = -Inf, ymax = Inf, alpha = .1) +
  geom_vline(xintercept = c(0), color='black', 
             size=rel(0.5), linetype = 3) +
  geom_hline(yintercept = c(0), color='black', 
             size=rel(0.5), linetype = 3) +
  
  geom_line(size = 1) +
  
  coord_cartesian(ylim = c(-7.5, 5)) +
  scale_y_reverse(breaks=c(-7.5, -5, -2.5, 0, 2.5, 5)) + 
  scale_x_continuous(breaks=c(-500, -250, 0, 250, 500, 750, 1000)) + 
  
  scale_colour_manual(breaks = c('Cz','FCz','CPz', 'Pz', 'Fz', 'GFP'),
                      values = c(viridis(option = 'D', 
                                         n = 5, 
                                         direction = -1, 
                                         begin = .1, 
                                         end =.90), 'black')) +
  
  scale_size_manual(values = c(rep(0.8, 5), 1), guide = F) +
  
  labs(y = expression(bold(paste("Amplitude (", mu, "V)"))), 
       x = expression(bold('Time (ms)')), 
       title = expression(bold(paste(Delta,"ERN across midline electrodes")))  ) +
  
  geom_segment(aes(x = -Inf, y = -7.5, xend = -Inf, yend = 5), 
               color = 'black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = -500, y = Inf, xend = 1000, yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme_classic() +
  theme(
    plot.title = element_text(color = "black", vjust = 1, size = 14, face = 'bold'),
    strip.background = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_text(color = 'black', size = 14, face = 'bold', 
                                margin = margin(t = 15)),
    axis.title.y = element_text(color = 'black', size = 14, face = 'bold', 
                                margin = margin(r = 15)),
    axis.text.x = element_text(color = 'black', size = 13),
    axis.text.y  = element_text(hjust = 1, color = 'black', size = 13),
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 11),
    legend.key.size = unit(0.7, 'cm'),
    legend.position = c(0.1, 0.75), 
    legend.box = "horizontal", 
    legend.background = element_blank()); elec_p 


# ----- Save the plot
# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_S1a.pdf', 
#                    elec_p, base_height = 5, base_width = 7)

# ----- Calculate charthesian coordiantes
# ----- for Mean ∆ERN topoplot
chanlocs$radianTheta <- pi/180*chanlocs$theta

chanlocs <- chanlocs %>%
  mutate(x = .$radius*sin(.$radianTheta),
         y = .$radius*cos(.$radianTheta))

# ----- Data for plot
to_p <- merge(Ave_Elect_long, select(chanlocs, x, y, Electrode), 'Electrode')
to_p <- select(to_p, Electrode, x, y, Time, Amplitude)
to_p <- filter(to_p, Time >= 0 & Time <= 100)
names(to_p) <- gsub(names(to_p),
                    pattern = '([[:upper:]])',
                    perl = TRUE,
                    replacement = '\\L\\1')

# ----- CREATE THE PLOT
t_plot  <- topoplot(to_p, contour = T, 
                    chan_marker = 'none', 
                    palette = 'D', limits = c(-6.5,6.5),
                    grid_res = 100); t_plot

# ----- ADD title
t_plot <- t_plot + 
  labs(title = 'Mean activity') + 
  theme(plot.title = element_text(hjust = .5, size = 19, face='bold'),
        legend.title = element_text(hjust = .5, size = 19, face='bold'),
        legend.text = element_text(size = 18))


# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_S1b.pdf', 
#                    t_plot, base_height = 5, base_width = 5)



# ----- 5) Plot ∆ERN Wave ---------------------------------
# ----- CREATE the plot 
ERN_p <- ggplot(filter(Ave_ERN, Time >= -500), 
       
       aes(Time, M_Amp, fill = Group, color = Group, linetype = Group)) + 
  
  theme_classic() + 
  
  annotate("rect", xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, alpha = .1) +
  geom_vline(xintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  geom_hline(yintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  
  geom_ribbon(aes(ymin = M_Amp - se, ymax = M_Amp + se), 
              alpha = .4, colour = NA) +
  
  geom_line( size = rel(1)) + 
  scale_y_reverse(breaks = c(-7.5 ,-5, -2.5, 0, 2.5, 5)) +
  scale_x_continuous(breaks = c(-500, -250, 0, 250, 500, 750, 1000)) +
  
  scale_color_viridis(option = 'A', discrete = T, 
                      direction = -1, begin = .05, end = .6) +
  scale_fill_viridis(option = 'A', discrete = T, 
                     direction = -1, begin = .05, end = .6) +
  
  labs(x = expression(bold("Time (ms)")), 
       y = expression(bold(paste("Amplitude (", mu, "V)"))), 
       title= expression(bold(paste("Grand averaged ", 
                                    Delta, 
                                    "ERN by social context")))) +
  
  geom_segment(aes(x = -Inf, y = -7.5, xend = -Inf, yend = 5), 
               color = 'black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = -500, y = Inf, xend = 1000, yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme(
    strip.background = element_blank(),
    axis.line = element_blank(),
    
    plot.title = element_text(color = "black", size = 14, face = 'bold'),
    axis.title.y = element_text(color = 'black', size = 14, face = 'bold',
                                margin = margin(r = 15)),
    axis.title.x = element_text(color = 'black', size = 14, face = 'bold', 
                                margin = margin(t = 15)),
    axis.text.x = element_text(color = 'black', size = 13),
    axis.text.y  = element_text(color = 'black', size = 13),
    legend.text = element_text(color = "black", size = 11),
    legend.key.size = unit(0.8, 'cm'),
    legend.key.width = unit(1, 'cm'),
    legend.title = element_blank(),
    
    legend.position = c(0.15, 0.9), 
    legend.box = "horizontal"); ERN_p


# SAVE PLOT
# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_4a.pdf', 
#                    ERN_p, base_height = 5, base_width = 7)


# ----- Calculate charthesian coordiantes
# ----- for ∆ERN by group topoplot
chanlocs$radianTheta <- pi/180*chanlocs$theta

chanlocs <- chanlocs %>%
  mutate(x = .$radius*sin(.$radianTheta),
         y = .$radius*cos(.$radianTheta))

# *****
# ----- Competiton data
to_p <- merge(Topo_ERN, select(chanlocs, x, y, Electrode), 'Electrode')
to_p <- select(to_p, Electrode, x, y, Time, M_Amp, Group)
to_p <- filter(to_p, Time >= 0 & Time <= 100, Group == 'Competition')
names(to_p) <- gsub(names(to_p),
                    pattern = '([[:upper:]])',
                    perl = TRUE,
                    replacement = '\\L\\1')
names(to_p)[5] <- 'amplitude'

# ----- CREATE the plot 
t_plot  <- topoplot(to_p, contour = T, 
                    chan_marker = 'none', 
                    palette = 'A', limits = c(-6.5,6.5),
                    grid_res = 100); t_plot

# ----- ADD title
t_plot <- t_plot + 
  labs(title = 'Competition') + 
  theme(plot.title = element_text(hjust = .5, size = 17, face = 'bold'),
        legend.title = element_text(hjust = .5, size = 17, face = 'bold'),
        legend.text = element_text(size = 17)); t_plot

# ----- SAVE the plot
# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_4b1.pdf', 
#                    t_plot, base_height = 4, base_width = 6)

# *****
# ----- Cooperation data
to_p <- merge(Topo_ERN, select(chanlocs, x, y, Electrode), 'Electrode')
to_p <- select(to_p, Electrode, x, y, Time, M_Amp, Group)
to_p <- filter(to_p, Time >= 0 & Time <= 100, Group == 'Cooperation')
names(to_p) <- gsub(names(to_p),
                    pattern = '([[:upper:]])',
                    perl = TRUE,
                    replacement = '\\L\\1')
names(to_p)[5] <- 'amplitude'

# ----- CREATE the plot 
t_plot  <- topoplot(to_p, contour = T, 
                    chan_marker = 'none', 
                    palette = 'A', limits = c(-6.5,6.5),
                    grid_res = 100); t_plot

# ----- ADD title
t_plot <- t_plot + 
  labs(title = 'Cooperation') + 
  theme(plot.title = element_text(hjust = .5, size = 17, face = 'bold'),
        legend.title = element_text(hjust = .5, size = 17, face = 'bold'),
        legend.text = element_text(size = 17)); t_plot

# ----- SAVE the plot
# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_4b2.pdf', 
#                    t_plot, base_height = 4, base_width = 6)





# ----- 5) Plot ∆ERN ERP Image ----------------------------
# ----- Arrange data frame accordint to electrodes amplitude
Ave_ro <- Subj_ERN %>% filter(Time >= 0 & Time <= 100) %>% 
  group_by(Group, Subject) %>% 
  dplyr::summarise(Mean_Amp = mean(M_Amp)) %>% 
  arrange(Group, Mean_Amp)

Subj_ERN$Subject <- factor(Subj_ERN$Subject, 
                           levels = as.character(Ave_ro$Subject))

# ----- CREATE the plot
erp_ern <- ggplot(filter(Subj_ERN, Time >= -500), 
       aes(x = Time, y = as.numeric(Subject))) + 
  
  geom_raster(aes(fill = M_Amp), interpolate = T) +
  
  annotate("rect", xmin = -500, xmax = 1000, 
           ymin = 0.5, ymax = 37.5, 
           alpha = .1, fill = NA, color = 'black') +
  geom_vline(xintercept = c(0), color = 'white', 
             size = 0.5, linetype=1) +
  
  scale_fill_gradientn(colors = brewer.pal(11, 'RdBu'), 
                       limits = c(-19, 19), breaks = c(-16, -8, 0, 8 , 16)) +
  scale_y_continuous(breaks = c(1, 76)) +
  scale_x_continuous(breaks = c(-500, -250,  0, 250, 500, 750, 1000)) +
  labs(x = expression(bold('Time (ms)')),
       y = expression(bold('Subject')), 
       title = expression(bold(paste('ERP image of individual ', Delta, 'ERN averages'))), 
       fill = expression(bold(paste('Amplitude (', mu, 'V)')))) +
  theme_classic() +
  geom_segment(aes(x = -Inf, y = 1, xend = -Inf, yend = 76), 
               color = 'black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = -500, y = -Inf, xend = 1000, yend = -Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  theme(
    plot.title = element_text(color = "black", face = 'bold', size = 14),
    axis.title.y = element_text(color = "black", face = 'bold', size = 14, 
                                margin = margin(r = 10)),
    axis.title.x = element_text(color = "black", face = 'bold' , size = 14,
                                margin = margin(t = 15)),
    axis.text.y  = element_text(color='black', size = 13),
    axis.text.x = element_text(color='black', size = 13),
    strip.background = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(color = "black", face = 'bold', size = 12),
    legend.text = element_text(color = "black", size = 12),
    legend.key.size = unit(1, 'cm'), 
    legend.position = "bottom"); erp_ern

# --- Reorder the colorbar
erp_ern <- erp_ern +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = .5, title.vjust = 1,
                               barwidth = 9,
                               barheight = 1.1)); erp_ern

# ----- SAVE the plot
# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_4c.pdf', 
#                    erp_ern, base_height = 5, base_width = 7)



# ----- 8) Plot ERPs --------------------------------------
# ----- CREATE the plot
ern_crn_erp <- ggplot(filter(Ave_ERP, Electrode == 'FCz', # <- select electrode to plot
              Time >= -500), 
       
       aes(Time, M_Amp, fill = Reaction, linetype = Group, color = Reaction)) + 
  
  annotate('text', x = -510, y = -1,
           label = 'paste(bold(FCz))', parse = TRUE, 
           size = 3, hjust = 0) +
  
  facet_wrap(~ Group, scales = 'free_y', ncol = 4) + 
  theme_classic() + 
  
  annotate("rect", xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, alpha = .1) +
  geom_vline(xintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  geom_hline(yintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  
  geom_ribbon(aes(ymin = M_Amp - se, ymax = M_Amp + se), alpha = .4, colour = NA) +
  geom_line( size = rel(1)) + 
  
  scale_y_reverse(breaks = c(-8 ,-4, 0,4, 8)) +
  scale_x_continuous(breaks = c(-500, -250, 0, 250, 500, 750, 1000)) +
  
  scale_color_viridis(option = 'C', discrete = T, direction = 1, end = .6) +
  scale_fill_viridis(option = 'C', discrete = T, direction = 1, end = .6) +
  
  labs(x = expression(bold("Time (ms)")), 
       y = expression(bold(paste("Amplitude (", mu, "V)"))), 
       title="Grand average ERPs") +
  
  geom_segment(aes(x = -Inf, y = -8, xend = -Inf, yend = 8), 
               color = 'black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = -500, y = Inf, xend = 1000, yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme(
    strip.background = element_blank(),
    axis.line = element_blank(),
    
    plot.title = element_text(color = "black", 
                              hjust = 0.5, vjust = 1, size = 15, face = 'bold'),
    strip.text = element_text(color = 'black', size = 13, face = 'bold'),
    axis.title.y = element_text(color = 'black',size = 13, face = 'bold'),
    axis.title.x = element_text(color = 'black', size = 13, face = 'bold',
                                margin = margin(t = 15), hjust = .2),
    axis.text.x = element_text(color = 'black', 
                               size = 12),
    axis.text.y  = element_text(color = 'black', 
                                size = 12, hjust = 1),
    legend.text = element_text(color = "black", 
                               size = 12),
    legend.key.size = unit(0.9, 'cm'),
    legend.title = element_blank(),
    
    legend.position = "right"); ern_crn_erp


# ----- SAVE the plot
# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_6a.pdf', 
#                    ern_crn_erp, base_height = 5, base_width = 12)


ern_crn_erp <- ggplot(filter(Ave_ERP, Electrode == 'Cz', # <- select electrode to plot
                             Time >= -500), 
                      
                      aes(Time, M_Amp, fill = Reaction, linetype = Group, color = Reaction)) + 
  
  annotate('text', x = -510, y = -1,
           label = 'paste(bold(Cz))', parse = TRUE, 
           size = 3, hjust = 0) +
  
  facet_wrap(~ Group, scales = 'free_y', ncol = 4) + 
  theme_classic() + 
  
  annotate("rect", xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, alpha = .1) +
  geom_vline(xintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  geom_hline(yintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  
  geom_ribbon(aes(ymin = M_Amp - se, ymax = M_Amp + se), alpha = .4, colour = NA) +
  geom_line( size = rel(1)) + 
  
  scale_y_reverse(breaks = c(-8 ,-4, 0, 4, 8)) +
  scale_x_continuous(breaks = c(-500, -250, 0, 250, 500, 750, 1000)) +
  
  scale_color_viridis(option = 'C', discrete = T, direction = 1, end = .6) +
  scale_fill_viridis(option = 'C', discrete = T, direction = 1, end = .6) +
  
  labs(x = expression(bold("Time (ms)")), 
       y = expression(bold(paste("Amplitude (", mu, "V)"))), 
       title="Grand Average ERPs") +
  
  geom_segment(aes(x = -Inf, y = -8, xend = -Inf, yend = 8), 
               color = 'black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = -500, y = Inf, xend = 1000, yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme(
    strip.background = element_blank(),
    axis.line = element_blank(),
    
    plot.title = element_text(color = "black", 
                              hjust = 0.5, vjust = 1, size = 15, face = 'bold'),
    strip.text = element_text(color = 'black', size = 13, face = 'bold'),
    axis.title.y = element_text(color = 'black',size = 13, face = 'bold'),
    axis.title.x = element_text(color = 'black', size = 13, face = 'bold',
                                margin = margin(t = 15), hjust = .2),
    axis.text.x = element_text(color = 'black', 
                               size = 12),
    axis.text.y  = element_text(color = 'black', 
                                size = 12, hjust = 1),
    legend.text = element_text(color = "black", 
                               size = 12),
    legend.key.size = unit(0.9, 'cm'),
    legend.title = element_blank(),
    
    legend.position = "right"); ern_crn_erp


# ----- SAVE the plot
# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_6b.pdf', 
#                    ern_crn_erp, base_height = 5, base_width = 12)



# # ------ 5) Prepare data for plot -----------------------------------
# # Calculate charthesian coordiantes
# chanlocs$radianTheta <- pi/180*chanlocs$theta
# 
# chanlocs <- chanlocs %>%
#   mutate(x = .$radius*sin(.$radianTheta),
#          y = .$radius*cos(.$radianTheta))
# 
# to_p <- merge(Ave_ERP, select(chanlocs, x, y, Electrode), 'Electrode')
# to_p <- select(to_p, Electrode, x, y, Time, M_Amp, Reaction, Group)
# to_p <- filter(to_p, Time >= 0 & Time <= 100, Reaction == 'Correct', Group == 'Competition')
# names(to_p) <- gsub(names(to_p),
#                     pattern = '([[:upper:]])',
#                     perl = TRUE,
#                     replacement = '\\L\\1')
# names(to_p)[5] <- 'amplitude'
# 
# # ------ 5) Plot Topographical Plot -----------------------------------
# t_plot  <- topoplot(to_p, contour = T, 
#                     chan_marker = 'none', 
#                     palette = 'B', limits = c(-8, 8), 
#                     grid_res = 100); t_plot
# 
# 
# 
# # SAVE PLOT
# save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_7b.pdf', 
#           t_plot, base_height = 5, base_width = 6)
