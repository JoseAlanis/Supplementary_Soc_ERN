# --- author: jose c. garcia alanis
# --- encoding: utf-8
# --- r version: 3.5.1 (2018-07-02) -- "Feather Spray"
# --- script version: Dez 2018
# --- content: analysis of post-error slowing

# --- 1) set paths and get workflow functions ----------------------------------
# path to project
setwd('/Volumes/TOSHIBA/manuscrips_and_data/soc_ern/')

# workflow functions
source('./r_functions/getPacks.R')
source('./r_functions/stdResid.R')
source('./r_functions/overDisp.R')
source('./r_functions/dataSummary.R')

# load multiple packages necessary for analysis
getPacks(c('dplyr', 'ggplot2', 'viridis'))


# --- 2) Import the data -------------------------------------------------------
errors <- read.table('./data_for_r/errors_data.txt')
ID <- as.integer(unique(gsub(errors$id, 
                             pattern = 'data_', 
                             replacement = '')))
ID <- sort(ID)

# row-bind individial data frames
for (i in ID) {
  # Function start time
  if (i == first(ID)) {
    starttime <-  Sys.time()
  }
  # 
  file = paste('./revision/rev_data/rev_PES/data_', i, '_PES.txt', sep = '')
  dat <- read.table(file, header = T)
  
  # Data from all subjects to data.frame.
  if (!exists('PES')) {
    PES <- data.frame()
  }
  PES <- rbind(PES, data.frame(dat))
  
  # Print when done.
  print(paste(i, 'done'))
  # Print function run time
  if (i == last(ID)) {
    print( paste('Runtime:',
                 as.numeric(difftime(Sys.time(), starttime, units = 'min')),
                 'mins'))
  }
}

# check numbers
unique(PES$ID)
# save data frame
write.table(PES, file = './data_for_r/post_error_data.txt', row.names = F)

# --- 3) Get personality data --------------------------------------------------
PES <- read.table('./data_for_r/post_error_data.txt', header = T)

# variable names
names(PES) <- tolower(names(PES))
# arrange
PES <- PES %>% 
  arrange(id, positioninsubject, desc(trial_index))

# summarise
# correct rt
PES_sum <- PES %>% 
  group_by(id, flankers, group, trial_index) %>%
  summarise(mean_rt = mean(rt),
            mean_corr_rt = unique(m_corr_rt),
            mean_amp = mean(m_amp))
# error amp
errors <- PES_sum %>% 
  filter(trial_index == 'Error') %>%
  group_by(id) %>%
  summarise(err_amp = mean(mean_amp)) %>%
  select(id, err_amp)
# post error rt
corrects_post_error <- PES_sum %>% 
  filter(trial_index == 'Correct post-error') %>%
  mutate(pes = mean_rt - mean_corr_rt) %>%
  select(id, flankers, group, pes)
  
PES <- merge(corrects_post_error, errors, 'id')

# personality data
perso <- read.table('./data_for_r/all_perso.txt', 
                    header = T)
unique(perso$id)

# merge with personality data
PES <- merge(select(PES, -group), perso, 'id')


# --- 3) Effect code variables -------------------------------------------------
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn', 'emmeans'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum","contr.poly"))

# effect code
contrasts(PES$group) <- contr.sum(2); contrasts(PES$group)
contrasts(PES$flankers) <- contr.sum(4); contrasts(PES$flankers)

mean(PES$err_amp)
sd(PES$err_amp)
PES <- within( PES, {
  err_amp <- err_amp - mean(err_amp, na.rm = T) 
})

# ----- 3) Fit the initial model ------
mod_pes <- lmer(data = PES, pes ~ flankers +  agency + 
                  err_amp + group*affiliation + 
                  (1|id), REML = F)
anova(mod_pes)
summary(mod_pes)

# Indentify outliers
pes_rm <- stdResid(data = PES, model = mod_pes, plot = T, show.bound = T)

# Re-fit model without outliers
mod_pes_1 <- lmer(data = filter(pes_rm, Outlier == 0),  
                pes ~ flankers +  agency + 
                  err_amp + group*affiliation + 
                  (1|id), REML = F)
# anova table
anova(mod_pes_1)
summary(mod_pes_1)
# CIs
confint(mod_pes_1)

# summary
tab_model(mod_pes_1, show.std = T)

# ----- 4) Follow-up analyses ---------
emm_options(lmerTest.limit = 2000)

# quick plot
plot_err_amp <- emmip(mod_pes_1, ~ err_amp, 
      at = list(agency = 0, 
                affiliation = 0,
                err_amp = c(-3.5, 3.5)),
      CIs = T) +
  labs(y = 'Estimates post-error slowing (ms)',
       x = 'Levels of ERN-amplitude',
       title = 'Simple slope of ERN-amplitude') +
  theme_classic() +
  coord_cartesian(ylim = c(0, 35)) +
  scale_x_discrete(breaks = c(-3.5, 3.5), labels = c('high', 'low')) +
  theme(plot.title = element_text(face = 'bold', hjust = .5),
        axis.text = element_text(color = 'black', size = 13), 
        axis.title.x = element_text(color = 'black', face = 'bold', size = 14, 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(color = 'black', size = 14, face = 'bold',
                                    margin = margin(r = 15))); plot_err_amp

# save plot
ggsave(plot_err_amp, 
       filename = './paper_figs/Fig_S2.pdf', 
       device = 'pdf',  width = 4, height = 5)

# quick plot
plot_group <- emmip(mod_pes_1, ~ group, 
                      at = list(agency = 0, 
                                affiliation = 0,
                                err_amp = 0),
                      CIs = T) +
  labs(y = 'Estimates post-error slowing (ms)',
       x = 'Levels of ERN-amplitude',
       title = 'Main effect of social context') +
  theme_classic() +
  coord_cartesian(ylim = c(0, 35)) +
  #scale_x_discrete(breaks = c(-3.5, 3.5), labels = c('high', 'low')) +
  theme(plot.title = element_text(face = 'bold', hjust = .5),
        axis.text = element_text(color = 'black', size = 13), 
        axis.title.x = element_text(color = 'black', face = 'bold', size = 14, 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(color = 'black', size = 14, face = 'bold',
                                    margin = margin(r = 15))); plot_group

# save plot
ggsave(plot_group, 
       filename = './paper_figs/Fig_S3.pdf', 
       device = 'pdf',  width = 4, height = 5)


# quick plot
plot_aff <- emmip(mod_pes_1, group ~ affiliation, 
                    at = list(agency = 0, 
                              affiliation = c(-4, 4),
                              err_amp = 0),
                    CIs = T) +
  labs(y = 'Estimates post-error slowing (ms)',
       x = 'Levels of affiliation',
       title = 'Simple slopes of affiliation') +
  theme_classic() +
  coord_cartesian(ylim = c(0, 40)) +
  scale_color_viridis(option = 'A', end = .6, direction = -1, discrete = T) + 
  scale_x_discrete(breaks = c(-4, 4), labels = c('low', 'high')) +
  theme(plot.title = element_text(face = 'bold', hjust = .5),
        axis.text = element_text(color = 'black', size = 13), 
        axis.title.x = element_text(color = 'black', face = 'bold', size = 14, 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(color = 'black', size = 14, face = 'bold',
                                    margin = margin(r = 15)),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.direction = 'horizontal'); plot_aff

# save plot
ggsave(plot_aff, 
       filename = './paper_figs/Fig_S4.pdf', 
       device = 'pdf',  width = 4, height = 5)



# Effect of Group
group_est <- emmeans(mod_pes_1, pairwise ~ affiliation | group, 
                     adjust = 'bonferroni',
                     at = list(agency = 0, 
                               affiliation = c(-4, 4), 
                               err_amp = 0)); group_est
# Compute CIs
confint(group_est)
# Plot results
emmip(mod_pes_1, ~ group, cov.reduce = FALSE, CIs = T, 
      at = list(affiliation = 0, err_amp = 0, agency = 0))

emmip(mod_pes_1, ~ err_amp, cov.reduce = FALSE, CIs = T, 
      at = list(affiliation = 0, err_amp = c(-4, 4), agency = 0))

emmip(mod_pes_1, affiliation ~ group, cov.reduce = FALSE, CIs = T, 
      at = list(affiliation = c(-4, 4), err_amp = 0, agency = 0))

# Effect of previous error amplitude
emtrends(mod_pes, var = 'Err_Amp', ~ 1, lmer.df = 'satterthwaite', 
         at=list(Agency = 0, Affiliation = 0))
# Effect of affiliation
emtrends(mod_pes, var = 'Affiliation', ~ 1, lmer.df = 'satterthwaite', 
         at=list(Err_Amp = 0, Agency = 0))
# Effect of agency
emtrends(mod_pes, var = 'Agency', ~ 1, lmer.df = 'satterthwaite', 
         at=list(Err_Amp = 0, Affiliation = 0))

# Effect of Affilaiton * ERN by Group 
aff_group <- emmeans(mod_pes, pairwise ~ Err_Amp | Affiliation + Group , 
                     adjust = 'fdr',  
                     at = list(Agency = 0, 
                               Affiliation = c(-4, 4), 
                               Err_Amp = c(-6.5, 2.5))); aff_group
# CIs
confint(aff_group)


# ----- 5) Plot effects ---------------

# PES by ERN amplitude
plot_pes1<- emmip(mod_pes, ~ Err_Amp, mult.name = "Group", cov.reduce = FALSE,
                  at = list(Err_Amp = c(-6.5, 2.5)),
                  CIs = T); plot_pes1

plot_pes1 <- plot_pes1 + coord_cartesian(ylim = c(0, 40)) +
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 40), 
               color = 'black', size = rel(1)) +
  geom_segment(aes(x = '-6.5', y = -Inf, xend = '2.5', yend = -Inf),
               color = 'black', size = rel(1)) +
  labs(y = 'Post-error slowing (ms)', x = 'Amplitude of the ERN (previous error)')+
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15))) +
  annotate('text', x = '-6.5', y = 0,
           label = expression(paste(beta, ' = -0.67*')), 
           parse = TRUE, 
           size = 5, hjust = 0); plot_pes1

# Save plot
cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_S1.pdf', 
                   plot_pes1,  base_height = 5, base_width = 4.5)


# PES by ERN amplitude, affiliation and group
plot_pes2 <- emmip(mod_pes, Affiliation ~ Err_Amp | Group, cov.reduce = FALSE, 
                   at = list(Affiliation = c(-4, 4), Err_Amp = c(-6.5, 2.5)), 
                   CIs = T); plot_pes2


plot_pes2 <- plot_pes2 + coord_cartesian(ylim = c(-20, 60)) +
  geom_segment(aes(x = -Inf, y = -20, xend = -Inf, yend = 60), 
               color = 'black', size = rel(1)) +
  geom_segment(aes(x = '-6.5', y = -Inf, xend = '2.5', yend = -Inf),
               color = 'black', size = rel(1)) +
  labs(y = 'Post-error slowing (ms)', x = 'Amplitude of the ERN (previous error)')+
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(color = 'black', size=14, face='bold'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15)),
        legend.text = element_text(color = 'black', size=12),
        legend.title = element_text(color = 'black', size = 12, face='bold')); plot_pes2

# Save plot
cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_S2.pdf', 
                   plot_pes2,  base_height = 5, base_width = 9)


