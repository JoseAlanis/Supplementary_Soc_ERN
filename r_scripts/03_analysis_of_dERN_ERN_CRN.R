# --- author: jose c. garcia alanis
# --- encoding: utf-8
# --- r version: 3.5.1 (2018-07-02) -- "Feather Spray"
# --- script version: Dez 2018
# --- content: analysis of behavioural correct reactions

# --- 1) Set paths and get workflow functions ----------------------------------
# path to project
setwd('/Volumes/TOSHIBA/manuscripts_and_data/soc_ern/')

# workflow functions
source('./r_functions/getPacks.R')
source('./r_functions/stdResid.R')
source('./r_functions/overDisp.R')
source('./r_functions/spR2.R')
source('./r_functions/dataSummary.R')
source('./r_functions/topoplot.R')

# install and load multiple R packages necessary for analysis.
getPacks(c('tidyr', 'dplyr', 'ggplot2', 'viridis'))


# get electrode locations
chanlocs <- read.table('./meta_dat/chanlocs_ftask.txt',
                       header = T, 
                       sep = ',')
# Rename column
names(chanlocs)[1] <- 'electrode'

# Calculate charthesian coordiantes topoplot
chanlocs$radianTheta <- pi/180*chanlocs$theta
chanlocs <- chanlocs %>%
  mutate(x = .$radius*sin(.$radianTheta),
         y = .$radius*cos(.$radianTheta))


# --- 3) Import data ----------------------------------------------------------
# personality data
perso <- read.table('./data_for_r/all_perso.txt', 
                    header = T)
unique(perso$id)
# Behavioural data
errors <- read.table('./data_for_r/errors_data.txt')

# physiological data
# dat_ave <- read.table('./revision/rev_data/phys_midline_incomp_ave.txt')
dat_ave <- read.table('./data_for_r/phys_all_incomp_ave.txt')

# -- bring physiological data in format for analysis ---
names(dat_ave) <- tolower(names(dat_ave))
dat_ave$flankers <- factor(dat_ave$flankers)
dat_ave$electrode <- factor(dat_ave$electrode)

# capitalize first letter in flankers variable
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
dat_ave$flankers <- firstup(as.character(dat_ave$flankers))
dat_ave$group <- firstup(as.character(dat_ave$group))
# to factor
dat_ave$flankers <- as.factor(dat_ave$flankers)
dat_ave$group <- as.factor(dat_ave$group)
# -- done --


# --- 3) compute delta-ERN (dERN) ----------------------------------------------
# (i.e., difference wave comparing correct and incorrect reponses)
response_ave <-  dat_ave %>% 
  filter(time >= 0 & time <= 100 & electrode %in% c('FCz', 'Cz'))

# spread reaction variable (to wide format)
dat_ave <- tidyr::spread(dat_ave, reaction, m_amp)
# compute difference wave
dat_ave$diff <- dat_ave$Incorrect - dat_ave$Correct


# --- 4) analyse dERN cross midline electrodes ---------------------------------
# Electrodes in question
midline <- c('Cz', 'Fz', 'FCz', 'CPz', 'Pz')

# Summarise by electrode and id
elect_id <- dat_ave %>% 
  group_by(id, electrode) %>% 
  filter((time >= 0 & time <= 100) & electrode %in% midline) %>%
  summarise(amplitude = mean(diff))

# -- effect code cathegorical predictors --
elect_id$electrode <- factor(elect_id$electrode, levels = midline)
contrasts(elect_id$electrode) <- contr.treatment(5, base = 1); contrasts(elect_id$electrode)

# -- change default contrasts options ! --
options(contrasts = c("contr.sum","contr.poly"))

# load Packages
getPacks(c('lme4', 'lmerTest', 'car', 'MuMIn', 'sjPlot', 'emmeans'))

# -- ** set up and fit by electrode model ** --
mod_elect_0 <- lmer(data = elect_id, 
                    amplitude ~ electrode + (1|id), 
                    REML = F)
anova(mod_elect_0)

# indetify outliers
elect_rm <- stdResid(data = data.frame(elect_id), model = mod_elect_0, 
                     return.data = T, plot = T,
                     show.loess = T, show.bound = T)

# re-fit without outliers
mod_elect <- lmer(data = filter(elect_rm, Outlier == 0), 
                  amplitude ~ electrode + (1|id), 
                  REML = F)
# anova table and model summary
anova(mod_elect)
summary(mod_elect)
# residuals ok?
qqPlot(resid(mod_elect))

# plot results
plot_model(mod_elect, 'std', auto.label = T)
# pairwise comparissons
emmeans(mod_elect, ~ electrode, 
        adjust = 'bonferroni', 
        contr = 'trt.vs.ctrl')
# R-squared
r.squaredGLMM(mod_elect)

# create summary table for fitte model
tab_model(file = './revision/rev_tables/mod_electrode_ERN.html',
          mod_elect, digits.p = 3, show.aic = T, show.std = T,
          title = 'Table 1: Results of linear mixed effects regression analysis of delta-ERN by electrode.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = c('Amplitude of delta-ERN'),
          pred.labels = c('(Intercept)',
                          'Fz', 
                          'FCz',
                          'CPz',
                          'Pz'))


# --- 5) Plot dERN cross midline electrodes -----------------------------------
# *****
# ***** RUN section 4 first *****
# *****
# summarise by electrode
elect_wide <- dat_ave %>% 
  group_by(time, electrode) %>% 
  summarise(amplitude = mean(diff))
# to wide format
elect_wide <- tidyr::spread(elect_wide, electrode, value = c('amplitude'))
# compute GFP
elect_wide$GFP <- apply(elect_wide[, 2:65], 1, sd)
# Back to long
elect_long <- tidyr::gather(elect_wide, electrode, amplitude, AF3:GFP)

# summarise by electrode
Ave_ro <- elect_long %>% filter(time >= 0 & time <= 100) %>% 
  group_by(electrode) %>% 
  dplyr::summarise(mean_amp = mean(amplitude)) %>% 
  arrange(mean_amp)

# order electrodes according to magnitude of dERN
elect_long$electrode <- factor(elect_long$electrode, 
                                   levels = as.character(Ave_ro$electrode))

# plot dERN-wave across midline electrodes
elec_p <- ggplot(filter(elect_long, 
                        electrode %in% c('Fz', 'FCz', 'Cz','CPz', 'Pz', 'GFP'),
                        time >= -500), 
                 
                 aes(time, amplitude, color = electrode, size = electrode)) +
  
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
                                         direction = 1, 
                                         begin = .1, 
                                         end =.90), 'black')) +
  
  scale_size_manual(values = c(rep(0.75, 5), 1), guide = F) +
  
  labs(y = expression(bold(paste("Amplitude of ", Delta, "ERN (", mu, "V)"))), 
       x = expression(bold('Time (ms)')), 
       title = "Mean evoked activity across midline electrodes" ) +
  
  geom_segment(aes(x = -Inf, y = -7.5, xend = -Inf, yend = 5), 
               color = 'black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = -500, y = Inf, xend = 1000, yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme_classic() +
  theme(
    plot.title = element_text(color = "black", size = 13, face = 'bold', 
                              vjust = .5),
    strip.background = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_text(color = 'black', size = 12, face = 'bold', 
                                margin = margin(t = 15)),
    axis.title.y = element_text(color = 'black', size = 12, face = 'bold', 
                                margin = margin(r = 15)),
    axis.text.x = element_text(color = 'black', size = 12),
    axis.text.y  = element_text(hjust = 1, color = 'black', size = 12),
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_text(color = "black", size = 12),
    legend.key.size = unit(0.5, 'cm'),
    legend.position = 'bottom', 
    legend.box = "horizontal", 
    legend.background = element_blank()); elec_p 

ggsave(elec_p, 
       filename = './paper_figs/Fig_S1a.pdf', 
       device = 'pdf',  width = 8, height = 6)

# -- ** plot topoplot 0-100 ms following motor response ** --
# data for topoplot
to_p <- merge(filter(elect_long, !electrode == 'GFP'), 
               select(chanlocs, x, y, electrode), 'electrode')
# choose time window
to_p <- filter(to_p, time >= 0 & time <= 100)
# name columns
names(to_p) <- gsub(names(to_p),
                    pattern = '([[:upper:]])',
                    perl = TRUE,
                    replacement = '\\L\\1')
# create topoplot
t_plot  <- topoplot(to_p, contour = T, 
                    chan_marker = 'none', 
                    palette = 'D', limits = c(-6,6),
                    grid_res = 100); t_plot
# add title
t_plot <- t_plot + 
  labs(title = 'Mean evoked activity') + 
  theme(plot.title = element_text(hjust = .5, size = 17, face = 'bold'),
        legend.title = element_text(hjust = .5, size = 17, face = 'bold'),
        legend.text = element_text(size = 17)); t_plot

ggsave(t_plot, 
       filename = './paper_figs/Fig_S1b.pdf', 
       device = 'pdf',  width = 5, height = 5)


# --- 6) Analyis of delta-ERN -------------------------------------------------
# only keep time window 0-100 ms after reponse
id_ave <-  dat_ave %>% 
  filter(time >= 0 & time <= 100 & electrode %in% c('FCz', 'Cz'))
# merge with behavioral data
id_ave <-  merge(filter(errors, flankers == 'Incompatible'), id_ave, c('id', 'flankers'))
# delete duplicates
id_ave <-  select(id_ave, -group)
# merge with personality
id_ave <- merge(id_ave, perso, 'id')


# ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn', 'emmeans'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum","contr.poly"))

# effect code cathegorical predictors
contrasts(id_ave$group) <- contr.sum(2); contrasts(id_ave$group)

# -- ** set up and fit "full model" (cointains all relevant variables) ** --
mod_ern_full_0 <- lmer(data = id_ave, 
                     diff ~ error_rate + d_motivation +
                       group*affiliation + 
                       group*agency + 
                       (1|id), 
                     REML = F)
# anova table
anova(mod_ern_full_0)

# identifie outlying observations
full_ave_rm <- stdResid(data = id_ave, model = mod_ern_full_0, 
                      return.data = T, plot = T, 
                      show.loess = T, show.bound = T)

# re-fit model without outliers
mod_ern_full <- lmer(data = filter(full_ave_rm, Outlier == 0), 
                     diff ~ error_rate + d_motivation +
                       group*affiliation + 
                       group*agency + 
                       (1|id), REML = F)
# anova table and model summatry
anova(mod_ern_full)
summary(mod_ern_full)
# resisuals ok?
qqPlot(resid(mod_ern_full))

# coefficent of deteminations
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_ern_full)

# effect code predictors
contrasts(id_ave$group) <- contr.sum(2); contrasts(id_ave$group)
# -- UNCOMMENT FOR DETAILED SIMPLE SLOPES --
# contrasts(id_ave$group) <- contr.treatment(2, base = 1); contrasts(id_ave$group)

# -- ** set up and fit a more parsimonious model ** --
# (cointains variables with significant contribution)
mod_ern_0 <- lmer(data = id_ave, 
                diff ~ group*affiliation + group*agency + 
                  (1|id), 
                REML = F)
# anova table and summary
anova(mod_ern_0)
summary(mod_ern_0)

# Identify outlying observations
ern_ave_rm <- stdResid(data = id_ave, model = mod_ern_0, 
                      return.data = T, plot = T, 
                      show.loess = T, show.bound = T)

# Re-fit model without outliers
mod_ern <- lmer(data = filter(ern_ave_rm, Outlier == 0), 
                diff ~ group*affiliation + group*agency + 
                  (1|id), REML = F)
# Anova table and model summary
anova(mod_ern)
summary(mod_ern)
# Residuals ok?
qqPlot(resid(mod_ern))
plot_model(mod_ern, type = 'diag')
tab_model(mod_ern, show.std = T)

# Compare models with and with-out 
# scores in control variables
anova(mod_ern_full_0, mod_ern_0)

# Coefficent of detemination
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_ern)

# compute effect sizes (semi partial R2) from anova table
amod <- anova(mod_ern); amod
amod <-  as.data.frame(amod); amod
amod$sp.R2 <- spR2(amod); amod

# descriptives
as.data.frame(id_ave %>% dplyr::group_by(group) %>% 
                dplyr::summarise(M = mean(diff), 
                                 SD = sd(diff), 
                                 SE = sd(diff) / sqrt(sum(!is.na(diff))), 
                                 Min = min(diff), Max = max(diff)))


# Create summary table for the fitted models
tab_model(file = './revision/rev_tables/delta_ERN.html',
          mod_ern_full, mod_ern,
          digits = 2,
          show.aic = T, 
          collapse.ci = T, 
          title = 'Table S5: Results of linear mixed effects regression analysis of delta-ERN by social context.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = c('Full model', 'Final model'), 
          pred.labels = c('(Intercept)',
                          'Error rate', 
                          'Engagement',
                          'Competition',
                          'Affiliation (aff)',
                          'Agency (ag)',
                          'Competition x aff',
                          'Competititon x ag'))


# --- 7) Plot simple slopes of affiliation and agency --------------------------
# allow higher number of observations
emm_options(lmerTest.limit = 4000)

# simple slopes of affiliation
emtrends(mod_ern, var = 'affiliation', lmer.df = 'satterthwaite', ~ 1, at = list(agecy = 0))
emtrends(mod_ern, var = 'affiliation', lmer.df = 'satterthwaite', pairwise ~ group, at = list(agecy = 0))

# simple slopes of agency
emtrends(mod_ern, var = 'agency', lmer.df = 'satterthwaite', ~ 1,  at = list(affiliation = 0))
emtrends(mod_ern, var = 'agency', lmer.df = 'satterthwaite', pairwise ~ group, at = list(affiliation = 0))

# save means
gr_means <- emmeans::emmeans(mod_ern, pairwise ~ group,
                                      lmer.df = 'satterthwaite',
                                      at = list(affiliation = 0, 
                                                agency = 0)); gr_means
# CIs
confint(gr_means)

# -- compute simple slopes for affiliation (i.e., low vs. high) ---
affiliation_means <- emmeans::emmeans(mod_ern, ~ group | affiliation,
                                      lmer.df = 'satterthwaite',
                                      at = list(affiliation = c(-4, 4), 
                                                agency = 0))

# contrast at low and high affiliation
pairs(affiliation_means, adjust = 'bonferroni')
confint(pairs(affiliation_means, adjust = 'bonferroni'))

# plot results ** affiliation **
pd = position_dodge(.2)
plt_aff <- ggplot(data = data.frame(affiliation_means), 
                    aes(x = as.factor(affiliation), 
                        y = emmean, 
                        color = as.factor(group),
                        group = group)) +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = .1, 
                size = .8, 
                position = pd, color='black') +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2, position = pd) + 
  scale_y_reverse() + 
  scale_x_discrete(breaks = c(-4, 4), labels = c('low', 'high')) +
  coord_cartesian(ylim = c(0, -10)) +
  scale_color_viridis(option = 'A', end = .6, discrete = T, direction = -1) +
  
  labs(title = 'Simple slopes of affiliation',
       y = expression(bold(paste("Amplitude of ", 
                                 Delta, "ERN (", mu, "V)"))),
       x = 'Levels of affiliation') +
  
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = -10), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = '-4', y = Inf, xend = '4', yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.text = element_text(color = 'black', size = 13), 
        axis.title.x = element_text(color = 'black', face = 'bold', size = 14, 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(color = 'black', size = 14),
        plot.title = element_text(face = 'bold', hjust = .5),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.direction = 'horizontal'); plt_aff

# save plot
ggsave(plt_aff, 
       filename = './paper_figs/Fig_5a.pdf', 
       device = 'pdf',  width = 4, height = 5)



# -- compute simple slopes for agency (i.e., low vs. high) --
agency_means <- emmeans::emmeans(mod_ern, ~ group | agency,
                                 lmer.df = 'satterthwaite',
                                 at = list(affiliation = 0, 
                                           agency = c(-14, 14)))

# Contrast at low and high affiliation
pairs(agency_means, adjust = 'bonferroni')
confint(pairs(agency_means, adjust = 'bonferroni'))

# plot results ** agency **
pd = position_dodge(.2)
plt_agency <- ggplot(data = data.frame(agency_means), 
                    aes(x = as.factor(agency), 
                        y = emmean, 
                        color = as.factor(group),
                        group = group)) +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = .1, 
                size = .8, 
                position = pd, color='black') +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2, position = pd) + 
  scale_y_reverse() + 
  scale_x_discrete(breaks = c(-14, 14), labels = c('low', 'high')) +
  coord_cartesian(ylim = c(0, -10)) +
  scale_color_viridis(option = 'A', end = .6, discrete = T, direction = -1) +
  
  labs(title = 'Simple slopes of agency',
       y = expression(bold(paste("Amplitude of ", 
                                 Delta, "ERN (", mu, "V)"))),
       x = 'Levels of agency') +
  
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = -10), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = '-14', y = Inf, xend = '14', yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.text = element_text(color = 'black', size = 13), 
        axis.title.x = element_text(color = 'black', face = 'bold', size = 14, 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(color = 'black', size = 14),
        plot.title = element_text(face = 'bold', hjust = .5),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.direction = 'horizontal'); plt_agency

# save plot
ggsave(plt_agency, 
       filename = './paper_figs/Fig_5b.pdf', 
       device = 'pdf',  width = 4, height = 5)



# --- 7) Plot dERN waveform ----------------------------------------------------
# dERN waveform
for_plot <- dat_ave %>% 
  filter(electrode %in% c('FCz', 'Cz')) %>% 
  group_by(group, time) %>% 
  summarise(Amplitude = mean(diff), se = sd(diff)/sqrt(length(diff))) 

# plot waveform
ERN_p <- ggplot(filter(for_plot), 
                
                aes(time, Amplitude, fill = group, color = group, linetype = group)) + 
  
  theme_classic() + 
  
  annotate("rect", xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, alpha = .1) +
  geom_vline(xintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  geom_hline(yintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  
  geom_ribbon(aes(ymin = Amplitude - se, ymax = Amplitude + se), 
              alpha = .3, colour = NA) +
  
  geom_line( size = rel(1)) + 
  scale_y_reverse(breaks = c(-7.5 ,-5, -2.5, 0, 2.5, 5)) +
  scale_x_continuous(breaks = c(-500, -250, 0, 250, 500, 750, 1000)) +
  
  scale_color_viridis(option = 'A', discrete = T, 
                      direction = -1, end = .6) +
  scale_fill_viridis(option = 'A', discrete = T, 
                     direction = -1, end = .6) +
  
  scale_linetype_manual(values = c(1, 6)) + 
  
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

# save plot
ggsave(ERN_p, 
       filename = './paper_figs/Fig_4a.pdf', 
       device = 'pdf',  width = 8, height = 5)


# --- 7) Analyse response ERPs -------------------------------------------------
# delete duplicates
response_ave <- select(response_ave, -c('group', 'flankers'))
response_ave <- merge(filter(errors, flankers == 'Incompatible'), response_ave, 'id')
response_ave <- merge(response_ave, perso, 'id')

# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn', 'emmeans'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum","contr.poly"))

# change order of levels
response_ave$reaction <-  factor(response_ave$reaction, levels = c('Incorrect', 'Correct'))

# Effect code cathegorical predictors
contrasts(response_ave$group) <- contr.sum(2); contrasts(response_ave$group)
contrasts(response_ave$reaction) <- contr.sum(2); contrasts(response_ave$reaction)

# # -- UNCOMMENT FOR DETAILED SIMPLE SLOPES --
# contrasts(response_ave$group) <- contr.treatment(2, base = 2); contrasts(response_ave$group)
# contrasts(response_ave$reaction) <- contr.treatment(2, base = 1); contrasts(response_ave$reaction)

# Set up and fit "full model" (cointains als relevant variables)
mod_response_full_0 <- lmer(data = response_ave, 
                       m_amp ~ error_rate + d_motivation +
                         reaction*group*affiliation + 
                         reaction*group*agency + 
                         (1|id/reaction), 
                       REML = F)
# Anova table
anova(mod_response_full_0)

# Identifie outlying observations
full_ave_rm <- stdResid(data = response_ave, model = mod_response_full_0, 
                        return.data = T, plot = T, 
                        show.loess = T, show.bound = T)

# Re-fit model without outliers
mod_response_full <- lmer(data = filter(full_ave_rm, Outlier == 0), 
                          m_amp ~ error_rate + d_motivation +
                            reaction*group*affiliation + 
                            reaction*group*agency + 
                            (1|id/reaction), REML = F)
# Anova table and model summatry
anova(mod_response_full)
summary(mod_response_full)
# Resisuals ok?
qqPlot(resid(mod_response_full))

tab_model(mod_response_full, show.std = T)
# Create summary table for the fitted models
tab_model(file = './revision/rev_tables/mod_ERN_CRN.html',
          mod_response_full,
          digits = 2,
          show.aic = T, 
          collapse.ci = T, 
          title = 'Table S6: Results of linear mixed effects regression analysis of ERN and CRN.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = c('Final model'), 
          pred.labels = c('(Intercept)',
                          'Error rate', 
                          'Engagement',
                          'Erroneous reaction (error)',
                          'Competition (comp)',
                          'Affiliation (aff)',
                          'Agency (ag)',
                          'Error x comp',
                          'Error x aff',
                          'Comp x aff',
                          'Error x ag',
                          'Comp x ag',
                          'Error x comp x aff',
                          'Error x comp x ag'))

# compute effect sizes (semi partial R2) from anova table
amod <- anova(mod_response_full); amod
amod <-  as.data.frame(amod); amod
amod$sp.R2 <- spR2(amod); amod

# descriptives
as.data.frame(response_ave %>% dplyr::group_by(reaction) %>% 
                dplyr::summarise(M = mean(m_amp), 
                                 SD = sd(m_amp), 
                                 SE = sd(m_amp) / sqrt(sum(!is.na(m_amp))), 
                                 Min = min(m_amp), Max = max(m_amp)))
# descriptives
as.data.frame(response_ave %>% dplyr::group_by(reaction, group) %>% 
                dplyr::summarise(M = mean(m_amp), 
                                 SD = sd(m_amp), 
                                 SE = sd(m_amp) / sqrt(sum(!is.na(m_amp))), 
                                 Min = min(m_amp), Max = max(m_amp)))

# -- pairwise comparissons --
emm_options(lmerTest.limit = 8000, infer = T) 

# reaction means
react_means <- emmeans(mod_response_full, pairwise ~ reaction,
                       adjust = 'bonferroni', lmer.df = 'satterthwaite',
                       at = list(affiliation = 0, error_rate = 0, 
                                 agency = 0, d_motivation = 0)); react_means
# CIs
confint(react_means)

# reaction by group
react_means <- emmeans(mod_response_full, pairwise ~ reaction | group,
                       adjust = 'bonferroni', lmer.df = 'satterthwaite',
                       at=list(affiliation = 0, error_rate = 0, 
                               agency = 0, d_motivation = 0)); react_means
as.data.frame(react_means$contrasts)
# CIs
confint(react_means)

react_means <- emmeans(mod_response_full, pairwise ~ group | reaction,
                       adjust = 'bonferroni', lmer.df = 'satterthwaite',
                       at=list(affiliation = 0, error_rate = 0, 
                               agency = 0, d_motivation = 0)); react_means
as.data.frame(react_means$contrasts)
# CIs
confint(react_means)


# affiliation means
aff_means <- emmeans(mod_response_full, pairwise ~ reaction | affiliation,
                     adjust = 'bonferroni', lmer.df = 'satterthwaite',
                     at = list(affiliation = c(-4, 4), 
                             agency = 0, d_motivation = 0)); aff_means
# CIs
confint(aff_means)

# affiliation means
aff_means <- emmeans(mod_response_full, pairwise ~ affiliation | reaction,
                     adjust = 'bonferroni', lmer.df = 'satterthwaite',
                     at = list(affiliation = c(-4, 4), 
                               agency = 0, d_motivation = 0)); aff_means
# CIs
confint(aff_means)


# agency means
ag_means <- emmeans(mod_response_full, pairwise ~ reaction | agency + group,
                    adjust = 'bonferroni', lmer.df = 'satterthwaite',
                    at = list(agency = c(-14, 14), 
                              affiliation = 0, d_motivation = 0)); ag_means
# CIs
confint(ag_means)

# agency means
ag_means <- emmeans(mod_response_full, pairwise ~ agency | reaction + group,
                    adjust = 'bonferroni', lmer.df = 'satterthwaite',
                    at = list(agency = c(-14, 14), 
                              affiliation = 0, d_motivation = 0)); ag_means
# CIs
confint(ag_means)



# # -- ** UNCOMMENT TO GET P- VALUES FOR SIMPLE SLOPES ** --
# # modify contrasts to get p-values 
# # in each level of reaction and group
# contrasts(full_ave_rm$group) <- contr.treatment(2, base = 1); contrasts(full_ave_rm$group)
# contrasts(full_ave_rm$group) <- contr.treatment(2, base = 2); contrasts(full_ave_rm$group)
# 
# contrasts(full_ave_rm$group) <- contr.sum(2); contrasts(full_ave_rm$group)
# contrasts(full_ave_rm$Reaction) <- contr.treatment(2, base = 2); contrasts(full_ave_rm$Reaction)
# 
# mod_response_full <- lmer(data = filter(full_ave_rm, Outlier == 0),
#                           M_Amp ~ error_rate + d_motivation +
#                             Reaction*group*affiliation +
#                             Reaction*group*agency +
#                             (1|id/Reaction), REML = F)
# # Anova table and model summatry
# anova(mod_response_full)
# summary(mod_response_full)



# --- 7) create interaction plota -------------------------------------------------
getPacks('plyr')
# -- compute simple slopes for affiliation (i.e., low vs. high) ---
affiliation_means <- emmeans::emmeans(mod_response_full, ~ affiliation | reaction + group,
                                      lmer.df = 'satterthwaite',
                                      at = list(affiliation = c(-4, 4), 
                                                d_motivation = 0,
                                                agency = 0))

affiliation_means <- data.frame(affiliation_means)
affiliation_means$reaction <- plyr::revalue(affiliation_means$reaction, c('Incorrect' = 'ERN',
                                                                          'Correct' = 'CRN'))

# plot results ** affiliation **
pd = position_dodge(.2)
plt_aff <- ggplot(data = data.frame(affiliation_means), 
                  aes(x = as.factor(affiliation), 
                      y = emmean, 
                      color = as.factor(group),
                      group = group)) +
  
  facet_wrap(~ reaction, scales = 'free') +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = .1, 
                size = .8, 
                position = pd, color='black') +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2, position = pd) + 
  scale_y_reverse(breaks=c(-7, -3.5, 0, 3.5, 7)) + 
  scale_x_discrete(breaks = c(-4, 4), labels = c('low', 'high')) +
  coord_cartesian(ylim = c(-7, 7)) +
  scale_color_viridis(option = 'A', end = .6, discrete = T, direction = -1) +
  
  labs(title = 'Simple slopes of affiliation',
       y = expression(bold(paste("Amplitude of ", 
                                 Delta, "ERN (", mu, "V)"))),
       x = 'Levels of affiliation') +
  
  geom_segment(aes(x = -Inf, y = 7, xend = -Inf, yend = -7), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = '-4', y = Inf, xend = '4', yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(color = 'black', face = 'bold', size = 14),
        axis.line = element_blank(),
        axis.text = element_text(color = 'black', size = 13), 
        axis.title.x = element_text(color = 'black', face = 'bold', size = 14, 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(color = 'black', size = 14),
        plot.title = element_text(face = 'bold', hjust = .5),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.direction = 'horizontal'); plt_aff

# save plot
ggsave(plt_aff, 
       filename = './paper_figs/Fig_7.pdf', 
       device = 'pdf',  width = 8.5, height = 5.5)



getPacks('plyr')
# -- compute simple slopes for affiliation (i.e., low vs. high) ---
agency_means <- emmeans::emmeans(mod_response_full, ~ agency | reaction + group,
                                      lmer.df = 'satterthwaite',
                                      at = list(agency = c(-14, 14), 
                                                d_motivation = 0,
                                                affiliation = 0))

agency_means <- data.frame(agency_means)
agency_means$reaction <- plyr::revalue(agency_means$reaction, c('Incorrect' = 'ERN',
                                                                          'Correct' = 'CRN'))

# plot results ** agency **
pd = position_dodge(.2)
plt_agency <- ggplot(data = data.frame(agency_means), 
                  aes(x = as.factor(agency), 
                      y = emmean, 
                      color = as.factor(group),
                      group = group)) +
  
  facet_wrap(~ reaction, scales = 'free') +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = .1, 
                size = .8, 
                position = pd, color='black') +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2, position = pd) + 
  scale_y_reverse(breaks=c(-7, -3.5, 0, 3.5, 7)) + 
  scale_x_discrete(breaks = c(-14, 14), labels = c('low', 'high')) +
  coord_cartesian(ylim = c(-7, 7)) +
  scale_color_viridis(option = 'A', end = .6, discrete = T, direction = -1) +
  
  labs(title = 'Simple slopes of agency',
       y = expression(bold(paste("Amplitude of ", 
                                 Delta, "ERN (", mu, "V)"))),
       x = 'Levels of agency') +
  
  geom_segment(aes(x = -Inf, y = 7, xend = -Inf, yend = -7), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = '-14', y = Inf, xend = '14', yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(color = 'black', face = 'bold', size = 14),
        axis.line = element_blank(),
        axis.text = element_text(color = 'black', size = 13), 
        axis.title.x = element_text(color = 'black', face = 'bold', size = 14, 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(color = 'black', size = 14),
        plot.title = element_text(face = 'bold', hjust = .5),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.direction = 'horizontal'); plt_agency

# save plot
ggsave(plt_agency, 
       filename = './paper_figs/Fig_8.pdf', 
       device = 'pdf',  width = 8.5, height = 5.5)