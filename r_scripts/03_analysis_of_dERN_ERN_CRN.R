# --- author: jose c. garcia alanis
# --- encoding: utf-8
# --- r version: 4.0.1 (2020-06-06) -- "See Things Now"
# --- script version: Jun 2020
# --- content: analysis of physiological data


# --- 1) load workflow functions and necessary packages ------------------------
# workflow functions
source('./r_functions/getPacks.R')
source('./r_functions/stdResid.R')
source('./r_functions/spR2.R')
source('./r_functions/dataSummary.R')
source('./r_functions/topoplot.R')

# load packages
getPacks(c('tidyr', 'dplyr', 'ggplot2', 'viridis'))

# load eeg channel locations
chanlocs <- read.table('../data/eeg/chanlocs_ftask.txt',
                       header = T,
                       sep = ',')
# rename column
names(chanlocs)[1] <- 'electrode'

# calculate charthesian coordiantes for topoplot
chanlocs$radianTheta <- pi/180*chanlocs$theta
chanlocs <- chanlocs %>%
  mutate(x = .$radius*sin(.$radianTheta),
         y = .$radius*cos(.$radianTheta))


# --- 2) Import data -----------------------------------------------------------
# load eeg data
dat_ave <- read.table('../data/eeg/all_electrodes_incompatible_average.txt')
# load behavioral data
errors <- read.table('../data/behavioral/behavioral_data.txt',
                    header = T)
errors <- na.omit(errors)
errors <- errors %>% select(-group)

# keep this for later
response_ave <-  dat_ave %>%
  filter(time >= 0 & time <= 100 & electrode %in% c('FCz', 'Cz'))

# --- 3) Compute delta-ERN (dERN) ----------------------------------------------
# to wide format
dat_ave <- dat_ave %>%
  pivot_wider(names_from = reaction, values_from = mean_amp)

# compute difference wave (i.e., comparing correct and incorrect reponses)
dat_ave <- dat_ave %>% mutate(diff = Incorrect - Correct)


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
tab_model(file = './results/tables/TABLE_electrode_effects_ERN.html',
          mod_elect, digits.p = 3, show.aic = T, show.std = T,
          title = 'Table: Results of linear mixed effects regression analysis of delta-ERN by electrode.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = 'Amplitude of delta-ERN',
          pred.labels = c('Intercept (electr. Cz)',
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
elect_wide <- elect_wide %>%
  pivot_wider(names_from = electrode, values_from = amplitude)
#elect_wide <- tidyr::spread(elect_wide, electrode, value = c('amplitude'))
# compute GFP
elect_wide$GFP <- apply(elect_wide[, 2:65], 1, sd)
# Back to long
elect_long <- elect_wide %>%
  pivot_longer(-time, names_to = 'electrode', values_to = 'amplitude')
# elect_long <- tidyr::gather(elect_wide, electrode, amplitude, AF3:GFP)

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
  geom_vline(xintercept = 0, color='black',
             size=rel(0.5), linetype = 3) +
  geom_hline(yintercept = 0, color='black',
             size=rel(0.5), linetype = 3) +
  
  geom_line(size = 1) +

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
       filename = './results/figures/Fig_S1a.pdf',
       device = 'pdf',  width = 8, height = 6)

# *** plot topoplot 0-100 ms following motor response ***
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
                    palette = 'RdBu', limits = c(-6,6),
                    grid_res = 100); t_plot
# add title
t_plot <- t_plot + 
  labs(title = 'Mean evoked activity') + 
  theme(plot.title = element_text(hjust = .5, size = 17, face = 'bold'),
        legend.title = element_text(hjust = .5, size = 17, face = 'bold'),
        legend.text = element_text(size = 17)); t_plot

ggsave(t_plot, 
       filename = './results/figures/Fig_S1b.pdf',
       device = 'pdf',  width = 5, height = 5)


# --- 6) Analyis of delta-ERN -------------------------------------------------
# only keep time window 0-100 ms after reponse
id_ave <-  dat_ave %>% 
  filter(time >= 0 & time <= 100 & electrode %in% c('FCz', 'Cz'))
# merge with behavioral data
id_ave <-  merge(filter(errors, flankers == 'Incompatible'), id_ave, c('id', 'flankers'))


# ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn', 'emmeans'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum", "contr.poly"))

# effect code cathegorical predictors
id_ave$group <- as.factor(id_ave$group)
contrasts(id_ave$group) <- contr.sum(2); contrasts(id_ave$group)

# grand-mean center the error rate of incompatible trials
id_ave <- id_ave %>%
  mutate(error_rate_p_conditon_centred =  error_rate_p_conditon - mean(error_rate_p_conditon))

# -- ** set up and fit "full model" (cointains all relevant variables) ** --
mod_ern_full_0 <- lmer(data = id_ave,
                       diff ~
                         error_rate_p_conditon_centred + diff_in_motivation_centred +
                           group*SC_centred +
                           group*MAE_centred + (1|id),
                       REML = F)
# anova table
anova(mod_ern_full_0)

# Compute resiudals and detect outliers
full_ave_rm <- stdResid(data = id_ave, model = mod_ern_full_0,
                      return.data = T, plot = T, 
                      show.loess = T, show.bound = T)

# In order to test the effect of outliers, re-fit the model without outliers
mod_ern_full <- lmer(data = filter(full_ave_rm, Outlier == 0),
                     diff ~
                       error_rate_p_conditon_centred + diff_in_motivation_centred +
                         group*SC_centred +
                         group*MAE_centred + (1|id),
                     REML = F)
# anova table and model summatry
anova(mod_ern_full)
summary(mod_ern_full)
# resisuals ok? They look ok.
qqPlot(resid(mod_ern_full))

# coefficent of deteminations
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_ern_full)

# effect code predictors
contrasts(id_ave$group) <- contr.sum(2); contrasts(id_ave$group)


# -- UNCOMMENT FOR DETAILED SIMPLE SLOPES --
# use contr.treatment() instead of contr.sum to look at estimates compared to
# a reference value, such as in:
# contrasts(id_ave$group) <- contr.treatment(2, base = 1); contrasts(id_ave$group)

# -- ** set up and fit a more parsimonious model ** --
# (cointains variables with significant contribution)
mod_ern_0 <- lmer(data = id_ave,
                  diff ~ group*SC_centred + group*MAE_centred + (1|id),
                  REML = F)
# anova table and summary
anova(mod_ern_0)
summary(mod_ern_0)

# Compute resiudals and detect outliers
ern_ave_rm <- stdResid(data = id_ave, model = mod_ern_0,
                      return.data = T, plot = T, 
                      show.loess = T, show.bound = T)

# In order to test the effect of outliers, re-fit the model without outliers
mod_ern <- lmer(data = filter(ern_ave_rm, Outlier == 0),
                diff ~ group*SC_centred + group*MAE_centred + (1|id),
                REML = F)
# Anova table and model summary
anova(mod_ern)
summary(mod_ern)
# Residuals ok?
qqPlot(resid(mod_ern))
# further model diagnostics
plot_model(mod_ern, type = 'diag')

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
tab_model(file = './results/tables/TABLE_S5_delta_ERN_analyis.html',
          mod_ern_full, mod_ern,
          digits = 2,
          show.aic = T, 
          collapse.ci = T,
          show.std = T,
          title = 'Table S5: Results of linear mixed effects regression analysis of delta-ERN by social context.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = c('Full model', 'Final model'), 
          pred.labels = c('(Intercept)',
                          'Error rate', 
                          'delta-Engagement',
                          'Competition',
                          'Affiliation (aff)',
                          'Agency (ag)',
                          'Competition x aff',
                          'Competititon x ag'))


# --- 7) Plot simple slopes of affiliation and agency --------------------------
# allow higher number of observations
emm_options(lmerTest.limit = 4000)

# slopes of affiliation
emtrends(mod_ern, var = 'SC_centred',
         lmer.df = 'satterthwaite', pairwise ~ group,
         at = list(MAE_centred = 0))

# slopes of agency
emtrends(mod_ern, var = 'MAE_centred',
         lmer.df = 'satterthwaite', pairwise ~ group,
         at = list(SC_centred = 0))

# effect of group (save means)
gr_means <- emmeans::emmeans(mod_ern, pairwise ~ group,
                                      lmer.df = 'satterthwaite',
                                      at = list(SC_centred = 0,
                                                MAE_centred = 0)); gr_means
# CIs
confint(gr_means)

sd(id_ave$SC_centred)
# *** compute simple slopes for affiliation (i.e., low vs. high) ***
affiliation_means <- emmeans::emmeans(mod_ern, ~ group | SC_centred,
                                      lmer.df = 'satterthwaite',
                                      at = list(SC_centred = c(-4, 4),
                                                MAE_centred = 0)); affiliation_means

# contrast at low and high affiliation
pairs(affiliation_means, adjust = 'bonferroni')
confint(pairs(affiliation_means, adjust = 'bonferroni'))

# plot results ** affiliation **
pd <- position_dodge(.2)
plt_aff <- ggplot(data = data.frame(affiliation_means),
                    aes(x = as.factor(SC_centred),
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
       filename = './results/figures/Fig_5a.pdf',
       device = 'pdf',  width = 4, height = 5)


sd(id_ave$MAE_centred)
# -- compute simple slopes for agency (i.e., low vs. high) --
agency_means <- emmeans::emmeans(mod_ern, ~ group | MAE_centred,
                                 lmer.df = 'satterthwaite',
                                 at = list(SC_centred = 0,
                                           MAE_centred = c(-14, 14)))

# Contrast at low and high affiliation
pairs(agency_means, adjust = 'bonferroni')
confint(pairs(agency_means, adjust = 'bonferroni'))

# plot results ** agency **
pd <- position_dodge(.2)
plt_agency <- ggplot(data = data.frame(agency_means), 
                    aes(x = as.factor(MAE_centred),
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
       filename = './results/figures/Fig_5b.pdf',
       device = 'pdf',  width = 4, height = 5)



# --- 8) Plot dERN waveform ----------------------------------------------------
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
       filename = './results/figures/Fig_4a.pdf',
       device = 'pdf',  width = 8, height = 5)


# --- 9) Analyse response ERPs -------------------------------------------------
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn', 'emmeans'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum","contr.poly"))

# merge with errors data and personality variables
response_ave <-  merge(filter(errors, flankers == 'Incompatible'), response_ave, c('id', 'flankers'))

# Effect code cathegorical predictors
response_ave$group <- as.factor(response_ave$group)
contrasts(response_ave$group) <- contr.sum(2); contrasts(response_ave$group)
# change order of levels
response_ave$reaction <-  factor(response_ave$reaction, levels = c('Incorrect', 'Correct'))
contrasts(response_ave$reaction) <- contr.sum(2); contrasts(response_ave$reaction)

# grand-mean center the error rate of incompatible trials
response_ave <- response_ave %>%
  mutate(error_rate_p_conditon_centred =  error_rate_p_conditon - mean(error_rate_p_conditon))

# # -- UNCOMMENT FOR DETAILED SIMPLE SLOPES --
# contrasts(response_ave$group) <- contr.treatment(2, base = 2); contrasts(response_ave$group)
# contrasts(response_ave$reaction) <- contr.treatment(2, base = 1); contrasts(response_ave$reaction)

# Set up and fit "full model" (cointains als relevant variables)
mod_response_full_0 <- lmer(data = response_ave,
                            mean_amp ~ error_rate_p_conditon_centred + diff_in_motivation_centred +
                              reaction*group*SC_centred +
                              reaction*group*MAE_centred + (1|id/reaction),
                            REML = F)
# Anova table
anova(mod_response_full_0)

# Identifie outlying observations
full_ave_rm <- stdResid(data = response_ave, model = mod_response_full_0, 
                        return.data = T, plot = T, 
                        show.loess = T, show.bound = T)

# Re-fit model without outliers
mod_response_full <- lmer(data = filter(full_ave_rm, Outlier == 0),
                          mean_amp ~ error_rate_p_conditon_centred + diff_in_motivation_centred +
                            reaction*group*SC_centred +
                            reaction*group*MAE_centred + (1|id/reaction),
                          REML = F)
# Anova table and model summatry
anova(mod_response_full)
summary(mod_response_full)
# Resisuals ok? A little heavy on the tails, but should be ok.
qqPlot(resid(mod_response_full))

tab_model(mod_response_full, show.std = T)
# Create summary table for the fitted models
tab_model(file = './results/tables/TABLE_model_ERN_vs_CRN.html',
          mod_response_full,
          digits = 2,
          show.aic = T, 
          collapse.ci = T,
          show.std = T,
          title = 'Table S6: Results of linear mixed effects regression analysis of ERN and CRN.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = 'Final model',
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
                dplyr::summarise(M = mean(mean_amp),
                                 SD = sd(mean_amp),
                                 SE = sd(mean_amp) / sqrt(sum(!is.na(mean_amp))),
                                 Min = min(mean_amp), Max = max(mean_amp)))
# descriptives
as.data.frame(response_ave %>% dplyr::group_by(reaction, group) %>% 
                dplyr::summarise(M = mean(m_amp), 
                                 SD = sd(m_amp), 
                                 SE = sd(m_amp) / sqrt(sum(!is.na(m_amp))), 
                                 Min = min(m_amp), Max = max(m_amp)))


# *** pairwise comparissons ***
emm_options(lmerTest.limit = 8000, infer = T) 


# reaction means
react_means <- emmeans(mod_response_full, pairwise ~ reaction,
                       adjust = 'bonferroni', lmer.df = 'satterthwaite',
                       at = list(SC_centred = 0, error_rate_p_conditon_centred = 0,
                                 MAE_centred = 0, diff_in_motivation_centred = 0)); react_means
# CIs
confint(react_means)

# reaction by group
react_means <- emmeans(mod_response_full, pairwise ~ reaction | group,
                       adjust = 'bonferroni', lmer.df = 'satterthwaite',
                       at = list(SC_centred = 0, error_rate_p_conditon_centred = 0,
                                 MAE_centred = 0, diff_in_motivation_centred = 0)); react_means
as.data.frame(react_means$contrasts)
# CIs
confint(react_means)

# group by reaction
react_means <- emmeans(mod_response_full, pairwise ~ group | reaction,
                       adjust = 'bonferroni', lmer.df = 'satterthwaite',
                       at = list(SC_centred = 0, error_rate_p_conditon_centred = 0,
                                 MAE_centred = 0, diff_in_motivation_centred = 0)); react_means
as.data.frame(react_means$contrasts)
# CIs
confint(react_means)


# reaction by affiliation (high vs low)
aff_means <- emmeans(mod_response_full, pairwise ~ reaction | SC_centred,
                     adjust = 'bonferroni', lmer.df = 'satterthwaite',
                     at = list(SC_centred = c(-4, 4),
                             MAE_centred = 0, diff_in_motivation_centred = 0)); aff_means
# CIs
confint(aff_means)

# affiliation (high vs low) by reaction
aff_means <- emmeans(mod_response_full, pairwise ~ SC_centred | reaction,
                     adjust = 'bonferroni', lmer.df = 'satterthwaite',
                     at = list(SC_centred = c(-4, 4),
                               MAE_centred = 0, diff_in_motivation_centred = 0)); aff_means
# CIs
confint(aff_means)


# reaction means by agency (high vs low)
ag_means <- emmeans(mod_response_full, pairwise ~ reaction | MAE_centred + group,
                    adjust = 'bonferroni', lmer.df = 'satterthwaite',
                    at = list(MAE_centred = c(-14, 14),
                              SC_centred = 0, diff_in_motivation_centred = 0)); ag_means
# CIs
confint(ag_means)

# agency means (high vs low) by reaction
ag_means <- emmeans(mod_response_full, pairwise ~ MAE_centred | reaction + group,
                    adjust = 'bonferroni', lmer.df = 'satterthwaite',
                    at = list(MAE_centred = c(-14, 14),
                              SC_centred = 0, diff_in_motivation_centred = 0)); ag_means
# CIs
confint(ag_means)



# # -- ** UNCOMMENT TO GET P-VALUES FOR SIMPLE SLOPES ** --
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



# --- 10) create interaction plots -------------------------------------------------
getPacks('plyr')
# -- compute simple slopes for affiliation (i.e., low vs. high) ---
affiliation_means <- emmeans::emmeans(mod_response_full, ~ SC_centred | reaction + group,
                                      lmer.df = 'satterthwaite',
                                      at = list(SC_centred = c(-4, 4),
                                                diff_in_motivation_centred = 0,
                                                MAE_centred = 0))

affiliation_means <- data.frame(affiliation_means)
affiliation_means$reaction <- plyr::revalue(affiliation_means$reaction, c('Incorrect' = 'ERN',
                                                                          'Correct' = 'CRN'))

# plot results ** affiliation **
pd <- position_dodge(.2)
plt_aff <- ggplot(data = data.frame(affiliation_means), 
                  aes(x = as.factor(SC_centred),
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
       filename = './results/figures/Fig_7.pdf',
       device = 'pdf',  width = 8.5, height = 5.5)



getPacks('plyr')
# -- compute simple slopes for affiliation (i.e., low vs. high) ---
agency_means <- emmeans::emmeans(mod_response_full, ~ MAE_centred | reaction + group,
                                      lmer.df = 'satterthwaite',
                                      at = list(MAE_centred = c(-14, 14),
                                                diff_in_motivation_centred = 0,
                                                SC_centred = 0))

agency_means <- data.frame(agency_means)
agency_means$reaction <- plyr::revalue(agency_means$reaction, c('Incorrect' = 'ERN',
                                                                          'Correct' = 'CRN'))

# plot results ** agency **
pd <- position_dodge(.2)
plt_agency <- ggplot(data = data.frame(agency_means), 
                  aes(x = as.factor(MAE_centred),
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
       filename = './results/figures/Fig_8.pdf',
       device = 'pdf',  width = 8.5, height = 5.5)
