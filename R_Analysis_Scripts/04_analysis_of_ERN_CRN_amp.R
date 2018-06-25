##### ##### #####     Analysis scrips for Alanis et al., 2018   ##### ##### #####
#                           LMER models for ERN and CRN
#                                     data

# Get helper functions
source('./R_Functions/getPacks.R')
source('./R_Functions/stdResid.R')

# Install and load multiple R packages necessary for analysis.
pkgs <- c('dplyr', 'reshape2',
          'lme4', 'lmerTest',
          'effects', 'emmeans', 'car', 'MuMIn',
          'ggplot2', 'viridis', 'sjPlot')

getPacks(pkgs)
rm(pkgs)

# ------ 1) READ IN THE DATA -------------------------------------------
load('~/Documents/Experiments/soc_ftask/data_for_r/Incomp_ERPs.RData')

# ------ 2) Select relevant observations ----------------------------
# --- Keep electrode CZ and FCz;
# --- time window 0 to 100 ms after response
ERN_CRN <- filter(ERP, Electrode == 'Cz' | Electrode == 'FCz')
ERN_CRN <- filter(ERN_CRN, Time >= 0 & Time <= 100)

# --- Descriptive statistics ---
# Amp by reaction
ERN_CRN %>% group_by(Reaction) %>% 
  dplyr::summarise(M = mean(Amplitude), 
                   SD = sd(Amplitude), 
                   SE = sd(Amplitude) / sqrt(sum(!is.na(Amplitude))))
# Amp by reaction by group
ERN_CRN %>% group_by(Reaction, Group) %>%  
  dplyr::summarise(M = mean(Amplitude), 
                   SD = sd(Amplitude), 
                   SE = sd(Amplitude) / sqrt(sum(!is.na(Amplitude))))

# Mean and SD of affialition - for simple slopes analyses
unique(select(ERN_CRN, Subject, Affiliation)) %>% 
  dplyr::summarise(M = mean(Affiliation), 
                   SD = sd(Affiliation))
# Mean and SD of agency - for simple slopes analyses
unique(select(ERN_CRN, Subject, Agency)) %>% 
  dplyr::summarise(M = mean(Agency), 
                   SD = sd(Agency))

# Mean and SD of errors - for simple slopes analyses
ERN_CRN %>% dplyr::summarise(M = mean(Total_Errors), 
                             SD = sd(Total_Errors))


# Effect code electrode variable
ERN_CRN$Electrode <- factor(ERN_CRN$Electrode)
contrasts(ERN_CRN$Electrode) <- contr.sum(2); contrasts(ERN_CRN$Electrode)
contrasts(ERN_CRN$Group) <- contr.sum(2); contrasts(ERN_CRN$Group)

# Relevel reaction variable
ERN_CRN$Reaction <- factor(ERN_CRN$Reaction, 
                           levels = c('Incorrect', 'Correct'))
# Effect code raction variable
contrasts(ERN_CRN$Reaction) <- contr.sum(2); contrasts(ERN_CRN$Reaction)


# --- Mean center number of errors ---
ERN_CRN <- within( ERN_CRN, {
  Tot_Errors <- Total_Errors - mean(Total_Errors, na.rm = T)
})


# --- SET UP AND FIT full model with
# --- interactions and control variables
mod_rns_full <- lmer(data = ERN_CRN, 
                     Amplitude ~ Tot_Errors + Motivation + 
                       Reaction*Group*Affiliation + Reaction*Group*Agency + 
                       (1|Subject/Reaction), REML = F)
# Anova table
anova(mod_rns_full)

# --- Find outliers
ERN_CRN_rm <- stdResid(data = ERN_CRN, mod_rns_full, 
                       main = expression('Residuals ' ['LMER model for phys. data']), 
                       xlab = expression('Fitted Values ' ['Mean Amp (µV)']),
                       ylab = 'Std. Pearson Residuals',
                       return.data = T, plot = T, 
                       show.bound = T)

# --- Re-fit without outliers
mod_rns_full_1 <- lmer(data = filter(ERN_CRN_rm, Outlier == 0), 
                     Amplitude ~ Tot_Errors + Motivation + 
                       Reaction*Group*Affiliation + Reaction*Group*Agency + 
                       (1|Subject/Reaction), REML = F)
# Anova table and coefficients
anova(mod_rns_full_1)
summary(mod_rns_full_1)
# Residuals ok?
qqPlot(resid(mod_rns_full_1)) # yes


# Coefficent of deteminations
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_rns_full_1)


# --- SET UP AND FIT parsimonious model
mod_ern <- lmer(data = ERN_CRN, 
                Amplitude ~ Motivation + 
                  Reaction*Group*Affiliation + Reaction*Group*Agency + 
                  (1|Subject/Reaction), REML = F)
# Anova table
anova(mod_ern)
# Residuals ok?
qqPlot(resid(mod_ern, 'pearson')) # a bit heavy on the tails.


# --- Find outliers
ERN_CRN_rm <- stdResid(data = ERN_CRN, mod_ern, 
                       main = expression('Residuals ' ['LMER model for phys. data']), 
                       xlab = expression('Fitted Values ' ['Mean Amp (µV)']),
                       ylab = 'Std. Pearson Residuals',
                       return.data = T, plot = T, 
                       show.bound = T)

# --- Re-fit without outliers
mod_ern_1 <- lmer(data = filter(ERN_CRN_rm, Outlier == 0), 
                  Amplitude ~ Motivation + 
                    Reaction*Group*Affiliation + Reaction*Group*Agency + 
                    (1|Subject/Reaction), REML = F)
# Anova table and coefficients
anova(mod_ern_1)
summary(mod_ern_1)
# Residuals ok?
qqPlot(resid(mod_ern_1, 'pearson')) # yes, better.

# Semi-Partial R2 (Edwards, et al., 2008)
# ((df numerator / df denominatot) x F) / 1 + ((df numerator / df denominatot) x F)
((1/75.996)*8.3930)/(1+((1/75.996)*8.3930)) # Motivation
((1/75.996)*272.6423)/(1+((1/75.996)*272.6423)) # Reaction
((1/75.855)*6.3513)/(1+((1/75.855)*6.3513)) # Reaction by Group
((1/75.913)*7.5238)/(1+((1/75.913)*7.5238)) # Agency by Group
((1/76.009)*6.1250)/(1+((1/76.009)*6.1250)) # Affiliation by reaction
((1/75.933)*8.2481)/(1+((1/75.933)*8.2481)) # Agency by reaction by group



# Compare models with and with-out 
# scores in control variables
anova(mod_rns_full, mod_ern)

# Coefficent of deteminations
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_ern_1)


# Build table
sjPlot::sjt.lmer(mod_rns_full_1, mod_ern_1, cell.spacing = 0.1,
                 show.aic = TRUE, p.numeric = FALSE,
                 string.est = 'Estimate',
                 string.ci = 'Conf. Int.',
                 string.p = 'p-value',
                 depvar.labels = c('Amplitude', 
                                   'Amplitude'),
                 pred.labels = c('N Errors', '∆Motivation',
                                 'Error', 'Competition', 
                                 'Affiliation', 'Agency',
                                 'Error x Competition',
                                 'Error x Affiliation',
                                 'Competition x Affiliation',
                                 'Error x Agency',
                                 'Competition x Agency',
                                 'Error x Competition x Affiliation',
                                 'Error x Competition x Agency') )


# ------ 3) Follow-up analyses model ∆ERN by group ------------------
emm_options(lmerTest.limit = 8000)

# ----- GROUP ESTIMATES --- ************
# Save group estimates
group_means <- emmeans(mod_ern_1, 
                       pairwise ~ Reaction | Group, 
                       lmer.df = 'satterthwaite', 
                       adjust='fdr', 
                       at = list(Affiliation = 0, 
                                 Agency =  0, 
                                 Motivation = 0))

# Group by reaction estimates
group_means <- emmeans(mod_ern_1,
                       pairwise ~ Group | Reaction,
                       lmer.df = 'satterthwaite', 
                       adjust='fdr', 
                       at = list(Affiliation = 0, 
                                 Agency =  0, 
                                 Motivation = 0))

# Effect of group
as.data.frame(group_means$contrasts)
group_means
# Compute CIs
confint(group_means)

# Save reaction estimates
reaction_means <- emmeans(mod_ern_1, 
                          pairwise ~ Reaction, 
                          lmer.df = 'satterthwaite',
                          adjust = 'fdr',
                          at = list(Affiliation = 0,
                                    Agency =  0, 
                                    Motivation = 0))

# Effect of reaction
reaction_means
as.data.frame(reaction_means$contrasts)
# Compute CIs
confint(reaction_means)




# ----- REACTION ESTIMATES --- ************
# Save reaction estimates by group
reaction_means <- emmeans(mod_ern_1, 
                          pairwise ~ Reaction | Group, 
                          lmer.df = 'satterthwaite', 
                          adjust = 'fdr',
                          at = list(Affiliation = 0,
                                    Agency =  0, 
                                    Motivation = 0)); reaction_means

reaction_means <- emmeans(mod_ern_1, 
                          pairwise ~ Group | Reaction, 
                          lmer.df = 'satterthwaite',
                          adjust = 'fdr',
                          at = list(Affiliation = 0,
                          Agency =  0, 
                          Motivation = 0)); reaction_means

# Effect of reaction
as.data.frame(reaction_means$emmeans)
reaction_means$contrasts
# Compute CIs
confint(reaction_means)



# ----- AFFILIATION ESTIMATES --- ************
# overall
emtrends(mod_ern_1, var = 'Affiliation', ~ 1,          
         at = list(Agency =  0,
                   Motivation = 0))

# by reaction
emtrends(mod_ern_1, var = 'Affiliation', 
         pairwise ~ Reaction, 
         lmer.df = 'satterthwaite', 
         adjust = 'fdr', 
         at = list(Agency =  0, 
                   Motivation = 0))

# by reaction by group
emtrends(mod_ern_1, var = 'Affiliation', 
         pairwise ~ Reaction | Group, 
         lmer.df = 'satterthwaite', 
         adjust = 'fdr',
         at = list(Agency =  0, 
                   Motivation = 0))

# emmeans at +1 SD aff and -1 SD aff
aff_react_group <- emmeans(mod_ern_1, 
                           pairwise ~ Reaction | Affiliation, 
                           at = list(Affiliation = c(-4, 4), 
                                     Agency =  0, 
                                     Motivation = 0), 
                           lmer.df = 'satterthwaite', 
                           adjust = 'fdr')

# Summary
aff_react_group
as.data.frame(aff_react_group$contrasts)
# CIs
confint(aff_react_group)



# ----- AGENCY ESTIMATES --- ************
# --- overall
emtrends(mod_ern_1, var = 'Agency', ~ 1,
         lmer.df = 'satterthwaite', 
         adjust = 'fdr',
         at = list(Affiliation = 0,
                   Motivation = 0))

# by reaction by group
emtrends(mod_ern_1, var = 'Agency', 
         pairwise ~ Reaction | Group, 
         lmer.df = 'satterthwaite', 
         adjust = 'fdr',
         at = list(Affiliation = 0,
                   Motivation = 0))

# emmeans at +1 SD aff and -1 SD agency
ag_react_group <- emmeans(mod_ern_1, 
                          pairwise ~ Reaction | Group * Agency, 
                          at = list(Agency = c(-14, 14), 
                                    Affiliation = 0,
                                    Motivation = 0), 
                          lmer.df = 'satterthwaite',
                          adjust = 'fdr')

# Summary
ag_react_group
as.data.frame(ag_react_group$contrasts)
# --- CIs
confint(ag_react_group)



# # ----- Refit for simple slopes of agency ***********
# --- ERN in Competiton
contrasts(ERN_CRN_rm$Reaction) <- contr.treatment(2, base = 1)
contrasts(ERN_CRN_rm$Reaction)

contrasts(ERN_CRN_rm$Group) <- contr.treatment(2, base = 1)
contrasts(ERN_CRN_rm$Group)

mod_ern_1 <- lmer(data = filter(ERN_CRN_rm, Outlier == 0),
                  Amplitude ~ Motivation +
                    Reaction*Group*Affiliation + Reaction*Group*Agency +
                    (1|Subject/Reaction), REML = F)
summary(mod_ern_1)

# --- CRN in Competiton
contrasts(ERN_CRN_rm$Reaction) <- contr.treatment(2, base = 2)
contrasts(ERN_CRN_rm$Reaction)

contrasts(ERN_CRN_rm$Group) <- contr.treatment(2, base = 1)
contrasts(ERN_CRN_rm$Group)

mod_ern_1 <- lmer(data = filter(ERN_CRN_rm, Outlier == 0),
                  Amplitude ~ Motivation +
                    Reaction*Group*Affiliation + Reaction*Group*Agency +
                    (1|Subject/Reaction), REML = F)
summary(mod_ern_1)

# --- ERN in Cooperation
contrasts(ERN_CRN_rm$Reaction) <- contr.treatment(2, base = 1)
contrasts(ERN_CRN_rm$Reaction)

contrasts(ERN_CRN_rm$Group) <- contr.treatment(2, base = 2)
contrasts(ERN_CRN_rm$Group)

mod_ern_1 <- lmer(data = filter(ERN_CRN_rm, Outlier == 0),
                  Amplitude ~ Motivation +
                    Reaction*Group*Affiliation + Reaction*Group*Agency +
                    (1|Subject/Reaction), REML = F)
summary(mod_ern_1)


# --- CRN in Cooperation
contrasts(ERN_CRN_rm$Reaction) <- contr.treatment(2, base = 2)
contrasts(ERN_CRN_rm$Reaction)

contrasts(ERN_CRN_rm$Group) <- contr.treatment(2, base = 2)
contrasts(ERN_CRN_rm$Group)

mod_ern_1 <- lmer(data = filter(ERN_CRN_rm, Outlier == 0),
                  Amplitude ~ Motivation +
                    Reaction*Group*Affiliation + Reaction*Group*Agency +
                    (1|Subject/Reaction), REML = F)
summary(mod_ern_1)



# ------ 4) PLOT Simple Slopes of Affiliation -----------------------

# --- Get effects of affiliation
dat_I <- allEffects(mod_ern_1, xlevels = 20)
# --- Quick plot
plot(dat_I, ylim = c(10, -10))

# --- Prepare data for plot
dat_I <- as.data.frame(dat_I[[2]])
dat_p <-filter(ERN_CRN_rm, Outlier == 0)
levels(dat_I$Reaction) <- c('CRN', 'ERN', NA)
levels(dat_p$Reaction) <- c('ERN', 'CRN', NA)


# --- Create Plot - ERN as a function of Affilaition by Group
ern_aff <- ggplot(filter(dat_p, Reaction == 'ERN'), 
                  aes(x = Affiliation, y = Amplitude, 
                      group = Subject, color = Group)) +
  
  stat_summary(fun.y = mean, geom = 'point', 
               size = 1, shape = 16, position = position_dodge(1)) +
  
  facet_wrap(~ Reaction) +
  
  geom_ribbon(data = filter(dat_I, Reaction == 'ERN'),
              aes(ymin = lower, ymax = upper, x = Affiliation, fill = Group), 
              alpha = .2, inherit.aes = F) +
  geom_line(data = filter(dat_I, Reaction == 'ERN'), 
            aes(x = Affiliation, y = fit, color = Group, linetype = Group), 
            inherit.aes = F, size = 0.8) +
  
  scale_y_reverse(breaks = c(10, 5, 0, -5, -10)) + 
  coord_cartesian( ylim = c(10, -10)) +
  
  scale_color_viridis(option = 'A', discrete = T, 
                      direction = -1, begin = .05, end = .6) +
  scale_fill_viridis(option = 'A', discrete = T, 
                     direction = -1, begin = .05, end = .6) +
  
  theme_classic() + 
  scale_x_continuous(breaks = c(-12, -6, 0, 6)) + 
  
  geom_segment(aes(x = -Inf, y = 10, xend = -Inf, yend = -10), 
               color = 'black', size = rel(1)) +
  geom_segment(aes(x = -12, y = Inf, xend = 6, yend = Inf), 
               color = 'black', size = rel(1)) +
  
  labs(x = expression(bold('Affiliation' [' centred'])), 
       y = expression(bold(paste('Estimated amplitdue (', mu, 'V)')))) +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, color = 'black', face = 'bold'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15)),
        legend.position = 'bottom', legend.key.width = unit(1, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)); ern_aff

temp <- expression(beta ['comp']== -0.26~' ' ~ beta ['coop']== -0.03) 

ern_aff <- ern_aff + annotate('text', x = -10.5, y = 10, 
                            label = as.character(temp), parse = T, size = 5, hjust = 0); ern_aff


# --- SAVE PLOT
# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_8a.pdf', 
#                    ern_aff, base_height = 5, base_width = 4.5)




# --- Create Plot - CRN as a function of Affilaition by Group
crn_aff <- ggplot(filter(dat_p, Reaction == 'CRN'), 
                  aes(x = Affiliation, y = Amplitude, group = Subject, color = Group)) +
  
  stat_summary(fun.y = mean, geom = 'point', size = 1, shape = 16, position = position_dodge(1)) +
  
  facet_wrap(~ Reaction) +
  
  geom_ribbon(data = filter(dat_I, Reaction == 'CRN'),
              aes(ymin = lower, ymax = upper, x = Affiliation, fill = Group), 
              alpha = .2, inherit.aes = F) +
  geom_line(data = filter(dat_I, Reaction == 'CRN'), 
            aes(x = Affiliation, y = fit, color = Group, linetype = Group), 
            inherit.aes = F, size = 0.8) +
  
  scale_y_reverse(breaks = c(10, 5, 0, -5, -10)) + 
  coord_cartesian( ylim = c(10, -10)) +
  
  scale_color_viridis(option = 'A', discrete = T, 
                      direction = -1, begin = .05, end = .6) +
  scale_fill_viridis(option = 'A', discrete = T, 
                     direction = -1, begin = .05, end = .6) +
  
  theme_classic() + 
  
  scale_x_continuous(breaks = c(-12, -6, 0, 6)) + 
  
  geom_segment(aes(x = -Inf, y = 10, xend = -Inf, yend = -10), 
               color = 'black', size = rel(1)) +
  geom_segment(aes(x = -12, y = Inf, xend = 6, yend = Inf), 
               color = 'black', size = rel(1)) +
  
  labs(x = expression(bold('Affiliation' [' centred'])), 
       y = expression(bold(paste('Estimated amplitdue (', mu, 'V)')))) +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, color = 'black', face = 'bold'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15)),
        legend.position = 'bottom', legend.key.width = unit(1, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)); crn_aff

temp <- expression(beta ['comp']== -0.12~' ' ~ beta ['coop']== -.03) 

crn_aff <- crn_aff + annotate('text', x = -10.5, y = 10, 
                              label = as.character(temp), parse = T, size = 5, hjust = 0); crn_aff


# --- SAVE PLOT
# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_8b.pdf', 
#                    crn_aff, base_height = 5, base_width = 4.5)



# ------ 4) PLOT Simple Slopes of Agency ----------------------------
# --- Get effects of affiliation
dat_I <- allEffects(mod_ern_1, xlevels = 20)
# --- Quick plot
plot(dat_I, ylim = c(10, -10))

# --- Prepare data for plot
dat_I <- as.data.frame(dat_I[[3]])
dat_p <-filter(ERN_CRN_rm, Outlier == 0)
levels(dat_I$Reaction) <- c('CRN', 'ERN', NA)
levels(dat_p$Reaction) <- c('ERN', 'CRN', NA)


# --- Create Plot - ERN as a function of Agency by Group
ern_ag <- ggplot(filter(dat_p, Reaction == 'ERN'), 
                  aes(x = Agency, y = Amplitude, 
                      group = Subject, color = Group)) +
  
  stat_summary(fun.y = mean, geom = 'point', 
               size = 1, shape = 16, position = position_dodge(1)) +
  
  facet_wrap(~ Reaction) +
  
  geom_ribbon(data = filter(dat_I, Reaction == 'ERN'),
              aes(ymin = lower, ymax = upper, x = Agency, fill = Group), 
              alpha = .2, inherit.aes = F) +
  geom_line(data = filter(dat_I, Reaction == 'ERN'), 
            aes(x = Agency, y = fit, color = Group, linetype = Group), 
            inherit.aes = F, size = 0.8) +
  
  scale_y_reverse(breaks = c(10, 5, 0, -5, -10)) + 
  coord_cartesian( ylim = c(10, -10)) +
  
  scale_color_viridis(option = 'A', discrete = T, 
                      direction = -1, begin = .05, end = .6) +
  scale_fill_viridis(option = 'A', discrete = T, 
                     direction = -1, begin = .05, end = .6) +
  
  theme_classic() + 
  scale_x_continuous(breaks = c(-45, -30, -15, 0, 15, 30)) + 
  
  geom_segment(aes(x = -Inf, y = 10, xend = -Inf, yend = -10), 
               color = 'black', size = rel(1)) +
  geom_segment(aes(x = -45, y = Inf, xend = 30, yend = Inf), 
               color = 'black', size = rel(1)) +
  
  labs(x = expression(bold('Agency' [' centred'])), 
       y = expression(bold(paste('Estimated amplitdue (', mu, 'V)')))) +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, color = 'black', face = 'bold'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15)),
        legend.position = 'bottom', legend.key.width = unit(1, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)); ern_ag



ern_ag <- ern_ag + annotate('text', x = -40, y = 10,
                              label = expression(paste(beta [comp], ' = 0.04        ', 
                                                       beta [coop], ' = -0.06')), 
                              parse = TRUE, 
                              size = 5, hjust = 0)


# --- SAVE PLOT
cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_9a.pdf', 
                   ern_ag, base_height = 5, base_width = 4.5)




# --- Create Plot - CRN as a function of Affilaition by Group
crn_ag <- ggplot(filter(dat_p, Reaction == 'CRN'), 
                  aes(x = Agency, y = Amplitude, 
                      group = Subject, color = Group)) +
  
  stat_summary(fun.y = mean, geom = 'point', size = 1, shape = 16, 
               position = position_dodge(1)) +
  
  facet_wrap(~ Reaction) +
  
  geom_ribbon(data = filter(dat_I, Reaction == 'CRN'),
              aes(ymin = lower, ymax = upper, x = Agency, fill = Group), 
              alpha = .2, inherit.aes = F) +
  geom_line(data = filter(dat_I, Reaction == 'CRN'), 
            aes(x = Agency, y = fit, color = Group, linetype = Group), 
            inherit.aes = F, size = 0.8) +
  
  scale_y_reverse(breaks = c(10, 5, 0, -5, -10)) + 
  coord_cartesian( ylim = c(10, -10)) +
  
  scale_color_viridis(option = 'A', discrete = T, 
                      direction = -1, begin = .05, end = .6) +
  scale_fill_viridis(option = 'A', discrete = T, 
                     direction = -1, begin = .05, end = .6) +
  
  theme_classic() + 
  
  scale_x_continuous(breaks = c(-45, -30, -15, 0, 15, 30)) + 
  
  geom_segment(aes(x = -Inf, y = 10, xend = -Inf, yend = -10), 
               color = 'black', size = rel(1)) +
  geom_segment(aes(x = -45, y = Inf, xend = 30, yend = Inf), 
               color = 'black', size = rel(1)) +
  
  labs(x = expression(bold('Agency' [' centred'])), 
       y = expression(bold(paste('Estimated amplitdue (', mu, 'V)')))) +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, color = 'black', face = 'bold'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15)),
        legend.position = 'bottom', legend.key.width = unit(1, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

crn_ag <- crn_ag + annotate('text', x = -40, y = 10,
                            label = expression(paste(beta [comp], ' = -0.03          ', 
                                                     beta [coop], ' = 0')), 
                            parse = TRUE, 
                            size = 5, hjust = 0)


# --- SAVE PLOT
cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_9b.pdf', 
                   crn_ag, base_height = 5, base_width = 4.5)
