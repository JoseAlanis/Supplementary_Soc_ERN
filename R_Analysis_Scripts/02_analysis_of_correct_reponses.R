##### ##### #####     Analysis scrips for Alanis et al., 2018   ##### ##### #####
#                           LMER models for behavioral 
#                             correct responses data

# Get helper functions
source('~/Documents/GitHub/Supplementary_Soc_ERN/R_Functions/getPacks.R')
source('/Documents/GitHub/Supplementary_Soc_ERN/R_Functions/stdResid.R')
source('~/Documents/GitHub/Supplementary_Soc_ERN/R_Functions/data_summary.R')


# Install and load multiple R packages necessary for analysis.
pkgs <- c('dplyr', 
          'lme4', 'lmerTest', 'sjstats',
          'effects', 'emmeans', 'car', 'MuMIn',
          'ggplot2', 'viridis', 'sjPlot')

getPacks(pkgs)
rm(pkgs)


# ------ READ in the data  --------------------------------------
load('~/Documents/Experiments/soc_ftask/data_for_r/Corrects_Data.RData')


# ------ 1) PLOT RT by flankers ---------------------------------
# ------ Distribution 
ggplot(Corrects, aes(M_RT, fill = Flankers)) +
  geom_histogram(color = 'black', bins = 15) + 
  facet_wrap(~ Flankers) +
  labs(x = '\n Reaction Time', y = 'Frequency \n') + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12, 
                                  colour = 'black', 
                                  face = 'bold'),
        axis.title.x = element_text(size = 14, 
                                    color = 'black',
                                    face = 'bold'),
        axis.title.y = element_text(size = 14, 
                                    color = 'black', 
                                    face = 'bold'),
        axis.text = element_text(size = 11, 
                                 color = 'black'),
        legend.position = 'none') + 
  geom_rug()


# ----- BOX-PLOT for Manuscript
corr_box <- ggplot(Corrects, 
                   aes(x = Flankers, y = M_RT, color = Flankers)) +
  
  geom_violin(fill = NA, 
              color='black', trim = T) + 
  
  geom_jitter(width = 0.35, 
              alpha = 0.7,
              size = 0.7) +
  
  stat_summary(fun.data = data_summary, color = 'black', shape = 23, 
               fill='black', size = .3) +
  
  labs(x = 'Trial Type', 
       y = expression(bold('Mean RT (ms)'))) + 
  
  theme_classic() +

  scale_y_continuous(breaks = c(100, 200, 300, 400, 500)) +
  
  geom_segment(aes(x = -Inf, y = 100, xend = -Inf, yend = 500), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = 'Compatible', y = -Inf, xend = 'Neutral', yend = -Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  scale_color_viridis(option = 'B', discrete = T, end = .9) + 
  
  theme(strip.background = element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_text(size = 14, 
                                    color = 'black',
                                    face = 'bold', margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, 
                                    color = 'black', 
                                    face = 'bold',  margin = margin(r = 15)),
        axis.text.x = element_text(size = 13, 
                                   color = 'black', angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 13, 
                                   color = 'black'),
        legend.position = 'none'); corr_box


# --- SAVE PLOT
save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_3a.pdf', 
          corr_box, base_height = 5, base_width = 3.5)


# ------ 2) EFFECT CODE and center predictors -------------------
# --- Effect code
contrasts(Corrects$Group) <- contr.sum(2)
contrasts(Corrects$Flankers) <- contr.sum(4)

# --- Center around zero
Corrects <- within(Corrects, {
  I_Errors <- Incomp_Errors - mean(Incomp_Errors, na.rm=T )
  Tot_Errors <- Total_Errors - mean(Total_Errors, na.rm=T )
})

# ------ 3) COMPUTE descriptive statistics ----------------------

# --- RT overall
Corrects %>% dplyr::summarise(M = mean(M_RT), 
                              SD = sd(M_RT), 
                              SE = sd(M_RT) / sqrt(sum(!is.na(M_RT))))

# --- RT by group
Corrects %>% group_by(Group) %>% 
  dplyr::summarise(M = mean(M_RT), 
                   SD = sd(M_RT), 
                   SE = sd(M_RT) / sqrt(sum(!is.na(M_RT))))

# --- RT by flankers
as.data.frame(Corrects %>% group_by(Flankers) %>% 
                dplyr::summarise(M = mean(M_RT), 
                                 SD = sd(M_RT), 
                                 SE = sd(M_RT) / sqrt(sum(!is.na(M_RT)))))

# --- PLOT distribution
hist(Corrects$M_RT)
rug(Corrects$M_RT)




# ------ 4) SET UP and FIT full model ---------------------------

# --- FULL MODEL with interactions, controlling for
# --- overall number of errors
mod_full <- lmer(data = Corrects, 
                 M_RT ~ Interest + Tot_Errors + 
                   Group*Flankers*Affiliation + 
                   Group*Flankers*Agency + (1|ID), 
                 REML = F)
anova(mod_full)

# --- Identify outlying observations
c_rm <- stdResid(data = Corrects, model = mod_full, plot = T, 
                 main = expression('Residuals ' ['LMER model for RT data']), 
                 xlab = expression('Fitted Values ' ['Mean RT']),
                 ylab = 'Std. Pearson Residuals', show.bound = T)

# --- Re-fit without outliers
mod_full_1 <- lmer(data = filter(c_rm, Outlier == 0), 
                   M_RT ~ Interest + Tot_Errors + 
                     Group*Flankers*Affiliation + 
                     Group*Flankers*Agency + (1|ID), 
                   REML = F)
anova(mod_full_1)


# ------ 5) SET UP and FIT the reported models ------------------
# --- MODEL including personality variables
mod_corrects <- lmer(data = Corrects, 
                     M_RT ~ Tot_Errors + Interest + 
                       Group + Flankers + Affiliation + Agency + (1|ID), 
                     REML = F)
anova(mod_corrects)
qqPlot(resid(mod_corrects, 'pearson'))

# --- Remove outliers
c_rm <- stdResid(data = Corrects, model = mod_corrects, plot = T, 
                 main = expression('Residuals ' ['LMER model for RT data']), 
                 xlab = expression('Fitted Values ' ['Mean RT']),
                 ylab = 'Std. Pearson Residuals', show.bound = T)

# --- Re-fit without outliers
mod_corrects_1 <- lmer(data = filter(c_rm, Outlier == 0), M_RT ~ Tot_Errors + Interest + 
                         Group + Flankers + Affiliation + Agency + (1|ID), REML = F)
anova(mod_corrects_1)
summary(mod_corrects_1)
car::qqPlot(resid(mod_corrects_1))

# --- Compare models with and without interactions
anova(mod_full, mod_corrects) # Interaction doesnt improve model


# --- Build tables
sjPlot::sjt.lmer(mod_full_1, mod_corrects_1,
                 show.aic = TRUE, p.numeric = FALSE,
                 string.est = 'Estimate',
                 string.ci = 'Conf. Int.',
                 string.p = 'p-value',
                 depvar.labels = c('Reaction Time', 'Reaction Time'))

sjPlot::sjt.lmer(mod_corrects_1, cell.spacing = 0.1,
                  show.aic = TRUE, p.numeric = FALSE,
                  string.est = 'Estimate',
                  string.ci = 'Conf. Int.',
                  string.p = 'p-value',
                  depvar.labels = c('Reaction Time'),
                  pred.labels = c('N Errors', '∆Motivation',
                                 'Competition', 'Compatible', 'Incompatible',
                                 'Neutral', 'Affiliation', 'Agency')
)

# Coefficient of determination
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_corrects_1)


# ------ 6) Follow UP analyses for RT model -------------------------
# --- Summary of simple slopes
rt_grid <- ref_grid(mod_corrects_1)
summary(rt_grid, infer = T)

# -- Save trial type estimates
rt_flank <- emmeans(mod_corrects_1, pairwise ~ Flankers, 
                    adjust = 'bonferroni', 
                    lmer.df = 'satterthwaite')

# --- Estimated marginal means for trail type
as.data.frame(rt_flank$contrasts)

# --- Compute CIs
confint(rt_flank)


# --- Save group estimates
rt_group <- emmeans(mod_corrects_1, pairwise ~ Group, 
                    type = 'response', 
                    adjust = 'bonferroni')

# --- Estimated marginal means for group
rt_group
# --- Compute CIs
confint(rt_group)


# --- Save affiliation estimates
rt_sc <- emtrends(mod_corrects_1, ~ 1, var='Affiliation')
rt_sc

# --- Save agency estimates
rt_ae <- emtrends(mod_corrects_1, ~ 1, var='Agency')
rt_ae



# ------ 7) Create FIGURES ------------------------------------------

# ------ 8) Plot trial type estimates -------------------------------
# --- PLOT estimates
corr_p <- ggplot(data = as.data.frame(rt_flank$emmeans), 
                aes(y = emmean, x = Flankers)) + 
  
  geom_hline(yintercept = mean(as.data.frame(rt_flank$emmeans)[, 2]), 
             linetype = 2) +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Flankers), 
                width = 0.2, size = 1, alpha = 0.7) + 
  geom_linerange(aes(ymin = emmean-SE, ymax = emmean+SE), color = 'gray',
                 size = 3) + 
  geom_point(shape = 18, size = 3.5, color = 'black') +
  labs(y = expression(bold('Estimated Mean RT (ms)')), x = 'Trial Type') +
  
  theme_classic() + 
  
  scale_y_continuous(breaks = c(250, 300, 350)) +
  
  scale_color_viridis(option = 'B', discrete = T, end = .90) +
  
  geom_segment(aes(x = -Inf, y = 250, xend = -Inf, yend = 350), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = 'Compatible', y = -Inf, xend = 'Neutral', yend = -Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme(strip.background = element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_text(size = 14, 
                                    color = 'black',
                                    face = 'bold', margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, 
                                    color = 'black', 
                                    face = 'bold',  margin = margin(r = 15)),
        axis.text.x = element_text(size = 13, 
                                   color = 'black', angle = 90, vjust = .5, hjust = 1),
        axis.text.y = element_text(size = 13, 
                                   color = 'black'),
        legend.position = 'none'); corr_p

# --- SAVE PLOT
save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_3b.pdf', 
          corr_p, base_height = 5, base_width = 2.5)


# --- Effect of Number of Errors 
dat_I <- filter(c_rm, Outlier == 0)
dat_I$pred <- predict(mod_corrects_1)

corr_err <- ggplot(dat_I, aes(x = Tot_Errors, y = pred, color=Flankers)) +
  
  annotate('text', x = -10, y = 450,
           label = 'paste(italic(ß), \' = -1.4***\')', parse = TRUE, size = 5) +
  
  geom_point(size = 1, alpha = .7) +
  
  geom_smooth(method='glm', color = 'black', fill='black', alpha = .2) +
  
  coord_cartesian(ylim = c(150, 450), xlim = c(-25, 50))  +

  theme_classic() +
  scale_color_viridis(option = 'B', end = .90,
                     discrete = T) +

  scale_x_continuous(breaks = c(-25, 0, 25, 50)) +
  scale_y_continuous(breaks = c(150, 250, 350, 450)) +

  geom_segment(aes(x = -Inf, y = 150, xend = -Inf, yend = 450),
               color='black', size=rel(1)) +
  geom_segment(aes(x = -25, y = -Inf, xend = 50, yend = -Inf),
               color='black', size=rel(1)) +
  
  labs(x = expression(bold('N Errors ' ['centred'])), 
       y = expression(bold('Mean RT (ms)'))) +
  
  theme(axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15)),
        legend.position = 'none'); corr_err

# --- SAVE PLOT
save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_3c.pdf', 
          corr_err, base_height = 5, base_width = 4)


# --- Effect of ∆Motivation
int_eff <- ggplot(dat_I, aes(x = Interest, y = pred, color=Flankers)) +
  
  annotate('text', x = -7, y = 450,
           label = 'paste(italic(ß), \' = 2.1*\')', parse = TRUE, size = 5) +
  
  geom_point(size = 1, alpha = .7) +
  
  geom_smooth(method='glm', color = 'black', fill='black', alpha = .2) +
  
  coord_cartesian(ylim = c(150, 450), xlim = c(-10, 10))  +

  theme_classic() +
  scale_color_viridis(option = 'B', end = .90,
                     discrete = T) +

  scale_x_continuous(breaks = c(-10, -5, 0, 5, 10)) +
  scale_y_continuous(breaks = c(150, 250, 350, 450)) +

  geom_segment(aes(x = -Inf, y = 150, xend = -Inf, yend = 450),
               color='black', size=rel(1)) +
  geom_segment(aes(x = -10, y = -Inf, xend = 10, yend = -Inf),
               color='black', size=rel(1)) +
  
  labs(x = expression(bold(paste(Delta,'Motivation ' ['centred']))), 
       y = expression(bold('Mean RT ' ['ms']))) +
  
  theme(axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15)),
        legend.position = 'none'); int_eff


# --- SAVE PLOT
save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_3d.pdf', 
          int_eff, base_height = 5, base_width = 4)

