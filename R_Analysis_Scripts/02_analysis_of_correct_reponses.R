##### ##### #####     Analysis scrips for Alanis et al., 2018   ##### ##### #####
#                           LMER models for behavioral 
#                             correct responses data

# Get helper functions
source('./R_Functions/getPacks.R')
source('./R_Functions/stdResid.R')
source('./R_Functions/data_summary.R')

# Install and load multiple R packages necessary for analysis.
pkgs <- c('dplyr', 'plyr',
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
  
  labs(x = 'Trial type', 
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
# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_3a.pdf', 
#           corr_box, base_height = 5, base_width = 3.5)


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
                 M_RT ~ Motivation + Tot_Errors + 
                   Group*Flankers*Affiliation + 
                   Group*Flankers*Agency + (1|ID), 
                 REML = F)
# Anova table
anova(mod_full)

# --- Identify outlying observations
c_rm <- stdResid(data = Corrects, model = mod_full, plot = T, 
                 main = expression('Residuals ' ['LMER model for RT data']), 
                 xlab = expression('Fitted Values ' ['Mean RT']),
                 ylab = 'Std. Pearson Residuals', show.bound = T)

# --- Re-fit without outliers
mod_full_1 <- lmer(data = filter(c_rm, Outlier == 0), 
                   M_RT ~ Motivation + Tot_Errors + 
                     Group*Flankers*Affiliation + 
                     Group*Flankers*Agency + (1|ID), 
                   REML = F)
# Anova table
anova(mod_full_1)


# ------ 5) SET UP and FIT the reported models ------------------
# --- MODEL including personality variables
mod_corrects <- lmer(data = Corrects, 
                     M_RT ~ Tot_Errors + Motivation + 
                       Group + Flankers + Affiliation + Agency + (1|ID), 
                     REML = F)
# Anova table
anova(mod_corrects)
# Residuals ok?
qqPlot(resid(mod_corrects, 'pearson')) # yes

# --- Remove outliers
c_rm <- stdResid(data = Corrects, model = mod_corrects, plot = T, 
                 main = expression('Residuals ' ['LMER model for RT data']), 
                 xlab = expression('Fitted Values ' ['Mean RT']),
                 ylab = 'Std. Pearson Residuals', show.bound = T)

# --- Re-fit without outliers
mod_corrects_1 <- lmer(data = filter(c_rm, Outlier == 0), 
                       M_RT ~ Tot_Errors + Motivation + 
                         Group + Flankers + Affiliation + Agency + (1|ID), 
                       REML = F)
# Anova table and coefficients
anova(mod_corrects_1)
summary(mod_corrects_1)
# Residuals ok?
qqPlot(resid(mod_corrects_1)) # yes

# Semi-Partial R2 (Edwards, et al., 2008)
# ((df numerator / df denominatot) x F) / 1 + ((df numerator / df denominatot) x F)
((1/75.831)*0.7702)/(1+((1/75.831)*0.7702)) # Group
((1/75.793)*0.7707)/(1+((1/75.793)*0.7707)) # Affiliation
((1/75.852)*2.8693)/(1+((1/75.852)*2.8693)) # Agency
((3/223.891)*162.6877)/(1+((3/223.891)*162.6877)) # Flankers
((1/75.819)*48.4859)/(1+((1/75.819)*48.4859)) # Errors
((1/75.819)*5.3411)/(1+((1/75.819)*5.3411)) # Motivation

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
                    adjust = 'fdr', 
                    lmer.df = 'satterthwaite')

# --- Estimated marginal means for trail type
rt_flank
as.data.frame(rt_flank$contrasts)

# --- Compute CIs
confint(rt_flank)


# --- Save group estimates
rt_group <- emmeans(mod_corrects_1, pairwise ~ Group, 
                    type = 'response', 
                    adjust = 'fdr')

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

# --- Save motivation estimates
rt_mo <- emtrends(mod_corrects_1, ~ 1, var='Motivation')
rt_mo

# --- Save motivation estimates
rt_err <- emtrends(mod_corrects_1, ~ 1, var='Tot_Errors')
rt_err



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
  labs(y = expression(bold('Estimated mean RT (ms)')), x = 'Trial type') +
  
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


# --- SAVE PLOT ---
# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_3b.pdf', 
#           corr_p, base_height = 5, base_width = 2.5)


# --- Effect of Number of Errors 
dat_p <- filter(c_rm, Outlier == 0)
dat_p$pred <- predict(mod_corrects_1)

dat_I <- allEffects(mod_corrects_1, xlevels = 20)
plot(dat_I, ylim = c(100, 400))
dat_I <- as.data.frame(dat_I[[1]])

# Create plot
rt_err <- ggplot(dat_p, 
                  aes(x = Tot_Errors, y = M_RT, group = interaction(Subject, Flankers), color = Flankers)) +

  stat_summary(fun.y = mean, geom = 'point', size = 1, shape = 16,
               position = position_dodge(.5)) +
  
  geom_ribbon(data = dat_I,
              aes(ymin = lower, ymax = upper, x = Tot_Errors), 
              alpha = .2, inherit.aes = F, fill = 'black') +
  geom_line(data = dat_I, 
            aes(x = Tot_Errors, y = fit), 
            inherit.aes = F, size = 0.8, color = 'black') +
  
  coord_cartesian( ylim = c(150, 450)) +
  
  scale_color_viridis(option = 'B', end = .90,
                      discrete = T) +
  
  theme_classic() + 
  
  scale_x_continuous(breaks = c(-20, 0, 20, 40, 60)) +
  scale_y_continuous(breaks = c(150, 250, 350, 450)) + 
  
  geom_segment(aes(x = -Inf, y = 150, xend = -Inf, yend = 450),
               color='black', size=rel(1)) +
  geom_segment(aes(x = -20, y = -Inf, xend = 60, yend = -Inf),
               color='black', size=rel(1)) +
  
  labs(x = expression(bold('Num. of errors ' ['centred'])), 
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
        legend.position = 'none'); rt_err

temp <- expression(beta == -1.4~'***')

rt_err <- rt_err + annotate('text', x = 10, y = 150, 
                              label = as.character(temp), parse = T, size = 5, hjust = 0); rt_err

# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_3c.pdf', 
#                    rt_err, base_height = 4.5, base_width = 4)


# --- Effect of ∆Motivation
dat_I <- allEffects(mod_corrects_1, xlevels = 20)
dat_I <- as.data.frame(dat_I[[2]])

# Create plot
rt_mot <- ggplot(dat_p, 
                 aes(x = Motivation, y = M_RT, group = interaction(Subject, Flankers), 
                     color = Flankers)) +
  
  stat_summary(fun.y = mean, geom = 'point', size = 1, shape = 16,
               position = position_dodge(.5)) +
  
  geom_ribbon(data = dat_I,
              aes(ymin = lower, ymax = upper, x = Motivation), 
              alpha = .2, inherit.aes = F, fill = 'black') +
  geom_line(data = dat_I, 
            aes(x = Motivation, y = fit), 
            inherit.aes = F, size = 0.8, color = 'black') +
  
  coord_cartesian( ylim = c(150, 450)) +
  
  scale_color_viridis(option = 'B', end = .90,
                      discrete = T) +
  
  theme_classic() + 
  
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
        legend.position = 'none'); rt_mot


temp <- expression(beta == 2.1~'*')

rt_mot <- rt_mot + annotate('text', x = -2.5, y = 150, 
                            label = as.character(temp), parse = T, size = 5, hjust = 0); rt_mot

# cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_3d.pdf', 
#                    rt_mot, base_height = 4.5, base_width = 4)