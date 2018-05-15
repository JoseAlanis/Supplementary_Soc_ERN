##### ##### #####     Analysis scrips for Alanis et al., 2018   ##### ##### #####
#                          GLMER models for behavioral
#                                 error data

# Get helper functions
source('~/Documents/GitHub/Supplementary_Soc_ERN/R_Functions/getPacks.R')
source('~/Documents/GitHub/Supplementary_Soc_ERN/R_Functions/stdResid.R')
source('~/Documents/GitHub/Supplementary_Soc_ERN/R_Functions/overDisp.R')
source('~/Documents/GitHub/Supplementary_Soc_ERN/R_Functions/data_summary.R')

# Install and load multiple R packages necessary for analysis.
pkgs <- c('dplyr', 
          'lme4', 'lmerTest',
          'effects', 'emmeans', 'car', 'MuMIn',
          'ggplot2', 'viridis')

getPacks(pkgs)
rm(pkgs)

# ------ READ in the data  --------------------------------------
load('~/Documents/Experiments/soc_ftask/data_for_r/Errors_Data.RData')


# ------ 1) Create PLOTs with descriptive data ------------------

# ------ Distribution
ggplot(Errors, aes(N_Errors, fill = Flankers)) +
  geom_histogram(color = 'black', bins = 9) + 
  facet_wrap(~ Flankers, scales = 'free_x') +
  labs(x = '\n Number of Errors', y = 'Frequency \n') + 
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


# ------ BOX-PLOT for Manuscript 
err_box <-  ggplot(Errors, 
                   aes(x = Flankers, y = N_Errors, color = Flankers)) +
  
  geom_violin(fill = NA, 
              color='black', trim = T) + 
  
  geom_jitter(width = 0.35, 
              alpha = 0.7,
              size = 0.7) +
  
  stat_summary(fun.data = data_summary, color = 'black', shape = 23, fill='black', size = .3) +

  labs(x = 'Trial Type', 
       y = 'Number of Errors') + 
  theme_classic() +
  
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 50), 
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
        legend.position = 'none'); err_box 


# SAVE PLOT
save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_2a.pdf', 
          err_box, base_height = 5, base_width = 3.5)



# ------ 2) COMPUTE and descriptive statistics  --------------------

# Errors overall
Errors %>% summarise(M = mean(Total_Errors), 
                     SD = sd(Total_Errors), 
                     SE = sd(Total_Errors) / sqrt(sum(!is.na(Total_Errors))))

# Errors by group
Errors %>% group_by(Group) %>% 
  summarise(M = mean(Total_Errors), 
            SD = sd(Total_Errors), 
            SE = sd(Total_Errors) / sqrt(sum(!is.na(Total_Errors))))

# Errors by flankers
as.data.frame(Errors %>% group_by(Flankers) %>% 
                summarise(M = mean(N_Errors), 
                          SD = sd(N_Errors), 
                          SE = sd(N_Errors) / sqrt(sum(!is.na(N_Errors)))))

# PLOT distribution
hist(Errors$Total_Errors)
rug(Errors$Total_Errors)


# ------ 3) SET UP and FIT full model ------------------------------

# # We document the results of a full model including all variables 
# # and interactions for the sake of completness. However, interpreting 
# # the estimates of the full model might be complicated. The more 
# # parsimonious models in section 4) and 7) yiel the same results
# # and their interpretation is straight foward.
# # Uncomment and run this section to fit the full model.

# # (UNCOMENT TO RUN):
# # DUMMY CODE cathegorical predictors
# contrasts(Errors$Group) <- contr.sum(2); contrasts(Errors$Group)
# contrasts(Errors$Flankers) <- contr.sum(4); contrasts(Errors$Flankers)
# 
# # FIT full model including all variables and interactions
# mod_err_full <-  glmer(data = Errors,
#                         N_Errors ~ Interest + Group*Flankers*Agency +
#                           Group*Flankers*Affiliation + (1|ID),
#                         family = poisson(link = 'log'), nAGQ = 20,
#                         control = glmerControl(optimizer="bobyqa"))
# Anova(mod_err_full, type = 'III')
# qqPlot(resid(mod_err_full))
# 
# # Coefficient of determination
# # R2m = only fixed effects, R2c = with random effects
# r.squaredGLMM(update(mod_err_full, nAGQ = 1))


# ------ 4) SET UP and FIT reported models  ------------------------

# ---- DUMMY CODE cathegorical predictors 
contrasts(Errors$Group) <- contr.sum(2); contrasts(Errors$Group)
contrasts(Errors$Flankers) <- contr.sum(4); contrasts(Errors$Flankers)

# ---- FIT model with interaction and controlling for interst
mod_err_inter <- glmer(data = Errors,
                    N_Errors ~ Interest + Flankers*Group + (1|ID),
                    family = poisson(link = 'log'),
                    control = glmerControl(optimizer="bobyqa"), nAGQ = 20)
Anova(mod_err_inter, type = 'III')

# ---- Compute resiudals and detect outliers
er_rm <- stdResid(data = Errors, model = mod_err_inter, plot = T, 
                 main = expression('Residuals ' ['Poisson model for error data']), 
                 xlab = expression('Fitted Values ' ['N Errors']),
                 ylab = 'Std. Pearson Residuals', show.bound = T, show.loess = T)

# ---- Re-fit without outliers
mod_err_inter_1 <- glmer(data = filter(er_rm, Outlier == 0),
                       N_Errors ~ Interest + Flankers*Group + (1|ID),
                       family = poisson(link = 'log'),
                       control = glmerControl(optimizer="bobyqa"), nAGQ = 20)
Anova(mod_err_inter_1, type = 'III')



# ----- FIT model with-out interaction
mod_errors <- glmer(data = Errors, 
                    N_Errors ~ Flankers + Group + (1|ID), 
                    family = poisson(link = 'log'),
                    control = glmerControl(optimizer="bobyqa"), nAGQ = 20)
Anova(mod_errors, type='III')
summary(mod_errors)
qqPlot(resid(mod_errors))

# Compute resiudals and detect outliers
e_rm <- stdResid(data = Errors, model = mod_errors, plot = T, 
         main = expression('Residuals ' ['Poisson model for error data']), 
         xlab = expression('Fitted Values ' ['N Errors']),
         ylab = 'Std. Pearson Residuals', show.bound = T, show.loess = T)

# Re-fit without outliers
mod_errors_1 <- glmer(data = filter(e_rm, Outlier == 0),  
                      N_Errors ~ Flankers + Group + (1|ID), 
                      family = poisson(link = 'log'), 
                      control = glmerControl(optimizer="bobyqa"), nAGQ = 20)
Anova(mod_errors_1, type='III')
summary(mod_errors_1)
qqPlot(resid(mod_errors_1))


# Coefficent of deteminations
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(update(mod_err_inter_1, nAGQ = 1)) # fit model by Laplace approximation
r.squaredGLMM(update(mod_errors_1, nAGQ = 1)) # fit model by Laplace approximation

# Compare models with and without interaction
anova(mod_err_inter, mod_errors) # Interaction doesn't improve the model

# Build table
sjPlot::sjt.glmer(mod_err_inter_1, mod_errors_1, exp.coef = F,
                  show.aic = TRUE, p.numeric = FALSE,
                  string.est = "Estimate",
                  string.ci = "Conf. Int.",
                  string.p = "p-value",
                  depvar.labels = c("Number of Errors", "Number of Errors"),
                  pred.labels = c("∆Motivation", "Compatible",
                                  "Incompatible", "Identical", 'Competition',
                                  "Compatible:Competition", "Incompatible:Competition", 
                                  "Identical:Competition" ))

# Check for overdisperion
overDisp(mod_errors_1)



# ------ 5) PAIRWISE CONTRASTS for error model ---------------------

# Summary of simple slopes
err_grid <- ref_grid(mod_errors_1)
summary(err_grid, infer=T)

# Save trial type estimates
est_err <- emmeans(mod_errors_1, pairwise ~ Flankers,
                   adjust = 'bonferroni')
as.data.frame(est_err$contrasts)
# Effect of trial type
# and CIs
confint(est_err)


# Save group estimates
est_group <- emmeans(mod_errors_1, pairwise ~ Group,
                   adjust = 'bonferroni')

# Effect of group
as.data.frame(est_group$contrasts)

# Compute CIs
confint(est_group)



# ------ 6) Plot log estimates for trial type ----------------------

# PLOT log-estimates
err_p <- ggplot(data = as.data.frame(emmeans(mod_errors_1, ~ Flankers)), 
       aes(y = emmean, x = Flankers)) + 
  
  geom_hline(yintercept = mean(as.data.frame(emmeans(mod_errors_1, ~ Flankers))[, 2]), 
             linetype=2) +
  
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL, color = Flankers), 
                width = 0.2, size = 1, alpha=0.7) + 
  geom_linerange(aes(ymin = emmean-SE, ymax = emmean+SE), color = 'gray',
                 size = 3) + 
  geom_point(shape = 18, size = 3.5, color = 'black') +
  
  labs(y = expression(bold('Estimated Incidence ' ['log scaled'])),
       x = 'Trial Type') +

  theme_classic() + 
  
  scale_color_viridis(option = 'B', discrete = T, end = .90) +
  
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 3), 
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
        legend.position = 'none'); err_p


# SAVE PLOT
save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_2b.pdf', 
          err_p, base_height = 5, base_width = 2.5)






# ------ 7) FOLLOW-UP analyses - Incompatible trials ---------------

# ----- Subset of data - only incompatible triasl
Errors_In <- filter(Errors, Flankers == 'Incompatible')

# ----- Effect code predictors
contrasts(Errors_In$Group) <- contr.sum(2); contrasts(Errors_In$Group)

# ----- FIT a full model controlling for ∆Motivation
mod_full_err_in <- glmer(data = Errors_In, N_Errors ~
                      Interest + Group*Affiliation + Group*Agency + (1|ID),
                    family = poisson(link = 'log'),
                    control = glmerControl(optimizer="bobyqa"), nAGQ = 20)
# ANOVA
car::Anova(mod_full_err_in, type='III')

# ----- Detect outlying observations
Ini_rm <- stdResid(data = Errors_In, mod_full_err_in, 
                  return.data = T, plot = T, 
                  show.loess = T, show.bound = T,
                  main = expression('Residuals ' ['Poisson model for in. error data']), 
                  xlab = expression('Fitted Values ' ['N Errors']),
                  ylab = 'Std. Pearson Residuals')

# ----- RE-FIT without outliers
mod_full_err_in_1 <- glmer(data = filter(Ini_rm, Outlier == 0), N_Errors ~
                           Interest + Group*Affiliation + Group*Agency + (1|ID),
                         family = poisson(link = 'log'),
                         control = glmerControl(optimizer="bobyqa"), nAGQ = 20)
# ANOVA
car::Anova(mod_full_err_in_1, type='III')


# FIT the model with-out ∆Motivation
mod_err_in <- glmer(data = Errors_In, 
                        N_Errors ~ Group + Affiliation + Agency + (1|ID), 
                    family = poisson(link = 'log'), 
                    control = glmerControl(optimizer = "bobyqa"), 
                    nAGQ = 20)
# ANOVA
Anova(mod_err_in, type='III')
qqPlot(resid(mod_err_in))

# Detect outlying observations
In_rm <- stdResid(data = Errors_In, mod_err_in, 
                  return.data = T, plot = T, 
                  show.loess = T, show.bound = T,
                  main = expression('Residuals ' ['Poisson model for in. error data']), 
                  xlab = expression('Fitted Values ' ['N Errors']),
                  ylab = 'Std. Pearson Residuals')

# RE-FIT without outliers
mod_err_in_1 <- glmer(data = filter(In_rm, Outlier == 0), 
                      N_Errors ~ Group + Affiliation + Agency + (1|ID), 
                      family = poisson(link = 'log'), 
                      control = glmerControl(optimizer = "bobyqa"), 
                      nAGQ = 20)
Anova(mod_err_in_1, type='III')
summary(mod_err_in_1)
qqPlot(resid(mod_err_in_1))


# Coefficent of deteminations
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(update(mod_full_err_in_1, nAGQ = 1)) # fit model by Laplace approximation
r.squaredGLMM(update(mod_err_in_1, nAGQ = 1)) # fit model by Laplace approximation

# Compare models with and without interaction
anova(mod_full_err_in, mod_err_in) # Interactions doesn't improve the model

# Build table
sjPlot::sjt.glmer(mod_full_err_in_1, mod_err_in_1, exp.coef = F,
                  show.aic = TRUE, p.numeric = FALSE,
                  string.est = "Estimate",
                  string.ci = "Conf. Int.",
                  string.p = "p-value",
                  depvar.labels = c("Number of Incomp. Errors", "Number of Incomp. Errors"),
                  pred.labels = c("∆Motivation", "Competition",
                                  "Affiliation", "Agency", 'Competition x Affiliation',
                                  'Competition x Agency') )

# Check for overdisperion
overDisp(mod_err_in_1)

# Save simple slopes of Affiliation
emm_trend_SC <- emtrends(mod_err_in_1, 
         var= 'Affiliation', 
         ~ 1)

# Effects of Affiliation
summary(emm_trend_SC, type = 'response')

# Save Simple slopes of Agency
emm_trend_MAE <- emtrends(mod_err_in_1, 
                      var= 'Agency', ~ 1,
                      transform ='response')

# Effect of Agency
summary(emm_trend_MAE, infer=T)



# ------ 8) CREATE Follow-Up FIGURES --------------------------------
# ----- Number of Incomp. Errors x Affiliation ------

dat_I <- filter(In_rm, Outlier == 0)
dat_I$pred <- predict(mod_err_in_1)


err_aff <- ggplot(dat_I, aes(x = Affiliation, y = pred)) +
  
  annotate("text", x = -9, y = 4,
           label = "paste(italic(ß), \" = 0.05*\")", parse = TRUE, 
           size = 5) +
  
  geom_point(colour="#6B186EFF", size = 1, alpha = .7) +
  
  geom_smooth(method="glm", color = 'black', fill="#6B186EFF", alpha = .2) +
  
  coord_cartesian(ylim = c(1, 4))  +
  
  theme_classic() + 
  
  scale_x_continuous(breaks = c(-12, -6, 0, 6)) + 
  scale_y_continuous(breaks = c(1, 2, 3, 4)) +
  
  geom_segment(aes(x = -Inf, y = 1, xend = -Inf, yend = 4), 
               color='black', size=rel(1)) +
  geom_segment(aes(x = -12, y = -Inf, xend = 6, yend = -Inf), 
               color='black', size=rel(1)) +
  
  labs(x =expression(bold('Affiliation' [' centred'])), 
       y = expression(bold('Est. Incidence of Incom. Errors ' ['log scaled']))) +
  
  theme(axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15)),
        legend.position = 'none'); err_aff

# SAVE PLOT
save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_2c.pdf', 
          err_aff, base_height = 5, base_width = 4)


# ----- Number of Incomp. Errors x Agency -----
err_ag <- ggplot(dat_I, aes(x = Agency, y = pred)) +
  
  annotate("text", x = -40, y = 45,
           label = "paste(italic(ß), \" = 0\")", parse = TRUE, size = 5) +
  
  geom_point(colour="#6B186EFF", size = 1, alpha = .7) +
  
  geom_smooth(method="glm", color = 'black', fill="#6B186EFF", alpha = .2) +
  
  coord_cartesian(ylim = c(1, 4))  +
  
  theme_classic() + 
  
  scale_x_continuous(breaks = c(-50, -25, 0, 25)) + 
  scale_y_continuous(breaks = c(1, 2, 3, 4)) +
  
  geom_segment(aes(x = -Inf, y = 1, xend = -Inf, yend = 4), 
               color='black', size=rel(1)) +
  geom_segment(aes(x = -50, y = -Inf, xend = 25, yend = -Inf), 
               color='black', size=rel(1)) +
  
  labs(x =expression(bold('Agency' [' centred'])), 
       y = expression(bold('Est. Incidence of Incom. Errors ' ['log scaled']))) +
  
  theme(axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15)),
        legend.position = 'none'); err_ag

# SAVE PLOT
save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_2d.pdf', 
          err_ag, base_height = 4, base_width = 3.5)



# ----- 10) ANALYSE RT - Incompatible trials -----------------------
Errors_In <- as.data.frame(Errors_In, row.names = 1:76)

## Check distribution
hist(Errors_In$M_RT, breaks=10)
rug(Errors_In$M_RT)

## Fit Model 
mod_err_RT <- lm(data = Errors_In, 
                 M_RT ~ Interest + Group + Agency + Affiliation)
anova(mod_err_RT)
summary(mod_err_RT)
plot(mod_err_RT)

## Remove Outliers? --> No effects
mod_err_RT <- lm(data = Errors_In[-c(18, 37), ],
                 M_RT ~ Interest + Group + Agency + Affiliation)
anova(mod_err_RT)
summary(mod_err_RT)
plot(mod_err_RT)

## COMPUTE CONFIDENCE INTERVALLS
confint(mod_err_RT)


