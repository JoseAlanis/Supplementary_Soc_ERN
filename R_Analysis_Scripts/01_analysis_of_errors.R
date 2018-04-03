##### ##### #####     Analysis scrips for Alanis et al., 2018   ##### ##### #####
#                          GLMER models for behavioral
#                                 error data

# Get helper functions
source('./Documents/GitHub/social_ERN_2018/R_Functions/getPacks.R')
source('./Documents/GitHub/social_ERN_2018/R_Functions/stdResid.R')
source('./Documents/GitHub/social_ERN_2018/R_Functions/overDisp.R')

# Install and load multiple R packages necessary for analysis.
pkgs <- c('dplyr', 
          'lme4', 'lmerTest',
          'effects', 'emmeans', 'car', 'MuMIn',
          'ggplot2', 'cowplot', 'viridis')

getPacks(pkgs)
rm(pkgs)

# ------ READ in the data  --------------------------------------
load('./Desktop/RefTask_2018_Final/for_upload/DATA/Errors_Data.RData')


# ------ 1) PLOT number of errors by flankers ----------------------
pdf('./Desktop/number_of_errors.pdf', height = 3.5, width = 6)

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

dev.off()



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
Errors %>% group_by(Flankers) %>% 
  summarise(M = mean(N_Errors), 
            SD = sd(N_Errors), 
            SE = sd(N_Errors) / sqrt(sum(!is.na(N_Errors))))

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

# DUMMY CODE cathegorical predictors 
contrasts(Errors$Group) <- contr.sum(2); contrasts(Errors$Group)
contrasts(Errors$Flankers) <- contr.sum(4); contrasts(Errors$Flankers)

## FIT model with interaction and controlling for interst
mod_err_inter <- glmer(data = Errors,
                    N_Errors ~ Interest + Flankers*Group + (1|ID),
                    family = poisson(link = 'log'),
                    control = glmerControl(optimizer="bobyqa"), nAGQ = 20)
Anova(mod_err_inter, type = 'III')


# FIT the reported models
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


# Compare models with and without interaction
anova(mod_err_inter, mod_errors) ## Interaction does not improve the model

# Coefficent of deteminations
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(update(mod_errors_1, nAGQ = 1)) # fit model by Laplace approximation

# Check overdisperion
overDisp(mod_errors_1)

# ------ 5) PAIRWISE CONTRASTS for error model ---------------------

# Summary of simple slopes
err_grid <- ref_grid(mod_errors_1)
summary(err_grid, infer=T, type = 'response')

# Save trial type estimates
est_err <- emmeans(mod_errors_1, pairwise ~ Flankers,
                   transform = 'response',
                   adjust = 'bonferroni')

# Effect of trial type
as.data.frame(est_err$emmeans)
est_err$contrasts
# Compute CIs
mutate(as.data.frame(est_err$contrasts), 
       LCL = estimate - SE * 1.96, 
       UCL = estimate + SE * 1.96)

# Save group estimates
est_group <- emmeans(mod_errors_1, pairwise ~ Group,
                   transform = 'response',
                   adjust = 'bonferroni')

# Effect of group
as.data.frame(est_group$emmeans)
est_group$contrasts
# Compute CIs
mutate(as.data.frame(est_group$contrasts), 
       LCL = estimate - SE * 1.96, 
       UCL = estimate + SE * 1.96)



# ------ 6) Plot log estimates for trial type ----------------------
ggplot(data = as.data.frame(emmeans(mod_errors_1, ~ Flankers)), 
       aes(y = emmean, x = Flankers)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, color = 'black', size = 0.8) + 
  geom_linerange(aes(ymin = emmean-SE, ymax = emmean+SE, color = Flankers), 
                 size = 3) + 
  geom_point(shape = 18, size = 3, color = 'black') +
  labs(y = expression( 'estimated incidence ' ['log scaled'])) +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color = 'black', size = 13),
        legend.position = 'none') +
  coord_flip(ylim = c(- 0.15, 3))



# ------ 7) FOLLOW-UP analyses - Incompatible trials ---------------

# Subset of data - only incompatible triasl
Errors_In <- filter(Errors, Flankers == 'Incompatible')

# Effect code predictors
contrasts(Errors_In$Group) <- contr.sum(2); contrasts(Errors_In$Group)

## MODEL without random structure (uncomment to run)
# mod_err_int <- glm(data = Errors_In, N_Errors ~ Group*SC_Centred + Group*MAE_Centred,
#                      family = poisson(link = 'log'))
# car::Anova(mod_err_int, type='III')
# plot(mod_err_int)

# FIT a full model controlling for interest
mod_full_err_in <- glmer(data = Errors_In, N_Errors ~
                      Interest + Group*Affiliation + Group*Agency + (1|ID),
                    family = poisson(link = 'log'),
                    control = glmerControl(optimizer="bobyqa"), nAGQ = 20)
car::Anova(mod_full_err_in, type='III')


# FIT the reported models
mod_err_in <- glmer(data = Errors_In, 
                        N_Errors ~ Group + Affiliation + Agency + (1|ID), 
                    family = poisson(link = 'log'), 
                    control = glmerControl(optimizer="bobyqa"), nAGQ = 20)
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
                      control = glmerControl(optimizer="bobyqa"), nAGQ = 20)
Anova(mod_err_in_1, type='III')
summary(mod_err_in_1)
qqPlot(resid(mod_err_in_1))


# Compare models with and without angency interaction
anova(mod_full_err_in, mod_err_in) ## Interactions does not improve the model


# Coefficent of deteminations
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(update(mod_err_in_1, nAGQ = 1)) # fit model by Laplace approximation

# Check overdisperion
overDisp(mod_err_in_1)
sjPlot::sjp.glmer(mod_err_in_1, type = 'pred', vars = 'Affiliation', show.ci = T)
sjPlot::sjp.glmer(mod_err_in_1, type = 'pred', vars = 'Agency', show.ci = T)

# Save Simple slopes of Affiliation
emm_trend_SC <- emtrends(mod_err_in_1, 
         var= 'Affiliation', 
         ~ 1, 
         transform ='response')

# Effects of Affiliation
summary(emm_trend_SC, infer=T)

# Save Simple slopes of Agency
emm_trend_MAE <- emtrends(mod_err_in_1, 
                      var= 'Agency', ~ 1,
                      transform ='response')

# Effect of Agency
summary(emm_trend_MAE, infer=T)


# Save estimates
af <- Effect(mod = mod_err_in_1, 
             c("Affiliation"), 
             xlevels = list(Affiliation = 20),  
             partial.residuals=T) 

ag <- Effect(mod = mod_err_in_1, 
             c("Agency"), 
             xlevels = list(Agency = 20),  
             partial.residuals=T) 


# ------ 8) CREATE FIGURES -----------------------------------------
# NUMBER OF ERRORS by Trial Type
err_flank <- ggplot(data = as.data.frame(est_err$emmeans), 
                    aes(x = Flankers, y = rate, fill = Flankers)) + 
  
  geom_bar(size = .5, width = .7, color=NA,
           position = position_dodge(0.5), 
           stat = 'identity') +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.25, color = 'black', size = 0.8, 
                position = position_dodge(0.5)) + 
  geom_linerange(aes(ymin = rate-SE, ymax = rate+SE), 
                 size = 2, color = 'black',
                 position = position_dodge(0.5)) +
  
  theme_classic() +
  scale_fill_viridis(option = 'B', begin = .30, end = .95,
                     discrete = T) +
  
  scale_y_continuous(breaks = c(0, 5, 10, 15)) +
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 15), 
               color = 'black', size = rel(1)) +
  geom_segment(aes(x = 'Neutral', y = -Inf, xend = 'Compatible', yend = -Inf), 
               color = 'black', size = rel(1)) +
  
  labs(x ='Trial Type', 
       y = expression(bold('N ' ['errors']))) +
  
  theme(axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12, 
                                    face = 'bold', 
                                    margin = margin(t = 10)),
        axis.text.y = element_text(size = 11, color = 'black'),
        axis.text.x =  element_text(size = 11, color = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, color = 'black'),
        plot.margin=unit(c(.8,.5,.5,.5),"cm"),
        legend.position = 'none') + 
  
  coord_flip(); print(err_flank)


# Number of Incomp. Errors x Affiliation
err_aff <- ggplot(as.data.frame(af), 
                  aes(Affiliation, fit)) +
  
  annotate("text", x = -9, y = 45,
           label = "paste(italic(ß), \" = 0.62*\")", parse = TRUE, size = 3) +
  
  geom_point(data = In_rm, aes(x = Affiliation, y = N_Errors), 
             colour="#C13A50FF", size = .9, alpha = .8) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              fill='gray', color = NA, alpha = .5) + 
  geom_line(color='black', size = .8) +
  
  coord_cartesian(ylim = c(0, 45))  +
  
  theme_classic() + 
  
  scale_x_continuous(breaks = c(-12, -6, 0, 6)) + 
  scale_y_continuous(breaks = c(0, 15, 30, 45)) +
  
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 45), 
               color='black', size=rel(1)) +
  geom_segment(aes(x = -12, y = -Inf, xend = 6, yend = -Inf), 
               color='black', size=rel(1)) +
  
  labs(x =expression(bold('Affiliation' [' centred'])), 
       y = expression(bold('N ' ['incom. errors']))) +
  
  theme(axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 11, color = 'black'),
        axis.text.y = element_text(size = 11, color = 'black'),
        axis.title.x = element_text(size = 12, face = 'bold', 
                                    margin = margin(t = 10)),
        axis.title.y = element_text(size = 12, face = 'bold', 
                                    margin = margin(r = 10)),
        plot.margin=unit(c(.8,.5,.5,.5),"cm"),
        legend.position = 'none'); print(err_aff)


# Number of Incomp. Errors x Agency
err_ag <- ggplot(as.data.frame(ag), 
                  aes(Agency, fit)) +
  
  annotate("text", x = -40, y = 45,
           label = "paste(italic(ß), \" = 0\")", parse = TRUE, size = 3) +
  
  geom_point(data = In_rm, aes(x = Agency, y = N_Errors), 
             colour="#C13A50FF", size = .9, alpha = .8) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              fill='gray', color = NA, alpha = .5) + 
  geom_line(color='black', size = .8) +
  
  coord_cartesian(ylim = c(0, 45))  +
  
  theme_classic() + 
  
  scale_x_continuous(breaks = c(-50, -25, 0, 25)) + 
  scale_y_continuous(breaks = c(0, 15, 30, 45)) +
  
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 45), 
               color='black', size=rel(1)) +
  geom_segment(aes(x = -50, y = -Inf, xend = 25, yend = -Inf), 
               color='black', size=rel(1)) +
  
  labs(x =expression(bold('Agency' [' centred'])), 
       y = expression(bold('N ' ['incom. errors']))) +
  
  theme(axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 11, color = 'black'),
        axis.text.y = element_text(size = 11, color = 'black'),
        axis.title.x = element_text(size = 12, face = 'bold', 
                                    margin = margin(t = 10)),
        axis.title.y = element_text(size = 12, face = 'bold', 
                                    margin = margin(r = 10)),
        plot.margin=unit(c(.8,.5,.5,.5),"cm"),
        legend.position = 'none'); print(err_ag)


# ------ 9) SAVE FIGURES as PDF ------------------------------------
# Upper panel
first_row <- plot_grid(err_flank, nrow = 1, labels = c('a)'), 
                       vjust = 1, hjust = -.5)
# Lower panes
second_row = plot_grid(err_aff, err_ag, nrow = 1, labels = c('b)', 'c)'), 
                       vjust = 1, hjust = -.5)
# Full figure
figure_errors <- plot_grid(first_row, second_row, ncol = 1)
# Save figure
save_plot('./Desktop/Figure_Errors.pdf', 
          figure_errors, base_height = 5.5, base_width = 5.5)





# ----- 10) ANALYSE RT - Incompatible trials -----------------------
Errors_In <- as.data.frame(Errors_In, row.names = 1:76)

## Check distribution
hist(Errors_In$M_RT)

## Fit Model 
mod_err_RT <- lm(data = Errors_In, M_RT ~ Interest + Group + Agency + Affiliation)
anova(mod_err_RT)
summary(mod_err_RT)
plot(mod_err_RT)

## Remove Outliers? --> No effects
mod_err_RT <- lm(data = Errors_In[-c(18, 37), ], M_RT ~ Interest + Group + Agency + Affiliation)
anova(mod_err_RT)
summary(mod_err_RT)
plot(mod_err_RT)

## COMPUTE CONFIDENCE INTERVALLS
confint(mod_err_RT)


