# --- author: jose c. garcia alanis
# --- encoding: utf-8
# --- r version: 4.0.1 (2020-06-06) -- "See Things Now"
# --- script version: Jun 2020
# --- content: analysis of error rates

# --- 1) load workflow functions and necessary packages ------------------------
# workflow functions
source('./r_functions/getPacks.R')
source('./r_functions/stdResid.R')
source('./r_functions/spR2.R')
source('./r_functions/dataSummary.R')

# load packages
getPacks(c('dplyr', 'ggplot2', 'viridis', 'sjPlot'))


# --- 2) Import data -----------------------------------------------------------
errors <- read.table('../data/behavioral/behavioral_data.txt',
                    header = T)

errors <- na.omit(errors)


# --- 3) Descriptives ----------------------------------------------------------
# Plot distribution of number of errors
error_dist_plot <- ggplot(errors,
                          aes(x = error_rate_p_conditon,
                              fill = flankers)) +
  
  geom_histogram(color = 'black', bins = 20) +
  facet_wrap(~ flankers, scales = 'free') +

  labs(x = '\n Error rate', y = 'Incidence \n') + 
  
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 30), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = 0, y = -Inf, xend = 0.5, yend = -Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme_classic() + 
  theme(axis.line = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, 
                                  colour = 'black', 
                                  face = 'bold'),
        axis.title.x = element_text(size = 14, 
                                    color = 'black',
                                    face = 'bold', hjust = 0.5),
        axis.title.y = element_text(size = 14, 
                                    color = 'black', 
                                    face = 'bold', hjust = 0.5),
        axis.text = element_text(size = 11, 
                                 color = 'black'),
        legend.position = 'none') + 
  geom_rug(); error_dist_plot
ggsave(error_dist_plot,
       filename = './results/figures/distribution_of_errors.pdf',
       device = 'pdf',  width = 10, height = 8)

# violin plot for manuscript
err_box <-  ggplot(errors,
                   aes(x = flankers,
                       y = error_rate_p_conditon,
                       color = flankers)) +
  
  geom_violin(fill = NA, 
              color='black', trim = T) + 
  
  geom_jitter(width = 0.35, 
              alpha = 0.7,
              size = 0.7) +
  
  stat_summary(fun.data = dataSummary, color = 'black', shape = 23, 
               fill='black', size = .3) +
  
  labs(x = 'Levels of flankers', 
       y = 'Error rate') + 
  theme_classic() +
  
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 0.5), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = 'Compatible', y = -Inf, xend = 'Neutral', yend = -Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  scale_color_viridis(option = 'B', discrete = T, end = .9) + 
  
  theme(strip.background = element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_text(size = 14, color = 'black', face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, color = 'black', face = 'bold', 
                                    margin = margin(r = 15)),
        axis.text.x = element_text(size = 13, color = 'black', 
                                   angle = 45, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 13, color = 'black'),
        legend.position = 'none'); err_box 

# save plot
ggsave(err_box, 
       filename = './results/figures/Fig_2a.pdf',
       device = 'pdf',  width = 4, height = 5)


# Compute descriptive statistics
# means and sd for overall number of errors
errors %>% dplyr::summarise(M_errors = mean(total_nr_of_errors),
                            SD_errors  = sd(total_nr_of_errors),
                            
                            M_error_rate = mean(error_rate_overall),
                            SD_error_rate  = sd(error_rate_overall)
                            )  %>%
  tab_df(title = 'Means and sd for errors (obverall)',
         file = './results/tables/errors_overall.html')

# means and sd for errors by group
errors %>% dplyr::group_by(group) %>%
  dplyr::summarise(M = mean(total_nr_of_errors),
            SD = sd(total_nr_of_errors),
            SE = sd(total_nr_of_errors) / sqrt(sum(!is.na(total_nr_of_errors)))) %>%
  tab_df(title = 'Means and sd for errors by group',
         file = './results/tables/errors_by_group.html')

# means and sd for error rates by flankers
errors %>% dplyr::group_by(flankers) %>%
  dplyr::summarise(M = mean(error_rate_p_conditon),
                   SD = sd(error_rate_p_conditon),
                   SE = sd(error_rate_p_conditon) / sqrt(sum(!is.na(error_rate_p_conditon))),
                   Min = min(error_rate_p_conditon),
                   Max = max(error_rate_p_conditon)) %>%
  tab_df(title = 'Means and sd for errors by flankers',
         file = './results/tables/errors_by_flankers.html')


# --- 4) Statistical analysis --------------------------------------------------
# *** Set up and fit the "full model" ***

# # For completness, we document the results of a "full model", which 
# # included all higher order interactions.
# # However, interpreting the estimates of the full model might be complicated.
# # The more parsimonious models described in the sections below show the same
# # results and their interpretation is straight foward.

# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum", "contr.poly"))

# DUMMY CODE cathegorical predictors
errors$group <- as.factor(errors$group)
contrasts(errors$group) <- contr.sum(2); contrasts(errors$group)

errors$flankers <- as.factor(errors$flankers)
contrasts(errors$flankers) <- contr.sum(4); contrasts(errors$flankers)

# ** fit full model including all variables and interactions **
mod_err_full <-  lmer(data = errors,
                      error_rate_p_conditon ~
                        diff_in_motivation_centred + group*flankers*MAE_centred +
                        group*flankers*SC_centred + (1|id),
                      REML = F)
# Anova table
anova_full <- anova(mod_err_full)
tab_df(cbind(row.names(anova_full), as.data.frame(anova_full)), digits = 5,
       title = 'ANOVA table for "full-model" of error rates (not reported)',
       file = './results/tables/ANOVA_full_mod_errors.html')

# residuals ok?
qqPlot(resid(mod_err_full))

# coefficients of determination
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(update(mod_err_full))


# *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
# *** Set up and fit reported models ***
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum", "contr.poly"))

# fffect code cathegorical predictors 
contrasts(errors$group) <- contr.sum(2); contrasts(errors$group)
contrasts(errors$flankers) <- contr.sum(4); contrasts(errors$flankers)

# Test social context by flakers interaction
# fit model controllig for motivation
mod_err_inter <- lmer(data = errors,
                      error_rate_p_conditon ~
                        diff_in_motivation_centred + flankers*group + (1|id),
                      REML = F)
# Anova table
anova_inter <- anova(mod_err_inter)
tab_df(cbind(row.names(anova_inter), as.data.frame(anova_inter)), digits = 5,
       title = 'ANOVA table for "interaction-model" of error rates',
       file = './results/tables/ANOVA_intreaction_mod_errors.html')

# Compute resiudals and detect outliers
er_rm <- stdResid(data = errors, 
                  model = mod_err_inter, 
                  plot = T, show.bound = T)

# In order to test the effect of outliers, re-fit the model without outliers
# (model does not change)
mod_err_inter_1 <- lmer(data = filter(er_rm, Outlier == 0),
                        error_rate_p_conditon ~
                        diff_in_motivation_centred + flankers*group + (1|id),
                        REML = F)

# Anova table and model summary
anova(mod_err_inter_1)
# Residuals ok? They look ok
qqPlot(resid(mod_err_inter_1))
# further model diagnostics (not too bad)
plot_model(mod_err_inter_1, 'diag')


# *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
# Since interaction and motivation were not significant:
# dro those terms, re-fit the model and see
# if effects change.
mod_errors <- lmer(data = errors,
                   error_rate_p_conditon ~ flankers + group + (1|id),
                   REML = F)
# anova table and model summary
anova(mod_errors)

# Compute resiudals and detect outliers
e_rm <- stdResid(data = errors, model = mod_errors, 
                 plot = T, show.bound = T)

# In order to test the effect of outliers, re-fit the model without outliers
# (model does not change)
mod_errors_1 <- lmer(data = filter(e_rm, Outlier == 0),
                     error_rate_p_conditon ~ flankers + group + (1|id),
                     REML = F)
# anova table and model summary
anova(mod_errors_1)
summary(mod_errors_1)
# residuals ok? They look ok
qqPlot(resid(mod_errors_1))

# compute effect sizes (semi partial R2) from anova table
amod <- anova(mod_errors_1); amod
amod <-  as.data.frame(amod); amod
amod$sp.R2 <- spR2(amod); amod

# coefficents of detemination
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_err_inter_1)
r.squaredGLMM(mod_errors_1)

# compare models with and without interaction
anova(mod_err_inter, mod_errors) # Interaction doesn't improve the model

# build table
# create summary table for the fitted models
tab_model(file = './results/tables/TABLE_S2_final_model_errors.html',
          mod_err_inter_1, mod_errors_1, digits = 2,
          show.aic = T, 
          collapse.ci = T,
          title = 'Table S2: Results of linear mixed effects regression analysis of error rates.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = c('Error rate', 'Error rate'),
          pred.labels = c('(Intercept)',
                          'Engagement',
                          'Compatible (C)',
                          'Identical (I)',
                          'Incompatible (In)',
                          'Competition',
                          'C x Competition',
                          'I x Competition',
                          'In x Competition'))


# --- 5) Pairwise contrasts for error model ------------------------------------
# load Packages
getPacks('emmeans')

# quick Plot
err_eff <- emmip(mod_errors_1, ~ flankers, CIs = T, 
                 lmer.df = 'satterthwaite')
err_eff + coord_cartesian(ylim = c(0, .2))

# group estimates
group_err <- emmeans(mod_errors_1, pairwise ~ group,
                   adjust = 'bonferroni', lmer.df = 'satterthwaite')
group_err
confint(group_err)

# save trial type estimates
est_err <- emmeans(mod_errors_1, pairwise ~ flankers,
                   adjust = 'bonferroni')
# effect of trial type
est_err
as.data.frame(est_err$contrasts)
# and CIs
confint(est_err)

# save group estimates
est_group <- emmeans(mod_errors_1, pairwise ~ group,
                     adjust = 'fdr')
# effect of group
est_group
as.data.frame(est_group$contrasts)
# CIs
confint(est_group)

# *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
# Plot model estimates for trial type
# plot estimates
err_p <- ggplot(data = as.data.frame(emmeans(mod_errors, ~ flankers)), 
                aes(y = emmean, x = flankers)) + 
  
  geom_hline(yintercept = mean(as.data.frame(emmeans(mod_errors, ~ flankers))[, 2]), 
             linetype=2) +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = flankers), 
                width = 0.2, size = 1, alpha = 0.7) + 
  geom_linerange(aes(ymin = emmean-SE, ymax = emmean+SE), color = 'gray',
                 size = 3) + 
  geom_point(shape = 18, size = 3.5, color = 'black') +
  
  labs(y = expression('Estimated error rate'),
       x = 'Trial type') +
  
  theme_classic() + 
  
  scale_color_viridis(option = 'B', discrete = T, end = .90) +
  
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = .25), 
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
        axis.text.x = element_text(size = 13, color = 'black', 
                                   angle = 45, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 13, 
                                   color = 'black'),
        legend.position = 'none'); err_p

# save plot
ggsave(err_p, 
       filename = './results/figures/Fig_2b.pdf',
       device = 'pdf',  width = 4, height = 5)


# --- 6) Follow-up analyses - Incompatible trials ------------------------------
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn'))

# -- Change default contrasts options ! --
options(contrasts = c("contr.sum", "contr.poly"))

# subset of data - only incompatible triasl
errors_in <- filter(errors, flankers == 'Incompatible')

# effect code predictors
contrasts(errors_in$group) <- contr.sum(2); contrasts(errors_in$group)

# # -- DON'T RUN --
# #  ** only use to dummy coding group factor to look at simple slopes **
# contrasts(errors_in$group) <- contr.treatment(2, base = 1); contrasts(errors_in$group)

# fit a full model controlling for motivation
mod_full_err_in <- lm(data = errors_in,
                      error_rate_p_conditon ~ diff_in_motivation_centred +
                        group*SC_centred + group*MAE_centred)
car::Anova(mod_full_err_in, type = 'III')

# Compute resiudals and detect outliers
Ini_rm <- stdResid(data = errors_in, mod_full_err_in, 
                   return.data = T, plot = T, show.bound = T)

# In order to test the effect of outliers, re-fit the model without outliers
mod_full_err_in_1 <- lm(data = filter(Ini_rm, Outlier == 0),
                        error_rate_p_conditon ~ diff_in_motivation_centred +
                          group*SC_centred + group*MAE_centred)
car::Anova(mod_full_err_in_1, type = 'III')
summary(mod_full_err_in_1)


# *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
# Since motivation has no effect, drop that variable and re-fit the model
mod_err_in <- lm(data = errors_in,
                 error_rate_p_conditon ~ group*SC_centred + group*MAE_centred)
car::Anova(mod_err_in, type = 'III')
summary(mod_err_in)

# Compute resiudals and detect outliers
In_rm <- stdResid(data = errors_in, mod_err_in, 
                  return.data = T, plot = T, 
                  show.loess = T, show.bound = T)

# In order to test the effect of outliers, re-fit the model without outliers
mod_err_in_1 <- lm(data = filter(In_rm, Outlier == 0), 
                   error_rate_p_conditon ~ group*SC_centred + group*MAE_centred)
car::Anova(mod_err_in_1, type='III')
summary(mod_err_in_1)
# residulas ok? The look ok
qqPlot(resid(mod_err_in_1))
tab_model(mod_err_in_1, show.std = T)

# semi-partial r2
sjstats::eta_sq(Anova(mod_err_in_1, type = 'III'), partial = F)

# coefficent of detemination
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_full_err_in_1) 
r.squaredGLMM(mod_err_in_1) 

# compare models with and without interaction
anova(mod_full_err_in, mod_err_in) # Interactions doesn't improve the model

# table for model summary
tab_model(file = './results/tables/TABLE_S2_final_model_incomp_errors.html',
          mod_full_err_in_1, mod_err_in_1, digits = 3,
          show.aic = T, 
          collapse.ci = T,
          show.std = T,
          title = 'Table S3: Results of linear regression analysis of error rates on incompatible trials.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = c('Error rate', 'Error rate'),
          pred.labels = c('(Intercept)',
                          'Engagement',
                          'Competition',
                          'Affiliation',
                          'Agency',
                          'Competition x affiliation',
                          'Competition x agency'))


# --- 7) Create interaction caontext by affiliation figure  --------------------
# affiliation high vs low means by group
affiliation_means <- emmeans::emmeans(mod_err_in_1, ~ SC_centred | group,
                                      CIs = T,
                                      lmer.df = 'satterthwaite',
                                      at = list(SC_centred = c(-4, 4),
                                                agency = 0))
# pairwise comparissons
pairs(affiliation_means, adjust = 'bonferroni')

# group means by affiliation high vs low
affiliation_means <- emmeans::emmeans(mod_err_in_1, ~ group | SC_centred,
                                 CIs = T, 
                                 lmer.df = 'satterthwaite',
                                 at = list(SC_centred = c(-4, 4),
                                           agency = 0))
# pairwise comparissons
pairs(affiliation_means, adjust = 'bonferroni')

# plot interaction affiliation x context
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
  scale_x_discrete(breaks = c(-4, 4), labels = c('low', 'high')) +
  coord_cartesian(ylim = c(0, .25)) +
  scale_color_manual(name="Social context",
                     values=c("#DE4968FF", "#000004FF")) +
  
  labs(title = 'Simple slopes of Affiliation',
       y = 'Estimated error rate (imcomp. trials)',
       x = 'Levels of Affiliation') +
  
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = .25), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = '-4', y = -Inf, xend = '4', yend = -Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.text = element_text(color = 'black', size = 13), 
        axis.title.x = element_text(color = 'black', face = 'bold', size = 14, 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(color = 'black', size = 14, face = 'bold',
                                    margin = margin(r = 15)),
        plot.title = element_text(face = 'bold', hjust = .5),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.direction = 'horizontal'); plt_aff

# save plot
ggsave(plt_aff, 
       filename = './results/figures/Fig_2c.pdf',
       device = 'pdf',  width = 4, height = 5)
