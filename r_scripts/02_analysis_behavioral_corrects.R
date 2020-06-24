# --- author: jose c. garcia alanis
# --- encoding: utf-8
# --- r version: 3.5.1 (2018-07-02) -- "Feather Spray"
# --- script version: Dez 2018
# --- content: analysis of behavioural correct reactions

# --- 1) load workflow functions and necessary packages ------------------------
# workflow functions
source('./r_functions/getPacks.R')
source('./r_functions/stdResid.R')
source('./r_functions/spR2.R')
source('./r_functions/dataSummary.R')

# load packages
getPacks(c('tidyr', 'dplyr', 'ggplot2', 'viridis'))

# --- 2) Import data -----------------------------------------------------------
corrects <- read.table('../data/behavioral/behavioral_data.txt',
                       header = T)

################################################################################
# 3) Descriptives --------------------------------------------------------------
# Plot dirtribution of correct reaction time
corr_box <- ggplot(corrects,
                   aes(x = flankers, y = mean_correct_rt, color = flankers)) +
  
  geom_violin(fill = NA, 
              color='black', trim = T) + 
  
  geom_jitter(width = 0.35, 
              alpha = 0.7,
              size = 0.7) +
  
  stat_summary(fun.data = dataSummary, color = 'black', shape = 23, 
               fill='black', size = .3) +
  
  labs(x = 'Trial type', 
       y = 'Mean RT (ms)') +
  
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
        axis.text.x = element_text(size = 13, color = 'black', 
                                   angle = 45, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 13, 
                                   color = 'black'),
        legend.position = 'none'); corr_box

# save plot
ggsave(corr_box, 
       filename = './results/figures/Fig_3a.pdf',
       device = 'pdf',  width = 4, height = 5)

################################################################################
# 4) Statistical analysis ------------------------------------------------------
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum", "contr.poly"))

# effect code
corrects$group <- as.factor(corrects$group)
contrasts(corrects$group) <- contr.sum(2); contrasts(corrects$group)
corrects$flankers <- as.factor(corrects$flankers)
contrasts(corrects$flankers) <- contr.sum(4); contrasts(corrects$flankers)

# *** Set up and fit the "full model" ***
# For completness. However, interpreting the estimates of the
# full model might be complicated. The more parsimonious models described
# in the sections below show the same results and their interpretation
# is straight foward.
mod_full <- lmer(data = corrects,
                 mean_correct_rt ~
                   diff_in_motivation_centred + error_rate_overall_centred +
                   group*flankers*SC_centred +
                   group*flankers*MAE_centred + (1|id),
                 REML = F)
# anova table
anova(mod_full)

# Compute resiudals and detect outliers
c_rm <- stdResid(data = corrects, model = mod_full, plot = T, 
                 main = expression('Residuals ' ['LMER model for RT data']), 
                 xlab = expression('Fitted Values ' ['Mean RT']),
                 ylab = 'Std. Pearson Residuals', show.bound = T)

# In order to test the effect of outliers, re-fit the model without outliers
# (model does not change much)
mod_full_1 <- lmer(data = filter(c_rm, Outlier == 0), 
                 mean_correct_rt ~
                   diff_in_motivation_centred + error_rate_overall_centred +
                   group*flankers*SC_centred +
                   group*flankers*MAE_centred + (1|id),
                   REML = F)
# anova table
anova(mod_full_1)
# model summary
tab_model(mod_full_1,
          title = 'Full model of correct RT.')


# ------------------------------------------------------------------------------
# *** Set up and fit reported models ***
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn', 'emmeans'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum", "contr.poly"))

# effect code
corrects$group <- as.factor(corrects$group)
contrasts(corrects$group) <- contr.sum(2); contrasts(corrects$group)
corrects$flankers <- as.factor(corrects$flankers)
contrasts(corrects$flankers) <- contr.sum(4); contrasts(corrects$flankers)

# more parsimonious model including personality variables
mod_corrects <- lmer(data = corrects, 
                     mean_correct_rt ~
                       diff_in_motivation_centred + error_rate_overall_centred +
                       group + flankers + SC_centred + MAE_centred + (1|id),
                     REML = F)
# anova table
anova(mod_corrects)

# Compute resiudals and detect outliers
c_rm <- stdResid(data = corrects, model = mod_corrects, plot = T, 
                 main = expression('Residuals ' ['LMER model for RT data']), 
                 xlab = expression('Fitted Values ' ['Mean RT']),
                 ylab = 'Std. Pearson Residuals', show.bound = T)

# In order to test the effect of outliers, re-fit the model without outliers
# (model does not change)
mod_corrects_1 <- lmer(data = filter(c_rm, Outlier == 0),
                       mean_correct_rt ~
                         diff_in_motivation_centred + error_rate_overall_centred +
                         group + flankers + SC_centred + MAE_centred + (1|id),
                       REML = F)
# anova table and model summary
anova(mod_corrects_1)
summary(mod_corrects_1)
# residuals ok? They look ok
qqPlot(resid(mod_corrects_1))

# compute effect sizes (semi partial R2) from anova table
amod <- anova(mod_corrects_1); amod
amod <-  as.data.frame(amod); amod
amod$sp.R2 <- spR2(amod); amod

# table for model summary
tab_model(file = './results/tables/TABLE_S4_final_model_corrects.html',
          mod_corrects_1, digits = 2,
          show.std = T,
          show.aic = T, 
          collapse.ci = T,
          title = 'Table S4: Results of linear mixed-effects regression analysis of correct reactions RT.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = 'Correct reactions RT',
          pred.labels = c('(Intercept)',
                          'delta-Engagement',
                          'Overall error rate',
                          'Competition',
                          'Compatible',
                          'Indetical',
                          'Incompatible',
                          'Affiliation',
                          'Agency'))

# descriptives
as.data.frame(corrects %>% dplyr::group_by(flankers) %>% 
                dplyr::summarise(M = mean(mean_correct_rt),
                                 SD = sd(mean_correct_rt),
                                 SE = sd(mean_correct_rt) / sqrt(sum(!is.na(mean_correct_rt))),
                                 Min = min(mean_correct_rt), Max = max(mean_correct_rt)))


# ------------------------------------------------------------------------------
# follow analyses and pairwise contrasts for correct RT model
group_corr <- emmeans(mod_corrects_1, pairwise ~ flankers,
                     adjust = 'bonferroni', lmer.df = 'satterthwaite')
as.data.frame(group_corr$contrasts)
confint(group_corr)

# plot model estimates for trial type
corr_p <- ggplot(data = data.frame(group_corr$emmeans), 
                aes(y = emmean, x = flankers)) + 
  
  geom_hline(yintercept = mean(as.data.frame(emmeans(mod_corrects_1, ~ flankers))[, 2]), 
             linetype = 2) +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = flankers), 
                width = 0.2, size = 1, alpha = 0.7) + 
  geom_linerange(aes(ymin = emmean-SE, ymax = emmean+SE), color = 'gray',
                 size = 3) + 
  geom_point(shape = 18, size = 3.5, color = 'black') +
  
  labs(y = expression(bold('Estimated error rate')),
       x = 'Trial type') +
  
  theme_classic() + 
  
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
        axis.text.x = element_text(size = 13, color = 'black', 
                                   angle = 45, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 13, 
                                   color = 'black'),
        legend.position = 'none'); corr_p

# save plot
ggsave(corr_p, 
       filename = './results/figures/Fig_3b.pdf',
       device = 'pdf',  width = 4, height = 5)


# create simple slope plots
# desciptives
corrects %>% summarise(sd_err = sd(error_rate_overall_centred),
                     sd_motiv = sd(diff_in_motivation_centred))

# slope of error rate
err_slope <- emmeans(mod_corrects_1, pairwise ~ error_rate_overall_centred,
                       at = list(error_rate_overall_centred = c(-0.05, 0.05),
                                 diff_in_motivation_centred = 0,
                                 MAE_centred = 0,
                                 SC_centred = 0),
                       adjust = 'bonferroni', lmer.df = 'satterthwaite'); err_slope


# plot slope of error rate
pd <- position_dodge(.2)
plt_err <- ggplot(data = data.frame(err_slope$emmeans), 
                  aes(x = as.factor(error_rate_overall_centred),
                      y = emmean, group = 1)) +

  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = .1, 
                size = .8, 
                position = pd, color = 'black') +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2, position = pd) + 
  scale_x_discrete(breaks = c('-0.05', '0.05'), labels = c('low', 'high')) +
  coord_cartesian(ylim = c(250, 375)) +
  
  labs(title = 'Simple slope of error rate',
       y = 'Estimated RT (ms)',
       x = 'Levels of error rate') +
  
  geom_segment(aes(x = -Inf, y = 250, xend = -Inf, yend = 375), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = '-0.05', y = -Inf, xend = '0.05', yend = -Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(color = 'black', face = 'bold', size = 14),
        axis.line = element_blank(),
        axis.text = element_text(color = 'black', size = 13), 
        axis.title.x = element_text(color = 'black', face = 'bold', size = 14, 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(color = 'black', size = 14, face = 'bold',
                                    margin = margin(r = 15)),
        plot.title = element_text(face = 'bold', hjust = .5),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.direction = 'horizontal'); plt_err

# save plot
ggsave(plt_err, 
       filename = './results/figures/Fig_3c.pdf',
       device = 'pdf',  width = 4, height = 5)

# slope of motivation
motiv_slope <- emmeans(mod_corrects_1, pairwise ~ diff_in_motivation_centred,
                       at = list(error_rate_overall_centred = 0,
                                 diff_in_motivation_centred = c(-2, 2),
                                 MAE_centred = 0,
                                 SC_centred = 0),
                       adjust = 'bonferroni', lmer.df = 'satterthwaite'); motiv_slope

# plot slope of motivation
pd <- position_dodge(.2)
plt_motiv <- ggplot(data = data.frame(motiv_slope$emmeans), 
                  aes(x = as.factor(diff_in_motivation_centred),
                      y = emmean, group = 1)) +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = .1, 
                size = .8, 
                position = pd, color = 'black') +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2, position = pd) + 
  scale_x_discrete(breaks = c('-2', '2'), labels = c('low', 'high')) +
  coord_cartesian(ylim = c(250, 375)) +
  
  labs(title = expression(bold(paste("Simple slope of ", 
                                     Delta, "Engagement"))),
       y = 'Estimated RT (ms)',
       x = expression(bold(paste("Levels of ", 
                                 Delta, "Engagement")))) +
  
  geom_segment(aes(x = -Inf, y = 250, xend = -Inf, yend = 375), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = '-2', y = -Inf, xend = '2', yend = -Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(color = 'black', face = 'bold', size = 14),
        axis.line = element_blank(),
        axis.text = element_text(color = 'black', size = 13), 
        axis.title.x = element_text(color = 'black', face = 'bold', size = 14, 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(color = 'black', size = 14, face = 'bold',
                                    margin = margin(r = 15)),
        plot.title = element_text(face = 'bold', hjust = .5),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.direction = 'horizontal'); plt_motiv

# save plot
ggsave(plt_motiv, 
       filename = './results/figures/Fig_3d.pdf',
       device = 'pdf',  width = 4, height = 5)
