# --- author: jose c. garcia alanis
# --- encoding: utf-8
# --- r version: 3.5.1 (2018-07-02) -- "Feather Spray"
# --- script version: Dez 2018
# --- content: analysis of behavioural correct reactions

# --- 1) Set paths and get workflow functions ----------------------------------
# path to project
setwd('/Volumes/TOSHIBA/manuscrips_and_data/soc_ern/')

# workflow functions
source('./r_functions/getPacks.R')
source('./r_functions/stdResid.R')
source('./r_functions/overDisp.R')
source('./r_functions/dataSummary.R')

# Install and load multiple R packages necessary for analysis.
getPacks(c('tidyr', 'dplyr',
           'lme4', 'lmerTest', 'car', 'MuMIn', 'emmeans', 
           'sjPlot', 'ggplot2', 'viridis'))


# --- 2) Import data -----------------------------------------------------------
# personality data
perso <- read.table('./data_for_r/all_perso.txt', 
                    header = T)
unique(perso$id)
# summarise (for later)
perso %>% summarise_if(is.numeric, sd)

# data frame containing behavioural errors
load('./data_for_r/errors_data.RData')
unique(errors_data$id)
# merge data frames
errors <- merge(errors_data, perso, 'id'); rm(errors_data, perso)
errors <-  dplyr::filter(errors, flankers == 'Incompatible')
errors$flankers <-  NULL

# corrects data
data_corrects <- read.table('./data_for_r/corrects.txt', header = T)
names(data_corrects)[1] <-  'id'
names(data_corrects)[2] <- 'flankers'

# capitalize first letter in flankers variable
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
data_corrects$flankers <- firstup(as.character(data_corrects$flankers))
# to factor
data_corrects$flankers <- as.factor(data_corrects$flankers)

# summarise for later
as.data.frame(corrects %>% group_by(flankers) %>% 
  summarise(RT = mean(M_RT), sd = sd(M_RT)))

# merge wich errors
corrects <-  merge(data_corrects, errors, c('id'))

# --- 3) Plot deistribution  ---------------------------------------------------
corr_box <- ggplot(corrects, 
                   aes(x = flankers, y = M_RT, color = flankers)) +
  
  geom_violin(fill = NA, 
              color='black', trim = T) + 
  
  geom_jitter(width = 0.35, 
              alpha = 0.7,
              size = 0.7) +
  
  stat_summary(fun.data = dataSummary, color = 'black', shape = 23, 
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
        axis.text.x = element_text(size = 13, color = 'black', 
                                   angle = 45, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 13, 
                                   color = 'black'),
        legend.position = 'none'); corr_box

# save plot
ggsave(corr_box, 
       filename = './paper_figs//Fig_3a.pdf', 
       device = 'pdf',  width = 4, height = 5)



# --- 4) effect code predictors ------------------------------------------------
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum","contr.poly"))

# effect code
contrasts(corrects$group) <- contr.sum(2); contrasts(corrects$group)
contrasts(corrects$flankers) <- contr.sum(4); contrasts(corrects$flankers)


# --- 5) full model  -----------------------------------------------------------
# fit full model
mod_full <- lmer(data = corrects, 
                 M_RT ~ d_motivation + total_error_rate + 
                  group*flankers*affiliation + 
                  group*flankers*agency + (1|id), 
                 REML = F)
# anova table
anova(mod_full)

# identify outlying observations
c_rm <- stdResid(data = corrects, model = mod_full, plot = T, 
                 main = expression('Residuals ' ['LMER model for RT data']), 
                 xlab = expression('Fitted Values ' ['Mean RT']),
                 ylab = 'Std. Pearson Residuals', show.bound = T)

# re-fit without outliers
mod_full_1 <- lmer(data = filter(c_rm, Outlier == 0), 
                   M_RT ~ d_motivation + total_error_rate + 
                     group*flankers*affiliation + 
                     group*flankers*agency + (1|id), 
                   REML = F)
# anova table
anova(mod_full_1)
tab_model(mod_full_1)


# --- 6) set up and fit the reported models ------------------------------------
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum","contr.poly"))

# effect code
contrasts(corrects$group) <- contr.sum(2); contrasts(corrects$group)
contrasts(corrects$flankers) <- contr.sum(4); contrasts(corrects$flankers)

# more parsimonious model including personality variables
mod_corrects <- lmer(data = corrects, 
                     M_RT ~ total_error_rate + d_motivation + 
                       group + flankers + affiliation + agency + 
                       (1|id), 
                     REML = F)
# anova table
anova(mod_corrects)
# residuals ok?
qqPlot(resid(mod_corrects, 'pearson'))

# identify outliers
c_rm <- stdResid(data = corrects, model = mod_corrects, plot = T, 
                 main = expression('Residuals ' ['LMER model for RT data']), 
                 xlab = expression('Fitted Values ' ['Mean RT']),
                 ylab = 'Std. Pearson Residuals', show.bound = T)

# re-fit without outliers
mod_corrects_1 <- lmer(data = filter(c_rm, Outlier == 0), 
                       M_RT ~ total_error_rate + d_motivation + 
                         group + flankers + affiliation + agency + 
                         (1|id), 
                       REML = F)
# anova table and model summary
anova(mod_corrects_1)
summary(mod_corrects_1)
# residuals ok?
qqPlot(resid(mod_corrects_1))
# model table
tab_model(mod_corrects_1, show.std = T)

# Semi-partial r2
((1/75.763)*0.5828)/(1+((1/75.763)*0.5828)) # context
((1/75.688)*0.9255)/(1+((1/75.688)*0.9255))
((1/75.704)*2.7196)/(1+((1/75.704)*2.7196))
((1/223.824)*190.0878)/(1+((1/223.824)*190.0878)) # flankers
((1/75.983)*44.7098)/(1+((1/75.983)*44.7098)) # error rate
((1/75.756)*5.8851)/(1+((1/75.756)*5.8851)) # d motivation

# table for model summary
tab_model(file = './revision/rev_tables/mod_corrects.html',
          mod_corrects_1, digits = 2,
          show.aic = T, 
          collapse.ci = T,
          title = 'Table S4: Results of linear mixed-effects regression analysis of correct reactions RT.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = c('Correct reactions RT'),
          pred.labels = c('(Intercept)',
                          'Overall error rate',
                          'delta-Engagement',
                          'Competition',
                          'Compatible',
                          'Indetical',
                          'Incompatible',
                          'Affiliation',
                          'Agency'))

# descriptives
as.data.frame(corrects %>% dplyr::group_by(flankers) %>% 
                dplyr::summarise(M = mean(M_RT), 
                                 SD = sd(M_RT), 
                                 SE = sd(M_RT) / sqrt(sum(!is.na(M_RT))), 
                                 Min = min(M_RT), Max = max(M_RT)))

# --- 6) Follow-up analyises ---------------------------------------------------
group_corr <- emmeans(mod_corrects_1, pairwise ~ flankers,
                     adjust = 'bonferroni', lmer.df = 'satterthwaite')
as.data.frame(group_corr$contrasts)
confint(group_corr)

# effects of overall error-rates
emtrends(mod_corrects_1, var = 'total_error_rate', 
         lmer.df = 'satterthwaite', ~ 1,
         at = list(affiliation = 0, 
                   agency = 0, d_motivation = 0))
# effects of overall delta-motivation
emtrends(mod_corrects_1, var = 'd_motivation', 
         lmer.df = 'satterthwaite', ~ 1,
         at = list(affiliation = 0, 
                   agency = 0, total_error_rate = 0))


# --- 7) Plot model estimates for trial type ------------------------------------
# plot estimates
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
       filename = './paper_figs/Fig_3b.pdf', 
       device = 'pdf',  width = 4, height = 5)


# --- 7) create simple slope plots ---------------------------------------------
# desciptives
errors %>% summarise(sd_err = sd(overall_e_rate), sd_motiv = sd(d_motivation))

err_slope <- emmeans(mod_corrects_1, pairwise ~ total_error_rate,
                       at = list(total_error_rate = c(-0.05, 0.05), 
                                 d_motivation = 0,
                                 agency = 0,
                                 affiliation = 0),
                       adjust = 'bonferroni', lmer.df = 'satterthwaite'); err_slope


# plot slope of error rate
pd = position_dodge(.2)
plt_err <- ggplot(data = data.frame(err_slope$emmeans), 
                  aes(x = as.factor(total_error_rate), 
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
       filename = './paper_figs/Fig_3c.pdf', 
       device = 'pdf',  width = 4, height = 5)



motiv_slope <- emmeans(mod_corrects_1, pairwise ~ d_motivation,
                       at = list(total_error_rate = 0, 
                                 d_motivation = c(-2, 2),
                                 agency = 0,
                                 affiliation = 0),
                       adjust = 'bonferroni', lmer.df = 'satterthwaite'); motiv_slope

# plot slope of motivation
pd = position_dodge(.2)
plt_motiv <- ggplot(data = data.frame(motiv_slope$emmeans), 
                  aes(x = as.factor(d_motivation), 
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
       filename = './paper_figs/Fig_3d.pdf', 
       device = 'pdf',  width = 4, height = 5)
