# --- author: jose c. garcia alanis
# --- encoding: utf-8
# --- r version: 3.5.1 (2018-07-02) -- "Feather Spray"
# --- script version: Dez 2018
# --- content: analysis of error rates

# --- 1) set paths and get workflow functions ----------------------------------
# path to project
setwd('/Volumes/TOSHIBA/manuscrips_and_data/soc_ern/')

# workflow functions
source('./r_functions/getPacks.R')
source('./r_functions/stdResid.R')
source('./r_functions/overDisp.R')
source('./r_functions/dataSummary.R')

# load multiple packages necessary for analysis
getPacks(c('dplyr', 'ggplot2', 'viridis'))


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

# --- 3) Create plots for descriptive data -------------------------------------
# distribution of number of errors
ggplot(errors, aes(e_rate, fill = flankers)) +
  
  geom_histogram(color = 'black', bins = 12) + 
  facet_wrap(~ flankers, scales = 'free') +
  
  scale_x_continuous(limits = c(0, .5)) +
  scale_y_continuous(limits = c(0, 33)) +
  
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
                                    face = 'bold', hjust = .2),
        axis.title.y = element_text(size = 14, 
                                    color = 'black', 
                                    face = 'bold', hjust = .15),
        axis.text = element_text(size = 11, 
                                 color = 'black'),
        legend.position = 'none') + 
  geom_rug()


# violin plot for manuscript
err_box <-  ggplot(errors, 
                   aes(x = flankers, y = e_rate, color = flankers)) +
  
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
       filename = './paper_figs/Fig_2a.pdf', 
       device = 'pdf',  width = 4, height = 5)


# --- 4) Compute descriptive statistics  --------------------------------------
# errors overall
errors %>% dplyr::summarise(M_errors = mean(t_errors), 
                            SD_errors  = sd(t_errors), 
                            
                            M_error_rate = mean(overall_e_rate), 
                            SD_error_rate  = sd(overall_e_rate)
                            )

# errors by group
errors %>% dplyr::group_by(group) %>% 
  dplyr::summarise(M = mean(t_errors), 
            SD = sd(t_errors), 
            SE = sd(t_errors) / sqrt(sum(!is.na(t_errors))))

# errors by flankers
as.data.frame(errors %>% dplyr::group_by(flankers) %>% 
                dplyr::summarise(M = mean(e_rate), 
                                 SD = sd(e_rate), 
                                 SE = sd(e_rate) / sqrt(sum(!is.na(e_rate))), 
                                 Min = min(e_rate), Max = max(e_rate)))

# --- 5) Set up and fit the "full model" ---------------------------------------

# # For completness, we document the results of a "full model", which 
# # included all higher order interactions. However, interpreting 
# # the estimates of the full model might be complicated. The more 
# # parsimonious models described in sections 7) and 7) show the same
# # results and their interpretation is straight foward.
# # Uncomment and run this section to fit the full model.

# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum","contr.poly"))

# (UNCOMENT TO RUN):
# DUMMY CODE cathegorical predictors
contrasts(errors$group) <- contr.sum(2); contrasts(errors$group)
contrasts(errors$flankers) <- contr.sum(4); contrasts(errors$flankers)

# ** fit full model including all variables and interactions **
mod_err_full <-  lmer(data = errors,
                      e_rate ~ d_motivation + group*flankers*agency +
                        group*flankers*affiliation + (1|id), 
                      REML = F)
# anova table
anova(mod_err_full)
# residuals ok?
qqPlot(resid(mod_err_full))

# coefficients of determination
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(update(mod_err_full))


# --- 6) Set up and fit reported models  ---------------------------------------
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn'))

# -- change default contrasts options ! --
options(contrasts = c("contr.sum","contr.poly"))

# fffect code cathegorical predictors 
contrasts(errors$group) <- contr.sum(2); contrasts(errors$group)
contrasts(errors$flankers) <- contr.sum(4); contrasts(errors$flankers)

# ***** ****** ***** ****** ***** ****** ***** ****** ***** ****** ***** ******
# fit model with interaction, controlling for motivation
# mod_err_inter <- glmer(data = errors,
#                        e_rate ~ d_motivation + flankers*group + (1|id),
#                        family = poisson(link = 'log'),
#                        control = glmerControl(optimizer='bobyqa'), 
#                        nAGQ = 20)

# fit model controllig for motivation
mod_err_inter <- lmer(data = errors,
                      e_rate ~ d_motivation + 
                        flankers*group + flankers*group + (1|id),
                      REML = F)
# nova table
anova(mod_err_inter)

# Compute resiudals and detect outliers
er_rm <- stdResid(data = errors, 
                  model = mod_err_inter, 
                  plot = T, show.bound = T)

# Re-fit without outliers
mod_err_inter_1 <- lmer(data = filter(er_rm, Outlier == 0),
                        e_rate ~ d_motivation + flankers*group + flankers*group + (1|id),
                        REML = F)

# Anova table and model summary
anova(mod_err_inter_1)
# Residuals ok?
qqPlot(resid(mod_err_inter_1))
plot_model(mod_err_inter_1, 'diag')

# ***** ****** ***** ****** ***** ****** ***** ****** ***** ****** ***** ******
# fit model without interaction
mod_errors <- lmer(data = errors, 
                    e_rate ~ flankers + group + (1|id), 
                   REML = F)
# anova table and model summary
anova(mod_errors)
summary(mod_errors)
qqPlot(resid(mod_errors))

# Compute resiudals and detect outliers
e_rm <- stdResid(data = errors, model = mod_errors, 
                 plot = T, show.bound = T)

# re-fit without outliers
mod_errors_1 <- lmer(data = filter(e_rm, Outlier == 0), 
                   e_rate ~ flankers + group + (1|id), 
                   REML = F)
# anova table and model summary
anova(mod_errors_1)
summary(mod_errors_1)
# residuals ok?
qqPlot(resid(mod_errors_1))
plot_model(mod_errors_1, 'diag')

# semi-partial r2
((3/149.161)*78.8424)/(1+((3/149.161)*78.8424))

# coefficents of detemination
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_err_inter_1)
r.squaredGLMM(mod_errors_1)

# compare models with and without interaction
anova(mod_err_inter, mod_errors) # Interaction doesn't improve the model

# build table
# create summary table for the fitted models
tab_model(file = './revision/rev_tables/mod_errors.html',
          mod_err_inter_1, mod_errors_1, digits.p = 3, show.aic = T, show.std = T, 
          title = 'Table 1: Results of linear mixed effects regression analysis of error rates.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = c('Error rate', 'Error rate'),
          pred.labels = c('(Intercept)',
                          'delta-Engagement',
                          'Compatible',
                          'Identical',
                          'Incompatible',
                          'Competition',
                          'Compatible x Competition',
                          'Identical x Competition',
                          'Incompatible x Competition'))

getPacks(c('emmeans'))
# quick Plot
err_eff <- emmip(mod_errors_1, ~ flankers, CIs = T, 
                 lmer.df = 'satterthwaite')
err_eff + coord_cartesian(ylim = c(0, .2))

# --- 7) Pairwise contrasts for error model ---------------------
# load Packages
getPacks(c('emmeans'))

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

# --- 8) Plot log estimates for trial type ------------------------------------
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
  
  labs(y = expression(bold('Estimated error rate')),
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
        axis.text.x = element_text(size = 13, 
                                   color = 'black', angle = 90, vjust = .5, hjust = 1),
        axis.text.y = element_text(size = 13, 
                                   color = 'black'),
        legend.position = 'none'); err_p

# save plot
ggsave(err_p, 
       filename = './paper_figs/Fig_2b.pdf', 
       device = 'pdf',  width = 4, height = 5)


# --- 9) FOLLOW-UP analyses - Incompatible trials ------------------------------
# load Packages
getPacks(c('lme4', 'lmerTest', 'sjPlot', 'car', 'MuMIn'))

# -- Change default contrasts options ! --
options(contrasts = c("contr.sum","contr.poly"))


# subset of data - only incompatible triasl
errors_in <- filter(errors, flankers == 'Incompatible')

# effect code predictors
contrasts(errors_in$group) <- contr.sum(2); contrasts(errors_in$group)

# # -- DON'T RUN --
# #  ** use dummy coding to look at simple slopes **
# contrasts(errors_in$group) <- contr.sum(2); contrasts(errors_in$group)

# fit a full model controlling for Motivation
mod_full_err_in <- lm(data = errors_in, e_rate ~
                          d_motivation + group*affiliation + group*agency)
car::Anova(mod_full_err_in, type = 'III')

# detect outlying observations
Ini_rm <- stdResid(data = errors_in, mod_full_err_in, 
                   return.data = T, plot = T, show.bound = T)

# re-fit without outliers
mod_full_err_in_1 <- lm(data = filter(Ini_rm, Outlier == 0), e_rate ~
                             d_motivation + group*affiliation + group*agency)
car::Anova(mod_full_err_in_1, type = 'III')
summary(mod_full_err_in_1)


# fit the model with-out Motivation
mod_err_in <- lm(data = errors_in, e_rate ~
                         group*affiliation + group*agency)
car::Anova(mod_err_in, type = 'III')
summary(mod_err_in)
qqPlot(resid(mod_err_in))

# detect outlying observations
In_rm <- stdResid(data = errors_in, mod_err_in, 
                  return.data = T, plot = T, 
                  show.loess = T, show.bound = T)

# re-fit without outliers
mod_err_in_1 <- lm(data = filter(In_rm, Outlier == 0), 
                   e_rate ~ group*affiliation + group*agency)
car::Anova(mod_err_in_1, type='III')
summary(mod_err_in_1)
# residulas ok?
qqPlot(resid(mod_err_in_1))
plot_model(mod_err_in_1, 'diag')

# semi-partial r2
sjstats::eta_sq(Anova(mod_err_in_1, type = 'III'), partial = F)

# coefficent of detemination
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(update(mod_full_err_in_1, nAGQ = 1)) # fit model by Laplace approximation
r.squaredGLMM(update(mod_err_in_1, nAGQ = 1)) # fit model by Laplace approximation

# compare models with and without interaction
anova(mod_full_err_in, mod_err_in) # Interactions doesn't improve the model

# table for model summary
tab_model(file = './revision/rev_tables/mod_incomp_errors.html',
          mod_full_err_in_1, mod_err_in_1, digits.p = 3, show.aic = T, show.std = T, 
          title = 'Table 1: Results of linear regression analysis of error rates on incompatible trials.',
          CSS = list(css.thead = 'padding:0.1cm;'),
          dv.labels = c('Error rate', 'Error rate'),
          pred.labels = c('(Intercept)',
                          'delta-Engagement',
                          'Competition',
                          'affiliation',
                          'agency',
                          'Competition x affiliation',
                          'Competition x agency'))

# --- 10) Interaction analysis -------------------------------------------------
# quick interaction plot
emmip(mod_err_in_1, affiliation ~ group, 
      at = list(affiliation = c(-4, 4), 
                agency = 0, d_motivation = 0), CIs = T)

# save simple slopes of Affiliation
emm_trend_SC <- emtrends(mod_err_in_1, 
                         var = 'affiliation', 
                         ~ 1)

# effects of Affiliation
summary(emm_trend_SC, infer = T)

# save Simple slopes of Agency
emm_trend_MAE <- emtrends(mod_err_in_1, 
                          var = 'agency', ~ 1,
                          transform ='response')

# effect of Agency
summary(emm_trend_MAE, infer = T)


# --- 11) CREATE Follow-Up FIGURES --------------------------------
# affiliation means
affiliation_means <- emmeans::emmeans(mod_err_in_1, ~ affiliation | group,
                                      CIs = T, 
                                      lmer.df = 'satterthwaite',
                                      at = list(end = 0, 
                                                affiliation = c(-4, 4), 
                                                agency = 0))
# pairwise comparissons
pairs(affiliation_means, adjust = 'bonferroni')

# affilaition means
affiliation_means <- emmeans::emmeans(mod_err_in_1, ~ group | affiliation,
                                 CIs = T, 
                                 lmer.df = 'satterthwaite',
                                 at = list(end = 0, 
                                           affiliation = c(-4, 4), 
                                           agency = 0))
# pairwise comparissons
pairs(affiliation_means, adjust = 'bonferroni')

# plot interaction affiliation x context
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
  scale_x_discrete(breaks = c(-4, 4), labels = c('low', 'high')) +
  coord_cartesian(ylim = c(0, .25)) +
  scale_color_manual(name="Social context",
                     values=c("#DE4968FF", "#000004FF")) +
  
  labs(title = 'Simple slopes of Affiliation',
       y = 'Estimated error rate',
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
        axis.title.y = element_text(color = 'black', size = 14),
        plot.title = element_text(face = 'bold', hjust = .5),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.direction = 'horizontal'); plt_aff

# save plot
ggsave(plt_aff, 
       filename = './paper_figs/Fig_2c.pdf', 
       device = 'pdf',  width = 4, height = 5)


# --- 12) ANALYSE RT - Incompatible trials -------------------------------------
Errors_In <- as.data.frame(Errors_In, row.names = 1:76)

# check distribution
hist(Errors_In$M_RT, breaks = 10)
rug(Errors_In$M_RT)

# fit Model 
mod_err_RT <- lm(data = Errors_In, 
                 M_RT ~ Interest + Group + Agency + Affiliation)
anova(mod_err_RT)
summary(mod_err_RT)
plot(mod_err_RT)

# remove Outliers? --> No effects
mod_err_RT <- lm(data = Errors_In[-c(18, 37), ],
                 M_RT ~ Interest + Group + Agency + Affiliation)
anova(mod_err_RT)
summary(mod_err_RT)
plot(mod_err_RT)

# compute CIs
confint(mod_err_RT)




