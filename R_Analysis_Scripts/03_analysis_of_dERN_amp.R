##### ##### #####     Analysis scrips for Alanis et al., 2018   ##### ##### #####
#                           LMER models for physiological
#                                     data

# Get helper functions
source('./R_Functions/getPacks.R')
source('./R_Functions/stdResid.R')

# Install and load multiple R packages necessary for analysis.
pkgs <- c('dplyr', 'reshape2',
          'lme4', 'lmerTest',
          'effects', 'emmeans', 'car', 'MuMIn',
          'ggplot2', 'viridis')

getPacks(pkgs)
rm(pkgs)



# ------ 1) Read in the data  -----------------------------------
load('~/Documents/Experiments/soc_ftask/data_for_r/Incomp_ERPs.RData')

# ------ 2) PLOT individual averages --------------------------------
ggplot(ERP[(ERP$Electrode == 'FCz') & ERP$Time > -500 , ], 
       aes(Time, Amplitude, color = Reaction, linetype = Group)) + 
  facet_wrap(~ Subject) +
  geom_vline(xintercept = c(0, 100), color = 'black', size = rel(0.5), linetype = 3) +
  geom_hline(yintercept = c(0), color = 'black', size = rel(0.5), linetype = 3) +
  geom_line(size = rel(0.7)) + 
  scale_color_viridis(option = 'D', discrete = T, 
                      begin = .15, end = .85, direction = -1) + 
  scale_y_reverse(breaks = c(-20, -10, 0, 10, 20)) +
  scale_linetype_manual(values = c(1,6)) +
  xlab("\n Time (ms)") + 
  ylab("Amplitude (µV) \n") +
  ggtitle("Average ERPs") +
  theme_classic() +
  geom_segment(aes(x = -Inf, y = -20, xend = -Inf, yend = 30), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = -500, y = Inf, xend = 1000, yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  theme(
    plot.title = element_text(color = "black", size = rel(1.5)),
    strip.background = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(color='black'),
    axis.text.y  = element_text(hjust = 1, color = 'black'),
    legend.text = element_text(color="black", size=rel(1)),
    legend.position = "bottom") 


# ------ 3) COMPUTE ∆ERN wave incorrect-correct ---------------------
# --- From long to wide
data_wide <- dcast(ERP, Group + Subject + Agency + Affiliation + Motivation + 
                     Total_Errors + Electrode + Time ~ Reaction,
                   value.var = c('Amplitude'))
# --- NAs created?
unique(is.na(data_wide))
# --- Compute ERN
data_wide <- data_wide %>% mutate(ERN = Incorrect - Correct)



# ------ 4) MODEL ∆ERN as function of midline electrodes --------
# --- AVERAGE for model
Ave_Elec <- plyr::ddply(filter(data_wide, (Time >= 0 & Time <= 100) & 
                                 (Electrode ==  'FCz' | Electrode == 'Cz' | Electrode == 'CPz')), 
                       c('Group', 'Electrode', 'Subject'), dplyr::summarise,
                       
                       N = sum(!is.na(ERN)),
                       M_Amp = mean(ERN, na.rm = T),
                       sd   = sd(ERN, na.rm = T),
                       se   = sd / sqrt(N),
                       # Compute t-staistic for confidence interval 
                       # (quantile t, with N-1 degrees of freedom)
                       # and use it to compute individual ci (multiplier for se):
                       ci   = se * qt(.95/2 + .5, N-1))

# --- DUMMY CODE electrode variable
Ave_Elec$Electrode <- factor(Ave_Elec$Electrode)
contrasts(Ave_Elec$Electrode) <- contr.sum(3); contrasts(Ave_Elec$Electrode)


# --- SET UP AND FIT the model
mod_elec <- lmer(data = Ave_Elec, 
                 M_Amp ~ Electrode + (1|Subject), 
                 REML = F)
anova(mod_elec)
summary(mod_elec)
qqPlot(resid(mod_elec, 'pearson'))

# --- Detect outlying observations
Elec_rm <- stdResid(data = Ave_Elec, mod_elec, 
                    plot = T, show.loess = T, 
                    show.bound = T)

# --- Re-fit without outliers
mod_elec_1 <- lmer(data = filter(Elec_rm, Outlier == 0), 
                   M_Amp ~ Electrode + (1|Subject), 
                   REML = F)
anova(mod_elec_1)
summary(mod_elec_1)
qqPlot(resid(mod_elec_1, 'pearson'))

# Coefficent of deteminations
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_elec_1)



# ------ 5) Follow-up analyses model electrodes --------------------
# --- Pairwise comparissons
emmeans(mod_elec_1, pairwise ~ Electrode, 
        adjust = 'bonferroni', lmer.df ='satterthwaite')

# --- Save electrode estimates
est_elec <- emmeans(mod_elec_1, pairwise ~ Electrode,
                   adjust = 'bonferroni', lmer.df ='satterthwaite')

# --- Effect of electrode and Plot results
est_elec;  emmip(mod_elec_1, ~ Electrode, CIs = T)
# --- Compute CIs
confint(est_elec)

# -----  6) PLOT estimates model electrodes ----------------------
# --- Create the plot
ggplot(data = as.data.frame(emmeans(mod_elec_1, ~ Electrode)), 
       aes(y = emmean, x = Electrode)) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, color = 'black', size = 0.8) + 
  geom_linerange(aes(ymin = emmean-SE, ymax = emmean+SE, color = Electrode), 
                 size = 3) + 
  geom_point(shape = 18, size = 3, color = 'black') +
  labs(y = expression(bold( 'estimated Amplitude ' ['µV'] ))) +
  theme(axis.title.y = element_text(color = 'black', size = 13,
                                    margin = margin(r=10)),
        axis.title.x = element_text(color = 'black', size = 13, face = 'bold',
                                    margin = margin(t=10)),
        axis.text = element_text(color = 'black', size = 11),
        legend.position = 'none') + 
  coord_cartesian(ylim = c(0, -6)) +
  scale_y_reverse()


# ------ 7) MODEL ∆ERN by context and personality ---------------
# ------ Only keep FCz and Cz
cz_fcz <- filter(data_wide, Electrode == 'Cz' | Electrode == 'FCz' )
cz_fcz <- filter(cz_fcz, Time >= 0 & Time <= 100 )

# --- DUMMY CODE variables
# Electrode
cz_fcz$Electrode <- factor(cz_fcz$Electrode)
contrasts(cz_fcz$Electrode) <- contr.sum(2); contrasts(cz_fcz$Electrode)
# Context
contrasts(cz_fcz$Group) <- contr.sum(2); contrasts(cz_fcz$Group)

# --- Mean center number of errors
cz_fcz <- within( cz_fcz, {
  Tot_Errors <- Total_Errors - mean(Total_Errors, na.rm = T)
})

# --- Descriptive statistics
cz_fcz %>% group_by(Group) %>% 
  dplyr::summarise(M = mean(ERN), 
            SD = sd(ERN))


# ***** ****** ******
# --- SET UP AND FIT full model 
# ---- controlling for ∆Motivation
# --- and number of errors
mod_ern_full <- lmer(data = cz_fcz, 
                     ERN ~ Tot_Errors + Motivation + 
                       Group*Affiliation + Group*Agency + 
                       (1|Subject), REML = F)
anova(mod_ern_full)

# --- Find outliers
cz_fcz_rm <- stdResid(data = cz_fcz, mod_ern_full, 
                      return.data = T, plot = T, 
                      show.loess = T, show.bound = T)

# --- Re-fit without outliers
mod_ern_full_1 <- lmer(data = filter(cz_fcz_rm, Outlier == 0), 
                       ERN ~ Tot_Errors + Motivation + 
                         Group*Affiliation + Group*Agency + 
                         (1|Subject), REML = F)
anova(mod_ern_full_1)

# Coefficent of deteminations
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_ern_full_1)



# ***** ****** ******
# --- SET UP AND FIT parsimonious model 
mod_ern <- lmer(data = cz_fcz, 
                ERN ~ Group*Affiliation + Group*Agency + 
                  (1|Subject), REML = F)
anova(mod_ern)
summary(mod_ern)
qqPlot(resid(mod_ern, 'pearson'))

# --- Find outliers
cz_fcz_rm <- stdResid(data = cz_fcz, mod_ern, 
                      return.data = T, plot = T, 
                      show.loess = T, show.bound = T)

# --- Re-fit without outliers
mod_ern_1 <- lmer(data = filter(cz_fcz_rm, Outlier == 0), 
                  ERN ~ Group*Affiliation + Group*Agency + 
                    (1|Subject), REML = F)
anova(mod_ern_1)
summary(mod_ern_1)
qqPlot(resid(mod_ern_1, 'pearson'))

# Compare models with and with-out 
# scores in control variables
anova(mod_ern_full, mod_ern)

# Coefficent of deteminations
# R2m = only fixed effects, R2c = with random effects
r.squaredGLMM(mod_ern_1)

# Build table
sjPlot::sjt.lmer(mod_ern_full_1, mod_ern_1, cell.spacing = 0.1,
                 show.aic = TRUE, p.numeric = FALSE,
                 string.est = "Estimate",
                 string.ci = "Conf. Int.",
                 string.p = "p-value",
                 depvar.labels = c("∆ERN Amplitude", 
                                   "∆ERN Amplitude"),
                 pred.labels = c("N Errors", "∆Motivation",
                                 "Competition", "Affiliation", 
                                 "Agency", "Competition x Affiliation",
                                 "Competition x Agency") )



# ------ 8) Follow-up analyses model ∆ERN by group ------------------
emm_options(lmerTest.limit = 4000)

# --- Save group estimates
group_means <- emmeans(mod_ern_1, pairwise ~ Group, 
                       lmer.df = 'satterthwaite')
# --- Effect of group
group_means
# --- Compute CIs
confint(group_means)


# --- Effect of Affiliation
emtrends(mod_ern_1, var = 'Affiliation', ~ 1, 
         lmer.df = 'satterthwaite')


# --- Effect of agency
emtrends(mod_ern_1, var = 'Agency', ~ 1, 
         lmer.df = 'satterthwaite')

# --- Effect of agency by group
emtrends(mod_ern_1, var = 'Agency', 
         pairwise ~ Group, lmer.df = 'satterthwaite')


# # ----- Refit for simple slopes
# contrasts(cz_fcz_rm$Group) <- contr.treatment(2, base = 1)
# contrasts(cz_fcz_rm$Group)
# mod_ern_1 <- lmer(data = filter(cz_fcz_rm, Outlier == 0),
#                   ERN ~ Group*Affiliation + Group*Agency +
#                     (1|Subject))
# summary(mod_ern_1)
# 
# contrasts(cz_fcz_rm$Group) <- contr.treatment(2, base = 2)
# contrasts(cz_fcz_rm$Group)
# mod_ern_1 <- lmer(data = filter(cz_fcz_rm, Outlier == 0),
#                   ERN ~ Group*Affiliation + Group*Agency +
#                     (1|Subject))
# summary(mod_ern_1)

# ------ 9) Plot estimates of ∆ERN by group analyses ----------------
# ------ Create the plot
ern_p <- ggplot(data = as.data.frame(group_means$emmeans), 
                aes(y = emmean, x = Group)) + 
  
  geom_hline(yintercept = mean(as.data.frame(group_means$emmeans)[, 2]), 
             linetype = 2) +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Group), 
                width = 0.2, size = 1, alpha = 1) + 
  geom_linerange(aes(ymin = emmean-SE, ymax = emmean+SE), 
                 color = 'gray',size = 3) + 
  geom_point(shape = 18, size = 3.5, color = 'black') +

  scale_color_viridis(option = 'D', discrete = T, 
                      direction = -1, end = .6) +
  
  scale_y_reverse(breaks = c(0, -2, -4, -6, -8)) + 
  
  geom_segment(aes(x = -Inf, y = -8, xend = -Inf, yend = 0), 
               color = 'black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = 'Competition', y = Inf, xend = 'Cooperation', yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  labs(y = expression(bold(paste("Estimated Amplitude (", mu, "V)"))), 
       x = 'Context') +
  
  theme_classic() + 
  
  theme(strip.background = element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_text(color = 'black', face = 'bold', size = 14,
                                    margin = margin(t = 15)),
        axis.title.y = element_text(color = 'black',face = 'bold', size = 14,  
                                    margin = margin(r = 15)),
        axis.text.x = element_text(color = 'black', size = 13, 
                                   angle = 90),
        axis.text.y = element_text(color = 'black', size = 13),
        legend.position = 'none'); ern_p






