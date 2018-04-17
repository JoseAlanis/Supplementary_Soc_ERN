##### ##### #####     Analysis scrips for Alanis et al., 2018   ##### ##### #####
#                           LMER models for physiological
#                                     data

# Get helper functions
source('./Documents/GitHub/Supplementary_Soc_ERN/R_Functions/getPacks.R')
source('./Documents/GitHub/Supplementary_Soc_ERN/R_Functions/stdResid.R')

# Install and load multiple R packages necessary for analysis.
pkgs <- c('dplyr', 'reshape2',
          'lme4', 'lmerTest',
          'effects', 'emmeans', 'car', 'MuMIn',
          'ggplot2', 'viridis')

getPacks(pkgs)
rm(pkgs)

# ------ READ in the data -------------------------------------------
load('~/Documents/Experiments/soc_ftask/data_for_r/Phys_Data.RData')



# ------ 1) PLOT individual averages --------------------------------
ggplot(Phys[ (Phys$Electrode == 'FCz') & Phys$Flankers == 'Incompatible' & Phys$Time > -500 , ], 
       aes(Time, Amp, color = Reaction, linetype = Group)) + 
  facet_wrap(~ Subject) +
  geom_vline(xintercept = c(0, 100), color = 'black', size = rel(0.5), linetype = 3) +
  geom_hline(yintercept = c(0), color = 'black', size = rel(0.5), linetype = 3) +
  geom_line(size = rel(0.7)) + 
  scale_color_viridis(option = 'D', discrete = T, begin = .15, end = .85, direction = -1) + 
  scale_y_reverse(breaks = c(-30, -15, 0, 15, 30)) +
  scale_linetype_manual(values = c(1,6)) +
  xlab("\n Time") + 
  ylab("Amplitude [µV] \n") +
  ggtitle("Grand Average ERP") +
  theme_classic() +
  geom_segment(aes(x = -Inf, y = -30, xend = -Inf, yend = 30), 
               color='black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = -500, y = Inf, xend = 1000, yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  theme(
    plot.title = element_text(color = "black", hjust = 0.5, vjust = 1, size = rel(1.5)),
    strip.background = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(color='black'),
    axis.text.y  = element_text(hjust = 1, color = 'black'),
    legend.text = element_text(color="black", size=rel(1)),
    legend.position = "bottom"
  ) 


# ------ 2) COMPUTE ∆ERN wave incorrect-correct ---------------------
# From long to wide
data_wide <- dcast(Phys, Group + Subject + Agency + Affiliation + Interest +
                     Incomp_Errors + Total_Errors + Flankers + Electrode + 
                     Time ~ Reaction,
                   value.var = c('Amp'))

data_wide <- data_wide[data_wide$Flankers == 'Incompatible', ]

# NAs created?
unique(is.na(data_wide))

# Compute ERN
data_wide <- data_wide %>% mutate(ERN = Incorrect - Correct)


# ------ PLOT ∆ERN ----------------------------------------------
# AVERAGE for plot
Ave_ERN <- plyr::ddply(data_wide, 
                       c('Group', 'Time', 'Electrode'), dplyr::summarise,
                       
                       N = sum(!is.na(ERN)),
                       M_Amp = mean(ERN, na.rm=T),
                       sd   = sd(ERN, na.rm=T),
                       se   = sd / sqrt(N),
                       # Compute t-staistic for confidence interval (quantile t, with N-1 degrees of freedom)
                       # and use it to compute individual ci (multiplier for se):
                       ci   = se * qt(.95/2 + .5, N-1)
                       )


# CREATE the plot
ggplot(filter(Ave_ERN, Electrode == 'FCz' | Electrode == 'Cz', Time >= -500), 
       
       aes(Time, M_Amp, fill = Group, color = Group)) + 
  
  facet_wrap(~ Electrode, scales = 'free_y', ncol = 4) + 
  theme_classic() + 
  
  annotate("rect", xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, alpha = .1) +
  geom_vline(xintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  geom_hline(yintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  
  geom_ribbon(aes(ymin = M_Amp - se, ymax = M_Amp + se), alpha = .3, colour = NA) +
  geom_line( size = rel(0.9)) + 
  
  scale_y_reverse(breaks = c(-7 ,-3.5, 0, 3.5, 7)) +
  scale_x_continuous(breaks = c(-500, 0, 500, 1000)) +
  
  scale_color_viridis(option = 'D', discrete = T, direction = 1, end = .6) +
  scale_fill_viridis(option = 'D', discrete = T, direction = 1, end = .6) +
  
  labs(x = expression(bold("Time " ['ms'])), 
       y = expression(bold("Amplitude " ['µV'])), 
       title="Grand Average ERN") +
  
  geom_segment(aes(x = -Inf, y = -7, xend = -Inf, yend = 7), color = 'black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = -500, y = Inf, xend = 1000, yend = Inf), color = 'black', size = rel(1), linetype = 1) +
  
  theme(
    strip.background = element_blank(),
    axis.line = element_blank(),
    
    plot.title = element_text(color = "black", hjust = 0.5, vjust = 1, size = 15, face = 'bold'),
    strip.text = element_text(color = 'black', size = 13, face = 'bold'),
    axis.title.y = element_text(color = 'black', size = 13, face = 'bold'),
    axis.title.x = element_text(color = 'black', size = 13, face = 'bold', 
                                margin = margin(t = 15), hjust = 0.25),
    axis.text.x = element_text(color = 'black', size = 12),
    axis.text.y  = element_text(hjust = 1, color = 'black', size = 12),
    legend.text = element_text(color = "black", size = 13),
    legend.key.size = unit(1, 'cm'),
    legend.title = element_blank(),
    
    legend.position = "bottom")



# ------ MODEL ∆ERN as function of electrodes -------------------
# AVERAGE for model
Ave_Elec <- plyr::ddply(filter(data_wide, (Time >= 0 & Time <= 100) & (Electrode == 'FCz' | Electrode == 'Cz' | Electrode == 'CPz')), 
                       c('Group', 'Electrode', 'Subject'), dplyr::summarise,
                       
                       N = sum(!is.na(ERN)),
                       M_Amp = mean(ERN, na.rm = T),
                       sd   = sd(ERN, na.rm = T),
                       se   = sd / sqrt(N),
                       # Compute t-staistic for confidence interval 
                       # (quantile t, with N-1 degrees of freedom)
                       # and use it to compute individual ci (multiplier for se):
                       ci   = se * qt(.95/2 + .5, N-1)
)

# DUMMY CODE electrode variable
Ave_Elec$Electrode <- factor(Ave_Elec$Electrode)
contrasts(Ave_Elec$Electrode) <- contr.sum(3)
contrasts(Ave_Elec$Electrode)

# SET UP AND FIT the model
mod_elec <- lmer(data = Ave_Elec, 
                 M_Amp ~ Electrode + (1|Subject))
anova(mod_elec)
summary(mod_elec)
qqPlot(resid(mod_elec, 'pearson'))

## Detect outlying observations
Elec_rm <- stdResid(data = Ave_Elec, mod_elec, 
                    plot = T, show.loess = T, show.bound = T)

## re-fit without outliers
mod_elec_1 <- lmer(data = filter(Elec_rm, Outlier == 0), 
                   M_Amp ~ Electrode + (1|Subject))
anova(mod_elec_1)
summary(mod_elec_1)

qqPlot(resid(mod_elec_1, 'pearson'))
stdResid(data = filter(Elec_rm, Outlier == 0), mod_elec_1, 
         return.data = F, plot = T, 
         show.loess = T, show.bound = T)

r.squaredGLMM(mod_elec_1)

# ------ Follow-up analyses model electrodes --------------------
# Pairwise comparissons
emmeans(mod_elec_1, pairwise ~ Electrode, 
        adjust = 'bonferroni', lmer.df ='satterthwaite')
# Plot results
emmip(mod_elec_1, ~ Electrode, CIs = T)

# Save electrode estimates
est_elec <- emmeans(mod_elec_1, pairwise ~ Electrode,
                   adjust = 'bonferroni', lmer.df ='satterthwaite')

# Effect of trial type
as.data.frame(est_elec$emmeans)
# Compute CIs
mutate(as.data.frame(est_elec$contrasts), 
       LCL = estimate - SE * 1.96, 
       UCL = estimate + SE * 1.96)

# Plot estimates
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


# ------ MODEL ∆ERN as function of context and personality ------
## Only keep FCz and Cz
cz_fcz <- filter(data_wide, Electrode == 'Cz' | Electrode == 'FCz' )
cz_fcz <- filter(cz_fcz, Time >= 0 & Time <= 100 )

# DUMMY CODE electrode variable
cz_fcz$Electrode <- factor(cz_fcz$Electrode)
contrasts(cz_fcz$Electrode) <- contr.sum(2)
contrasts(cz_fcz$Electrode)

contrasts(cz_fcz$Group) <- contr.sum(2)
contrasts(cz_fcz$Group)

cz_fcz <- within( cz_fcz, {
  Tot_Errors <- Total_Errors - mean( Total_Errors, na.rm=T )
  I_Errors <- Incomp_Errors - mean( Incomp_Errors, na.rm=T )
})

# ### Descriptive statistics
cz_fcz %>% group_by(Group) %>% summarise(M = mean(ERN), 
                     SD = sd(ERN), 
                     SE = sd(ERN) / sqrt(sum(!is.na(ERN))))



# SET UP AND FIT full model
mod_ern_full <- lmer(data = cz_fcz, 
                ERN ~ Tot_Errors + Interest + Group*Affiliation + Group*Agency + 
                  (1 |Subject))
anova(mod_ern_full)


# SET UP AND FIT the model
mod_ern <- lmer(data = cz_fcz, 
                ERN ~ Group*Affiliation + Group*Agency + 
                  (1|Subject))
anova(mod_ern)
summary(mod_ern)
qqPlot(resid(mod_ern, 'pearson'))

# visreg::visreg(mod_ern, 'Time', by='Time')
# visreg::visreg2d(mod_ern, 'Time', 'Affiliation')

# Find outliers
cz_fcz_rm <- stdResid(data = cz_fcz, mod_ern, 
                      return.data = T, plot = T, 
                      show.loess = T, show.bound = T)

# re-fit without outliers
mod_ern_1 <- lmer(data = filter(cz_fcz_rm, Outlier == 0), 
                  ERN ~ Group*Affiliation + Group*Agency + 
                    (1|Subject))
anova(mod_ern_1)
summary(mod_ern_1)
qqPlot(resid(mod_ern_1, 'pearson'))

stdResid(data = filter(cz_fcz_rm, Outlier == 0), mod_ern_1, 
         return.data = F, plot = T, 
         show.loess = T, show.bound = T)

r.squaredGLMM(mod_ern_1)



# ------ Follow-up analyses model ∆ERN by group ---------------------
emm_options(pbkrtest.limit = 4000)
summary(ref_grid(mod_ern_1), infer=T)

# Save group estimates
group_means <- emmeans(mod_ern_1, pairwise ~ Group)

# Effect of group
as.data.frame(group_means$emmeans)
group_means$contrasts
# Compute CIs
mutate(as.data.frame(group_means$contrasts), 
       LCL = estimate - SE * 1.96, 
       UCL = estimate + SE * 1.96)


# Effect of affiliation by group
emtrends(mod_ern_1, var = 'Affiliation', ~ 1)

# Effect of agency
emtrends(mod_ern_1, var = 'Agency', ~ 1)

# Effect of agency by group
emtrends(mod_ern_1, var = 'Agency', pairwise ~ Group)



# --- Refit for simple slopes---

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

# ------ Save ∆ERN averages for later analyses ------------------
Ave_ERN <- plyr::ddply(cz_fcz, 
                       c('Group', 'Subject'), dplyr::summarise,
                       ERN_Amp = mean(ERN, na.rm=T)
                       )
write.table(Ave_ERN, './Desktop/RefTask_2018_Data/DATA/Ave_ERN.txt', row.names = F)


# ------ PLOT ∆ERN by Group -------------------------------------
# AVERAGE for plot
Ave_Group <- plyr::ddply(filter(data_wide, Electrode == 'Cz' | 
                                  Electrode == 'FCz', Time >= -500), 
                        c('Group', 'Time'), dplyr::summarise,
                        
                        N = sum(!is.na(ERN)),
                        M_Amp = mean(ERN, na.rm = T),
                        sd   = sd(ERN, na.rm = T),
                        se   = sd / sqrt(N),
                        # Compute t-staistic for confidence interval 
                        # (quantile t, with N-1 degrees of freedom)
                        # and use it to compute individual ci (multiplier for se):
                        ci   = se * qt(.95/2 + .5, N-1)
)

# CREATE the plot
pdf('./Desktop/GA_ERN.pdf', width = 5.5, height = 4)

ggplot(Ave_Group, 
       
       aes(Time, M_Amp, fill = Group, color = Group)) + 
  
  theme_classic() + 
  
  annotate("rect", xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf, alpha = .1) +
  geom_vline(xintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  geom_hline(yintercept = c(0), color='black', size = rel(0.5), linetype = 3) +
  
  geom_ribbon(aes(ymin = M_Amp - ci, ymax = M_Amp + ci), alpha = .25, colour = NA) +
  
  geom_line( size = rel(0.8) ) + 
  
  scale_y_reverse(breaks = c(-7, -3.5, 0, 3.5, 7)) +
  scale_x_continuous(breaks = c(-500, -250, 0, 250,  500, 750, 1000)) +
  
  scale_color_viridis(option = 'B', discrete = T, 
                      direction = -1, end = .6) +
  scale_fill_viridis(option = 'B', discrete = T, 
                     direction = -1, end = .6) +
  
  labs(x = expression(bold("Time " ['ms'])), 
       y = expression(bold("Amplitude " ['µV'])), 
       title="a) Grand Average ERN") +
  
  geom_segment(aes(x = -Inf, y = -7, xend = -Inf, yend = 7), 
               color = 'black', size = rel(1), linetype = 1) +
  geom_segment(aes(x = -500, y = Inf, xend = 1000, yend = Inf), 
               color = 'black', size = rel(1), linetype = 1) +
  
  theme(
    strip.background = element_blank(),
    axis.line = element_blank(),
    
    plot.title = element_text(color = "black", hjust = 0, vjust = 1, 
                              size = 15, face = 'bold', margin = margin(b = 15)),
    strip.text = element_text(color = 'black', size = 13, face = 'bold'),
    axis.title.y = element_text(color = 'black', size = 13, face = 'bold'),
    axis.title.x = element_text(color = 'black', size = 13, face = 'bold', 
                                margin = margin(t = 15), hjust = .5),
    axis.text.x = element_text(color = 'black', size = 12),
    axis.text.y  = element_text(hjust = 1, color = 'black', size = 12),
    legend.text = element_text(color = "black", size = 12),
    legend.key.size = unit(.8, 'cm'),
    legend.title = element_blank(),
    
    legend.position = "bottom")

dev.off()










