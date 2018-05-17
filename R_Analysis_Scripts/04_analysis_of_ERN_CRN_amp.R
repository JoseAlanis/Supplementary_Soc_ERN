##### ##### #####     Analysis scrips for Alanis et al., 2018   ##### ##### #####
#                           LMER models for physiological
#                                     data

# Get helper functions
source('~/Documents/GitHub/Supplementary_Soc_ERN/R_Functions/getPacks.R')
source('~/Documents/GitHub/Supplementary_Soc_ERN/R_Functions/stdResid.R')

# Install and load multiple R packages necessary for analysis.
pkgs <- c('dplyr', 'reshape2',
          'lme4', 'lmerTest',
          'effects', 'emmeans', 'car', 'MuMIn',
          'ggplot2', 'viridis')

getPacks(pkgs)
rm(pkgs)

# ------ READ in the data -------------------------------------------
load('~/Documents/Experiments/soc_ftask/data_for_r/Phys_Data.RData')

# ------ 1) Select relevant observations ----------------------------
# Incompatible trials, electrode CZ and Fz and
# time window 0 to 100 ms after response
Phys <- filter(Phys, Flankers == 'Incompatible')
Phys <- filter(Phys, Electrode == 'Cz' | Electrode == 'FCz')

ERN_CRN <- filter(Phys, Time >= 0 & Time <= 100)

# ### Descriptive statistics
ERN_CRN %>% group_by(Reaction) %>% summarise(M = mean(Amp), 
                                         SD = sd(Amp), 
                                         SE = sd(Amp) / sqrt(sum(!is.na(Amp))))

ERN_CRN %>% summarise(M = mean(Affiliation), 
                          SD = sd(Affiliation))

ERN_CRN %>% summarise(M = mean(Agency), 
                      SD = sd(Agency))


# DUMMY CODE electrode variable
ERN_CRN$Electrode <- factor(ERN_CRN$Electrode)
contrasts(ERN_CRN$Electrode) <- contr.sum(2); contrasts(ERN_CRN$Electrode)
contrasts(ERN_CRN$Group) <- contr.sum(2); contrasts(ERN_CRN$Group)

ERN_CRN$Reaction <- factor(ERN_CRN$Reaction, levels = c("Incorrect", "Correct"))
contrasts(ERN_CRN$Reaction) <- contr.sum(2); contrasts(ERN_CRN$Reaction)

contrasts(ERN_CRN$Group) <- contr.treatment(2, base=2); contrasts(ERN_CRN$Group)
contrasts(ERN_CRN$Reaction) <- contr.treatment(2, base=1); contrasts(ERN_CRN$Reaction)


# SET UP AND FIT full model
mod_rns_full <- lmer(data = ERN_CRN, 
                     Amp ~ Total_Errors + Interest + 
                       Reaction*Group*Affiliation + Reaction*Group*Agency + 
                       (1+Reaction|Subject))
anova(mod_rns_full)


# SET UP AND FIT the reported model
mod_ern <- lmer(data = ERN_CRN, 
                Amp ~ Interest + 
                  Reaction*Affiliation + Reaction*Group*Agency + 
                  (1+Reaction|Subject))
anova(mod_ern)
summary(mod_ern)
qqPlot(resid(mod_ern, 'pearson'))


# Find outliers
ERN_CRN_rm <- stdResid(data = ERN_CRN, mod_ern, 
                      return.data = T, plot = T, 
                      show.loess = T, show.bound = T)

# re-fit without outliers
mod_ern_1 <- lmer(data = filter(ERN_CRN_rm, Outlier == 0), 
                  Amp ~ Interest + 
                    Reaction*Group*Affiliation + Reaction*Group*Agency + 
                    (1+Reaction|Subject))
anova(mod_ern_1)
summary(mod_ern_1)
qqPlot(resid(mod_ern_1, 'pearson'))

r.squaredGLMM(mod_ern_1)

visreg::visreg(mod_ern, 'Reaction', by='Group', partial = F,  ylim =c(8, -8))
visreg::visreg(mod_ern, 'Affiliation', by='Reaction', partial=F, ylim =c(8, -8), overlay=T)
visreg::visreg(mod_ern, 'Agency', by = 'Group', partial=F, ylim =c(8, -8), cond=list(Reaction = 'Incorrect'))
visreg::visreg(mod_ern, 'Agency', by = 'Group', partial=F, ylim =c(8, -8), cond=list(Reaction = 'Correct'))

# ------ Follow-up analyses model ∆ERN by group ---------------------
emm_options(pbkrtest.limit = 1000)
summary(ref_grid(mod_ern_1), infer=T)


# Save group estimates
group_means <- emmeans(mod_ern_1, 
                       pairwise ~ Group)

# Effect of group
as.data.frame(group_means$emmeans)
group_means$contrasts
# Compute CIs
mutate(as.data.frame(group_means$contrasts), 
       LCL = estimate - SE * 1.96, 
       UCL = estimate + SE * 1.96)


# Save reaction estimates
reaction_means <- emmeans(mod_ern_1, 
                          pairwise ~ Reaction)

# Effect of group
as.data.frame(reaction_means$emmeans)
reaction_means$contrasts
# Compute CIs
mutate(as.data.frame(reaction_means$contrasts), 
       LCL = estimate - SE * 1.96, 
       UCL = estimate + SE * 1.96)


# Save reaction estimates by group
reaction_means <- emmeans(mod_ern_1, 
                          pairwise ~ Reaction | Group); reaction_means

reaction_means <- emmeans(mod_ern_1, 
                          pairwise ~ Group | Reaction); reaction_means

# Effect of group
as.data.frame(reaction_means$emmeans)
reaction_means$contrasts
# Compute CIs
mutate(as.data.frame(reaction_means$contrasts), 
       LCL = estimate - SE * 1.96, 
       UCL = estimate + SE * 1.96)



# Effects of affiliation
# overall
emtrends(mod_ern_1, var = 'Affiliation', ~ 1)

# by reaction
emtrends(mod_ern_1, var = 'Affiliation', pairwise ~ Reaction)

# emmeans at +1 SD aff and -1 SD aff
slope_aff <- ref_grid(mod_ern_1, at = list(Affiliation = c(-3.7,3.7)))
aff_means <- emmeans(slope_aff, pairwise ~ Reaction | Affiliation)
# Compute CIs
mutate(as.data.frame(aff_means$contrasts), 
       LCL = estimate - SE * 1.96, 
       UCL = estimate + SE * 1.96)

# Effect of agency
# overall
emtrends(mod_ern_1, var = 'Agency', ~ 1)

emtrends(mod_ern_1, var = 'Agency', pairwise ~ Group | Reaction)

# emmeans at +1 SD aff and -1 SD aff
slope_ag <- ref_grid(mod_ern_1, at = list(Agency = c(-26, 26)))
ag_means <- emmeans(slope_ag, pairwise ~ Reaction | Agency)
# Compute CIs
mutate(as.data.frame(ag_means$contrasts), 
       LCL = estimate - SE * 1.96, 
       UCL = estimate + SE * 1.96)
