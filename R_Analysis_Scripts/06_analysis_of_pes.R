##### ##### #####     Analysis scrips for Alanis et al., 2018   ##### ##### #####
#                           Post-error slowing

# Get helper functions
source('~/Documents/GitHub/Supplementary_Soc_ERN/R_Functions/getPacks.R')
source('~/Documents/GitHub/Supplementary_Soc_ERN/R_Functions/stdResid.R')

# Install and load multiple R packages necessary for analysis.
pkgs <- c('dplyr', 'reshape2',
          'lme4', 'lmerTest',
          'effects','emmeans', 'car', 'MuMIn',
          'ggplot2', 'viridis', 
          'sjPlot', 'cowplot')

getPacks(pkgs)
rm(pkgs)

# ----- 1) Read in the data -----------
PES <- read.table('~/Documents/Experiments/soc_ftask/data_for_r/PES.txt', header=T)

# ----- 2) Dummy code and center ------
contrasts(PES$Group) <- contr.sum(2) 
contrasts(PES$Group)

contrasts(PES$PE_Flankers) <- contr.sum(4) 
contrasts(PES$PE_Flankers)

mean(PES$Err_Amp)
PES <- within( PES, {
  Err_Amp <- Err_Amp - mean(Err_Amp, na.rm = T) 
})

# FIT the model
mod_pes <- lmer(data = PES, pe_slowing ~ Motivation + PE_Flankers +  Agency + 
                  Err_Amp*Group*Affiliation +
                  (1|Subject), REML = F)
anova(mod_pes)
summary(mod_pes)

pes_rm <- stdResid(data = PES, model = mod_pes, plot = T, show.bound = T)

mod_pes <- lmer(data = filter(pes_rm, Outlier == 0), pe_slowing ~ Motivation + 
                  PE_Flankers + Agency +
                  Err_Amp*Group*Affiliation + 
                  (1|Subject), REML = F)
anova(mod_pes)
summary(mod_pes)

# ----- 3) Follow-up analyses ------
emm_options(lmerTest.limit = 2000)

# --- Effect of Group
group_est <- emmeans(mod_pes, pairwise ~ Group, 
                     adjust = 'fdr',
                     at = list(Agency = 0, Affiliation = 0, Err_Amp = 0)); group_est
# Compute CIs
confint(group_est)
# Plot results
emmip(mod_pes, ~ Group, cov.reduce = FALSE, CIs = T, 
      at = list(Affiliation = 0, Err_Amp = 0, Agency = 0))

# Effect of previous error amplitude
emtrends(mod_pes, var = 'Err_Amp', ~ 1, lmer.df = 'satterthwaite', 
         at=list(Agency = 0, Affiliation = 0))
# Effect of affiliation
emtrends(mod_pes, var = 'Affiliation', ~ 1, lmer.df = 'satterthwaite', 
         at=list(Err_Amp = 0, Agency = 0))
# Effect of agency
emtrends(mod_pes, var = 'Agency', ~ 1, lmer.df = 'satterthwaite', 
         at=list(Err_Amp = 0, Affiliation = 0))

# Effect of Affilaiton * ERN by Group 
aff_group <- emmeans(mod_pes, pairwise ~ Err_Amp | Affiliation + Group , 
                     adjust = 'fdr',  
                     at = list(Agency = 0, Affiliation = c(-4, 4), Err_Amp = c(-6.5, 2.5))); aff_group
confint(aff_group)





plot_pes1<- emmip(mod_pes, ~ Err_Amp, mult.name = "Group", cov.reduce = FALSE,
                  at = list(Err_Amp = c(-6.5, 2.5)),
                  CIs = T); plot_pes1

plot_pes1 <- plot_pes1 + coord_cartesian(ylim = c(0, 40)) +
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 40), 
               color = 'black', size = rel(1)) +
  geom_segment(aes(x = '-6.5', y = -Inf, xend = '2.5', yend = -Inf),
               color = 'black', size = rel(1)) +
  labs(y = 'Post-error slowing (ms)', x = 'Amplitude of the ERN (previous error)')+
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15))) +
  annotate('text', x = '-6.5', y = 0,
           label = expression(paste(beta, ' = -0.67*')), 
           parse = TRUE, 
           size = 5, hjust = 0); plot_pes1

cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_S1.pdf', 
                   plot_pes1,  base_height = 5, base_width = 4.5)



plot_pes2 <- emmip(mod_pes, Affiliation ~ Err_Amp | Group, cov.reduce = FALSE, 
                   at = list(Affiliation = c(-4, 4), Err_Amp = c(-6.5, 2.5)), 
                   CIs = T); plot_pes2


plot_pes2 <- plot_pes2 + coord_cartesian(ylim = c(-20, 60)) +
  geom_segment(aes(x = -Inf, y = -20, xend = -Inf, yend = 60), 
               color = 'black', size = rel(1)) +
  geom_segment(aes(x = '-6.5', y = -Inf, xend = '2.5', yend = -Inf),
               color = 'black', size = rel(1)) +
  labs(y = 'Post-error slowing (ms)', x = 'Amplitude of the ERN (previous error)')+
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(color = 'black', size=14, face='bold'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = rel(1.1)),
        axis.ticks.length = unit(.1, 'cm'),
        axis.text.x = element_text(size = 13, color = 'black'),
        axis.text.y = element_text(size = 13, color = 'black'),
        axis.title.x = element_text(size = 14, face = 'bold', 
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 14, face = 'bold', 
                                    margin = margin(r = 15)),
        legend.text = element_text(color = 'black', size=12),
        legend.title = element_text(color = 'black', size = 12, face='bold')); plot_pes2

cowplot::save_plot('~/Documents/Experiments/soc_ftask/paper_figs/Fig_S2.pdf', 
                   plot_pes2,  base_height = 5, base_width = 9)


