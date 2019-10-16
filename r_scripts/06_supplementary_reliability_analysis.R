# --- author: Jose C. Garcia Alanis
# --- encoding: utf-8
# --- r version: 3.6.0 (2019-04-26)
# --- content: reliability analysis
# --- version: Wed Oct 16 13:08:24 2019

# comments:
# SINGLE TRIAL DATA are needed for this analysis step.

# --- 1) Set paths  -----------------------------------------------------------
project_path <- '/Volumes/TOSHIBA/manuscripts_and_data/soc_ern/'
setwd(project_path); rm(list = ls()) # <- clean workspace

# --- 2) Get helper functions and packages ------------------------------------
source('./r_functions/getPacks.R')

# Install and load multiple R packages necessary for analysis.
getPacks(c('dplyr'))

# --- 3) Import data ----------------------------------------------------------
# get data frame containing behavioral data
errors <- read.table('./revision/rev_data/errors_data.txt', header = T)
# get demographic data
perso <- read.table('./data_for_r/all_perso.txt', header = T)

# subject identifier to iteger (for loop)
ID <- as.integer(unique(gsub(errors$id, 
                             pattern = 'data_', 
                             replacement = '')))

# also extract group variable and age
Group <- select(errors, id, group)
names(Group) <- c('ID', 'Group')

Age <- select(perso, id, Age)
names(Age) <- c('ID', 'Age')


# *** Loop through files and keep relevant observations ***
# (i.e., FCz and Cz, impopatible trials), then bind them together.

# set seed
set.seed(46)

# minumum number of trials to allow in analysis
min_trials <- 6

# start loop
for (i in ID) {
  # Function start time
  if (i == first(ID)) {
    starttime <-  Sys.time()
  }
  # path to file in question
  file = paste('./revision/rev_data/rev_Trial_Amps/data_', i, '_Trial_Amps.txt', 
               sep = '')
  # select relevant trials and time window
  dat <- read.table(file, header = T)
  dat <- dat %>% filter(Flankers == 'incompatible'
                        & Electrode %in% c('Cz', 'FCz') 
                        & Reaction == 'Incorrect'
                        & Range == '0 - 100 ms') %>%
    group_by(ID, Trial_Nr) %>% 
    # compute mean amplitude
    summarise(m_amp = mean(Mean_Amp))
  
  
  # check if combined data.frame already exists
  if (!exists('single_trials')) {
    single_trials <- data.frame()
  }
  single_trials <- rbind(single_trials, data.frame(dat))
  
  # number of trials
  n_t <- dat %>% 
    select(ID, Trial_Nr)
  
  # if subject has enough trials
  if (length(unique(n_t$Trial_Nr)) >= min_trials) {
    trials <- dat %>% 
      select(ID, Trial_Nr, m_amp) %>%
      filter(Trial_Nr %in% sample(Trial_Nr, min_trials)) %>%
      tidyr::spread(Trial_Nr, m_amp)
    
    # rename columns
    names(trials) <- c('ID', as.character(seq(min_trials)))
    
    # check if combined wide data.frame already exists
    if (!exists('sample_trials_wide')) {
      sample_trials_wide <- data.frame()
    }
    # data from all subjects to data.frame
    sample_trials_wide <- rbind(sample_trials_wide, data.frame(trials))
  }

  # Print when done.
  print(paste(i, 'done'))
  # Print function run time
  if (i == last(ID)) {
    print( paste('Runtime:',
                 as.numeric(difftime(Sys.time(), starttime, units = 'min')),
                 'mins'))
  }
}

# --- 4) Compute reliability measures -----------------------------------------
# check number of trials
n_trials <- single_trials %>% 
  select(ID, Trial_Nr) %>% 
  unique() %>% 
  group_by(ID) %>% 
  summarise(N = sum(!is.na(Trial_Nr))) %>% 
  arrange(N); summary(n_trials)

# merge with group information
sample_trials_wide <- merge(Group, sample_trials_wide, 'ID')
single_trials <- merge(Group, single_trials, 'ID')

# compute alpha for overall sample
psych::alpha(select(sample_trials_wide, contains('X')))

# compute alpha for coopration group
psych::alpha(select(filter(sample_trials_wide, Group == 'Cooperation'), contains('X')))
# compute alpha for competition group
psych::alpha(select(filter(sample_trials_wide, Group == 'Competition'), contains('X')))

# --- 4) Compute mean age in groups -------------------------------------------
# ids
coop <- filter(sample_trials_wide, Group == 'Cooperation')$ID
comp <- filter(sample_trials_wide, Group == 'Competition')$ID

# individual ages
age_in_cooperation_group <- filter(perso, id %in% coop)$Age
age_in_competition_group <- filter(perso, id %in% comp)$Age

# means
mean(age_in_cooperation_group)
mean(age_in_competition_group)

# sample sizes
table(filter(perso, id %in% sample_trials_wide$ID)$group)

# --- 5) Check experimental effects -------------------------------------------
# load necessary packages
require(lme4)
require(lmerTest)

# only keep individuals with enough trials
data_for_mod <- single_trials %>% 
  filter(ID %in% sample_trials_wide$ID) %>%
  mutate(ID = factor(ID))

# combine trials and extraversion meausres
extra <- select(perso, id, affiliation, agency)
names(extra) <- c('ID', 'affiliation', 'agency')
data_for_mod <- merge(data_for_mod, extra, 'ID')

# set up model
group_mod <- lmer(data = data_for_mod, 
                  m_amp ~ Group*affiliation + Group*agency + (1|ID),
                  contrasts = list(Group = 'contr.sum'))

# check if reported effects hold
anova(group_mod, ddf = 'Kenward-Roger') # <- interactions + affiliation effect hold
