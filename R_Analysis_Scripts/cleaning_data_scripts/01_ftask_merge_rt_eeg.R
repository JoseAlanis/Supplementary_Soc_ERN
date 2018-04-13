##### ##### #####     Analysis scrips for Alanis et al., 2018   ##### ##### #####
#                         Merge behavioural and eeg data
#


# Load necessary packages
require(dplyr)
#se <- function(x) sqrt(var(x)/length(x)) # uncomment if required

# ----- 1) READ in behavioural data ---------------------------------
load('~/Documents/Experiments/soc_ftask/data_for_r/all_behav.RData')

# Only keep relevent trials
# i.e, trials where subject perfomed a response
rt_data <- all_behav[!all_behav$RT == 88888 & !all_behav$RT == 99999, ] 
rt_data$Reaction <- factor(rt_data$Reaction)

rt_data <- rt_data %>% select(ID, Group, Trial_Nr, RT, Reaction, Flankers)

# Read in meta-data: channel locations
chans <- read.table('~/Documents/Experiments/soc_ftask/meta_dat/chanlocs_ftask.txt', 
                    sep=',', header=T)
chans$labels <- as.character(chans$labels)

for (i in 1:76) {
  
  # --- SAVE time for runtime summary ---
  start.time <- Sys.time()
  
  
  # ------ ***** 1) Load particpnat's data        ***** -------------
  # --- Paricipant in question ---
  ID <- paste('data_', i, sep = '')
  
  # --- CHECK for trials containing EEG artifacts ---
    if (file.size(paste('~/Documents/Experiments/soc_ftask/from_matlab/RefTask_Flagged_B500_300/', 
                        i, 
                        '_RefTask_Flagged.txt',sep='') ) > 0) {
      # create a vector with trial indetifier
      x <- read.table(paste('~/Documents/Experiments/soc_ftask/from_matlab/RefTask_Flagged_B500_300/', 
                            i, 
                            '_RefTask_Flagged.txt', sep=''))
      x <- as.numeric(x$V1)
    } else {
      # empty vector if no bad trials found
      x <- c() 
    }
  
  # --- CHECK if nr of trials matches the nr of EEG epochs ---
  EEG1 <- length(count.fields(paste('~/Documents/Experiments/soc_ftask/from_matlab/RefTask_RT_from_EEG/RT',
                                    i, 
                                    '_RefTask.txt', sep = ''))) - 1
  Pre1 <- nrow(rt_data[rt_data$ID %in% ID , ])  
  
  # --- Select paricipant data ---
  behav <- rt_data[rt_data$ID %in% ID , ]
  
  #  --- data_64 has less trials (problems with the EEG machine) ---
  if (ID == 'data_64') {
    behav <- behav[127:nrow(behav), ]
  }
  
  # --- Create trial variable to match epoch nr ---
  behav$Trial <- behav$Trial_Nr
  behav$Trial_Nr <- 1:nrow(behav)

  # --- GET EEG epochs, name columns, paste participant ID
  # --- and remove artifact distorted epochs
  phys <- read.table(paste('~/Documents/Experiments/soc_ftask/from_matlab/RefTask_for_R_B500_300/',
                           i,
                           '_bined_data.txt', sep = ''))
  # --- name columns and add subject id
  names(phys) <- c(chans$labels, 'Trial_Nr', 'Time')
  phys$ID <- ID
  
  
  # ------ ***** 2) Merge behav. and phys. data   ***** -------------
  temp <- merge(behav, phys, c('ID', 'Trial_Nr'))
  rm(phys, behav)
  
  temp <- temp %>% arrange(ID, Trial_Nr)
  
  
  # ------ ***** 3) Create PES data set for analysis ***** ----------
  # --- Indetify correct trials preceded by an error
  # --- (this is to control for post error processes 
  # --- in correct trials, e.g., post erros slowing)
  Index <- which((temp$Reaction == 'correct' ) & 
                   dplyr::lag(temp$Reaction) == 'incorrect')
  
  # --- Indetify correct trials preceded by correct trials
  Index_2 <- which((temp$Reaction == 'correct' ) & 
                   dplyr::lag(temp$Reaction) == 'correct')
  
  # --- Create PES data.frame ---
  PES <- filter(temp, Trial_Nr %in% temp[Index, 'Trial_Nr'])
  PES <- rbind(PES, filter(temp, Trial_Nr %in% 
                           (temp[Index, 'Trial_Nr']-1) ) )
  PES <- select(PES, ID:Trial, Fz, FCz, Cz, CPz, Pz, Time)
  
  PES <- mutate(PES, Trial_Index = ifelse(Trial_Nr %in% 
                                          temp[Index, 'Trial_Nr'], 'N+1', 'Error'))
  
  # --- Transform data frame and filter relevant obs. ---
  PES <- filter(PES, Time >= 0 & Time <= 100)
  PES <- tidyr::gather(PES, Electrode, Amplitude, Fz:Pz, 
                       factor_key = TRUE)
  PES <- filter(PES, Electrode == 'FCz' | Electrode == 'Cz')
  PES$Electrode <- factor(PES$Electrode)
  
  # --- Average data.frame  ---
  PES <- plyr::ddply(PES,  c('Group', 'ID', 'Flankers', 
                             'Reaction', 'Trial_Nr', 'RT', 'Trial_Index'), 
                     dplyr::summarise,
                     M_Amp = mean(Amplitude, na.rm = T)
                     )
  # --- Arrange data.frame ---
  PES <- PES %>% arrange(ID, Trial_Nr)
  
  # --- Chunking by trial index and trial no. ---
  PES <- PES %>%
    group_by(ID, Trial_Index) %>%
    mutate(positionInSubject = 1:n())
  
  # --- Remove trials distored by EEG artifacts ---
  PES <- filter(PES, !positionInSubject %in% 
                  PES[PES$Trial_Nr %in% x , ]$positionInSubject)
  
  
  # ******    SAME FOR POST-CORRECT TRIALS    *****
  PCT <- filter(temp, Trial_Nr %in% temp[Index_2, 'Trial_Nr'])
  PCT <- filter(PCT, Time >= 0 & Time <= 100)
  PCT <- select(PCT, ID:Trial, Fz, FCz, Cz, CPz, Pz, Time)
  PCT <- gather(PCT, Electrode, Amplitude, Fz:Pz, 
                       factor_key = TRUE)
  PCT <- filter(PCT, Electrode == 'FCz' | Electrode == 'Cz')
  PCT$Electrode <- factor(PCT$Electrode)
  # --- Average data.frame  ---
  PCT <- plyr::ddply(PCT,  c('ID', 'Flankers'), 
                     dplyr::summarise,
                     
                     M_RT = mean(RT, na.rm = T))
  # --- Arrange data.frame ---
  PCT <- reshape2::dcast(PCT, ID ~ Flankers, value.var = 'M_RT')
  PCT <- mutate(PCT, M_Corr_RT = mean(compatible:neutral))
  
  # --- Merge both data.frame
  PES <- merge(PES, PCT, 'ID')
  
  # --- Save the data.frame ---
  write.table(PES, file = paste('~/Documents/Experiments/soc_ftask/merge_output/PES/', 
                                ID, 
                                '_PES.txt', sep=''), 
              row.names = F, sep = '\t')
  
  # --- remove data.frame after saving ---
  rm(PES, PCT)
  
  
  
  # ------ ***** 4) Create trial amplitude averages ***** ----------
  # --- Remove trials distored by EEG artifacts
  temp <- temp[!temp$Trial_Nr %in% x , ]
  # --- REMOVE post error trials
  temp <- temp[!temp$Trial_Nr %in% temp[Index, 'Trial_Nr'],  ]
  # --- REMOVE trials with RT < 100 ms
  temp <- temp[!temp$RT < 100  , ]
  temp <- filter(temp, Reaction == 'correct' |  Reaction == 'incorrect')
  # --- arrange data.frame ---
  temp <- temp %>% arrange(ID, Trial_Nr)


  # --- From wide to long format ---
  temp_long <- tidyr::gather(temp, 
                             Electrode, 
                             Amplitude, 
                             Fp1:O2, 
                             factor_key = TRUE)

  # --- COMPUTE peak latencies and mean amplitudes for
  # --- "ERN" and "Pe" and store the in a data.frame
  Latencies_ERN <- plyr::ddply(filter( temp_long, Time >= 0 & Time <= 100,
                                       Electrode %in% c('Fz', 'FCz', 'Cz', 'CPz', 'Pz')), 
                           c('ID', 'Group', 'Trial_Nr', 'RT', 'Reaction', 'Flankers', 'Electrode'), 
                           
                           dplyr::summarise, 
                           
                           Mean_Amp = mean(Amplitude), 
                           lat = Time[which.min(Amplitude)],
                           Range = '0 - 100 ms')
  # --- Pe ---
  Latencies_Pe <- plyr::ddply(filter( temp_long, Time > 150 & Time < 300, 
                              Electrode %in% c('Fz', 'FCz', 'Cz', 'CPz', 'Pz')),
                              c('ID', 'Group', 'Trial_Nr', 'RT', 'Reaction', 'Flankers', 'Electrode'), 
                              
                              dplyr::summarise, 
                              
                              Mean_Amp = mean(Amplitude), 
                              lat = Time[which.max(Amplitude)],
                              Range = '150 - 300 ms')
  
  # --- Bin ERN and Pe latencies together ---
  Trial_Amps <- rbind(Latencies_ERN, Latencies_Pe)
  # --- Save data.frame ---
  write.table(Trial_Amps, 
              file = paste('~/Documents/Experiments/soc_ftask/merge_output/Trial_Amplitudes/', 
                           ID, 
                           '_Trial_Amps.txt', sep = ''), 
              row.names = F, sep = '\t')
  
  # --- Remove data.frame after saving ---
  rm(temp_long, Trial_Amps)
  
  
  # --- COMPUTE indivdual average --- 
  temp <- temp %>% 
    select(-Trial_Nr, -Trial, -RT) %>% 
    group_by(ID, Group, Reaction, Flankers, Time) %>% 
    summarise_all(funs(M = mean)) 
    # summarise_all(funs(M = mean, SE = se)) # alternative to add SE
    
  # --- Save data.frame ---
  write.table(temp, file = paste('~/Documents/Experiments/soc_ftask/merge_output/Averaged_ERPs/', 
                                 ID, 
                                 '_Ave.txt', sep = ''), 
              row.names = F, sep = '\t')
  # --- Remove data.frame after saving ---
  rm(temp)
  
  
  # ------ ***** 5) Function run-time ***** -------------------------
  time.taken <- Sys.time() - start.time
  
  
  # ------ ***** 6) Print summary  ***** ----------------------------
  print(paste(ID, 'has', EEG1, 'EEG Segments and', Pre1, 'behavioural trials;', 
              nrow(Latencies_ERN)/5, 'qualifiy for analysis.   Approx. function runtime = ', 
              as.character.POSIXt(round(time.taken, 2))))
  
}            


