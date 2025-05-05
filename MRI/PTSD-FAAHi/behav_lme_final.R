###############################################################################
# Title: behav_lme.R
# Author: Ryann Tansey
# Date: November 23, 2023
# Last updated: March 12, 2024
# 
# A script to do the group analysis on the fMRI emotional conflict task data 
#   for the PTSD/FAAHi study (PI: Dr. Leah Mayo). See: [publication here]
#
# For a detailed description of the task, see Etkin et al. (2006): "Resolving
#   emotional conflict: a role for the rostral anterior cingulate cortex in 
#   modulation activity in the amygdala."
###############################################################################

library(lme4)
library(lmerTest)
library(psych)
library(dplyr)
library(emmeans)
library(ggplot2)

# Read in the data from Matlab
alldata <- read.csv('/Users/ryann.tansey/Documents/projects/PTSD_FAAHi/emotion_conflict_task/behav_alldata.csv')

# Create variables to test for the interaction (previous trial, current trial)
alldata$prev_trial <- with(alldata, ifelse(grepl("con_1", category_str, fixed = TRUE) |
                                           grepl("high", category_str, fixed = TRUE),
                                           'incongruent', 'congruent'))

alldata$curr_trial <- with(alldata, ifelse(grepl("con_1", category_str, fixed = TRUE) |
                                           grepl("con_2", category_str, fixed = TRUE),
                                           'congruent', 'incongruent'))

alldata$emotion <- with(alldata, ifelse(grepl("happy", code, fixed = TRUE),
                                        'happy', 'fear'))

# Create a variable for sex
alldata$sex <- with(alldata, ifelse(grepl("K", subject, fixed = TRUE),
                                    'female', 'male'))

# Create a variable for sex of the face in the trial
alldata$sex_trial <- with(alldata, ifelse(grepl("female", code, fixed = TRUE),
                                          'female', 'male'))


hits <- alldata[which(alldata$stim_type == 'hit'),]


# Read in accuracy spreadsheets
# One with prev_trial and curr_trial (to test interaction)
accur_int <- read.csv('/Users/ryann.tansey/Documents/projects/PTSD_FAAHi/emotion_conflict_task/behav_accur_etkin.csv')
# One with emotion and curr_trial
accur_emo <- read.csv('/Users/ryann.tansey/Documents/projects/PTSD_FAAHi/emotion_conflict_task/behav_accur_emotion.csv')


# Run the LMEs

# REACTION TIME
# MODEL 1: PrevTrial*CurrTrial interaction in whole sample
# Testing the finding from the Etkin et al. (2006) paper
RT_lme1 <- lmer(time ~ prev_trial*curr_trial + sex + (1|subject), data = hits)
summary(RT_lme1)
# REML criterion at convergence: 176501.1
# Fixed effects:
#                                             Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                                  6481.31      86.36    86.98  75.052  < 2e-16 ***
# prev_trialincongruent                          98.29      30.73 10346.04   3.199  0.00138 ** 
# curr_trialincongruent                         524.87      31.21 10346.38  16.815  < 2e-16 ***
# sexmale                                        36.81     218.92    78.84   0.168  0.86688    
# prev_trialincongruent:curr_trialincongruent   -59.70      44.33 10346.19  -1.347  0.17809    



# MODEL 2: Emotion*CurrTrial interaction in whole sample
RT_lme2 <- lmer(time ~ emotion*curr_trial + sex + (1|subject), data = hits)
summary(RT_lme2)
# REML criterion at convergence: 176493.3
# Fixed effects:
#                                     Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                         6584.30      86.22    86.57  76.370  < 2e-16 ***
# emotionhappy                        -106.24      30.72 10346.42  -3.459 0.000545 ***
# curr_trialincongruent                483.64      31.20 10346.47  15.502  < 2e-16 ***
# sexmale                               37.44     218.82    78.84   0.171 0.864573    
# emotionhappy:curr_trialincongruent    18.73      44.29 10346.25   0.423 0.672438 

#--------------------------------------------------------------------------------------------

# MODEL 3: Grp + PrevTrial*CurrTrial
RT_lme3 <- lmer(time ~ group + prev_trial*curr_trial + sex + (1|subject), data = hits)
summary(RT_lme3)
# REML criterion at convergence: 176487.6
#Fixed effects:
#                                             Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                                  6379.55     118.56    81.83  53.807  < 2e-16 ***
# groupplacebo                                  195.11     156.37    77.75   1.248  0.21587    
# prev_trialincongruent                          98.29      30.73 10346.03   3.199  0.00138 ** 
# curr_trialincongruent                         524.88      31.21 10346.38  16.816  < 2e-16 ***
# sexmale                                        73.48     220.12    77.82   0.334  0.73943    
# prev_trialincongruent:curr_trialincongruent   -59.71      44.33 10346.18  -1.347  0.17805


# MODEL 4: Grp + Emotion*CurrTrial
RT_lme4 <- lmer(time ~ group + emotion*curr_trial + sex + (1|subject), data = hits)
summary(RT_lme4)
# REML criterion at convergence: 176479.8
# Fixed effects:
#                                     Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                         6482.44     118.43    81.62  54.737  < 2e-16 ***
# groupplacebo                         195.33     156.29    77.75   1.250 0.215115    
# emotionhappy                        -106.26      30.72 10346.41  -3.459 0.000544 ***
# curr_trialincongruent                483.64      31.20 10346.47  15.502  < 2e-16 ***
# sexmale                               74.15     220.01    77.82   0.337 0.737012    
# emotionhappy:curr_trialincongruent    18.74      44.29 10346.25   0.423 0.672178 

#--------------------------------------------------------------------------------------------

# MODEL 5: AEA + PrevTrial*CurrTrial
RT_lme5 <- lmer(time ~ AEA + prev_trial*curr_trial + sex + (1|subject), data = hits)
summary(RT_lme5)
# REML criterion at convergence: 174152.1
# Fixed effects:
#                                             Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                                  6533.13     125.05    80.52  52.243  < 2e-16 ***
# AEA                                           -15.34      26.37    76.77  -0.582  0.56245    
# prev_trialincongruent                          96.02      30.86 10211.04   3.111  0.00187 ** 
# curr_trialincongruent                         530.62      31.36 10211.37  16.921  < 2e-16 ***
# sexmale                                        66.11     227.39    76.82   0.291  0.77205    
# prev_trialincongruent:curr_trialincongruent   -64.40      44.54 10211.18  -1.446  0.14822    


# MODEL 6: AEA + Emotion*CurrTrial
RT_lme6 <- lmer(time ~ AEA + emotion*curr_trial + sex + (1|subject), data = hits)
summary(RT_lme6)
# REML criterion at convergence: 174143.7
# Fixed effects:
#                                     Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                         6633.92     124.92    80.35  53.103  < 2e-16 ***
# AEA                                  -15.40      26.36    76.77  -0.584 0.560788    
# emotionhappy                        -103.67      30.85 10211.39  -3.360 0.000781 ***
# curr_trialincongruent                489.52      31.35 10211.47  15.615  < 2e-16 ***
# sexmale                               66.82     227.28    76.82   0.294 0.769536    
# emotionhappy:curr_trialincongruent    13.97      44.49 10211.24   0.314 0.753565  

#--------------------------------------------------------------------------------------------

# MODEL 7: PCL + PrevTrial*CurrTrial
RT_lme7 <- lmer(time ~ PCL + prev_trial*curr_trial + sex + (1|subject), data = hits)
summary(RT_lme7)
# REML criterion at convergence: 176495.8
# Fixed effects:
#                                             Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)                                  6525.097    228.726    78.875  28.528  < 2e-16 ***
# PCL                                            -1.143      5.524    77.802  -0.207  0.83662    
# prev_trialincongruent                          98.294     30.728 10346.038   3.199  0.00138 ** 
# curr_trialincongruent                         524.875     31.214 10346.373  16.815  < 2e-16 ***
# sexmale                                        41.526    221.431    77.845   0.188  0.85173    
# prev_trialincongruent:curr_trialincongruent   -59.697     44.329 10346.181  -1.347  0.17812

# MODEL 8: PCL + Emotion*CurrTrial
RT_lme8 <- lmer(time ~ PCL + emotion*curr_trial + sex + (1|subject), data = hits)
summary(RT_lme8)
# REML criterion at convergence: 176488
# Fixed effects:
#                                     Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)                         6627.573    228.587    78.834  28.994  < 2e-16 ***
# PCL                                   -1.129      5.521    77.803  -0.205 0.838457    
# emotionhappy                        -106.244     30.717 10346.411  -3.459 0.000545 ***
# curr_trialincongruent                483.648     31.198 10346.461  15.502  < 2e-16 ***
# sexmale                               42.098    221.327    77.846   0.190 0.849642    
# emotionhappy:curr_trialincongruent    18.733     44.288 10346.237   0.423 0.672312 

#--------------------------------------------------------------------------------------------

# MODEL 9: CAPS_wk0 + CAPS_wk12 + PrevTrial*CurrTrial
RT_lme9 <- lmer(time ~ CAPS_wk0 + CAPS_wk12 + prev_trial*curr_trial + sex + (1|subject), data = hits)
summary(RT_lme9)
# REML criterion at convergence: 159612
# Fixed effects:
#                                             Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                                 6147.239    401.252   69.131  15.320  < 2e-16 ***
# CAPS_wk0                                      13.611     13.549   68.773   1.005  0.31859    
# CAPS_wk12                                     -5.337      8.076   68.779  -0.661  0.51091    
# prev_trialincongruent                         93.106     32.189 9362.003   2.893  0.00383 ** 
# curr_trialincongruent                        530.375     32.680 9362.328  16.229  < 2e-16 ***
# sexmale                                       50.723    223.159   68.867   0.227  0.82087    
# prev_trialincongruent:curr_trialincongruent  -45.925     46.395 9362.141  -0.990  0.32227    

# MODEL 10: CAPS_wk0 + CAPS_wk12 + Emotion*CurrTrial
RT_lme10 <- lmer(time ~ CAPS_wk0 + CAPS_wk12 + emotion*curr_trial + sex + (1|subject), data = hits)
summary(RT_lme10)
# REML criterion at convergence: 159590.5
# Fixed effects:
#                                    Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                        6267.501    401.046   69.119  15.628  < 2e-16 ***
# CAPS_wk0                             13.609     13.542   68.774   1.005    0.318    
# CAPS_wk12                            -5.322      8.072   68.780  -0.659    0.512    
# emotionhappy                       -146.469     32.153 9362.399  -4.555 5.29e-06 ***
# curr_trialincongruent               487.725     32.658 9362.455  14.934  < 2e-16 ***
# sexmale                              51.265    223.055   68.868   0.230    0.819    
# emotionhappy:curr_trialincongruent   34.977     46.316 9362.243   0.755    0.450      


#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

# ACCURACY

# MODEL 1: PrevTrial*CurrTrial interaction in the whole sample
# Testing the finding from the Etkin et al. (2006) paper
acc_lme1 <- lmer(accuracy ~ prev_trial*curr_trial + sex + (1|subject), data = accur_int)
summary(acc_lme1)
# REML criterion at convergence: -716.7
# Fixed effects:
#                                               Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                   0.932515   0.012222 121.697409  76.297  < 2e-16 ***
# prev_trialincongruent                        -0.024998   0.008987 240.000001  -2.782 0.005839 ** 
# curr_trialincongruent                        -0.099163   0.008987 240.000001 -11.034  < 2e-16 ***
# sexmale                                      -0.018283   0.028353  78.999997  -0.645 0.520905    
# prev_trialincongruent:curr_trialincongruent   0.043655   0.012709 240.000001   3.435 0.000698 ***

interaction.plot(x.factor = accur_int$prev_trial,
                 trace.factor = accur_int$curr_trial,
                 response = accur_int$accuracy,
                 trace.label = 'Current trial',
                 xlab = 'Previous trial',
                 ylab = 'Accuracy',
                 xtick = TRUE,
                 main = 'Interaction plot for accuracy')

pairs(emmeans(acc_lme1, ~ prev_trial*curr_trial), adjust = 'fdr')
# contrast                                        estimate      SE  df t.ratio p.value
# congruent congruent - incongruent congruent       0.0250 0.00899 240   2.782  0.0070
# congruent congruent - congruent incongruent       0.0992 0.00899 240  11.034  <.0001
# congruent congruent - incongruent incongruent     0.0805 0.00899 240   8.958  <.0001
# incongruent congruent - congruent incongruent     0.0742 0.00899 240   8.253  <.0001
# incongruent congruent - incongruent incongruent   0.0555 0.00899 240   6.177  <.0001
# congruent incongruent - incongruent incongruent  -0.0187 0.00899 240  -2.076  0.0390



# MODEL 2: Emotion*CurrTrial in whole sample
acc_lme2 <- lmer(accuracy ~ emotion*curr_trial + sex + (1|subject), data = accur_emo)
summary(acc_lme2)
# REML criterion at convergence: -657.9
# Fixed effects:
# Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                          0.916550   0.012570 133.875778  72.918  < 2e-16 ***
# emotionhappy                         0.005839   0.010155 240.000000   0.575    0.566    
# curr_trialincongruent               -0.066733   0.010155 240.000000  -6.571 3.07e-10 ***
# sexmale                             -0.017796   0.028380  79.000000  -0.627    0.532    
# emotionhappy:curr_trialincongruent  -0.021216   0.014361 240.000000  -1.477    0.141 

#--------------------------------------------------------------------------------------------

# MODEL 3: Grp + PrevTrial*CurrTrial
acc_lme3 <- lmer(accuracy ~ group + prev_trial*curr_trial + sex + (1|subject), data = accur_int)
summary(acc_lme3)
# REML criterion at convergence: -712.7
# Fixed effects:
#                                               Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                   0.947319   0.016095  99.451816  58.856  < 2e-16 ***
# groupplacebo                                 -0.028374   0.020205  77.999998  -1.404 0.164197    
# prev_trialincongruent                        -0.024998   0.008987 240.000001  -2.782 0.005839 ** 
# curr_trialincongruent                        -0.099163   0.008987 240.000001 -11.034  < 2e-16 ***
# sexmale                                      -0.023629   0.028436  77.999998  -0.831 0.408542    
# prev_trialincongruent:curr_trialincongruent   0.043655   0.012709 240.000001   3.435 0.000698 ***


# MODEL 4: Grp + Emotion*CurrTrial
acc_lme4 <- lmer(accuracy ~ group + emotion*curr_trial + sex + (1|subject), data = accur_emo)
summary(acc_lme4)
# REML criterion at convergence: -653.9
# Fixed effects:
#                                     Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                          0.931353   0.016368 105.561320  56.902  < 2e-16 ***
# groupplacebo                        -0.028371   0.020225  78.000000  -1.403    0.165    
# emotionhappy                         0.005839   0.010155 240.000000   0.575    0.566    
# curr_trialincongruent               -0.066733   0.010155 240.000000  -6.571 3.07e-10 ***
# sexmale                             -0.023141   0.028464  78.000000  -0.813    0.419    
# emotionhappy:curr_trialincongruent  -0.021216   0.014361 240.000000  -1.477    0.141 

#--------------------------------------------------------------------------------------------

# MODEL 5: AEA + PrevTrial*CurrTrial
acc_lme5 <- lmer(accuracy ~ AEA + prev_trial*curr_trial + sex + (1|subject), data = accur_int)
summary(acc_lme5)
# REML criterion at convergence: -696.1
# Fixed effects:
#                                               Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                   0.924561   0.016909  96.306800  54.679  < 2e-16 ***
# AEA                                           0.002290   0.003408  76.999996   0.672 0.503601    
# prev_trialincongruent                        -0.025038   0.009072 237.000001  -2.760 0.006237 ** 
# curr_trialincongruent                        -0.099835   0.009072 237.000001 -11.004  < 2e-16 ***
# sexmale                                      -0.022091   0.029383  76.999997  -0.752 0.454435    
# prev_trialincongruent:curr_trialincongruent   0.043360   0.012830 237.000001   3.379 0.000849 ***

# MODEL 6: AEA + Emotion*CurrTrial
acc_lme6 <- lmer(accuracy ~ AEA + emotion*curr_trial + sex + (1|subject), data = accur_emo)
summary(acc_lme6)
# REML criterion at convergence: -640.2
# Fixed effects:
#                                       Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                          0.907782   0.017166 101.557073  52.884  < 2e-16 ***
# AEA                                  0.002261   0.003412  76.999997   0.663    0.509    
# emotionhappy                         0.007648   0.010205 237.000001   0.749    0.454    
# curr_trialincongruent               -0.066554   0.010205 237.000001  -6.522 4.16e-10 ***
# sexmale                             -0.021559   0.029415  76.999997  -0.733    0.466    
# emotionhappy:curr_trialincongruent  -0.023217   0.014432 237.000001  -1.609    0.109

#--------------------------------------------------------------------------------------------

# MODEL 7: PCL + PrevTrial*CurrTrial
acc_lme7 <- lmer(accuracy ~ PCL + prev_trial*curr_trial + sex + (1|subject), data = accur_int)
summary(acc_lme7)
# REML criterion at convergence: -704
# Fixed effects:
#                                               Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                  9.296e-01  3.004e-02  8.348e+01  30.946  < 2e-16 ***
# PCL                                          7.576e-05  7.157e-04  7.800e+01   0.106 0.915962    
# prev_trialincongruent                       -2.500e-02  8.987e-03  2.400e+02  -2.782 0.005839 ** 
# curr_trialincongruent                       -9.916e-02  8.987e-03  2.400e+02 -11.034  < 2e-16 ***
# sexmale                                     -1.859e-02  2.868e-02  7.800e+01  -0.648 0.518719    
# prev_trialincongruent:curr_trialincongruent  4.365e-02  1.271e-02  2.400e+02   3.435 0.000698 ***


# MODEL 8: PCL + Emotion*CurrTrial
acc_lme8 <- lmer(accuracy ~ PCL + emotion*curr_trial + sex + (1|subject), data = accur_emo)
summary(acc_lme8)
# REML criterion at convergence: -645.2
# Fixed effects:
#                                     Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                         9.140e-01  3.021e-02  8.500e+01  30.257  < 2e-16 ***
# PCL                                 6.684e-05  7.164e-04  7.800e+01   0.093    0.926    
# emotionhappy                        5.839e-03  1.016e-02  2.400e+02   0.575    0.566    
# curr_trialincongruent              -6.673e-02  1.016e-02  2.400e+02  -6.571 3.07e-10 ***
# sexmale                            -1.807e-02  2.871e-02  7.800e+01  -0.629    0.531    
# emotionhappy:curr_trialincongruent -2.122e-02  1.436e-02  2.400e+02  -1.477    0.141 

#--------------------------------------------------------------------------------------------

# MODEL 9: CAPS_wk0 + CAPS_wk12 + PrevTrial*CurrTrial
acc_lme9 <- lmer(accuracy ~ CAPS_wk0 + CAPS_wk12 + prev_trial*curr_trial + sex + (1|subject), data = accur_int)
summary(acc_lme9)
# REML criterion at convergence: -624.1
# Fixed effects:
#                                               Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                  9.617e-01  4.984e-02  7.093e+01  19.296  < 2e-16 ***
# CAPS_wk0                                    -5.207e-04  1.674e-03  6.900e+01  -0.311  0.75662    
# CAPS_wk12                                   -4.708e-04  9.975e-04  6.900e+01  -0.472  0.63846    
# prev_trialincongruent                       -2.229e-02  9.546e-03  2.160e+02  -2.335  0.02047 *  
# curr_trialincongruent                       -9.476e-02  9.546e-03  2.160e+02  -9.926  < 2e-16 ***
# sexmale                                     -2.379e-02  2.755e-02  6.900e+01  -0.863  0.39088    
# prev_trialincongruent:curr_trialincongruent  4.024e-02  1.350e-02  2.160e+02   2.981  0.00321 ** 


# MODEL 10: CAPS_wk0 + CAPS_wk12 + Emotion*CurrTrial
acc_lme10 <- lmer(accuracy ~ CAPS_wk0 + CAPS_wk12 + emotion*curr_trial + sex + (1|subject), data = accur_emo)
summary(acc_lme10)
# REML criterion at convergence: -571.6
# Fixed effects:
#                                      Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                         9.459e-01  4.992e-02  7.147e+01  18.948  < 2e-16 ***
# CAPS_wk0                           -5.459e-04  1.673e-03  6.900e+01  -0.326    0.745    
# CAPS_wk12                          -4.714e-04  9.972e-04  6.900e+01  -0.473    0.638    
# emotionhappy                        1.017e-02  1.078e-02  2.160e+02   0.944    0.346    
# curr_trialincongruent              -6.294e-02  1.078e-02  2.160e+02  -5.839 1.91e-08 ***
# sexmale                            -2.337e-02  2.755e-02  6.900e+01  -0.848    0.399    
# emotionhappy:curr_trialincongruent -2.353e-02  1.524e-02  2.160e+02  -1.544    0.124  
