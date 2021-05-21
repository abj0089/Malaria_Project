
rm(list = ls(all = T))

pacman::p_load(readxl, tidyverse , rjags, runjags, HDInterval)

setwd('/Users/ball4364/Desktop/Malaria/Data/Our data/')

source('../../Scripts/Functions/SMA_RD_models.R')

##############
# ADMISSIONS #

read_csv('complete_data.csv', guess_max = 10000) %>%
  mutate(Site = ifelse(Site == 'Junju (Kilifi B)', 'Kilifi B', Site),
         Site = ifelse(Site == "Kadziununi-Mauveni (Kilifi C)", 'Kilifi C', Site),
         Site = ifelse(Site == 'Ngerenya (Kilifi A)', 'Kilifi A', Site),
         Study_ID = paste(Site, Period, sep = ' '),
         Study_Number = as.numeric(as.factor(Study_ID))) %>%
  dplyr::select(Site, Study_ID, Study_Number, Period, 
                Hb_level, Deep_breathing, Unconscious_Observed, APVU, BCS,
                Inter_costal_recession, 
                Pallor, Blood_transfusion, 
                Total_Population) -> all_data

###########
# SURVYES #

load(file = 'formatted_surveys.RData')

formatted_surveys %>%
  group_by(Study_Number, Number, Positives) %>%
  summarise(mn = Min_Age - 1,
            mx = Max_Age - 1 , 
            PR = sum(Positives / Number) / n()) %>%
  ungroup() -> converted_pr_determ

source('../../Scripts/Smith_algorithm.R')

smith_et_al_deterministic <- convert(converted_pr_determ)

converted_pr_determ$PR_Smith <- smith_et_al_deterministic

#######
# SMA #
#######

all_data %>%
  mutate(SMA = Hb_level < 5) %>%
  group_by(Study_ID, Study_Number) %>%
  summarise(N = sum(SMA, na.rm = T),
            PY = unique(Total_Population)) %>%
  mutate(N = ifelse(is.na(N), 0, N)) -> case_numbers

all_data %>%
  mutate(SMA = as.numeric(Hb_level < 5)) %>%
  filter(!(is.na(Hb_level) & is.na(Pallor) & is.na(Blood_transfusion))) %>% 
  mutate(Prob_class = ifelse((is.na(Blood_transfusion) & is.na(Pallor)) | (!Blood_transfusion & Pallor != 'Severe')  | (is.na(Blood_transfusion) & Pallor != 'Severe')  | (!Blood_transfusion & is.na(Pallor)), 1,
                             # if Pallor is NA and BT is true 
                             ifelse((is.na(Pallor) & Blood_transfusion) | ((Pallor != 'Severe' & Blood_transfusion)), 2, 
                                    ifelse((is.na(Blood_transfusion) & Pallor == 'Severe') | (!Blood_transfusion & Pallor == 'Severe'), 3, 4)))) %>%
  dplyr::select(Study_Number, SMA, Prob_class) -> statuses_pre

statuses_pre %>%
  mutate(NA_stat = is.na(SMA)) %>%
  filter(NA_stat) %>%
  group_by(Study_Number, Prob_class) %>%
  summarise(N = n()) -> not_diagnosed

not_diagnosed %>%
  full_join(distinct(statuses_pre, Study_Number)) %>%
  ungroup() %>% 
  spread(Prob_class, N, fill = 0)%>%
  dplyr::select(-6)-> to_add

statuses_pre %>%
  filter(!is.na(SMA)) %>%
  group_by(Study_Number, Prob_class) %>%
  summarise(N_positives = sum(SMA), N_patients = n()) %>%
  ungroup() %>%
  arrange(Study_Number) -> statuses

S_dist <-  c( 0.03916203, 0.04464933, 0.0464907, 0.04703639,
              0.04468962, 0.04430289, 0.04088593, 0.04048751, 0.0385978,
              0.03586601, 0.03665516, 0.02850493, 0.03038655, 0.02493861,
              0.02169105, 0.01506840, 0.01481482, 0.01519726, 0.01474023,
              0.01553814, 0.01193805, 0.01220409, 0.01202119, 0.01202534,
              0.01283806, 0.01083728, 0.01095784, 0.01092874, 0.01067932,
              0.01137331, 0.01035001, 0.01035001, 0.01035001, 0.01035001,
              0.01035001, 0.008342426, 0.008342426, 0.008342426, 0.008342426,
              0.008342426, 0.00660059, 0.00660059, 0.00660059, 0.00660059,
              0.00660059, 0.004597861, 0.004597861, 0.004597861, 0.004597861,
              0.004597861, 0.003960774, 0.003960774, 0.003960774, 0.003960774,
              0.003960774, 0.003045687, 0.003045687, 0.003045687, 0.003045687,
              0.003045687, 0.002761077, 0.002761077, 0.002761077, 0.002761077,
              0.002761077, 0.001275367, 0.001275367, 0.001275367, 0.001275367,
              0.001275367, 0.001275367, 0.001275367, 0.001275367, 0.001275367,
              0.001275367, 0.001275367, 0.001275367, 0.001275367, 0.001275367,
              0.001275367, 0.001275367, 0.001275367, 0.001275367, 0.001275367,
              0.001275367) 

###########
# FOR GLM #
###########

vars <- c('PR_corrected', 'PR_corrected_mic', 'PR', 'P_prime', 'rho_SMA', 'rho_RD', 'intercept', 'beta1', 'beta2', 'beta3', 'estimated_admissions', 'r')

jags_data <- list(
  
  SITES = 1:max(formatted_surveys$Study_Number),
  probs = 1:4, 
  patient_index = statuses$Prob_class, 
  patient_status = statuses$N_positives,
  N_patients = statuses$N_patients, 
  PATIENTS = 1:nrow(statuses), 
  #site_indexes = as.matrix(dplyr::select(status_indexes, Min_index, Max_index)),
  actual_admissions = case_numbers$N, 
  N_NA = as.matrix(dplyr::select(to_add, 2:5)), 
  
  admissions = case_numbers$N,
  person_years = case_numbers$PY,
  
  Pf_positives = formatted_surveys$Positives,
  Total_surveyed = formatted_surveys$Number,
  survey_type = as.numeric(formatted_surveys$Detection_Method == 'RDT'), 
  age_range = as.matrix(dplyr::select(formatted_surveys, Min_Age, Max_Age)),
  ages = 1:85,
  S = S_dist,
  b = 1.807675, 
  c = 0.070376, 
  alpha = 9.422033, 
  s = 1 - 0.446326

)

N_samples <- 20000
thin <- 10
adapt <- 2500
burnin <- 10000
n.chains <- 6
  
################################
# POISSON VS NEGATIVE BINOMIAL #
################################

# Poisson 
#########

posteriors_GLM_logistic_poisson <-  run.jags(model = poisson_GLM_logisitic,
                                        sample = N_samples,
                                        data = jags_data,
                                        adapt = adapt,
                                        burnin = burnin,
                                        thin = thin,
                                        n.chains = n.chains,
                                        monitor = vars,
                                        method = "parallel")

DIC_GLM_logistic_poisson <- extract(posteriors_GLM_logistic_poisson, what = 'DIC')

# Negative binomial  
###################

posteriors_GLM_logistic_nb <-  run.jags(model = neg_bin_GLM_logisitic,
                                        sample = N_samples,
                                        data = jags_data,
                                        adapt = adapt,
                                        burnin = burnin,
                                        thin = thin,
                                        n.chains = n.chains,
                                        monitor = vars,
                                        method = "parallel")

DIC_GLM_logistic_nb <- extract(posteriors_GLM_logistic_nb, what = 'DIC')

#########
# TERMS #
#########

posteriors_GLM_linear_nb <-  run.jags(model = neg_bin_GLM_linear,
                                      sample = N_samples,
                                      data = jags_data,
                                      adapt = adapt,
                                      burnin = burnin,
                                      thin = thin,
                                      n.chains = n.chains,
                                      monitor = vars,
                                      method = "parallel")

DIC_GLM_linear_nb <- extract(posteriors_GLM_linear_nb, what = 'DIC')

posteriors_GLM_intercept_nb <-  run.jags(model = neg_bin_GLM_intercept,
                                         sample = N_samples,
                                         data = jags_data,
                                         adapt = adapt,
                                         burnin = burnin,
                                         thin = thin,
                                         n.chains = n.chains,
                                         monitor = vars,
                                         method = "parallel")

DIC_GLM_intercept_nb <- extract(posteriors_GLM_intercept_nb, what = 'DIC')

# save for later

all_fits <- list(Poisson_Logistic = posteriors_GLM_logistic_poisson, 
                 NB_Logistic = posteriors_GLM_logistic_nb,
                 NB_linear = posteriors_GLM_linear_nb,
                 NB_intercept = posteriors_GLM_intercept_nb)

save(all_fits, file = '../../Output/PfPR_model_fits_SMA.RData')
  
all_DIC <- list(Poisson_logisitc = DIC_GLM_logistic_poisson, 
                NB_logisitic = DIC_GLM_logistic_nb,
                 NB_linear = DIC_GLM_linear_nb,
                 NB_intercept = DIC_GLM_intercept_nb)

save(all_DIC, file = '../../Output/PfPR_model_DIC_SMA.RData')



