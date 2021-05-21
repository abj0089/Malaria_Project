
rm(list = ls(all = T))

pacman::p_load(readxl, tidyverse, patchwork, wesanderson, gratia, rjags, runjags, HDInterval, lme4)

setwd('/Users/ball4364/Desktop/Malaria/Data/Our data/')

source('../Functions/composite_models.R')

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

read_csv('complete_data.csv') %>%
  group_by(Site, Period) %>%
  distinct(Detection_Method) %>%
  mutate(Site = ifelse(Site == "Junju (Kilifi B)", 'Kilifi B', Site),
         Site = ifelse(Site == "Kadziununi-Mauveni (Kilifi C)", 'Kilifi C', Site),
         Site = ifelse(Site == "Ngerenya (Kilifi A)", 'Kilifi A', Site),
         Detection_Method = ifelse(Detection_Method == 'Microsopy', 'Microscopy', Detection_Method)) -> sampling_data

read_xlsx('Parasite_surveys.xlsx') %>%
  separate(Site, into = c('Site', 'Country'), sep = ',') %>%
  mutate(Country = gsub(pattern = ' ', replacement = '', Country),
         Study_ID = paste(Site, Period, sep = ' '),
         Study_Number = as.numeric(as.factor(Study_ID)),
         Min_Age  = round(Min_Age) + 1, 
         Max_Age = round(Max_Age) + 1,
         Max_Age = ifelse(Max_Age > 85, 85, Max_Age)) %>%
  left_join(sampling_data) %>%
  arrange(Study_Number) -> formatted_surveys 

save(formatted_surveys, file = '../Output/formatted_surveys.RData')

formatted_surveys %>%
  group_by(Study_Number, Number, Positives, ) %>%
  summarise(mn = Min_Age - 1,
            mx = Max_Age - 1 , 
            PR = sum(Positives / Number) / n()) %>%
  ungroup() -> converted_pr_determ

source('../Scripts/Smith_algorithm.R')

smith_et_al_deterministic <- convert(converted_pr_determ)

converted_pr_determ$PR_Smith <- smith_et_al_deterministic

#######
# SMA #
#######

all_data %>%
  mutate(CM = ifelse(BCS < 3 | APVU == 'U' | Unconscious_Observed, T, F), 
         CM = ifelse(is.na(CM), F, CM), 
         C = Hb_level < 5 | Deep_breathing | CM) %>%
  group_by(Study_ID, Study_Number) %>%
  summarise(N = sum(C, na.rm = T),
            PY = unique(Total_Population)) %>%
  mutate(N = ifelse(is.na(N), 0, N)) -> case_numbers

all_data %>%
  mutate(SMA = as.numeric(Hb_level < 5),
         RD = as.numeric(Deep_breathing),
         CM = ifelse(BCS < 3 | APVU == 'U' | Unconscious_Observed, T, F),
         CM = ifelse(is.na(CM), F, CM)) %>% 
  filter(!(is.na(SMA) & is.na(Pallor) & is.na(Blood_transfusion) & is.na(RD) & is.na(Inter_costal_recession) & is.na(CM))) %>% 
  mutate(Prob_class_1 = ifelse((is.na(Blood_transfusion) & is.na(Pallor)) | (!Blood_transfusion & Pallor != 'Severe')  | (is.na(Blood_transfusion) & Pallor != 'Severe')  | (!Blood_transfusion & is.na(Pallor)), 1,
                               ifelse((is.na(Pallor) & Blood_transfusion) | ((Pallor != 'Severe' & Blood_transfusion)), 2, 
                                      ifelse((is.na(Blood_transfusion) & Pallor == 'Severe') | (!Blood_transfusion & Pallor == 'Severe'), 3, 4))),
         Prob_class_2 = ifelse(is.na(Inter_costal_recession) | Inter_costal_recession == F, 1, 2)) %>%
  dplyr::select(Study_Number, SMA, RD, CM, Prob_class_1, Prob_class_2) %>%
  arrange(Study_Number) -> statuses

statuses %>%
  mutate(RN = row_number()) %>%
  group_by(Study_Number) %>%
  summarise(Min_index = min(RN), 
            Max_index = max(RN)) -> status_indexes

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
  probs_SMA = 1:4, 
  probs_RD = 1:2, 
  prob_index = as.matrix(dplyr::select(statuses, Prob_class_1, Prob_class_2)), 
  patient_status = as.matrix(dplyr::select(statuses, SMA, RD, CM)),
  PATIENTS = 1:nrow(statuses), 
  site_indexes = as.matrix(dplyr::select(status_indexes, Min_index, Max_index)),

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
thin <- 2
adapt <- 2500
burnin <- 10000
n.chains <- 6

# Poisson
#########

posteriors_GLM_logistic_poisson <-  run.jags(model = poisson_GLM_logistic,
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

posteriors_GLM_logistic_nb <-  run.jags(model = neg_bin_GLM_logistic,
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

save(all_fits, file = '../Output/PfPR_model_fits_Composite.RData')

all_DIC <- list(Poisson_logisitc = DIC_GLM_logistic_poisson,
                NB_logisitic = DIC_GLM_logistic_nb,
                NB_linear = DIC_GLM_linear_nb,
                NB_intercept = DIC_GLM_intercept_nb)

save(all_DIC, file = '../Output/PfPR_composite_DIC.RData')

###########
# GET FIT #
###########

do.call(rbind, all_fits$NB_linear$mcmc) %>%
  as.data.frame() -> posteriors

posteriors %>%  
  gather(Parameter, Value, 1:ncol(.)) %>% 
  group_by(Parameter) %>%
  summarise(Lower95 = hdi(Value)[1],
            Upper95 = hdi(Value)[2],
            Median = median(Value)) %>%
  dplyr::select(Parameter, Median, Lower95, Upper95) -> processed_posteriors

processed_posteriors %>%
  filter(str_detect(Parameter, pattern = 'estimated_admissions')) %>% 
  mutate(Study_Number = parse_number(Parameter)) %>%
  left_join(case_numbers) %>%
  mutate(Median_rate = (Median / PY) * 1000,
         Lower_rate = (Lower95 / PY) * 1000, 
         Upper_rate = (Upper95  / PY) * 1000) %>%
  dplyr::select(-Parameter, -Median, -Lower95, -Upper95) -> estimated_admissions

processed_posteriors %>%
  filter(str_detect(Parameter, pattern = 'PR_corrected') & !str_detect(Parameter, pattern = 'PR_corrected_mic')) %>% 
  mutate(Study_Number = parse_number(Parameter)) %>%
  left_join(estimated_admissions) -> corrected_PRs

processed_posteriors %>%
  filter(str_detect(Parameter, pattern = 'PR_corrected_mic')) %>% 
  mutate(Study_Number = parse_number(Parameter))  %>% 
  rename(Median_PR_mic = Median, Lower95_mic = Lower95, Upper95_mic = Upper95) %>%
  dplyr::select(-Parameter) %>% 
  left_join(corrected_PRs) %>% 
  left_join(converted_pr_determ) %>% 
  left_join(mutate(sampling_data, Study_ID = paste(Site, Period))) %>%
  ggplot(aes(x = Median_PR_mic, y = Median, shape = Detection_Method))+
  geom_point()
  

posteriors[, c('intercept', 'beta1')] %>% 
  as.matrix()-> chains_parms

PR_range <- seq(0, 0.75, length.out = 1000)

model_matrix <- cbind(1, PR_range)

fits  <- model_matrix %*% t(chains_parms)

medians <- apply(X = fits , FUN = median, MARGIN = 1)
intervals <- apply(X = fits, FUN = quantile, MARGIN = 1, c(0.01, 0.025, 0.975, 0.99))

fit_data_frame <- data.frame(PR = PR_range * 100, 
                             Median_fit = exp(medians + log(1000)), 
                             Upper_95 = exp(intervals[3, ] + log(1000)), 
                             Lower_95 = exp(intervals[2, ] + log(1000)),
                             Upper_99 = exp(intervals[4, ] + log(1000)),
                             Lower_99 = exp(intervals[1, ] + log(1000)),
                             Type = 'Composite rate of severe malaria phenotypes')

corrected_PRs %>%
  left_join(converted_pr_determ) -> plot_C

plot_C %>%
  left_join(mutate(sampling_data, Study_ID = paste(Site, Period))) %>%
  mutate(Type = 'Composite rate of severe malaria phenotypes',
         Median = Median * 100,
         Lower95 = Lower95 * 100,
         Upper95 = Upper95 * 100) %>%
  ggplot(aes(x = Median, y =  Median_rate, pch = Detection_Method))+
  annotate('segment', x = 0, y = 8, xend = 25, yend = 8, 
           arrow = arrow(type = 'closed', length = unit(0.15, "cm")), col = 'grey50')+
  geom_vline(aes(xintercept = 0), col = 'grey50', alpha = 0.5, inherit.aes = F)+
  geom_vline(aes(xintercept = 25), col = 'grey50', alpha = 0.5, inherit.aes = F)+
  geom_vline(aes(xintercept = 50), col = 'grey50', alpha = 0.5, inherit.aes = F)+
  geom_vline(aes(xintercept = 75), col = 'grey50', alpha = 0.5, inherit.aes = F)+
  annotate(geom = 'label', x = 0, y = 8.9, label = '2.06 [1.58-2.73] fold\nincrease in admissions', inherit.aes = F,
           hjust = -0.05, vjust = 1.2, family = 'Helvetica', col = 'grey50', size = 2.5, fill = 'white', label.padding = unit(0.1, "lines"),
           label.r = unit(0, "lines"),
           label.size = 0)+
  geom_linerange(aes(ymin = Lower_rate, ymax = Upper_rate), col = 'grey33', lineend = 'round', size = 0.33)+
  geom_linerange(aes(xmin = Lower95, xmax = Upper95), col = 'grey33', lineend = 'round', size = 0.33)+
  geom_point(fill = 'white',colour = 'grey33', stroke = 0.5, size = 1.25)+
  scale_shape_manual('', values = c(21, 19))+
  geom_ribbon(aes(x = PR, ymin = Lower_99, ymax = Upper_99),
              inherit.aes = F, data = fit_data_frame, fill = wes_palette("Darjeeling1")[5], alpha = 0.2)+
  geom_ribbon(aes(x = PR, ymin = Lower_95, ymax = Upper_95),
              inherit.aes = F, data = fit_data_frame, fill = wes_palette("Darjeeling1")[5], alpha = 0.4)+
  geom_line(aes(x = PR, y = Lower_99),
            inherit.aes = F, data = fit_data_frame, col = 'black', lty = 1, size = 0.1, alpha = 0.75)+
  geom_line(aes(x = PR, y = Upper_99),
            inherit.aes = F, data = fit_data_frame, col = 'black', lty = 1, size = 0.1, alpha = 0.75)+
  geom_line(aes(x = PR, y = Lower_95),
            inherit.aes = F, data = fit_data_frame, col = 'black', lty = 1, size = 0.1, alpha = 0.75)+
  geom_line(aes(x = PR, y = Upper_95),
            inherit.aes = F, data = fit_data_frame, col = 'black', lty = 1, size = 0.1, alpha = 0.75)+
  geom_line(aes(x = PR, y = Median_fit),
            inherit.aes = F, data = fit_data_frame , size = 0.66, col = 'black', lineend = 'round')+
  scale_x_continuous(breaks = seq(0, 70, 10))+
  scale_y_continuous(breaks = seq(0, 100, 1))+
  ylab('Annual admissions per 1000 children (3 months - 9 years)\n')+
  xlab(expression(atop('', italic(Pf)~PR['2-10'])))+
  theme(plot.title = element_text(hjust = 0),
        legend.position = 'none',
        legend.background = element_blank()) 

