
rm(list = ls(all = T))

pacman::p_load(readxl, tidyverse, patchwork, wesanderson, gratia, rjags, runjags, HDInterval, lme4)

setwd('/Users/ball4364/Desktop/Malaria/Data/Our data/')

read_csv('complete_data.csv', guess_max = 10000) %>%
  mutate(Site = ifelse(Site == 'Junju (Kilifi B)', 'Kilifi B', Site),
         Site = ifelse(Site == "Kadziununi-Mauveni (Kilifi C)", 'Kilifi C', Site),
         Site = ifelse(Site == 'Ngerenya (Kilifi A)', 'Kilifi A', Site),
         Study_ID = paste(Site, Period, sep = ' '),
         Study_Number = as.numeric(as.factor(Study_ID))) %>%
  dplyr::select(Site, Study_ID, Study_Number, Period, Age_years,
                Hb_level, Deep_breathing, Unconscious_Observed, APVU, BCS,
                Inter_costal_recession, 
                Pallor, Blood_transfusion, 
                N_Age) -> all_data

read_csv('location_data.csv') %>%
  mutate(Site = ifelse(Site == 'Junju (Kilifi B)', 'Kilifi B', Site),
         Site = ifelse(Site == "Kadziununi-Mauveni (Kilifi C)", 'Kilifi C', Site),
         Site = ifelse(Site == 'Ngerenya (Kilifi A)', 'Kilifi A', Site),
         Study_ID = paste(Site, Period, sep = ' '),
         Study_Number = as.numeric(as.factor(Study_ID))) %>%
  arrange(Study_Number, Age_years) -> location_data

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

formatted_surveys %>%
  group_by(Study_Number, Number, Positives) %>%
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
  mutate(Study_ID = paste(Site, Period, sep = ' '),
         Study_Number = as.numeric(as.factor(Study_ID))) %>% 
  left_join(converted_pr_determ) %>% 
  mutate(CM = ifelse(BCS < 3 | APVU == 'U' | Unconscious_Observed, T, F), 
         CM = ifelse(is.na(CM), F, CM), 
         C = Hb_level < 5 | Deep_breathing | CM,
         C_alt = Hb_level < 5 | Pallor == 'Severe' | Blood_transfusion | Deep_breathing | Inter_costal_recession | CM) -> with_PR

with_PR %>%
  group_by(Age_years, PR_discrete) %>%
  summarise(N = sum(C, na.rm = T),
            N_alt = sum(C_alt, na.rm = T),
            PY = sum(unique(N_Age)), 
            Rate = (N / PY) * 1000,
            Rate_alt = (N_alt / PY) * 1000) %>%
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
         Prob_class_2 = ifelse(is.na(Inter_costal_recession) | Inter_costal_recession == F, 1, 2),
         Age_years = Age_years) %>%
  dplyr::select(Study_Number, Age_years, N_Age, SMA, RD, CM, Prob_class_1, Prob_class_2) %>%
  arrange(Study_Number, Age_years) -> statuses

statuses %>%
  mutate(Composite =  SMA | RD | CM,
         RN = row_number()) %>%
  group_by(Study_Number, Age_years, N_Age) %>% 
  summarise(Composite = sum(Composite, na.rm = T), 
            Min = min(RN), 
            Max = max(RN)) %>%
  ungroup() %>% 
  mutate(Age_years = ifelse(Age_years == 0, 0, Age_years - 3/12),
         Age_plus = ifelse(Age_years == 0, 1 - 3/12, Age_years + 1),
         Age_delta = ifelse(Age_years == 0, 1 - 3/12, 1)) -> actual_admissions

actual_admissions %>%
  mutate(RN = 1:n()) %>%
  group_by(Study_Number) %>%
  summarise(Min = min(RN), 
            Max = max(RN)) -> mins_max

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

jags_data <- list(
  
  SITE_AGES = 1:nrow(actual_admissions),
  actual_admissions = actual_admissions$Composite,
  Min = actual_admissions$Min,
  Max = actual_admissions$Max,
  site_indexes = actual_admissions$Study_Number,
  actual_ages = actual_admissions$Age_years,
  PY = actual_admissions$N_Age,
  lead_age = actual_admissions$Age_plus,
  delta_age = actual_admissions$Age_delta,
  Min2 = mins_max$Min,
  Max2 = mins_max$Max,
  
  SITES = 1:max(formatted_surveys$Study_Number),
  probs_SMA = 1:4, 
  probs_RD = 1:2, 
  prob_index = as.matrix(dplyr::select(statuses, Prob_class_1, Prob_class_2)), 
  patient_status = as.matrix(dplyr::select(statuses, SMA, RD, CM)),
  PATIENTS = 1:nrow(statuses), 
  
  survey_type = as.numeric(formatted_surveys$Detection_Method == 'RDT'), 
  Pf_positives = formatted_surveys$Positives,
  Total_surveyed = formatted_surveys$Number,
  age_range = as.matrix(dplyr::select(formatted_surveys, Min_Age, Max_Age)),
  ages = 1:85,
  S = S_dist,
  b = 1.807675, 
  c = 0.070376, 
  alpha = 9.422033, 
  s = 1 - 0.446326
  
)

N_samples <- 25000
thin <- 4
adapt <- 2500
burnin <- 25000
n.chains <- 4


age_model <-' model{
  
  ######################
  # X-AXIS UNCERTAINTY #
  ######################
  
  for(age in ages){
      
    P[age] <- 1 - exp(- (b * age))
      
    F[age] <- 1 - (age >= alpha) * s * (1 - exp(- c * (age - alpha)))
      
    PFS[age] <- P[age] * F[age] * S[age]
    
  }
    
  for(survey in SITES){

    PR[survey] <- P_prime[survey] * sum((PFS[age_range[survey, 1]:age_range[survey, 2]])) / sum(S[age_range[survey, 1]:age_range[survey, 2]])
    
    Pf_positives[survey] ~ dbinom(PR[survey], Total_surveyed[survey])

    PR_corrected_mic[survey] <- P_prime[survey] * (sum(PFS[3:10]) / sum(S[3:10]))
    
    PR_corrected_RDT[survey] <- phi(probit(PR_corrected_mic[survey]) * slope_PR + intercept_PR)

    PR_corrected[survey] <- (1 - survey_type[survey]) * PR_corrected_mic[survey] + survey_type[survey] * PR_corrected_RDT[survey]
    
  }
  
  ######################
  # Y-AXIS UNCERTAINTY #
  ######################
  
  # PREDICT PATIENT STATUS
  
  for(patient in PATIENTS){
    
    patient_status[patient, 1] ~ dbern(rho_SMA[prob_index[patient, 1]])
    patient_status[patient, 2] ~ dbern(rho_RD[prob_index[patient, 2]])
    
    all_status[patient] <- sum(patient_status[patient, 1:3])
    
    composite_status[patient] <- all_status[patient] > 0
    
  }
  
  # CUMULATIVE SUMS
  
  for(site_age in SITE_AGES){

      # Rate corrections
      
      estimated_admissions[site_age] <- sum(composite_status[Min[site_age]:Max[site_age]])
      
      correction[site_age] <- lambda[site_age] / (estimated_admissions[site_age] * (estimated_admissions[site_age] > 0) + lambda[site_age] * (estimated_admissions[site_age] == 0))
      
      function[site_age] <- ((ra[site_indexes[site_age]] ^ sh[site_indexes[site_age]]) / exp(loggam(sh[site_indexes[site_age]]))) * (actual_ages[site_age] ^ (sh[site_indexes[site_age]] - 1)) * exp(-ra[site_indexes[site_age]] * actual_ages[site_age])
      function_lead[site_age] <- ((ra[site_indexes[site_age]] ^ sh[site_indexes[site_age]]) / exp(loggam(sh[site_indexes[site_age]]))) * (lead_age[site_age] ^ (sh[site_indexes[site_age]] - 1)) * exp(-ra[site_indexes[site_age]] * lead_age[site_age])
        
      integrate[site_age] <- (function[site_age] + (function_lead[site_age] - function[site_age])) * delta_age[site_age]
      denominator[site_age] <- sum(integrate[Min2[site_indexes[site_age]]:Max2[site_indexes[site_age]]])
      integrator[site_age] <- (PR_corrected[site_indexes[site_age]] <= PR_crit) * 0.1 + (PR_corrected[site_indexes[site_age]] > PR_crit) * (integrate[site_age] / denominator[site_age])
      
      lambda[site_age] <- exp(log(PY[site_age]) + alpha_parm[1] + betas[1] * PR_corrected[site_indexes[site_age]] + RE[site_indexes[site_age]]) * integrator[site_age]
      
      actual_admissions[site_age] ~ dpois(lambda[site_age] * correction[site_age])

  }
  
  ##########
  # PRIORS #
  ##########
  
  sigma ~ dunif(0, 10)
  
  PR_crit ~ dunif(0, 1)
  
  sd ~ dunif(0, 10)
  
  # Admissions model 
  
  for(p in probs_SMA){
    
    rho_SMA[p] ~ dbeta(1, 1)
    
  }
  
  for(p in probs_RD){
    
    rho_RD[p] ~ dbeta(1, 1)
    
  }
  
  for(s in SITES){
  
      P_prime[s] ~ dbeta(1, 1)

      ######### REs ##########
      
      RE[s] ~ dnorm(0, 1 / (sd ^ 2))
        
      m[s] <- exp(alpha_parm[2] + betas[2] * PR_corrected[s])

      ra[s] <- (m[s] + sqrt(m[s] ^ 2 + 4 * sigma ^ 2)) / (2 * sigma ^ 2)
      sh[s] <- 1 + m[s] * ra[s]

  }
  
  for(i in 1:2){
  
      betas[i] ~ dnorm(0, 1 / 10)

  }
  
  for(i in 1:2){

    alpha_parm[i] ~ dnorm(0, 1 / 10)

  }
  
  slope_PR ~ dnorm(0.97, 1 / (0.01 ^ 2))
  intercept_PR ~ dnorm(-0.22, 1 / (0.005 ^ 2))
  
}'

age_mod <-  run.jags(model = age_model,
                     sample = N_samples,
                     data = jags_data,
                     adapt = adapt,
                     burnin = burnin,
                     thin = thin,
                     n.chains = n.chains,
                     monitor = c('betas', 'alpha_parm', 'r', 'sd', 'sigma',
                                 'RE', 'estimated_admissions', 'PR_crit',
                                 'rho_SMA', 'rho_RD', 'PR_corrected', 'sigma'),
                     method = "parallel")

save(age_mod, file = '../../Output/PfPR_gamma_age_model_2.RData')
load(file = '../../Output/PfPR_gamma_age_model_2.RData')

add.summary(age_mod, vars = c('betas', 'alpha_parm', 'sd', 'sigma', 'PR_crit'))

plot(age_mod, vars = c('betas', 'alpha_parm'))

do.call(rbind, age_mod$mcmc) %>%
  as.data.frame()  -> posteriors

posteriors[, str_detect(names(posteriors), 'beta') | str_detect(names(posteriors), 'alpha_parm') | str_detect(names(posteriors), 'sigma') | str_detect(names(posteriors), 'PR_crit')] %>% 
  sample_n(10000) %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(Iteration = 1:n()) -> chains_parms

names(chains_parms) <- gsub(pattern = '[', replacement = '', x = names(chains_parms), fixed = T)
names(chains_parms) <- gsub(pattern = ']', replacement = '', x = names(chains_parms), fixed = T)

preds <- expand.grid(Age =  seq((min(actual_admissions$Age_years)), (10), by = 1/12), 
                     PR =  seq(0, 0.75, by = 0.01), 
                     Iteration = 1:max(chains_parms$Iteration))

preds %>%
  left_join(chains_parms) -> preds_combined

log_logistic <- function(PR, alpha2, beta2, sigma, Age){
  
  m <- exp(alpha2 + beta2 * PR)
  
  ra <- ( m + sqrt( m^2 + 4*sigma^2 ) ) / ( 2 * sigma^2 )
  sh <- 1 + m * ra
  
  dgamma(Age, shape = sh, rate = ra)
  
}

get_preds <- function(PR, alpha1, beta1, RE1 = 0){
  
  func <- exp(log(1000) + alpha1 + beta1 * PR + RE1)
  
  return(func)
  
}

memory.limit(size = 100000)

preds_combined %>%
  group_by(PR, Iteration) %>%
  mutate(Age_plus = lead(Age)) %>%  
  filter(!is.na(Age_plus)) %>%
  ungroup() %>% 
  mutate(LLD = log_logistic(PR, alpha2 = alpha_parm2, beta2 = betas2, sigma = sigma, Age),
         LLDP = log_logistic(PR, alpha2 = alpha_parm2, beta2 = betas2, sigma = sigma, Age_plus),
         Integrate = (LLD + (LLDP - LLD)) * (Age_plus - Age), 
         Integrate = ifelse(PR <= PR_crit, 1 / length(unique(preds$Age)), Integrate),
         Fit = get_preds(PR, alpha1 = alpha_parm1, beta1 = betas1)) -> with_fits

preds_combined %>%
  mutate(Age_Fit = (alpha_parm2 + betas2 * PR)) %>%
  group_by(PR) %>%
  summarise(Age = exp(median(Age_Fit)),
            Age_Up = exp(quantile(Age_Fit, probs = c(0.025, 0.975))[2]),
            Age_low = exp(quantile(Age_Fit, probs = c(0.025, 0.975))[1])) %>%
  mutate(PR = PR * 100)-> summarised_median

with_fits %>%
  ungroup() %>%
  filter(!is.na(Fit)) %>%
  group_by(PR, Iteration) %>%
  mutate(Denominator = sum(Integrate),
         Integrate_2 = Integrate / Denominator) %>%
  group_by(PR, Age) %>%
  summarise(Median_fit = (median(Integrate_2 * Fit, na.rm = T)), 
            Lower_fit = (quantile(Integrate_2 * Fit, probs = c(0.025, 0.975), na.rm = T)[2]),
            Upper_fit = (quantile(Integrate_2 * Fit, probs = c(0.025, 0.975), na.rm = T)[1]),
            Median_PR_cut = (median(PR_crit, na.rm = T)) * 100,
            Lower_PR_cut = (quantile(PR_crit, probs = c(0.025, 0.975), na.rm = T)[2]) * 100,
            Upper_PR_cut = (quantile(PR_crit, probs = c(0.025, 0.975), na.rm = T)[1]) * 100,
            IQR = IQR(Integrate_2 * Fit)) %>%
  mutate(PR = PR * 100)-> summarised

summarised %>%
  ggplot(aes(x = Age, y = PR, fill = (Median_fit), z = (Median_fit)))+
  geom_raster(interpolate = T)+
  geom_contour(col = 'white', lineend = 'round', lty = 5, binwidth = 0.15, size = 0.2)+
  geom_ribbon(aes(xmin = Age_low , xmax = Age_Up , y = PR), alpha = 0.25,
              fill = 'white', data = filter(summarised_median, PR >= unique(summarised$Median_PR_cut) * 0.993), 
              inherit.aes = F, lty = 3, size = 0.5)+
  geom_path(aes(x = Age, y = PR),
            col = 'white', arrow = arrow(ends = "first", type = "closed", length = unit(0.1, "inches")),
            data = filter(summarised_median, PR >= unique(summarised$Median_PR_cut) * 0.993), inherit.aes = F, lineend = 'round')+
  annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = unique(summarised$Median_PR_cut), 
           fill = 'grey75', col = 'transparent', alpha = 0.65, inhearit.aes = F) +
  annotate(geom = 'text', x = 7.5, y = unique(summarised$Median_PR_cut),
           label = paste0('PR\' = ', round(unique(summarised$Median_PR_cut), 2),
                          '% [', round(unique(summarised$Upper_PR_cut), 2), '-',
                          round(unique(summarised$Lower_PR_cut), 2), ']'),
           family = 'Helvetica', col = 'white', vjust = 1.5, face = 'italic', size = 3)+
  scale_fill_gradientn(expression('Admission \n rate'), 
                       colours = wes_palette("Zissou1", 100, type = "continuous"))+
  scale_x_continuous(expand = c(0, 0), breaks = c(0, seq(1, 9, 1) - 3/12), labels = c('3m', seq(1, 9, 1)))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 100, 10))+
  ylab(expression(atop(italic(Pf)~PR['2-10'], ''))) +
  xlab('\nAge (years)')+ 
  ggtitle('a.')+
  theme(
    legend.title = element_text(colour = 'white', family = 'Helvetica', size = 8),
    legend.text = element_text(colour = 'white', family = 'Helvetica', size = 7),
    legend.title.align = -1,
    legend.position = c(0.925, 0.775),
    legend.key.size = unit(1, "lines"),
    legend.background = element_rect(fill = 'transparent', colour = 'transparent')) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0),
         size = guide_legend(title.position="top", title.hjust = 0))

summarised %>%
  filter(PR %in% c(0, 25, 50, 75)) %>%
  mutate(PR = factor(PR, levels = c(0, 25, 50, 75), labels = c('<1%', '25%', '50%', '75%'))) %>%
  ggplot(aes(x = Age, y = Median_fit, col = PR, fill = PR, group = PR))+
  geom_ribbon(aes(x = Age, ymin = Lower_fit , ymax = Upper_fit), alpha = 0.5, size = 0.2, col = 'transparent')+
  geom_line(lineend = 'round')+
  scale_fill_manual(expression(atop('', italic(Pf)~PR['2-10'])), values = wes_palette("Zissou1", 4, type = "continuous"))+
  scale_colour_manual(expression(atop('', italic(Pf)~PR['2-10'])), values = wes_palette("Zissou1", 4, type = "continuous"))+
  scale_x_continuous(breaks = c(0, seq(1, 9, 1) - 3/12), labels = c('3m', seq(1, 9, 1)))+
  ylab('Annual admissions per 1000 children (3 months - 9 years)\n')+
  xlab('\nAge (years)')+
  ggtitle('b.')+
  theme(    legend.title = element_text(colour = 'grey33', family = 'Helvetica', size = 8),
            legend.text = element_text(colour = 'grey33', family = 'Helvetica', size = 7),
            legend.title.align = -1,
            legend.position = c(0.9, 0.85),
            legend.key.size = unit(1, "lines"),
            legend.background = element_rect(fill = 'transparent', colour = 'transparent')) 

######################
# Fits for each site #
######################

posteriors[, str_detect(names(posteriors), 'RE')] %>%
  as.matrix(ncol = 1) %>%
  as.data.frame() %>%
  mutate(Iteration = 1:n()) %>%
  sample_n(5000) %>%
  gather(Parm, Value, -Iteration) %>% 
  mutate(RE = 'RE_1',
         Study_number = parse_number(Parm)) %>% 
  dplyr::select(-Parm) %>% 
  spread(RE, Value) -> REs

posteriors[, str_detect(names(posteriors), 'PR_corrected')] %>% 
  mutate(Iteration = 1:n()) %>%
  filter(Iteration %in% REs$Iteration) %>% 
  gather(Study_number, PR, - Iteration) %>% 
  mutate(Study_number = parse_number(Study_number)) -> the_PRs

posteriors[, str_detect(names(posteriors), 'beta') | str_detect(names(posteriors), 'alpha_parm') | str_detect(names(posteriors), 'PR_crit') | str_detect(names(posteriors), 'sigma')] %>%
  mutate(Iteration = 1:n()) %>%
  filter(Iteration %in% REs$Iteration) %>% 
  as.matrix() %>%
  as.data.frame() -> parms 

names(parms) <- gsub(pattern = '[', replacement = '', x = names(parms), fixed = T)
names(parms) <- gsub(pattern = ']', replacement = '', x = names(parms), fixed = T)

preds_RE <- expand.grid(Age_years =  c(min(actual_admissions$Age_years), seq(1, 10, 1) - 3/12),
                        Iteration = unique(parms$Iteration))

REs %>%
  left_join(the_PRs) %>% 
  left_join(parms) %>%
  left_join(preds_RE) %>%
  group_by(Study_number, Iteration) %>%
  mutate(Age_plus = lead(Age_years)) %>%  
  filter(!is.na(Age_plus)) %>% 
  ungroup() %>%  
  mutate(LLD = log_logistic(PR, alpha2 = alpha_parm2, beta2 = betas2, sigma = sigma, Age_years),
         LLDP = log_logistic(PR, alpha2 = alpha_parm2, beta2 = betas2, sigma = sigma, Age_plus),
         Integrate = (LLD + (LLDP - LLD)) * (Age_plus - Age_years), 
         Integrate = ifelse(PR <= PR_crit, 1 / length(unique(preds_RE$Age_years)), Integrate),
         Fit = get_preds(PR, alpha1 = alpha_parm1, beta1 = betas1)) -> with_fits_sites

with_fits_sites %>%
  ungroup() %>%
  filter(!is.na(Fit)) %>%
  group_by(Study_number, Iteration) %>%
  mutate(Denominator = sum(Integrate),
         Integrate_2 = Integrate / Denominator) %>%
  group_by(Study_number, Age_years) %>%
  summarise(Median_fit = (median(Integrate_2 * Fit, na.rm = T)), 
            Lower_fit = (quantile(Integrate_2 * Fit, probs = c(0.025, 0.975), na.rm = T)[2]),
            Upper_fit = (quantile(Integrate_2 * Fit, probs = c(0.025, 0.975), na.rm = T)[1]),
            Median_PR = (median(PR, na.rm = T)) * 100, 
            Lower_PR = (quantile(PR, probs = c(0.025, 0.975), na.rm = T)[2]) * 100,
            Upper_PR = (quantile(PR, probs = c(0.025, 0.975), na.rm = T)[1]) * 100) %>%
  rename(Study_Number = Study_number)-> summarised_sites

all_data %>%
  mutate(Study_ID = paste(Site, Period, sep = ' '),
         Study_Number = as.numeric(as.factor(Study_ID))) %>% 
  mutate(CM = ifelse(BCS < 3 | APVU == 'U' | Unconscious_Observed, T, F), 
         CM = ifelse(is.na(CM), F, CM), 
         C = Hb_level < 5 | Deep_breathing | CM,
         C_alt = Hb_level < 5 | Pallor == 'Severe' | Blood_transfusion | Deep_breathing | Inter_costal_recession | CM) %>%
  group_by(Study_ID, Study_Number, Age_years) %>%
  summarise(N = sum(C, na.rm = T),
            N_alt = sum(C_alt, na.rm = T),
            PY = sum(unique(N_Age)), 
            Rate = (N / PY) * 1000,
            Rate_alt = (N_alt / PY) * 1000) %>%
  ungroup() %>%
  mutate(RN = 1:n())-> the_cases

the_cases %>%
  distinct(Study_ID) %>%
  mutate(RN = 1:n()) -> just_sites

summarised_sites %>%
  left_join(just_sites, by = c('Study_Number' = 'RN')) -> summarised_sites

posteriors[, str_detect(names(posteriors), 'estimated_admissions')] %>% 
  sample_n(10000) %>%
  as.matrix() %>%
  as.data.frame() %>%
  gather(Parm, Value) %>% 
  mutate(RN = parse_number(Parm)) %>%
  dplyr::select(-Parm) %>% 
  group_by(RN) %>%
  summarise(Median_site_rate = median(Value),
            Upper_site_rate = (quantile(Value, probs = c(0.025, 0.975))[2]),
            Lower_site_rate = (quantile(Value, probs = c(0.025, 0.975))[1])
  ) %>% 
  ungroup() %>%
  left_join(the_cases) %>% 
  mutate(Age_years = ifelse(Age_years > 0, Age_years - 0.25, Age_years)) %>% 
  mutate(Median_site_rate = Median_site_rate / PY * 1000, 
         Upper_site_rate = Upper_site_rate / PY * 1000, 
         Lower_site_rate = Lower_site_rate / PY * 1000) -> site_data


summarised_sites %>%
  full_join(site_data) %>%
  ggplot(aes(x = Age_years))+
  geom_bar(aes(y = Rate, group = factor(Age_years)),
           stat = 'identity', position ="identity", fill = 'grey75', alpha = 1)+
  geom_linerange(aes(ymin = Lower_site_rate, ymax = Upper_site_rate), col = 'black')+
  geom_point(aes(y = Median_site_rate), col = 'black', pch = 21, fill = 'white')+
  geom_ribbon(aes(ymin = Lower_fit, ymax = Upper_fit, group = Study_ID), fill = wes_palette("Darjeeling1")[5], alpha = 0.5)+
  geom_line(aes(y = Median_fit, group = Study_ID), col = 'blue')+
  scale_x_continuous(breaks = c(0, seq(1, 9, 1) - 3 / 12), labels = c('3m-', paste0(seq(1, 9, 1), '-')))+
  ylab('Annual admissions per 1000 children (3 months - 9 years)\n')+
  xlab('\n Age (years)') +
  facet_wrap(.~ fct_reorder(Study_ID, Median_PR), ncol = 5)+
  theme(legend.position = 'none', 
        strip.text = element_text(size = 9, face = 'bold')) 


