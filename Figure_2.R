
rm(list = ls(all = T))

pacman::p_load(readxl, tidyvers, patchwork, wesanderson, rjags, runjags, HDInterval)

setwd('/Users/ball4364/Desktop/Malaria/Data/Our data/')


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

all_data %>%
  distinct(Study_Number, Total_Population) %>%
  rename(PY = Total_Population)-> total_pops

all_data %>% 
  filter(!(is.na(APVU) & is.na(Unconscious_Observed) & is.na(BCS))) %>%
  mutate(CM = ifelse(BCS < 3 | APVU == 'U' | Unconscious_Observed, T, F)) %>%
  group_by(Study_Number) %>%
  summarise(N = sum(CM, na.rm = T),
            PY = unique(Total_Population)) %>%
  mutate(N = ifelse(is.na(N), 0, N),
         Median_rate = (N / PY) * 1000) %>%
  dplyr::select(Study_Number, Median_rate)-> case_numbers

##################################################
# Load model fits 
##################################################

load(file = '../../Output/PfPR_model_fits_CM.RData')

do.call(rbind, all_fits$NB_intercept$mcmc) %>%
  as.data.frame() -> posteriors_CM

posteriors_CM %>%  
  gather(Parameter, Value, 1:ncol(.)) %>% 
  group_by(Parameter) %>%
  summarise(Lower95 = quantile(Value, 0.025),
            Upper95 = quantile(Value, 0.975),
            Lower99 = quantile(Value, 0.01),
            Upper99 = quantile(Value, 0.99),
            Median = median(Value)) %>%
  dplyr::select(Parameter, Median, Lower95, Upper95, Lower99, Upper99) -> processed_posteriors_CM 

rm(all_fits)

load(file = '../../Output/PfPR_model_fits_SMA.RData')

do.call(rbind, all_fits$NB_Logistic$mcmc) %>%
  as.data.frame() -> posteriors_SMA

posteriors_SMA %>%  
  gather(Parameter, Value, 1:ncol(.)) %>% 
  group_by(Parameter) %>%
  summarise(Lower95 = quantile(Value, 0.025),
            Upper95 = quantile(Value, 0.975),
            Lower99 = quantile(Value, 0.01),
            Upper99 = quantile(Value, 0.99),
            Median = median(Value)) %>%
  dplyr::select(Parameter, Median, Lower95, Upper95, Lower99, Upper99) -> processed_posteriors_SMA

rm(all_fits)

load(file = '../../Output/PfPR_model_fits_RD.RData')

do.call(rbind, all_fits$NB_linear$mcmc) %>%
  as.data.frame() -> posteriors_RD

posteriors_RD %>%  
  gather(Parameter, Value, 1:ncol(.)) %>% 
  group_by(Parameter) %>%
  summarise(Lower95 = quantile(Value, 0.025),
            Upper95 = quantile(Value, 0.975),
            Lower99 = quantile(Value, 0.01),
            Upper99 = quantile(Value, 0.99),
            Median = median(Value)) %>%
  dplyr::select(Parameter, Median, Lower95, Upper95, Lower99, Upper99) -> processed_posteriors_RD

rm(all_fits)

#####################################################
# PREDS
#####################################################

########################
# SMA
########################

processed_posteriors_SMA %>%
  filter(str_detect(Parameter, pattern = 'estimated_admissions')) %>% 
  mutate(Study_Number = parse_number(Parameter)) %>%
  left_join(total_pops) %>%
  mutate(Median_rate = (Median / PY) * 1000,
         Lower_rate95 = (Lower95 / PY) * 1000, 
         Upper_rate95 = (Upper95  / PY) * 1000,
         Lower_rate99 = (Lower99 / PY) * 1000, 
         Upper_rate99 = (Upper99  / PY) * 1000) %>%
  dplyr::select(-Parameter, -Median, -Lower95, -Upper95, -Lower99, -Upper99) -> estimated_admissions_SMA

processed_posteriors_SMA %>%
  filter(str_detect(Parameter, pattern = 'PR_corrected') & !str_detect(Parameter, pattern = 'PR_corrected_mic')) %>% 
  mutate(Study_Number = parse_number(Parameter)) %>%
  left_join(estimated_admissions_SMA) %>%
  mutate(Pathology = 'a. severe malaria anemia')-> corrected_PRs_SMA

posteriors_SMA[, c('intercept', 'beta1', 'beta2')] %>%
  #sample_n(20000) %>%
  as.matrix()-> chains_parms_SMA

PR_range_SMA <- seq(0, max(corrected_PRs_SMA$Upper95), length.out = 500)

fits_SMA <- matrix(NA, nrow = length(PR_range_SMA), ncol = nrow(posteriors_SMA))

for(i in 1:length(PR_range_SMA)){
  
  fits_SMA[i, ] <- posteriors_SMA[, 'intercept'] + posteriors_SMA[, 'beta1'] / (1 + exp(-posteriors_SMA[, 'beta2'] * PR_range_SMA[i] ))
  
}

medians_SMA <- apply(X = fits_SMA , FUN = median, MARGIN = 1)
#intervals_SMA <- apply(X = fits_SMA, FUN = hdi, MARGIN = 1)
intervals_SMA <- apply(X = fits_SMA, FUN = quantile, MARGIN = 1, c(0.01, 0.025, 0.975, 0.99))

fit_data_frame_SMA <- data.frame(PR = PR_range_SMA * 100, 
                                 Median_fit = exp(medians_SMA + log(1000)), 
                                 Upper_95 = exp(intervals_SMA[3, ] + log(1000)), 
                                 Lower_95 = exp(intervals_SMA[2, ] + log(1000)),
                                 Upper_99 = exp(intervals_SMA[4, ] + log(1000)),
                                 Lower_99 = exp(intervals_SMA[1, ] + log(1000)),
                                 Pathology = 'a. severe malaria anemia')

########################
# RD
########################

processed_posteriors_RD %>%
  filter(str_detect(Parameter, pattern = 'estimated_admissions')) %>% 
  mutate(Study_Number = parse_number(Parameter)) %>%
  left_join(total_pops) %>%
  mutate(Median_rate = (Median / PY) * 1000,
         Lower_rate95 = (Lower95 / PY) * 1000, 
         Upper_rate95 = (Upper95  / PY) * 1000,
         Lower_rate99 = (Lower99 / PY) * 1000, 
         Upper_rate99 = (Upper99  / PY) * 1000) %>%
  dplyr::select(-Parameter, -Median, -Lower95, -Upper95, -Lower99, -Upper99) -> estimated_admissions_RD

processed_posteriors_RD %>%
  filter(str_detect(Parameter, pattern = 'PR_corrected') & !str_detect(Parameter, pattern = 'PR_corrected_mic')) %>% 
  mutate(Study_Number = parse_number(Parameter)) %>%
  left_join(estimated_admissions_RD) %>%
  mutate(Pathology = 'b. respiratory distress')-> corrected_PRs_RD

posteriors_RD[, c('intercept', 'beta1')] %>% #, 'beta2', 'beta3')] %>%
  #sample_n(20000) %>%
  as.matrix()-> chains_parms_RD

PR_range_RD <- seq(0, max(corrected_PRs_RD$Upper95), length.out = 100)

model_matrix_RD <- cbind(1, PR_range_RD)#, PR_range ^ 3)

fits_RD  <- model_matrix_RD %*% t(chains_parms_RD)

medians_RD <- apply(X = fits_RD , FUN = median, MARGIN = 1)
#intervals_RD <- apply(X = fits_RD, FUN = hdi, MARGIN = 1)
intervals_RD <- apply(X = fits_RD, FUN = quantile, MARGIN = 1, c(0.01, 0.025, 0.975, 0.99))

fit_data_frame_RD <- data.frame(PR = PR_range_RD * 100, 
                                Median_fit = exp(medians_RD + log(1000)), 
                                Upper_95 = exp(intervals_RD[3, ] + log(1000)), 
                                Lower_95 = exp(intervals_RD[2, ] + log(1000)),
                                Upper_99 = exp(intervals_RD[4, ] + log(1000)),
                                Lower_99 = exp(intervals_RD[1, ] + log(1000)),
                                Pathology = 'b. respiratory distress')

########################
# CM
########################

processed_posteriors_CM %>%
  filter(str_detect(Parameter, pattern = 'PR_corrected') & !str_detect(Parameter, pattern = 'PR_corrected_mic')) %>% 
  mutate(Study_Number = parse_number(Parameter)) %>%
  left_join(case_numbers) %>%
  mutate(Pathology = 'c. cerebral malaria')-> corrected_PRs_CM

PR_range_CM <- seq(0, max(corrected_PRs_CM$Upper95), length.out = 100)

intercept <- filter(processed_posteriors_CM, str_detect(Parameter, 'intercept'))

fit_data_frame_CM <- data.frame(PR = PR_range_CM * 100, 
                                Median_fit = exp(as.numeric(intercept[, 2]) + log(1000)), 
                                Upper_95 = exp(as.numeric(intercept[, 4]) + log(1000)), 
                                Lower_95 = exp(as.numeric(intercept[, 3]) + log(1000)),
                                Upper_99 = exp(as.numeric(intercept[, 6]) + log(1000)),
                                Lower_99 = exp(as.numeric(intercept[, 5]) + log(1000)),
                                Pathology = 'c. cerebral malaria') %>%
  mutate(Pathology = factor(Pathology, levels = c('a. severe malaria anemia', 'b. respiratory distress', 'c. cerebral malaria')))


all_fits <- bind_rows(fit_data_frame_RD, fit_data_frame_SMA, fit_data_frame_CM) %>%
  mutate(Pathology = factor(Pathology, levels = c('a. severe malaria anemia', 'b. respiratory distress', 'c. cerebral malaria')))

corrected_PRs_all <- bind_rows(corrected_PRs_RD, corrected_PRs_SMA, corrected_PRs_CM)%>%
  mutate(Pathology = factor(Pathology, levels = c('a. severe malaria anemia', 'b. respiratory distress', 'c. cerebral malaria')))

all_data %>%
  distinct(Site, Period, Study_Number) -> SP 

read_csv('/complete_data.csv') %>%
  group_by(Site, Period) %>%
  distinct(Detection_Method) %>%
  mutate(Site = ifelse(Site == "Junju (Kilifi B)", 'Kilifi B', Site),
         Site = ifelse(Site == "Kadziununi-Mauveni (Kilifi C)", 'Kilifi C', Site),
         Site = ifelse(Site == "Ngerenya (Kilifi A)", 'Kilifi A', Site),
         Detection_Method = ifelse(Detection_Method == 'Microsopy', 'Microscopy', Detection_Method)) %>%
  left_join(SP) %>%
  distinct(Detection_Method, Study_Number)-> sampling_data

all_fits %>%
  mutate(RATE = max(Upper_99)) %>%
  distinct(Pathology, RATE) %>%
  mutate(PR  = 0, RATE = RATE) -> labs

corrected_PRs_all %>%
  left_join(sampling_data) %>%
  mutate(Median = Median * 100,
         Lower95 = Lower95 * 100,
         Upper95 = Upper95 * 100,
         Lower99 = Lower99 * 100,
         Upper99 = Upper99 * 100) %>% 
  ggplot(aes(x = Median, y =  Median_rate))+
  geom_linerange(aes(ymin = Lower_rate95, ymax = Upper_rate95), col = 'grey33', lineend = 'round', size = 0.25)+
  geom_linerange(aes(xmin = Lower95, xmax = Upper95), col = 'grey33', lineend = 'round', size = 0.25)+
  geom_point(aes(shape = Detection_Method), fill = 'white', colour = 'grey33', stroke = 0.5, size = 1.5)+
  scale_shape_manual('', values = c(21, 19))+
  geom_ribbon(aes(x = PR, ymin = Lower_99, ymax = Upper_99),
              inherit.aes = F, data = all_fits, fill = wes_palette("Darjeeling1")[5], alpha = 0.2)+
  geom_ribbon(aes(x = PR, ymin = Lower_95, ymax = Upper_95),
              inherit.aes = F, data = all_fits, fill = wes_palette("Darjeeling1")[5], alpha = 0.4)+
  geom_line(aes(x = PR, y = Lower_99),
            inherit.aes = F, data = all_fits, col = 'black', lty = 1, size = 0.1, alpha = 0.75)+
  geom_line(aes(x = PR, y = Upper_99),
            inherit.aes = F, data = all_fits, col = 'black', lty = 1, size = 0.1, alpha = 0.75)+
  geom_line(aes(x = PR, y = Lower_95),
            inherit.aes = F, data = all_fits, col = 'black', lty = 1, size = 0.1, alpha = 0.75)+
  geom_line(aes(x = PR, y = Upper_95),
            inherit.aes = F, data = all_fits, col = 'black', lty = 1, size = 0.1, alpha = 0.75)+
  geom_line(aes(x = PR, y = Median_fit),
            inherit.aes = F, data = all_fits , size = 0.66, col = 'black', lineend = 'round')+
  geom_label(aes(x = PR, y = RATE, label = Pathology), data = labs, hjust = 0.05, family = 'Helvetica', col = 'grey33',
             fill = 'white', label.padding = unit(0.1, "lines"),
             label.r = unit(0, "lines"),
             label.size = 0) +
  scale_x_continuous(breaks = seq(0, 70, 10))+
  scale_y_continuous(breaks = seq(0, 100, 1))+
  facet_wrap(~ Pathology, ncol = 3)+
  ylab('Annual admissions per 1000 children (3 months - 9 years)\n')+
  xlab(expression(atop('', italic(Pf)~PR['2-10']))) +
  theme(strip.text = element_blank(),
        legend.position = 'none') 

