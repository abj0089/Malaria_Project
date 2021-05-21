
rm(list = ls(all = T))

setwd('/Users/ball4364/iCloudDrive/EEID/Malaria/Data/')

pacman::p_load(tidyverse, readxl)

##########################
### ADMISSIONS DATA ######
##########################

read_excel('Malaria admissions - 3m-9yr_all (070520).xlsx', guess_max = 10000) %>% 
  mutate(Site = ifelse(Site == "Kadzinuni-Mauveni", "Kadziununi-Mauveni (Kilifi C)", 
                       ifelse(Site == "Junju", "Junju (Kilifi B)",
                              ifelse(Site == "Ngerenya", "Ngerenya (Kilifi A)", Site))),
         Cerebral = as.numeric(BCS < 3 | `unconscious observed` == 'Y' | APVU == 'U'),
         Cerebral = ifelse(is.na(Cerebral), 0, Cerebral),
         Anemia = ifelse(is.na(SMA), 0, 1),
         Respiratory = ifelse(`deep breathing` == 'Y', 1, 0),
         Respiratory = ifelse(is.na(Respiratory), 0, Respiratory),
         Period = ifelse(Period == '2006-2007', "2006-07", Period),
         Period = ifelse(Site == 'Kadziununi-Mauveni (Kilifi C)', '2018-19', Period)) %>%
  dplyr::select(Site, 
                Period, 
                Admission_month = `month adm`, 
                Admission_year = `year adm`,
                Sex = gender,
                Age_years = `Age years`,
                Age_months = `Age mths`,
                Temperature = temp,
                Fever = fever, 
                Fever_days = `Fever days`,
                Convulsions = convulsions,
                Fits_per_24hrs = `Fits/24HRS`,
                Nasal_flaring = `Nasal Flaring`,
                Inter_costal_recession = `Inter costal recession`,
                Deep_breathing = `deep breathing`,
                Unconscious_Observed = `unconscious observed`,
                APVU = `APVU`,
                BCS = `BCS`,
                Pallor = `pallor`,                
                Hb_level = `hblevel`, 
                Blood_transfusion = `bloodtransfusion`,
                Absconded_or_discharged = `AD`,
                SMA = SMA) %>%
  mutate(Sex = tolower(Sex),
         Fever = ifelse(Fever == "Y", T, ifelse(Fever == 'N', F, Fever)),
         Convulsions = ifelse(Convulsions == "Y", T, ifelse(Convulsions == 'N', F, NA)),
         Nasal_flaring = ifelse(Nasal_flaring == "Y", T, ifelse(Nasal_flaring == 'N', F, NA)),
         Inter_costal_recession = ifelse(Inter_costal_recession == "Y", T, ifelse(Inter_costal_recession == 'N', F, NA)),
         Deep_breathing = ifelse(Deep_breathing == "Y", T, ifelse(Deep_breathing == 'N', F, NA)),
         Unconscious_Observed = ifelse(Unconscious_Observed == "Y", T, ifelse(Unconscious_Observed == 'N', F, NA)),
         Blood_transfusion = ifelse(Blood_transfusion == "Y", T, ifelse(Blood_transfusion == 'N', F, NA)),
         SMA = ifelse(SMA == 'Y', T, SMA),
         Absconded_or_discharged = ifelse(Absconded_or_discharged == 'A', 
                                          "Absconded",
                                          ifelse(Absconded_or_discharged == 'D', 
                                                 'Discharged',
                                                 ifelse(Absconded_or_discharged == 'Unknown', NA, Absconded_or_discharged))))-> admissions_data

##########################
##### LOCATION DATA ######
##########################

b <- 1.8
A <- seq(0, 20, 0.1)
c <- 0.07

P_A <- function(A, b){
  
  1 - exp(- b * A)
  
}

F_A <- function(A, c, alpha, s){
  
  (A < alpha) * 1 + (A >= alpha) * (1 - s * (1 - exp(- c * (A - alpha))))
  
}

plot(F_A(A, 0.07, 9.5, 0.64))


read_excel('Malaria High Summary 3m-9y 2006-20 (140520).xlsx', 
           sheet = 'Background features') %>% 
  dplyr::select(Country, 
                Site, 
                Period = Dates, 
                PR_dates = `Dates of PR data`,
                Detection_Method = `Method used`, 
                latitude = `Lat community`, 
                longitude = `Long community`,
                HbAS = HbAS,
                Hospital_min_dist = `Min dist to hospital`,
                Hospital_max_dist = `Max distance to hospital`, 
                ITN = `ITN coverage`, 
                IRS = IRS, 
                Rainfall = `Av an. rainfall in period`, 
                Altitude = `Altitude MaSL`,
                PR_c = `Age corrected PR 2-10`,
                N_examined = `Number examined`,
                N_positive = `Number PF positive`) %>%
  mutate(IRS = ifelse(IRS %in% c('[N]', 'N'), F, T),
         PR = (N_positive / N_examined) * 100,
         R = PR / PR_c,
         R_me = ) -> locations_formatted

##########################
######## AGE DATA ######## 
##########################

uganda_age_sheets <- excel_sheets('Reveiw Age Summaries Uganda 2012-19 (080520) 2.xlsx')

uganda_age_profiles <- data.frame()

for(i in uganda_age_sheets){
  
  read_excel('Reveiw Age Summaries Uganda 2012-19 (080520) 2.xlsx', sheet = i) %>% 
    mutate(file_id = i) %>%
    bind_rows(uganda_age_profiles, .) -> uganda_age_profiles
  
}

uganda_age_profiles %>%
  mutate(Period = str_trunc(file_id, width = 7, side = 'left', ellipsis = ''),
         Site = str_remove(file_id, paste0(' ', Period)), 
         Age = ifelse(Age == '3m-11m', '-0', Age),
         Age = as.numeric(parse_number(Age))) %>%
  dplyr::select(Site, Period, Age_years = Age, N_Age = PYOR) -> uganda_age_profiles_formatted

kenya_tanzania_files <- excel_sheets('Review Age Summaries Kenya & TZ 2006-20 (080520).xlsx')

kenya_tanzania_age_profiles <- data.frame()

for(i in kenya_tanzania_files){
  
  read_excel('Review Age Summaries Kenya & TZ 2006-20 (080520).xlsx', sheet = i) %>% 
    mutate(file_id = i) %>%
    bind_rows(kenya_tanzania_age_profiles, .) -> kenya_tanzania_age_profiles
  
}

kenya_tanzania_age_profiles %>%
  mutate(Period = str_trunc(file_id, width = 7, side = 'left', ellipsis = ''),
         Site = str_remove(file_id, paste0(' ', Period)), 
         Age = ifelse(Age == '3m-11m', '-0', Age),
         Age = as.numeric(parse_number(Age))) %>%
  dplyr::select(Site, Period, Age_years = Age, N_Age = PYOR) -> kenya_tanzania_age_profiles_formatted

kenya_tanzania_age_profiles_formatted %>%
  bind_rows(uganda_age_profiles_formatted) %>%
  mutate(Site = ifelse(Site == "Kadzinuni-Mauveni", "Kadziununi-Mauveni (Kilifi C)", 
                       ifelse(Site == "Junju", "Junju (Kilifi B)",
                              ifelse(Site == "Ngerenya", "Ngerenya (Kilifi A)", Site))),
         Period = ifelse(Period == '2017-19' & Site == "Kadziununi-Mauveni (Kilifi C)", '2018-19', Period)) %>%
  group_by(Site, Period) %>%
  mutate(Total_Population = sum(N_Age)) %>%
  ungroup()-> all_age_data

############
# COMBINED #
############

admissions_data %>%
  left_join(locations_formatted, by = c('Site', 'Period')) %>%
  left_join(all_age_data, by = c('Site', 'Period', 'Age_years')) -> combined

write_csv(combined, 'complete_data.csv')

















