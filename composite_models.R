##########################
# Poisson vs NB 
##########################

poisson_GLM_logistic <- 'model{
  
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
  
  for(patient in PATIENTS){
    
    patient_status[patient, 1] ~ dbern(rho_SMA[prob_index[patient, 1]])
    patient_status[patient, 2] ~ dbern(rho_RD[prob_index[patient, 2]])

    both_status[patient] <- sum(patient_status[patient, 1:3])

    composite_status[patient] <- both_status[patient] > 0

  }

  for(hospital in SITES){

    estimated_admissions[hospital] <- sum(composite_status[site_indexes[hospital, 1]:site_indexes[hospital, 2]])

    correction[hospital] <- lambda[hospital] / (estimated_admissions[hospital] * (estimated_admissions[hospital] > 0) + lambda[hospital] * (estimated_admissions[hospital] == 0))

  }
  
  ####################
  # REGRESSION MODEL #
  ####################
  
  for(hospital in SITES){

    lambda[hospital] <- exp(intercept + log(person_years[hospital]) + beta1 / (1 + exp(-beta2 * PR_corrected[hospital])))

    admissions[hospital] ~ dpois(lambda[hospital] * correction[hospital])

  }
  
  ##########
  # PRIORS #
  ##########
  
  # PR Model
  
  for(survey in SITES){
    
    P_prime[survey] ~ dbeta(1, 1)
    
  }
  
  # Admissions model 
  
  for(p in probs_SMA){
  
    rho_SMA[p] ~ dbeta(1, 1)
  
  }
  
  for(p in probs_RD){
  
    rho_RD[p] ~ dbeta(1, 1)
  
  }
  
  # Regression
  
  intercept ~ dnorm(0, 1 / 100)
  
  beta1 ~ dnorm(0, 1 / 100)
  beta2_log ~ dnorm(0, 1 / 100)
  beta2 <- exp(beta2_log)
  
  slope_PR ~ dnorm(0.97, 1 / (0.01 ^ 2))
  intercept_PR ~ dnorm(-0.22, 1 / (0.005 ^ 2))
  
}'


neg_bin_GLM_logistic <- 'model{
  
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
  
  for(patient in PATIENTS){
    
    patient_status[patient, 1] ~ dbern(rho_SMA[prob_index[patient, 1]])
    patient_status[patient, 2] ~ dbern(rho_RD[prob_index[patient, 2]])

    all_status[patient] <- sum(patient_status[patient, 1:3])

    composite_status[patient] <- all_status[patient] > 0

  }

  for(hospital in SITES){

    estimated_admissions[hospital] <- sum(composite_status[site_indexes[hospital, 1]:site_indexes[hospital, 2]])

    correction[hospital] <- lambda[hospital] / (estimated_admissions[hospital] * (estimated_admissions[hospital] > 0) + lambda[hospital] * (estimated_admissions[hospital] == 0))

  }
  
  ####################
  # REGRESSION MODEL #
  ####################
  
  for(hospital in SITES){

    lambda[hospital] <- exp(intercept + log(person_years[hospital]) + beta1 / (1 + exp(-beta2 * PR_corrected[hospital])))
    
    pnb[hospital] <- r / (r + lambda[hospital] * correction[hospital])
    
    admissions[hospital] ~ dnegbin(pnb[hospital], r) 

  }
  
  ##########
  # PRIORS #
  ##########
  
  # PR Model
  
  for(survey in SITES){
    
    P_prime[survey] ~ dbeta(1, 1)
    
  }
  
  # Admissions model 
  
  for(p in probs_SMA){
  
    rho_SMA[p] ~ dbeta(1, 1)
  
  }
  
  for(p in probs_RD){
  
    rho_RD[p] ~ dbeta(1, 1)
  
  }
    
  # Regression
  
  intercept ~ dnorm(0, 1 / 10)
  
  beta1 ~ dnorm(0, 1 / 100)
  beta2_log ~ dnorm(0, 1 / 100)
  beta2 <- exp(beta2_log)
  
  r ~ dunif(0, 100)
  
  slope_PR ~ dnorm(0.97, 1 / (0.01 ^ 2))
  intercept_PR ~ dnorm(-0.22, 1 / (0.005 ^ 2))

}'

#######################################
# NB GLM cubic vs Linear vs Intercept #
#######################################

neg_bin_GLM_intercept <- 'model{
  
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
  
  for(patient in PATIENTS){
    
    patient_status[patient, 1] ~ dbern(rho_SMA[prob_index[patient, 1]])
    patient_status[patient, 2] ~ dbern(rho_RD[prob_index[patient, 2]])

    both_status[patient] <- sum(patient_status[patient, 1:3])

    composite_status[patient] <- both_status[patient] > 0

  }

  for(hospital in SITES){

    estimated_admissions[hospital] <- sum(composite_status[site_indexes[hospital, 1]:site_indexes[hospital, 2]])

    correction[hospital] <- lambda[hospital] / (estimated_admissions[hospital] * (estimated_admissions[hospital] > 0) + lambda[hospital] * (estimated_admissions[hospital] == 0))

  }
  
  ####################
  # REGRESSION MODEL #
  ####################
  
  for(hospital in SITES){

    lambda[hospital] <- exp(intercept + log(person_years[hospital]))
    
    pnb[hospital] <- r / (r + lambda[hospital] * correction[hospital])
    
    admissions[hospital] ~ dnegbin(pnb[hospital], r) 

  }
  
  ##########
  # PRIORS #
  ##########
  
  # PR Model
  
  for(survey in SITES){
    
    P_prime[survey] ~ dbeta(1, 1)
    
  }
  
  # Admissions model 
  
  for(p in probs_SMA){
  
    rho_SMA[p] ~ dbeta(1, 1)
  
  }
  
  for(p in probs_RD){
  
    rho_RD[p] ~ dbeta(1, 1)
  
  }
  
  # Regression
  
  intercept ~ dnorm(0, 1 / 100)

  r ~ dunif(0, 100)
  
  slope_PR ~ dnorm(0.97, 1 / (0.01 ^ 2))
  intercept_PR ~ dnorm(-0.22, 1 / (0.005 ^ 2))

}'




neg_bin_GLM_linear <- 'model{
  
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
  
  for(patient in PATIENTS){
    
    patient_status[patient, 1] ~ dbern(rho_SMA[prob_index[patient, 1]])
    patient_status[patient, 2] ~ dbern(rho_RD[prob_index[patient, 2]])

    both_status[patient] <- sum(patient_status[patient, 1:3])

    composite_status[patient] <- both_status[patient] > 0

  }

  for(hospital in SITES){

    estimated_admissions[hospital] <- sum(composite_status[site_indexes[hospital, 1]:site_indexes[hospital, 2]])

    correction[hospital] <- lambda[hospital] / (estimated_admissions[hospital] * (estimated_admissions[hospital] > 0) + lambda[hospital] * (estimated_admissions[hospital] == 0))

  }
  
  ####################
  # REGRESSION MODEL #
  ####################
  
  for(hospital in SITES){

    lambda[hospital] <- exp(intercept + beta1 * PR_corrected[hospital] + log(person_years[hospital]))
    
    pnb[hospital] <- r / (r + lambda[hospital] * correction[hospital])
    
    admissions[hospital] ~ dnegbin(pnb[hospital], r) 

  }
  
  ##########
  # PRIORS #
  ##########
  
  # PR Model
  
  for(survey in SITES){
    
    P_prime[survey] ~ dbeta(1, 1)
    
  }
  
  # Admissions model 
  
  for(p in probs_SMA){
  
    rho_SMA[p] ~ dbeta(1, 1)
  
  }
  
  for(p in probs_RD){
  
    rho_RD[p] ~ dbeta(1, 1)
  
  }
  
  # Regression
  
  intercept ~ dnorm(0, 1 / 100)
  
  beta1 ~ dnorm(0, 1 / 100)

  r ~ dunif(0, 100)

  slope_PR ~ dnorm(0.97, 1 / (0.01 ^ 2))
  intercept_PR ~ dnorm(-0.22, 1 / (0.005 ^ 2))
  
}'

