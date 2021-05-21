#######
# GLM #
#######

##########################
# Poisson vs NB 
##########################

poisson_GLM_logisitic <- 'model{
  
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

    patient_status[patient] ~ dbinom(rho[patient_index[patient]], N_patients[patient])

  }

  for(hospital in SITES){

    estimated_admissions[hospital] <- actual_admissions[hospital] + sum(inprod(rho[probs], N_NA[hospital, probs]))

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
  
  for(p in probs){
  
    rho[p] ~ dbeta(1, 1)
  
  }
  
  # Regression
  
  intercept ~ dnorm(0, 1 / 100)
  
  beta1 ~ dnorm(0, 1 / 100)
  beta_log ~ dnorm(0, 1 / 100)
  beta2 <- exp(beta_log)
  
  slope_PR ~ dnorm(0.97, 1 / (0.01 ^ 2))
  intercept_PR ~ dnorm(-0.22, 1 / (0.005 ^ 2))

}'

neg_bin_GLM_logisitic <- 'model{
  
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

    patient_status[patient] ~ dbinom(rho[patient_index[patient]], N_patients[patient])

  }

  for(hospital in SITES){

    estimated_admissions[hospital] <- actual_admissions[hospital] + sum(inprod(rho[probs], N_NA[hospital, probs]))

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
  
  for(p in probs){
  
    rho[p] ~ dbeta(1, 1)
  
  }
  
  # Regression
  
  intercept ~ dnorm(0, 1 / 100)
  
  beta1 ~ dnorm(0, 1 / 100)
  beta_log ~ dnorm(0, 1 / 100)
  beta2 <- exp(beta_log)

  r ~ dunif(0, 100)
  
  slope_PR ~ dnorm(0.97, 1 / (0.01 ^ 2))
  intercept_PR ~ dnorm(-0.22, 1 / (0.005 ^ 2))

}'


  
#######################################
# NB GLM cubic vs Linear vs Intercept #
#######################################

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

    patient_status[patient] ~ dbinom(rho[patient_index[patient]], N_patients[patient])

  }

  for(hospital in SITES){

    estimated_admissions[hospital] <- actual_admissions[hospital] + sum(inprod(rho[probs], N_NA[hospital, probs]))

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
  
  for(p in probs){
  
    rho[p] ~ dbeta(1, 1)
  
  }
  
  # Regression
  
  intercept ~ dnorm(0, 1 / 100)
  
  beta1 ~ dnorm(0, 1 / 100)

  r ~ dunif(0, 100)
  
  slope_PR ~ dnorm(0.97, 1 / (0.01 ^ 2))
  intercept_PR ~ dnorm(-0.22, 1 / (0.005 ^ 2))

}'


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

    patient_status[patient] ~ dbinom(rho[patient_index[patient]], N_patients[patient])

  }

  for(hospital in SITES){

    estimated_admissions[hospital] <- actual_admissions[hospital] + sum(inprod(rho[probs], N_NA[hospital, probs]))

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
  
  for(p in probs){
  
    rho[p] ~ dbeta(1, 1)
  
  }
  
  # Regression
  
  intercept ~ dnorm(0, 1 / 100)

  r ~ dunif(0, 100)
  
  slope_PR ~ dnorm(0.97, 1 / (0.01 ^ 2))
  intercept_PR ~ dnorm(-0.22, 1 / (0.005 ^ 2))

}'



