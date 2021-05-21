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

  ####################
  # REGRESSION MODEL #
  ####################
  
  for(hospital in SITES){

    lambda[hospital] <- exp(intercept + log(person_years[hospital]) + beta1 / (1 + exp(-beta2 * PR_corrected[hospital])))

    admissions[hospital] ~ dpois(lambda[hospital])

  }
  
  ##########
  # PRIORS #
  ##########
  
  # PR Model
  
  for(survey in SITES){
    
    P_prime[survey] ~ dbeta(1, 1)
    
  }
  
  # Regression
  
  intercept ~ dnorm(0, 1 / 100)
  
  beta1 ~ dnorm(0, 1 / 100)
  beta2 ~ dnorm(0, 1 / 100)
  
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

  ####################
  # REGRESSION MODEL #
  ####################
  
  for(hospital in SITES){

    lambda[hospital] <- exp(intercept + log(person_years[hospital]) + beta1 / (1 + exp(-beta2 * PR_corrected[hospital])))
    
    pnb[hospital] <- r / (r + lambda[hospital])
    
    admissions[hospital] ~ dnegbin(pnb[hospital], r) 

  }
  
  ##########
  # PRIORS #
  ##########
  
  # PR Model
  
  for(survey in SITES){
    
    P_prime[survey] ~ dbeta(1, 1)
    
  }

  
  # Regression
  
  intercept ~ dnorm(0, 1 / 100)
  
  beta1 ~ dnorm(0, 1 / 100)
  beta2 ~ dnorm(0, 1 / 100)

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
  
  ####################
  # REGRESSION MODEL #
  ####################
  
  for(hospital in SITES){

    lambda[hospital] <- exp(intercept + beta1 * PR_corrected[hospital] + log(person_years[hospital]))
    
    pnb[hospital] <- r / (r + lambda[hospital])
    
    admissions[hospital] ~ dnegbin(pnb[hospital], r) 

  }
  
  ##########
  # PRIORS #
  ##########
  
  # PR Model
  
  for(survey in SITES){
    
    P_prime[survey] ~ dbeta(1, 1)
    
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
  
  ####################
  # REGRESSION MODEL #
  ####################
  
  for(hospital in SITES){

    lambda[hospital] <- exp(intercept + log(person_years[hospital]))
    
    pnb[hospital] <- r / (r + lambda[hospital])
    
    admissions[hospital] ~ dnegbin(pnb[hospital], r) 

  }
  
  ##########
  # PRIORS #
  ##########
  
  # PR Model
  
  for(survey in SITES){
    
    P_prime[survey] ~ dbeta(1, 1)
    
  }
  
  # Regression
  
  intercept ~ dnorm(0, 1 / 100)

  r ~ dunif(0, 100)
  
  slope_PR ~ dnorm(0.97, 1 / (0.01 ^ 2))
  intercept_PR ~ dnorm(-0.22, 1 / (0.005 ^ 2))

}'

