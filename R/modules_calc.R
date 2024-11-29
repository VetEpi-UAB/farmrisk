
#### Origin #### 


origin_expression<-quote({
  #Probability of a herd being free of disease
  h_free<-((1-h_prev )/(1-h_prev*plan_sensi))*(1-h_pos)
  
  #Probability of a herd being infected
  h_inf<-1-h_free
  
  #probability that an animal in a herd is infected
  a_inf<-h_inf*w_prev
  
  #Probability that an infectd animal is infectious
  a_inf_infectious<-a_inf*infectious_inf
  
  #Probability that an animal in a herd is a BVD PI
  a_inf_pi<-h_inf*w_prev_pi
  
  #Incidence rate (asuming steady-state)
  w_inc_rate<-(w_prev/(1-w_prev))/inf_time
  
  #Cumulative incidence during susceptibility window
  w_inc_cum<-1-(1-w_inc_rate)^sus_time
  
  #Within-herd prevalence of Trojan cows
  w_prev_tr<-w_inc_cum*pregnant_p
  
  #Probability that cow in a herd is a BVD TR
  a_inf_tr<-h_inf*w_prev_tr
  
  #Probability that an animal is infected in any form of the disease
  #(for BVD that means TI, PI, or TR)
  a_inf_all<-1-(1-a_inf)*(1-a_inf_pi)*(1-a_inf_tr)
  
  
  #Output for module analysis
  output<-a_inf_all
  
})


#### Test #### 


test_expression<-quote({
  
  #Test adjusted sensitivity for recent infection (if test time is over test optimun time, test sensi does not change)
  
  #Equivalent to: ifelse(test_time>test_opt_time, test_opt_time, test_time)
  test_time_adj <- test_time * (test_time <= test_opt_time) + test_opt_time * (test_time > test_opt_time)
  
  #NAs form 0/0 are 0
  test_sensi_adj<-mcnode_na_rm(test_sensi*(test_time_adj/test_opt_time))
  
  #Probability a tested infected animal is a false negative
  #NAs form 0/0 are 0
  a_false_neg<-mcnode_na_rm(
    test_useful_i*(a_inf*(1-test_sensi_adj))/((1-test_sensi_adj)*a_inf+(test_spec*(1-a_inf))))
  
  #Probability an infectd animal is not tested
  a_no_test<-(1-test_useful_i)*a_inf
  
  #Probability an infected animal is no detected (false negative or not tested)
  a_no_detect<-a_false_neg+a_no_test
  
  #Probability a not detected animal is infectious
  a_no_detect_infectious<-a_no_detect*infectious_inf
  
  #Probability an infected animal is detected (and removed)
  a_detect <- (1-a_false_neg)*a_inf
  
  #Probability a detected animal is infectious
  a_detect_infectious <- a_detect*infectious_inf

  #FOR BVD PI ANIMALS
  #Probability a tested PI animal is a false negative
  #NAs form 0/0 are 0
  a_false_neg_pi<-mcnode_na_rm(
    test_useful_pi*(a_inf_pi*(1-test_sensi))/((1-test_sensi)*a_inf_pi+(test_spec*(1-a_inf_pi))))
  
  #Probability a PI animal is not tested
  a_no_test_pi<-(1-test_useful_pi)*a_inf_pi
  
  #Probability a PI animal is no detected (false negative or not tested)
  a_no_detect_pi<-a_false_neg_pi+a_no_test_pi
  
  #Probability a PI animal is detected (and removed)
  a_detect_pi <- (1-a_false_neg_pi)*a_inf_pi
  
  
  #FOR BVD TR COWS
  #Probability a tested trojan cow is a false negative
  a_false_neg_tr<-mcnode_na_rm(
    test_useful_tr*(a_inf_tr*(1-test_sensi))/((1-test_sensi)*a_inf_tr+(test_spec*(1-a_inf_tr))))
  
  #Probability an trojan cow is not tested
  a_no_test_tr<-(1-test_useful_tr)*a_inf_tr
  
  #Probability an trojan cow is no detected (false negative or not tested)
  a_no_detect_tr<-a_false_neg_tr+a_no_test_tr
  
  #Probability an trojan cow is detected (and removed)
  a_detect_tr <- 1-a_no_detect_tr
  
  #Output for module analysis
  output<-1-(1-a_no_detect)*(1-a_no_detect_pi)*(1-a_no_detect_tr)
  
})


#### Direct contact #### 


dir_contact_expression<-quote({
  #Probability one animal is infected due to direct contact with an infected animal
  a_dir_contact<-direct*other_b_inf_infectious*inf_dc
  
  #Probability one animal is infected due to direct contact with a BVD PI animal
  a_dir_contact_pi<-direct*other_b_inf_pi*inf_dc_pi
  
  #Output for module analysis
  output<-1-(1-a_dir_contact)*(1-a_dir_contact_pi)
  
})


#### Indirect contact #### 


indir_contact_expression<-quote({
  #Probability of cleaning and disinfection efficacy
  hygiene_eff<-cleaning*cleaning_eff*disinfection*disinfection_eff
  
  #Probability of pathogen survival at time contact
  survival_contact_time<-log10(1 + (10^survival_init- 1) * exp(-survival_k * risk_contact_time))
  
  
  #FOR INFECTIOUS ANIMALS
  #Probability of surface contamination by an infectious animal
  inf_surface_init<-prev_b_inf_infectious*indir
  
  #Probability one animal is infected due to one contact with a contaminated surface
  a_indir_contact_one<-inf_surface_init*(1-hygiene_eff)*survival_contact_time*(inf_ic)^indir_level
  
  #Probability one animal is infected due to contact in any of times it has contact with the surface
  a_indir_contact<-1-(1-a_indir_contact_one)^n_times
  
  
  
  #FOR BVD PI ANIMALS
  #Probability of surface contamination by a pi animal
  inf_surface_init_pi<-prev_b_inf_pi*indir
  
  #Probability one animal is infected due to one contact with a contaminated surface
  a_indir_contact_one_pi<-inf_surface_init_pi*(1-hygiene_eff)*survival_contact_time*(inf_ic_pi)^indir_level
  
  #Probability one animal is infected due to contact in any of times it has contact with the surface
  a_indir_contact_pi<-1-(1-a_indir_contact_one_pi)^n_times
  
  #Output for module analysis
  output<-1-(1-a_indir_contact)*(1-a_indir_contact_pi)

})



#### Time #### 


time_expression<-quote({
  
  #Rate at which infected animals contact with the risk point surface
  contamination_rate <- other_contact_rate*other_inf_a
  
  #Rate at which animals contact with the risk point surface
  risk_contact_rate <- contamination_rate+contact_rate
  
  #Time between contamination and risk contact
  #Assuming random exponential variable (Poisson process)
  #risk_contact_time<-mcstoc(rexp, rate=risk_contact_rate)
  
  #Using only mean contact time to reduce variability
  risk_contact_time<-1/risk_contact_rate
  
})



#### Area effect infection #### 

area_inf_expression<-quote({
  #Probability a farm gets infected if it has n_times close herds, 
  #with a h_inf probability of being infected
  area_inf<-1-((1-(1-(1-area_effect*h_inf*inf_ic^indir_level)))^exponent)^(risk_days/360)
  
  #Output for module analysis
  output<-area_inf
}) 