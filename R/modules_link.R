
#### Farm link expression ####

farm_link_expression<-quote({
  #Assuming that the probability of the animal being infected just before testing is negligible
  test_time<-test_opt_time
})



#### Unknown Farm link expression ####

unk_link_expression<-quote({
  
  #OTHER FARMS (in current journeys)
  
  #PROBABILITY MULTIPLE ORIGINS - Did your animals and the others come from the same farm? (veh_epiunit)
  #if the vehicle veh_epiunit="Yes", probability of sharing is 0 (veh_mult_origin=1)
  #if the vehicle veh_epiunit="Unknown", probability of sharing is is taken from literature (veh_mult_origin=veh_mult_origin_p)
  #if the vehicle veh_epiunit="No", probability of sharing is 1 (veh_mult_origin=1)
  
  #Equivalent to: ifelse(veh_epiunit_yes,0, ifelse(veh_epiunit_unk,veh_mult_origin_p,1))
  veh_mult_origin <- (1 - veh_epiunit_yes) * (veh_epiunit_unk * veh_mult_origin_p + (1 - veh_epiunit_unk))
  
  #NUMBER OF ORIGINS IN A SHARED TRANSPORT
  #Value is taken from the literature if the probability of multiple origins is not zero.
  
  #Equivalent to:ifelse(veh_mult_origin>0,veh_farms_from_n,0)
  other_farms_n <- veh_farms_from_n * (veh_mult_origin>0)
  
  #Other animals (by farm) in the same vehicle
  #NAs form 0/0 are 0
  
  other_animals_n<-mcnode_na_rm(other_all_animals_n/other_farms_n)
  

  #PREVIOUS FARMS (in previous journeys)
  
  #NUMBER OF PREVIOUS FARMS - Is the vehicle shared with other ruminant farms? (veh_otherfarms)
  #If veh_own=TRUE, number of previous farms is given by user bsg survey (prev_farms_n=veh_otherfarms_cattle_n)
  #If veh_own=FALSE, number of previous farms is taken from literature (prev_farms_n=veh_farms_from_n)
  #If the vehicle veh_otherfarms=FALSE ("No"), number of previous farms is 0 (prev_farms_n=0)
  #Probability the previous journey had animals from several farms (veh_mult_origin_p)
  
  #Equivalent to: ifelse(veh_own,ifelse(veh_otherfarms,veh_otherfarms_cattle_n,0),veh_farms_from_n)
  #Expected number of farms if the vehicle had more than one origin
  prev_farms_mult_n <- veh_own * veh_otherfarms * veh_otherfarms_cattle_n +
    (1 - veh_own) * veh_farms_from_n
  
  #Number of farms
  #If prev_farms_mult=0 it means there are no previous farms
  #If prev_farms_mult=>0 
    # the previous movement only carried animals from one (1) farm with a probability of 1-veh_mult_origin_p
    # previous movement carried animals from more than one (prev_farms_mult_n) holding with a probability of veh_mult_origin_p
  prev_farms_n<-(prev_farms_mult_n>0) *veh_mult_origin_p + prev_farms_mult_n*(1-veh_mult_origin_p)

  #Previous animals (by farm) in the same vehicle
  prev_animals_n<-mcnode_na_rm(prev_all_animals_n/prev_farms_n)
  
  
})


#### Transport link expression ####

transport_link_expression<-quote({
  
  #Probabilty that the vehicle carries at least one infected animal
  prev_b_inf_infectious<-prev_b_inf_infectious_agg
  prev_b_inf_pi<-prev_b_inf_pi_agg
  
  #Probabilty that the vehicle carries at least one infected animal
  other_b_inf_infectious<-other_b_inf_infectious_agg
  other_b_inf_pi<-other_b_inf_pi_agg
  
  #PROBABILITY SHARED TRANSPORT - Were only the animal for this farm transported in this movement? (veh_exclusive)
  #if the vehicle  veh_exclusive="Yes", probability of sharing is 0 (veh_share=0)
  #if the vehicle veh_exclusive="Unknown", probability of sharing is taken from literature (veh_share=veh_share_p)
  #if the vehicle veh_exclusive="No", probability of sharing is 1 (veh_share=1)
  
  #Equivalent to: ifelse(veh_exclusive_yes,0, ifelse(veh_exclusive_unk,veh_share_p,1))
  veh_share <- (1 - veh_exclusive_yes) * (1 - veh_exclusive_unk * (1 - veh_share_p))
  
  
  #PROBABILITY DIRECT CONTACT: other animals in the vehicle and animals come from diferent farms
  direct<-veh_share*veh_mult_origin*(1-veh_none)

  #PROBABILITY INDIRECT CONTACT:
  #If veh_own=TRUE and veh_otherfarms=FALSE then the probability of indirect contact is 0
  #Else the probability of indirect contact is 1
  indir<-1-(veh_own*(1-veh_otherfarms))*(1-veh_none)
  

  #PROBABILITY CLEANING AND DISINFECTION
  #If veh_cleaning = "Each time the vehicle is used"  then veh_cleaning=TRUE
  cleaning <- veh_cleaning
  disinfection <- veh_cleaning
  
  #Time between unloading animals in a farm and loading new animals
  risk_contact_time<-veh_time_between
  
  #Level of indirect contact risk (1=High, 2=Moderate, 3=Low, 4=Very low)
  indir_level<-1
  
  #Number of times indirect contact happens
  n_times<-1
  
})

#### Quarantine link I expression ####

quarantine_I_link_expression<-quote({

  #Time when diagnostic tests are performed during quarantine
  test_time<-quarantine_test_time
  
})

#### Quarantine II link expression ####
quarantine_II_link_expression<-quote({
  #Probability direct contact
  direct<-(1-quarantine)*quarantine_direct

  #Probability at least one animal is infectious during quarantine but is detected and will be removed
  prev_b_inf_infectious<-b_detect_infectious
  prev_b_inf_pi<-b_detect_pi
  
  #Probability indirect contact
  indir<-quarantine*(1-quarantine_direct)*quarantine_equipment
  
  cleaning<-quarantine_cleaning
  disinfection<-quarantine_disinfection
  
  #Time between quarantine visit and herd visit
  risk_contact_time<-(1-quarantine_visits_never)*quarantine_visits_time
  
  #Level of indirect contact risk (1=High, 2=Moderate, 3=Low, 4=Very low)
  indir_level<-2
  
  #Number of times indirect contact happens
  n_times<-mcnode_na_rm(quarantine_time/quarantine_frequency)
})


#### Weighted risk #### 

quarantine_weight_expression<-quote({
  a_entry_all<-1-((1-a_no_detect)*(1-a_no_detect_pi)*(1-a_no_detect_tr)*(1-a_indir_contact)*(1-a_indir_contact_pi))
  transport_a_contact_all_weight<-a_entry_all*(transport_a_contact_all/(transport_a_contact_all+farm_a_no_detect_all))
  farm_a_no_detect_all_weight<-a_entry_all*(farm_a_no_detect_all/(transport_a_contact_all+farm_a_no_detect_all))
}) 

#### Fattening link expression ####

fattening_link_expression<-quote({
  #Probability indirect contact
  indir<-(1-fattening_direct)*fattening_indirect*fattening_animal
  
  
  cleaning<-1-indir
  disinfection<-1-indir
  
  #Probability at least one animal in infected
  prev_b_inf<-1-((1-b_no_detect)*(1-b_no_detect_pi)*(1-b_contact))
  
  #Time between quarantine visit and herd visit
  risk_contact_time<- fattening_frequency
  
  #Level of indirect contact risk (1=High, 2=Moderate, 3=Low, 4=Very low)
  indir_level<-2
  
  #Number of times indirect contact happens
  n_times<-mcnode_na_rm(fattening_time/fattening_frequency)
  
})


#### Visits link expression ####

visit_link_expression<-quote({
  
  #Probability fomites entering the farm (no excusive boots/equipment/vehicle_enters)
  #Equivalent to: visit_livestock AND ifelse(visit_enter>=0,visit_enter,visit_enter_p)
  indir<-visit_livestock*(visit_enter * (visit_enter >= 0) + visit_enter_p * (visit_enter < 0))
  
  #Disease control plan sensitivity in visited farms (assumed to be unknown)
  plan_sensi<-0
  h_pos<-0
  
  #Proportion of pregnant cows
  #not relevant for direct or indirect contact, but needed for origin calc
  pregnant_p<-0
  
  #Probability cleaning and disinfecting fomites in visited farms
  #visit_cleaning=1 ("Always")
  #visit_cleaning=0.5 ("Sometimes"),
  #visit_cleaning=0 ("Never"),
  #If visit_cleaning=-1 ("Unknown"), then visit_cleaning_p from bibliography
  #Equivalent to: ifelse(visit_cleaning<0,visit_cleaning_p,visit_cleaning)
  cleaning <- visit_cleaning * (visit_cleaning >= 0) + visit_cleaning_p * (visit_cleaning < 0)
  
  disinfection<-cleaning
  
  #Previous farms visited by fomite
  farms_n<-visit_farms_n*visit_enter_p
  
  #Previous animals (by farm) contacted by fomite
  animals_n<-visit_animals_n*animal_category_p
  
  #Time between quarantine visit and herd visit
  risk_contact_time<-visit_time_between
  
  #Level of indirect contact risk (1=High, 2=Moderate, 3=Low, 4=Very low)
  #if they enter in direct contact with herd animals: Moderate,
  #if they are close to farm animals (eg. they help loading other animals): Low, 
  #if they don't contact any animal: Very low
  indir_level<-(visit_direct*2)+(visit_close*(1-visit_direct)*3)+((1-visit_close)*(1-visit_direct)*4)

  #Number of times indirect contact happens
  n_times<-risk_days/visit_frequency
  
})


#### Neighbour link expression ####

neighbour_link_expression<-quote({
  
  #We assume probability of direct or indirect contact with neighbour farms 
  #is dependent on the number of farms
  farms_n<-neighbour_n
  
  #PROBABILITY DIRECT CONTACT (one neighbour farm)
  direct<-neighbour_direct
  other_b_inf_infectious<-h_prev
  other_b_inf_pi<-h_prev

  #ZONE EFFECT (PROBABILITY INDIRECT CONTACT one neighbour farm)
  area_effect<-neighbour_farms
  exponent<-farms_n
  indir_level<-3
  h_inf<-h_prev
  
})

#### Wildlife area link expression ####

wildlife_area_link_expression<-quote({
  
  #We assume probability of contact with
  #is dependent on the number of wildlife visits per day 
  #and wildlife pathogen prevalence
  #ZONE EFFECT (PROBABILITY INDIRECT CONTACT infected wildlife)
  area_effect<-access*fencing_wildlife+access*(1-fencing_perimeter)
  exponent<-contact_rate_wildlife
  indir_level<-3
  h_inf<-wl_prev
  
})

#### Wildlife contact water link expression ####

wildlife_water_link_expression<-quote({
  farms_n<-1

  contact_point_n<-livestock_units*contact_point_ratio*contact_point_p
  
  #1 Time link
  other_contact_rate<-mcnode_na_rm(contact_rate_wildlife/contact_point_n)
  other_inf_a<-wl_prev
  contact_rate<-mcnode_na_rm(contact_rate_livestock/contact_point_n)
  
  #2 Indirect contact with mud at waterer
  indir<-(access*fencing_wildlife+access*(1-fencing_perimeter))*mud_p*contact_point_access_wildlife
  
  #Probability cleaning and disinfecting contact_points
  cleaning<-0
  disinfection<-0
  
  #Probability of infected wild host visit
  prev_b_inf_infectious<-other_contact_rate*other_inf_a
  prev_b_inf_pi<-0
  
  #Level of indirect contact risk (1=High, 2=Moderate, 3=Low, 4=Very low)
  #if they enter in direct contact with animal: medium, if not: minimun
  indir_level<-2
  
  #Number of times indirect contact happens
  #surface_p represents the probability the indirect contact happens in a certain surface (eg. mud 50% vs pastic 50%)
  n_times<-risk_days
  
})


#### Animal mix link expression ####

mix_link_expression<-quote({
  
  #Direct is the probability of each animal of having direct contact with animals from the other herd during pasture
  #If pasture_cattlebull=true, then we assume  75%  animals of our farm will have direct contact with the other farm animals
  #If pasture_cattlebull=false, then we assume  50%  animals of our farm will have direct contact with the other farm animals
  
  #PROBABILITY DIRECT CONTACT
  direct<-pasture_share
  pasture_mix_animals_n<-pasture_other_animals_n*((pasture_cattlebull*0.75)+((1-pasture_cattlebull)*0.5))
  other_b_inf_infectious<-1-(1-a_no_detect_infectious)^(pasture_other_farm_n*pasture_mix_animals_n)
  other_b_inf_pi<-1-(1-a_no_detect_pi)^(pasture_other_farm_n*pasture_mix_animals_n)
  
})

