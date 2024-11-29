#' Calculate risk of disease entry in cattle farms
#' 
#' @description
#' This function calculates the probability of disease entry through different pathways in cattle farms
#' 
#' @param farm_id Character. Unique identifier for the farm
#' @param risk_days Numeric. Time period (in days) for risk assessment calculation. Default is 360
#' @param calc_summary Logical. Whether to calculate summary for all nodes. Default: TRUE
#' @param admin_wif Logical. Whether to run all admin what-if scenarios. Default: TRUE
#' @param model_analysis Logical. Whether to run model analysis. Default: FALSE
#' @param alternative_params List. Optional alternative parameters for simulation
#' @param alternative_formula Formula. Optional formula for parameter modifications
#' @param inputs_path Character. Path to the directory containing input files. Default is "input_files"
#' @param outputs_path Character. Path to the directory where output files will be saved. Default is "output_files"
#' 
#' @details
#' The function evaluates multiple disease entry pathways including:
#' * Purchase of animals
#' * Pasture movements
#' * Visits
#' * Neighboring farms
#' * Wildlife contact
#' 
#' For each pathway, the function calculates:
#' * Probability of disease entry
#' * Relative risk compared to baseline
#' * Pathway-specific risks
#' 
#' @return
#' Returns a mcmodule containing:
#' * Risk probability summaries
#' * Pathway-specific risk calculations
#' * Comparative risk assessments
#' * Detailed statistics for each disease pathway
#' 
#' @examples
#' \dontrun{
#' farm_risk <- farmrisk_cattle("FARM001", risk_days = 360)
#' }
#' 
#' @export

farmrisk_cattle<-function(farm_id=NULL,
                   risk_days=360,
                   calc_summary = TRUE,
                   admin_wif = TRUE,
                   model_analysis = FALSE,
                   alternative_params = NULL,
                   alternative_formula = NULL,
                   input_path=paste0(system.file("input_files/",package="farmrisk"),"/"),
                   output_path=paste0(system.file("output_files/",package="farmrisk"),"/"),
                   forms_path=paste0(system.file("forms/",package="farmrisk"),"/")){

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

set_up_simulation(farm_id = farm_id, envir = parent.frame(n = 2))
if(exists('mov')&&'pasture'%in%mov$mov_type){
# ---
# title: "Animal pasture pathogen introduction pathway"
# author: "Natalia Ciria"
# editor: visual
# bibliography: references.bib
# execute:
#   output: false
# ---
# 

# 
# ## Description
# 
# This pathway analyses the probability of pathogen introduction when animals go to pastured. The pathway consists of several modules. Each module has its own set of data preparation, evaluation and combination steps. The risk of introduction is also assessed under different what-if biosecurity scenarios.
# 
## ------------------------------------------------------------------------------------------------------------
message("\nPASTURE PATHWAY: ")

# 
# ## Risk days
# 
## ------------------------------------------------------------------------------------------------------------
from_date<-mov$pasture_from_date[!is.na(mov$pasture_from_date)]
to_date<-mov$pasture_to_date[!is.na(mov$pasture_to_date)]
mov$pasture_days<-mov$pasture_to_date-mov$pasture_from_date

from_date_s<-c(min(to_date)-1, sort(from_date))
to_date_s<-c(sort(to_date),max(from_date)+1)

diff_date<-to_date_s-from_date_s
cont_date<-diff_date>5 #Assume pasture continuity if gap is under 5 days

cont_to_date<-ifelse(cont_date,NA, to_date_s)
cont_from_date<-ifelse(cont_date,NA, from_date_s)

pasture_risk_days<-data.frame(cont_to_date,cont_from_date)%>%
  fill(cont_to_date, .direction ="up")%>%
  fill(cont_from_date, .direction ="down")%>%
  mutate(cont_diff_date=cont_to_date-cont_from_date)%>%
  na.omit()%>%
  unique()%>%
  summarise(sum=sum(cont_diff_date))

mov$pasture_risk_days<-mov$pasture_to_date-mov$pasture_from_date

risk_days<-risk_days-pasture_risk_days$sum

# 
# ## Pasture data preparation
# 
# ### Prepare data
# 
# Tidy data from biosecurity and movements surveys
# 
## ------------------------------------------------------------------------------------------------------------
#Tidy number of animals by animal category
pasture_animal<-mov%>%
  filter(mov_type=="pasture")%>%
  tidy_prefix("pasture_animals_p", rm_prefix="pasture")%>%
  left_join(get_columns(mov, c("cattle_n","cattlebull", "farms_n","share"), "pasture"))%>%
  left_join(tidy_census_table(),relationship = "many-to-many")%>%
  mutate(animals_p=as.numeric(sub("%", "", animals_p))/100, #ojo puede ser más de 1
         pasture_own_animals_n=animals_n*animals_p, #animals from user farm
         pasture_own_farm_n=1,
         pasture_other_animals_n=pasture_cattle_n/pasture_farms_n,#animals form other farms
         pasture_other_animals_n=ifelse(is.na(pasture_other_animals_n),0,pasture_other_animals_n),
         pasture_other_farm_n=ifelse(is.na(pasture_farms_n),0,pasture_farms_n))%>%
  left_join(animal)%>%
  mutate(livestock_units=animals_n*livestock_units) 



#Tidy health status by disease and farm
pasture_status_origin<-tidy_status_table(pathway="pasture", module="origin")

#Tidy diagnostic tests before pasture_transport by disease, farm and type of test
pasture_test_origin<-tidy_test_table(pathway="pasture", module="origin")


#Find homogeneous grups
pasture_health_origin<-pasture_animal%>%
  left_join(pasture_status_origin,relationship = "many-to-many")%>%
  left_join(pasture_test_origin)%>%
  tidy_group()


#Get the smallest region code
pasture_region<-tidy_panel(panel="pasture")%>%
  filter(!is.na(pasture_id))%>%
  mutate(pasture_gid0=bsg$farm_gid0)%>%
  get_region_code()

#Get mov vehicles table
pasture_to_veh<-filter(mov, mov_type=="pasture")%>%
  tidy_prefix("pasture_to_veh", rm_prefix="pasture_to")%>%
  mutate(veh_direction="pasture_to")

pasture_from_veh<-filter(mov, mov_type=="pasture")%>%
  tidy_prefix("pasture_from_veh", rm_prefix="pasture_from")%>%
  mutate(veh_direction="pasture_from")

#provisional
pasture_from_veh<-pasture_to_veh%>%
  mutate(veh_direction="pasture_from")

pasture_veh<-bind_rows(pasture_to_veh, pasture_from_veh)%>%
  left_join(tidy_prefix(tidy_panel(panel="veh"),"veh", rm_prefix=FALSE), relationship = "many-to-many")%>%
  filter(!is.na(veh_id))

# 
# Merge tidy user data with admin inputs
# 
## ------------------------------------------------------------------------------------------------------------
#Origin
pasture_origin_data<-pasture_animal%>%
  left_join(pasture_region)%>%
  left_join(pathogen_animal,relationship = "many-to-many")%>%
  hjoin_region(pathogen_region)%>%
  left_join(pathogen)%>%
  left_join(pasture_health_origin)%>%
  left_join(pathogen_status)

# 
## ------------------------------------------------------------------------------------------------------------
#Test
pasture_test_data<-pasture_test_origin%>%
  left_join(pathogen_test)

#Used m2 in each movement
pasture_m2_used<-pasture_animal%>%
  left_join(pasture_region, relationship = "many-to-many")%>%
  hjoin_region(animal_region)%>%
  mutate(m2_used=veh_m2animal*animals_n)%>%
  group_by(farm_id, mov_id)%>%
  summarize(m2_used=sum(m2_used),
            animals_total=sum(animals_n))


#Estimate the number of animals from other farms that could fit in the vehicle
pasture_veh_animal<-pasture_veh%>%
  left_join(pasture_m2_used)%>%
  left_join(pasture_animal, relationship = "many-to-many")%>%
  left_join(pasture_region)%>%
  hjoin_region(animal_region)%>%
  mutate(veh_size=ifelse(is.na(veh_size),veh_size_ext,veh_size),
         m2_free=veh_size-m2_used)%>%
  mutate(animals_p=animals_n/animals_total,
         other_all_animals_n=(animals_p*m2_free)/veh_m2animal,
         other_all_animals_n=ifelse(other_all_animals_n>0,other_all_animals_n,2),
         prev_all_animals_n=veh_size/veh_m2animal,
         surface=ifelse(veh_cleaning=="each","metal","soil"))

# 
# #### Merge all
# 
## ------------------------------------------------------------------------------------------------------------
#All data
pasture_data<-pasture_origin_data%>%
  left_join(pasture_test_data)%>%
  mutate(animals_p=NULL)%>% #desambiguate
  left_join(pasture_veh_animal, relationship = "many-to-many")%>%
  left_join(pathogen_surface)


# 
# #### Separate data by vehicle (direction)
# 
## ------------------------------------------------------------------------------------------------------------
#Separate data by vehicle (direction)
pasture_data_to<-filter(pasture_data, veh_direction=="pasture_to")
message("\npasture_data_to (",paste(dim(pasture_data_to), collapse=", "),") created")

pasture_data_from<-filter(pasture_data, veh_direction=="pasture_from")
message("\npasture_data_from (",paste(dim(pasture_data_from), collapse=", "),") created")

# 
# #### Admin what-if
# 
## ------------------------------------------------------------------------------------------------------------
if(admin_wif){
  #There are no admin what-ifs in this module, but we set up wif for compatibility
  pasture_data_to<-set_up_wif(pasture_data_to)
  pasture_data_from<-set_up_wif(pasture_data_from)

}

# 
# ## Unknown farm module
# 
# ### Prepare data
# 
## ------------------------------------------------------------------------------------------------------------
#Other origin keys
unk_origin_by<-c("mov_id", "farm_id","animal_category", "pathogen")
#Other origin
pasture_unk_origin_data<-pasture_animal%>%
  left_join(pasture_region)%>%
  left_join(pathogen_animal, relationship = "many-to-many")%>%
  left_join(pathogen)%>%
  hjoin_region(pathogen_region)%>%
  #left_join(pasture_health_origin)%>%
  mutate(status="unk")%>%
  left_join(pathogen_status)%>%
  mutate(animals_p=NULL)%>% #desambiguate
  left_join(pasture_veh_animal, relationship = "many-to-many")%>%
  add_group_id(by=unk_origin_by)%>%
  mutate(pregnant_p=0.5) #Does not affect risk but needed to calc origin_expression

#Separate data by vehicle (direction)
#Vehicle to pasture
pasture_to_unk_origin_data<-filter(pasture_unk_origin_data, veh_direction=="pasture_to")
message("\npasture_to_unk_origin_data (",paste(dim(pasture_to_unk_origin_data), collapse=", "),") created")

#Vehicle from pasture
pasture_from_unk_origin_data<-filter(pasture_unk_origin_data, veh_direction=="pasture_from")
message("\npasture_from_unk_origin_data (",paste(dim(pasture_from_unk_origin_data), collapse=", "),") created")

# 
## ------------------------------------------------------------------------------------------------------------
pasture_to_unk_origin_data<-group_match(pasture_to_unk_origin_data,pasture_data_to,unk_origin_by)

pasture_from_unk_origin_data<-group_match(pasture_from_unk_origin_data,pasture_data_from,unk_origin_by)

# 
# ### Evaluate (to)
# 
## ------------------------------------------------------------------------------------------------------------
unk_origin_expression<-c(unk_link=unk_link_expression,
                         unk_origin=origin_expression)

pasture_to_unk_farm<-eval_model_expression(model_expression = unk_origin_expression,
                                     data=pasture_to_unk_origin_data)

pasture_to_unk_farm<-get_totals(mcmodule=pasture_to_unk_farm, mcnodes=c("a_inf_infectious","a_inf_pi"),
                     animals_n="other_animals_n", 
                     farms_n="veh_farms_from_n",
                     prefix="other",
                     all_mcnodes=FALSE)

pasture_to_unk_farm<-get_totals(mcmodule=pasture_to_unk_farm, mcnodes=c("a_inf_infectious","a_inf_pi"),
                     animals_n="prev_animals_n",
                     farms_n="veh_farms_from_n",
                     prefix="prev",
                     all_mcnodes=FALSE)

pasture_to_unk_farm<-get_agg_totals(mcmodule=pasture_to_unk_farm, 
                         mcnode=c("prev_b_inf_infectious"),
                         keys_names=c("pathogen", "farm_id", "mov_id"),
                         agg_variates=TRUE)

pasture_to_unk_farm<-get_agg_totals(mcmodule=pasture_to_unk_farm, 
                         mcnode=c("prev_b_inf_pi"),
                         keys_names=c("pathogen", "farm_id", "mov_id"),
                         agg_variates=TRUE)

pasture_to_unk_farm<-get_agg_totals(mcmodule=pasture_to_unk_farm, 
                         mcnode=c("other_b_inf_infectious"),
                         keys_names=c("pathogen", "farm_id", "mov_id"),
                         agg_variates=TRUE)

pasture_to_unk_farm<-get_agg_totals(mcmodule=pasture_to_unk_farm, 
                         mcnode=c("other_b_inf_pi"),
                         keys_names=c("pathogen", "farm_id", "mov_id"),
                         agg_variates=TRUE)

pasture_to_unk_farm<-add_prefix(pasture_to_unk_farm)

# 
# ### Evaluate (from)
# 
## ------------------------------------------------------------------------------------------------------------
#Check if trasnport to and from are equal
from_equal_to<-all(
  pasture_data_to[, !names(pasture_data_to) %in% c("veh_direction")]==
    pasture_data_from[, !names(pasture_data_from) %in% c("veh_direction")],
  na.rm=TRUE)

#If it they are not not equal calcualte unk farm "from" module
if(!from_equal_to){
  
  pasture_from_unk_farm<-eval_model_expression(model_expression = unk_origin_expression,
                                       data=pasture_from_unk_origin_data)
  
  pasture_from_unk_farm<-get_totals(mcmodule=pasture_from_unk_farm, mcnodes=c("a_inf_infectious","a_inf_pi"),
                       animals_n="other_animals_n", 
                       farms_n="veh_farms_from_n",
                       prefix="other",
                       all_mcnodes=FALSE)
  
  pasture_from_unk_farm<-get_totals(mcmodule=pasture_from_unk_farm, mcnodes=c("a_inf_infectious","a_inf_pi"),
                       animals_n="prev_animals_n",
                       farms_n="veh_farms_from_n",
                       prefix="prev",
                       all_mcnodes=FALSE)
  
  pasture_from_unk_farm<-get_agg_totals(mcmodule=pasture_from_unk_farm, 
                           mcnode=c("prev_b_inf_infectious"),
                           keys_names=c("pathogen", "farm_id", "mov_id"),
                           agg_variates=TRUE)
  
  pasture_from_unk_farm<-get_agg_totals(mcmodule=pasture_from_unk_farm, 
                           mcnode=c("prev_b_inf_pi"),
                           keys_names=c("pathogen", "farm_id", "mov_id"),
                           agg_variates=TRUE)
  
  pasture_from_unk_farm<-get_agg_totals(mcmodule=pasture_from_unk_farm, 
                           mcnode=c("other_b_inf_infectious"),
                           keys_names=c("pathogen", "farm_id", "mov_id"),
                           agg_variates=TRUE)
  
  pasture_from_unk_farm<-get_agg_totals(mcmodule=pasture_from_unk_farm, 
                           mcnode=c("other_b_inf_pi"),
                           keys_names=c("pathogen", "farm_id", "mov_id"),
                           agg_variates=TRUE)
  
  pasture_from_unk_farm<-add_prefix(pasture_from_unk_farm)
}

# 
# ## Pasture transport module
# 
# It is assumed that the probability of infection of an animal is the same as the **probability of infection of the animals currently in the vehicle** (per farm), but the probability of at least one animal being infected depends on the free space in the vehicle (`animals_other_n`) and the number of farms from which it is shared (or number of farms shared by default).
# 
# ### Evaluate (to)
# 
## ------------------------------------------------------------------------------------------------------------
transport_expression<-c(transport_link=transport_link_expression,
                        dir_contact_transport= dir_contact_expression,
                        indir_contact_transport=indir_contact_expression)

pasture_to_transport<-eval_model_expression(
  model_expression = transport_expression,
  prev_mcmodule = pasture_to_unk_farm,
  data = pasture_data_to)
  #create_nodes = FALSE)

pasture_to_transport<-get_totals(mcmodule=pasture_to_transport, mcnodes=c("a_dir_contact","a_dir_contact_pi"),
                     animals_n="pasture_own_animals_n",
                     farms_n="pasture_own_farm_n")

pasture_to_transport<-get_totals(mcmodule=pasture_to_transport, mcnodes=c("a_indir_contact","a_indir_contact_pi"),
                     animals_n="pasture_own_animals_n",
                     farms_n="pasture_own_farm_n")

pasture_to_transport<-get_totals(mcmodule=pasture_to_transport, mcnodes=c("a_dir_contact_all","a_indir_contact_all"),
                     animals_n="pasture_own_animals_n",
                     farms_n="pasture_own_farm_n",
                     name = "a_contact_all")

pasture_to_transport<-get_agg_totals(mcmodule=pasture_to_transport, mcnode=c("a_contact_all"))


pasture_to_transport<-add_prefix(pasture_to_transport)

# 
# ### Evaluate (from)
# 
## ------------------------------------------------------------------------------------------------------------
#If it is equal copy unk farm "to" module
if(from_equal_to){
  pasture_from_transport<-add_prefix(pasture_to_transport, 
                             prefix="pasture_from_transport",
                             rewrite_module="pasture_to_transport")
  message("\npasture_from_transport mcmodule copied from pasture_to_transport")
#Else calcualte unk farm "from" module
}else{
  pasture_from_transport<-eval_model_expression(
    model_expression = transport_expression,
    prev_mcmodule = pasture_from_unk_farm,
    data = pasture_data_from)
    #create_nodes = FALSE)

  pasture_from_transport<-get_totals(mcmodule=pasture_from_transport, mcnodes=c("a_dir_contact","a_dir_contact_pi"),
                       animals_n="pasture_own_animals_n",
                       farms_n="pasture_own_farm_n")
  
  pasture_from_transport<-get_totals(mcmodule=pasture_from_transport, mcnodes=c("a_indir_contact","a_indir_contact_pi"),
                       animals_n="pasture_own_animals_n",
                       farms_n="pasture_own_farm_n")
  
  pasture_from_transport<-get_totals(mcmodule=pasture_from_transport, mcnodes=c("a_dir_contact_all","a_indir_contact_all"),
                       animals_n="pasture_own_animals_n",
                       farms_n="pasture_own_farm_n",
                       name = "a_contact_all")
  
  pasture_from_transport<-get_agg_totals(mcmodule=pasture_from_transport, mcnode=c("a_contact_all"))
  
  pasture_from_transport<-get_agg_totals(mcmodule=pasture_from_transport, mcnode=c("b_contact_all"))

  pasture_from_transport<-add_prefix(pasture_from_transport)
}

# 
# ### Combine modules transport from and to pasture
# 
## ------------------------------------------------------------------------------------------------------------
pasture_transport<-combine_modules(pasture_to_unk_farm, pasture_to_transport)
pasture_transport<-combine_modules(pasture_transport, pasture_from_transport)

pasture_transport<-get_totals(mcmodule=pasture_transport, mcnodes=c("pasture_to_transport_a_dir_contact_all","pasture_from_transport_a_dir_contact_all"),
                       animals_n="pasture_own_animals_n",
                       farms_n="pasture_own_farm_n",
                       name = "a_dir_contact_all")

pasture_transport<-get_totals(mcmodule=pasture_transport, mcnodes=c("pasture_to_transport_a_indir_contact_all","pasture_from_transport_a_indir_contact_all"),
                       animals_n="pasture_own_animals_n",
                       farms_n="pasture_own_farm_n",
                       name = "a_indir_contact_all")

pasture_transport<-get_totals(mcmodule=pasture_transport, mcnodes=c("pasture_to_transport_a_contact_all","pasture_from_transport_a_contact_all"),
                       animals_n="pasture_own_animals_n",
                       farms_n="pasture_own_farm_n",
                       name = "a_contact_all")

pasture_transport<-get_agg_totals(mcmodule=pasture_transport, mcnode=c("a_contact_all"))

pasture_transport<-get_agg_totals(mcmodule=pasture_transport, mcnode=c("b_contact_all"))

pasture_transport<-add_prefix(pasture_transport)

# 
# ## Other farms module
# 
# ### Admin what-if
# 
## ------------------------------------------------------------------------------------------------------------
pasture_data_to<-pasture_data_to%>%
  mutate(pregnant_p=0.5) #Does not affect risk but needed to calc origin_expression

if(admin_wif){
  #Use, and not share, your own vehicle
  pasture_data_to<-wif_no_share_veh(pasture_data_to, scenario="No shared pasture transport", prev_wif="wif_own_veh")
  
  #All farms sharing pastures test and remove positive animals before going to pasture
  pasture_data_to<-wif_test(pasture_data_to, scenario="Screening before pasture")
  
  pasture_data_to<-wif_no_share_pasture(pasture_data_to, scenario="No share pasture")


}

# 
# ### Evaluate
# 
## ------------------------------------------------------------------------------------------------------------
pasture_expression<-c(farm_link=farm_link_expression,
                       origin = origin_expression,
                       test_origin = test_expression)

pasture_herd<-eval_model_expression(model_expression = pasture_expression,
                                     data=pasture_data_to)

pasture_herd<-get_totals(mcmodule=pasture_herd, mcnodes=c("a_inf","a_inf_pi","a_inf_tr"),
                 animals_n="pasture_other_animals_n",
                 farms_n="pasture_other_farm_n")

pasture_herd<-get_totals(mcmodule=pasture_herd, mcnodes=c("a_no_detect","a_no_detect_pi","a_no_detect_tr"),
                 animals_n="pasture_other_animals_n",
                 farms_n="pasture_other_farm_n")


pasture_herd<-get_agg_totals(mcmodule=pasture_herd, mcnode=c("a_inf_all"))

pasture_herd<-get_agg_totals(mcmodule=pasture_herd, mcnode=c("b_inf_all"))

pasture_herd<-get_agg_totals(mcmodule=pasture_herd, mcnode=c("a_no_detect_all"))

pasture_herd<-get_totals(mcmodule=pasture_herd, mcnode=c("a_no_detect_infectious"), all_mcnodes = FALSE,
                 animals_n="pasture_other_animals_n",
                 farms_n="pasture_other_farm_n")

pasture_herd<-add_prefix(pasture_herd)

# 
# ## Animal mix module
# 
# ### Evaluate
# 
## ------------------------------------------------------------------------------------------------------------
mix_expression<-c(mix_link=mix_link_expression,
                  dir_contact_mix= dir_contact_expression)

pasture_mix<-eval_model_expression(
  model_expression = mix_expression,
  prev_mcmodule = pasture_herd,
  data = pasture_data_to)
  #create_nodes = FALSE)

pasture_mix<-get_totals(mcmodule=pasture_mix, mcnodes=c("a_dir_contact","a_dir_contact_pi"),
                animals_n="pasture_own_animals_n",
                farms_n="pasture_own_farm_n")

pasture_mix<-get_agg_totals(mcmodule=pasture_mix, mcnode=c("a_dir_contact_all"))

pasture_mix<-get_agg_totals(mcmodule=pasture_mix, mcnode=c("b_dir_contact_all"))

pasture_mix<-add_prefix(pasture_mix)

# 
# ### Combine modules pasture herd, mix and transport
# 
## ------------------------------------------------------------------------------------------------------------
pasture<-combine_modules(pasture_herd, pasture_mix)

pasture<-combine_modules(pasture_transport, pasture)


#By lot
pasture<-get_agg_totals(mcmodule=pasture,mcnode="pasture_transport_a_contact_all", keys_names=c("pathogen", "farm_id","scenario_id", "mov_id", "animal_category"), suffix="lot")

pasture<-get_agg_totals(mcmodule=pasture,mcnode="pasture_mix_a_dir_contact_all", keys_names=c("pathogen", "farm_id","scenario_id", "mov_id", "animal_category"), suffix="lot")

pasture<-at_least_one(mcmodule=pasture, mcnodes = c("pasture_transport_a_contact_all_lot","pasture_mix_a_dir_contact_all_lot"), name="pasture_a_livestocK_all_lot")

#By mov
pasture<-get_agg_totals(mcmodule=pasture,mcnode="pasture_transport_b_contact_all", keys_names=c("pathogen", "farm_id","scenario_id", "mov_id"), suffix="mov")

pasture<-get_agg_totals(mcmodule=pasture,mcnode="pasture_mix_b_dir_contact_all", keys_names=c("pathogen", "farm_id","scenario_id", "mov_id"), suffix="mov", data_name = "pasture")

pasture<-at_least_one(mcmodule=pasture, mcnodes = c("pasture_transport_b_contact_all_mov","pasture_mix_b_dir_contact_all_mov"), name="pasture_b_livestocK_all_mov")

#By pathogen (all movements)
pasture<-at_least_one(mcmodule=pasture, mcnodes=c("pasture_transport_a_contact_all_agg","pasture_mix_a_dir_contact_all_agg"), name="pasture_a_livestocK_all_agg")

#By pathogen (all movements)
pasture<-at_least_one(mcmodule=pasture, mcnodes=c("pasture_transport_b_contact_all_agg","pasture_mix_b_dir_contact_all_agg"), name="pasture_b_livestock_all_agg")

# 
# ## Wildlife area effect module
# 
# ### Prepare data
# 
## ------------------------------------------------------------------------------------------------------------
pasture_wildlife_area_data<-pasture_region%>%
  left_join(tidy_panel(bsg, "pasture"))%>%
  left_join(tidy_prefix(mov, "pasture_days",rm_prefix = FALSE))%>%
  mutate(surface="soil",
         access=TRUE,
         fencing_wildlife=TRUE,
         fencing_perimeter=FALSE)%>%
  left_join(pathogen_surface)%>%
  hjoin_region(wildlife_region)%>%
  hjoin_region(wildlife_pathogen_region)%>%
  left_join(wildlife_point_density)%>%
  hjoin_region(pathogen_region)%>%
  left_join(pathogen)%>%
  #Filter those that have wildlife prevalence
  filter(pathogen%in%
           unique(wildlife_pathogen_region$pathogen[wildlife_pathogen_region$wl_prev_mode>0]))

message("\npasture_wildlife_area_data (",paste(dim(pasture_wildlife_area_data), collapse=", "),") created")

# 
# ### Prepare what-if
# 
## ------------------------------------------------------------------------------------------------------------
if(admin_wif){
  #There are no admin what-ifs in this module, but we set up wif for compatibility
  pasture_wildlife_area_data<-set_up_wif(pasture_wildlife_area_data)
}

# 
# ### Evaluate
# 
## ------------------------------------------------------------------------------------------------------------
wildlife_area_expression<-c(wildlife_area_link = wildlife_area_link_expression,                         area_inf_wildlife_area=area_inf_expression)  

pasture_wildlife<-eval_model_expression(model_expression = wildlife_area_expression,                              data=pasture_wildlife_area_data,
                             param_names = c(risk_days="pasture_days")) 
pasture_wildlife<-get_agg_totals(mcmodule=pasture_wildlife, mcnode=c("area_inf"))

#By mov
pasture_wildlife<-get_agg_totals(mcmodule=pasture_wildlife,mcnode="area_inf", keys_names=c("pathogen", "farm_id","scenario_id", "mov_id"), suffix="mov")

pasture_wildlife<-add_prefix(pasture_wildlife)

# 
# ## Wildlife contact points module
# 
# ### Prepare data
# 
## ------------------------------------------------------------------------------------------------------------
pasture_wildlife_water_data<-pasture_wildlife_area_data%>%
  pivot_longer(
    cols = ends_with("_p")&!pasture_mud_p&!pasture_sum_p,
    names_to = "waterpoint_type",
    values_to = "point_p",
    values_drop_na = TRUE)%>%
  left_join(pasture_animal)%>%
  mutate(
    contact_point_type=sub("_p", "", waterpoint_type),
    mud_p=ifelse(contact_point_type%in%c("low", "high"),ifelse(pasture_mud,as.numeric(gsub("%", "", pasture_mud_p)),0),100)/100,
    contact_point_p=as.numeric(sub("%", "", point_p))/100,
    contact_point_type=gsub("pasture_","",contact_point_type)
  )

pasture_wildlife_water_data<-pasture_wildlife_water_data%>%
  left_join(animal_point)%>%
  left_join(wildlife_point)

message("\npasture_wildlife_water_data (",paste(dim(pasture_wildlife_water_data), collapse=", "),") created")

# 
# ### Admin what-if
# 
## ------------------------------------------------------------------------------------------------------------
if(admin_wif){
  #Avoid mud on wateres
  pasture_wildlife_water_data<-wif_no_mud(pasture_wildlife_water_data, scenario="Avoid mud on wateres")
  pasture_wildlife_water_data<-wif_high_waterer(pasture_wildlife_water_data, scenario="Use high wateres")
}

# 
# ### Evaluate
# 
## ------------------------------------------------------------------------------------------------------------
wildlife_water_expression<-c(
  wildlife_link = wildlife_water_link_expression,
                            time_wildlife = time_expression,
                            indir_contact_wildlife=indir_contact_expression)

#debugonce(get_node_list)
pasture_wildlife_water<-eval_model_expression(
  model_expression = wildlife_water_expression,
  data=pasture_wildlife_water_data,
  param_names = c(risk_days="pasture_days"))

pasture_wildlife_water<-get_totals(mcmodule=pasture_wildlife_water, mcnode=c("a_indir_contact"))

pasture_wildlife_water<-get_agg_totals(mcmodule=pasture_wildlife_water, mcnode=c("b_indir_contact"))

#By mov
pasture_wildlife_water<-get_agg_totals(mcmodule=pasture_wildlife_water,mcnode="b_indir_contact", keys_names=c("pathogen", "farm_id","scenario_id", "mov_id"), suffix="mov")

pasture_wildlife_water<-add_prefix(pasture_wildlife_water)

# 
# ## Combine modules wildlife area, waterpoints and pasture
# 
## ------------------------------------------------------------------------------------------------------------
pasture_wildlife<-combine_modules(pasture_wildlife, pasture_wildlife_water)

pasture_wildlife<-at_least_one(mcmodule=pasture_wildlife, mcnodes=c("pasture_wildlife_area_inf_agg","pasture_wildlife_water_b_indir_contact_agg"), name="pasture_wildlife_inf_agg")

pasture<-combine_modules(pasture, pasture_wildlife)

#All pasture modules
pasture<-at_least_one(mcmodule=pasture, mcnodes=c("pasture_b_livestock_all_agg","pasture_wildlife_inf_agg"), name="pasture_inf_agg")

#By movement
pasture<-at_least_one(mcmodule=pasture, mcnodes = c("pasture_wildlife_area_inf_mov","pasture_wildlife_water_b_indir_contact_mov"), name="pasture_wildlife_all_mov")

pasture<-at_least_one(mcmodule=pasture, mcnodes = c("pasture_b_livestocK_all_mov","pasture_wildlife_all_mov"), name="pasture_all_mov")

message("\nPASTURE PATHWAY OK!")

# 
# ## Save results
# 
}
if(exists('mov')&&'purchase'%in%mov$mov_type){
# ---
# title: "Animal Purchase pathogen introduction pathway"
# author: "Natalia Ciria"
# editor: visual
# bibliography: references.bib
# execute:
#   output: false
# ---
# 

# 
# ## Description
# 
# This pathway analyses the probability of pathogen introduction when animals are purchased. The pathway consists of several modules such as farm, quarantine, transport and fattening. Each module has its own set of data preparation, evaluation and combination steps. The risk of introduction is also assessed under different what-if biosecurity scenarios.
# 
## ------------------------------------------------------------------------------------------------------------
message("\nPURCHASE PATHWAY: ")

# 
# ## Farm module
# 
# ### Prepare data
# 
# Tidy data from biosecurity and movements surveys
# 
## ------------------------------------------------------------------------------------------------------------
#Tidy number of animals by animal category
purchase_animal<-tidy_animal_table(pathway = "purchase")

#Tidy health status by disease and farm
purchase_status_origin<-tidy_status_table(pathway="purchase", module="origin")

#Tidy diagnostic tests before transport by disease, farm and type of test
purchase_test_origin<-tidy_test_table(pathway="purchase", module="origin")


#Find homogeneous grups
purchase_health_origin<-purchase_animal%>%
  left_join(purchase_status_origin, relationship = "many-to-many")%>%
  left_join(purchase_test_origin)%>%
  left_join(tidy_prefix(mov,"purchase_quarantine", rm_prefix="purchase"))%>%
  tidy_group()


#Get the smallest region code
purchase_region<-get_region_code(mov, pathway = "purchase")

veh<-tidy_prefix(tidy_panel(panel="veh"),"veh", rm_prefix=FALSE)%>%
  filter(!is.na(veh_id))

#Get mov vehicles table
purchase_veh<-filter(mov, mov_type=="purchase")%>%
  tidy_prefix("purchase_veh", rm_prefix="purchase")%>%
  left_join(veh)

# 
# Merge tidy user data with admin inputs
# 
## ------------------------------------------------------------------------------------------------------------
#Origin
purchase_origin_data<-purchase_animal%>%
  left_join(purchase_region)%>%
  left_join(pathogen_animal,relationship = "many-to-many")%>%
  left_join(pathogen)%>%
  hjoin_region(pathogen_region, add_keys=c("animal_category"))%>%
  left_join(purchase_health_origin)%>%
  mutate(status=ifelse(status=="unk_plan", "unk",as.character(status)))%>%
  left_join(pathogen_status)

# 
## ------------------------------------------------------------------------------------------------------------
#Test
purchase_test_data<-purchase_test_origin%>%
  left_join(pathogen_test)

#Used m2 in each movement
purchase_m2_used<-purchase_animal%>%
  left_join(purchase_region)%>%
  hjoin_region(animal_region)%>%
  mutate(m2_used=veh_m2animal*animals_n*farms_total)%>%
  group_by(farm_id, mov_id)%>%
  summarize(m2_used=sum(m2_used),
            animals_total=sum(animals_n*farms_total))


#Estimate the number of animals from other farms that could fit in the vehicle
purchase_veh_animal<-purchase_veh%>%
  left_join(purchase_m2_used)%>%
  left_join(purchase_animal)%>%
  left_join(purchase_region)%>%
  hjoin_region(animal_region)%>%
  mutate(veh_size=ifelse(is.na(veh_size),veh_size_ext,veh_size),
         m2_free=veh_size-m2_used)%>%
  mutate(animals_p=animals_n/animals_total,
         other_all_animals_n=(animals_p*m2_free)/veh_m2animal,
         other_all_animals_n=ifelse(other_all_animals_n>0,other_all_animals_n,2),
         prev_all_animals_n=veh_size/veh_m2animal,
         surface=ifelse(veh_cleaning=="each","metal","soil"))

# 
# #### Merge all
# 
## ------------------------------------------------------------------------------------------------------------
#All data
purchase_data<-purchase_origin_data%>%
  left_join(purchase_test_data)%>%
  left_join(purchase_veh_animal)%>%
  left_join(pathogen_surface)

message("\npurchase_data (",paste(dim(purchase_data), collapse=", "),") created")

# 
# #### Admin what-if
# 
## ------------------------------------------------------------------------------------------------------------
if(admin_wif){
  #Test  before transport
  purchase_data<-wif_test(purchase_data, scenario="Test before transport")
  
  #Use, and not share, your own vehicle
  #purchase_data<-wif_own_veh(purchase_data, scenario="Use, and not share, your own vehicle (purchase)")
  
  #Use, and not share, your own vehicle
  purchase_data<-wif_no_share_veh(purchase_data, scenario="No shared purchase transport", prev_wif="wif_own_veh")
  
  #Clean the vehicle
  #purchase_data<-wif_clean(purchase_data, scenario="Clean and disinfect the vehicle between transports")
  
}

# 
# ### Evaluate
# 
## ------------------------------------------------------------------------------------------------------------
purchase_expression<-c(farm_link=farm_link_expression,
                       origin = origin_expression,
                       test_origin = test_expression)

farm<-eval_model_expression(model_expression = purchase_expression,
                                     data=purchase_data)

farm<-get_totals(mcmodule=farm, mcnodes=c("a_inf","a_inf_pi","a_inf_tr"))

farm<-get_totals(mcmodule=farm, mcnodes=c("a_no_detect","a_no_detect_pi","a_no_detect_tr"))

farm<-get_agg_totals(mcmodule=farm, mcnode=c("a_inf_all"))

farm<-get_agg_totals(mcmodule=farm, mcnode=c("b_inf_all"))

farm<-get_agg_totals(mcmodule=farm, mcnode=c("a_no_detect_all"))

farm<-get_totals(mcmodule=farm, mcnode=c("a_no_detect_infectious"), all_mcnodes = FALSE)

farm<-add_prefix(farm)

# 
# ## Other farm module
# 
# ### Prepare data
# 
## ------------------------------------------------------------------------------------------------------------
#Other origin keys
unk_origin_by<-c("mov_id", "farm_id","animal_category", "pathogen")
#Other origin
purchase_unk_origin_data<-purchase_animal%>%
  left_join(purchase_region)%>%
  left_join(pathogen_animal)%>%
  left_join(pathogen)%>%
  hjoin_region(pathogen_region, add_keys ="animal_category")%>%
  #left_join(purchase_health_origin)%>%
  mutate(status="unk")%>%
  left_join(pathogen_status)%>%
  left_join(purchase_veh_animal)%>%
  add_group_id(by=unk_origin_by)

message("\npurchase_unk_origin_data (",paste(dim(purchase_unk_origin_data), collapse=", "),") created")

# 
## ------------------------------------------------------------------------------------------------------------
purchase_unk_origin_data<-group_match(purchase_unk_origin_data,purchase_data,unk_origin_by)

# 
# ### Evaluate
# 
## ------------------------------------------------------------------------------------------------------------
unk_origin_expression<-c(unk_link=unk_link_expression,
                         unk_origin=origin_expression)

unk_farm<-eval_model_expression(model_expression = unk_origin_expression,
                                     data=purchase_unk_origin_data)

unk_farm<-get_totals(mcmodule=unk_farm, mcnodes=c("a_inf_infectious","a_inf_pi"),
                     animals_n="other_animals_n", 
                     farms_n="veh_farms_from_n",
                     prefix="other",
                     all_mcnodes=FALSE)

unk_farm<-get_totals(mcmodule=unk_farm, mcnodes=c("a_inf_infectious","a_inf_pi"),
                     animals_n="prev_animals_n",
                     farms_n="veh_farms_from_n",
                     prefix="prev",
                     all_mcnodes=FALSE)

unk_farm<-get_agg_totals(mcmodule=unk_farm, 
                         mcnode=c("prev_b_inf_infectious"),
                         keys_names=c("pathogen", "farm_id", "mov_id"),
                         agg_variates=TRUE)

unk_farm<-get_agg_totals(mcmodule=unk_farm, 
                         mcnode=c("prev_b_inf_pi"),
                         keys_names=c("pathogen", "farm_id", "mov_id"),
                         agg_variates=TRUE)

unk_farm<-get_agg_totals(mcmodule=unk_farm, 
                         mcnode=c("other_b_inf_infectious"),
                         keys_names=c("pathogen", "farm_id", "mov_id"),
                         agg_variates=TRUE)

unk_farm<-get_agg_totals(mcmodule=unk_farm, 
                         mcnode=c("other_b_inf_pi"),
                         keys_names=c("pathogen", "farm_id", "mov_id"),
                         agg_variates=TRUE)

unk_farm<-add_prefix(unk_farm)

# 
# ### Combine modules farm and unknown farm
# 
## ------------------------------------------------------------------------------------------------------------
purchase<-combine_modules(farm, unk_farm)

# 
# ## Transport module
# 
# It is assumed that the probability of infection of an animal is the same as the **probability of infection of the animals currently in the vehicle** (per farm), but the probability of at least one animal being infected depends on the free space in the vehicle (`animals_other_n`) and the number of farms from which it is shared (or number of farms shared by default).
# 
# ### Evaluate
# 
## ------------------------------------------------------------------------------------------------------------
transport_expression<-c(transport_link=transport_link_expression,
                        dir_contact_transport= dir_contact_expression,
                        indir_contact_transport=indir_contact_expression)

transport<-eval_model_expression(
  model_expression = transport_expression,
  prev_mcmodule = purchase,
  data = purchase_data)
  #create_nodes = FALSE)


transport<-get_totals(mcmodule=transport, mcnodes=c("a_dir_contact","a_dir_contact_pi"))

transport<-get_totals(mcmodule=transport, mcnodes=c("a_indir_contact","a_indir_contact_pi"))

transport<-get_totals(mcmodule=transport, mcnodes=c("a_dir_contact_all","a_indir_contact_all"), name = "a_contact_all")

transport<-get_agg_totals(mcmodule=transport, mcnode=c("a_contact_all"))

transport<-add_prefix(transport)

# 
# ### Combine modules purchase and transport
# 
## ------------------------------------------------------------------------------------------------------------
purchase<-combine_modules(purchase, transport)

purchase<-get_totals(mcmodule=purchase, mcnodes=c("transport_a_contact_all","farm_a_no_detect_all"), animals_n="farm_animals_n", name="purchase_a_inf_all")

purchase<-get_totals(mcmodule=purchase, mcnodes=c("transport_a_contact_all","farm_a_no_detect"), animals_n="farm_animals_n", name="purchase_a_inf")

# 
# ## Quarantine module
# 
# ### Prepare data
# 
## ------------------------------------------------------------------------------------------------------------
#Tidy diagnostic tests before transport by disease, farm and type of test
purchase_test_quarantine<-tidy_test_table(pathway="purchase", module="quarantine")

#Tidy quarantine parameters
quarantine_bsg_data<-tidy_prefix(bsg, module="quarantine", rm_prefix = FALSE)

#Add movement and tests info
quarantine_data<-quarantine_bsg_data%>%
  left_join(purchase_origin_data[,!names(purchase_origin_data)%in%"test"])%>%
  left_join(purchase_test_quarantine)%>%
  left_join(pathogen_test)%>%
  mutate(quarantine_test_time="end",
         surface=ifelse(quarantine_cleaning=="each","metal","soil"))%>%
  left_join(pathogen_surface)

message("\nquarantine_data (",paste(dim(quarantine_data), collapse=", "),") created")

# 
# ### Admin what-if
# 
## ------------------------------------------------------------------------------------------------------------
if(admin_wif){
  #Match index
  quarantine_by<-c("mov_id", "farm_id","animal_category", "pathogen","status")

  quarantine_data<-group_match(quarantine_data,purchase_data,quarantine_by)
  
  #Clean quarantine equipement
  quarantine_data<-wif_quarantine(quarantine_data, scenario="Quarantine new animals")
  
  #Test  during quarantine
  quarantine_data<-wif_test(quarantine_data, scenario="Test during quarantine", prev_wif = "wif_quarantine")
  
  #Clean quarantine equipement
  quarantine_data<-wif_clean(quarantine_data, scenario="Clean and disinfect quarantine equipment", prev_wif = "wif_quarantine")
  
  #Exclusive quarantine material
  quarantine_data<-wif_exclusive_material(quarantine_data, scenario="Exclusive quarantine material", prev_wif="wif_quarantine")
  
  #Test during quarantine with exclusive material
  quarantine_data<-combine_wif(data=quarantine_data, scenario="Test during quarantine with exclusive material", wif_list=c("quarantine","exclusive_material","test"))
  
}

# 
# ### Evaluate I
# 
## ------------------------------------------------------------------------------------------------------------
quarantine_I_expression<-c(quarantine_I_link=quarantine_I_link_expression,
                     test_quarantine=test_expression)

quarantine_I<-eval_model_expression(
  model_expression = quarantine_I_expression,
  prev_mcmodule = purchase,
  param_names=c(a_inf ="purchase_a_inf",
                a_inf_pi="farm_a_no_detect_pi",
                a_inf_tr="farm_a_no_detect_tr",
                b_contact="transport_b_contact_all"),
  data=quarantine_data,
  match_prev = TRUE)

quarantine_I<-get_totals(mcmodule=quarantine_I, mcnodes=c("a_no_detect","a_no_detect_pi","a_no_detect_tr"))

quarantine_I<-get_agg_totals(mcmodule=quarantine_I, mcnode=c("a_no_detect_all"))

quarantine_I<-get_totals(mcmodule=quarantine_I, mcnode=c("a_detect_infectious","a_detect_pi"), name="a_detect_infectious_all")

quarantine_I<-get_agg_totals(mcmodule=quarantine_I, mcnode=c("a_detect_infectious_all"))

quarantine_I<-add_prefix(quarantine_I)

# 
# ### Combine modules purchase and quarantine I
# 
## ------------------------------------------------------------------------------------------------------------
purchase<-combine_modules(purchase, quarantine_I)

# 
# ### Evaluate II
# 
## ------------------------------------------------------------------------------------------------------------
quarantine_II_expression<-c(quarantine_II_link=quarantine_II_link_expression,
                            indir_contact_quarantine=indir_contact_expression,
                     quarantine_weight=quarantine_weight_expression)

quarantine_II<-eval_model_expression(
  model_expression = quarantine_II_expression,
  prev_mcmodule = purchase,
  param_names=c(b_detect_infectious="quarantine_I_b_detect_infectious",
                b_detect_pi="quarantine_I_b_detect_pi",
                b_no_detect_infectious="quarantine_b_no_detect_infectious",
                b_no_detect_pi="quarantine_b_no_detect_pi",
                b_contact="transport_b_contact_all"),
  data=quarantine_data,
  match_prev = TRUE)


quarantine_II<-get_totals(mcmodule=quarantine_II, mcnodes=c("a_indir_contact", "a_indir_contact_pi"))

quarantine_II<-get_agg_totals(mcmodule=quarantine_II, mcnode=c("b_indir_contact_all"))

quarantine_II<-get_totals(mcmodule=quarantine_II, mcnodes=c("farm_a_no_detect_all_weight"))

quarantine_II<-get_totals(mcmodule=quarantine_II, mcnodes=c("transport_a_contact_all_weight"))

quarantine_II<-add_prefix(quarantine_II)

# 
# ### Combine modules quarantine I and quarantine II
# 
## ------------------------------------------------------------------------------------------------------------
quarantine<-combine_modules(quarantine_I, quarantine_II)


quarantine<-get_totals(mcmodule=quarantine, mcnodes=c("quarantine_I_a_no_detect_all","quarantine_II_a_indir_contact_all"), name = "a_entry_all")

quarantine<-get_totals(mcmodule=quarantine, mcnodes=c("quarantine_I_f_no_detect_all","quarantine_II_f_indir_contact_all"), name = "f_entry_all")


quarantine<-get_totals(mcmodule=quarantine, mcnode=c("a_entry_all"))

quarantine<-get_agg_totals(mcmodule=quarantine, mcnode=c("a_entry_all"))

quarantine<-get_agg_totals(mcmodule=quarantine, mcnode=c("b_entry_all"))

quarantine<-add_prefix(quarantine)

# 
# ### Combine modules purchase and quarantine
# 
## ------------------------------------------------------------------------------------------------------------
purchase<-combine_modules(purchase, quarantine)

# 
# ## Fattening module
# 
# ### Prepare data
# 
## ------------------------------------------------------------------------------------------------------------
#Tidy quarantine parameters
fattening_bsg_data<-tidy_prefix(bsg, module="fattening", rm_prefix = FALSE)

#Add movement and tests info
fattening_data<-fattening_bsg_data%>%
  left_join(purchase_origin_data)%>%
  mutate(surface="soil")%>%
  left_join(pathogen_surface)%>%
  mutate(test=NA)

message("\nfattening_data (",paste(dim(fattening_data), collapse=", "),") created")

# 
# ### Admin what-if
# 
## ------------------------------------------------------------------------------------------------------------
if(admin_wif){
  #Test  before transport
  fattening_by<-quarantine_by
  fattening_data<-group_match(fattening_data, purchase_data, by=fattening_by)
}

# 
# ### Evaluate
# 
## ------------------------------------------------------------------------------------------------------------
fattening_expression<-c(fattening_link=fattening_link_expression,
                     indir_contact_fattening=indir_contact_expression)

fattening<-eval_model_expression(
  model_expression = fattening_expression,
  prev_mcmodule = purchase,
  param_names=c(b_no_detect="farm_b_no_detect",
                b_no_detect_pi="farm_b_no_detect_pi",
                b_contact="farm_b_no_detect_all"),
  data=fattening_data)


fattening<-get_totals(mcmodule=fattening, mcnodes=c("a_indir_contact"))

fattening<-get_agg_totals(mcmodule=fattening, mcnode=c("a_indir_contact"))


fattening<-add_prefix(fattening)

# 
# ### Combine modules purchase and fattening
# 
## ------------------------------------------------------------------------------------------------------------
purchase<-combine_modules(purchase, fattening)

purchase<-get_totals(mcmodule=purchase, mcnodes=c("quarantine_a_entry_all_all","fattening_a_indir_contact_all"), name = "a_entry", data_name = "quarantine")

purchase<-get_agg_totals(mcmodule=purchase, mcnode=c("farm_b_no_detect_all"), data_name = "quarantine")

purchase<-get_agg_totals(mcmodule=purchase, mcnode=c("quarantine_II_farm_b_no_detect_all_weight"), data_name = "quarantine")

purchase<-get_agg_totals(mcmodule=purchase, mcnode=c("transport_b_contact_all"), data_name = "quarantine")

purchase<-get_agg_totals(mcmodule=purchase, mcnode=c("quarantine_II_transport_b_contact_all_weight"), data_name = "quarantine")

purchase<-get_agg_totals(mcmodule=purchase, mcnode=c("purchase_b_inf_all"), data_name = "quarantine")


purchase<-get_agg_totals(mcmodule=purchase, mcnode=c("b_entry"), data_name = "quarantine")

purchase<-get_agg_totals(mcmodule=purchase, mcnode=c("b_entry"), keys_names= c("pathogen", "farm_id","scenario_id", "mov_id"), suffix="mov", data_name = "quarantine")

purchase<-get_agg_totals(mcmodule=purchase, mcnode=c("b_entry"), keys_names= c("pathogen", "farm_id","scenario_id", "mov_id" ,"animal_category"), suffix="mov_animal", data_name = "quarantine")

message("\nPURCHASE PATHWAY OK!")

# 
# ## Save results
# 
}
if(sum(bsg$visit_type_n,bsg$visit_veh_type_n)>0){
# ---
# title: "Farm visits pathogen introduction pathway"
# author: "Natalia Ciria"
# editor: visual
# bibliography: references.bib
# execute:
#   output: false
# ---
# 

# 
# ## Description
# 
## ------------------------------------------------------------------------------------------------------------
message("\nVISITS PATHWAY: ")

# 
# ## Visits module
# 
# ### Prepare data
# 
## ------------------------------------------------------------------------------------------------------------
visit_data<-get_region_code(bsg, pathway = "farm")

visit_veh_data<-visit_data%>%
  left_join(tidy_prefix(tidy_panel(panel="visit_veh"),"visit_veh", rm_prefix = "veh"))%>%
  filter(!is.na(visit_type))%>%
  mutate(visit_type=as.factor(visit_type),
         visit_frequency=as.factor(visit_frequency),
         visit_boots_enter=as.factor(visit_boots_enter),
         visit_boots_cleaning=as.factor(visit_boots_cleaning))


visit_people_data<-visit_data%>%
  left_join(tidy_prefix(tidy_panel(select(bsg,!contains("_veh")), panel="visit"),"visit", rm_prefix = FALSE))%>%
  filter(!is.na(visit_type))%>%
  mutate(visit_type=as.factor(visit_type),
         visit_frequency=as.factor(visit_frequency),
         visit_boots_enter=as.factor(visit_boots_enter),
         visit_boots_cleaning=as.factor(visit_boots_cleaning))



visit_data<-visit_people_data%>%
  bind_rows(visit_veh_data)%>%
  mutate(visit_id=ifelse(is.na(visit_id),visit_veh_id,visit_id),
         #If is people, direct contact is TRUE (implicit in question)
         visit_direct=ifelse(is.na(visit_veh_id),TRUE,visit_direct),
         visit_veh_id=NULL,
         #provisional step to TEST visit link
         prev_b_inf_TEST=0.6)


visit_data<-visit_data%>%
  hjoin_region(visit_region)%>%
  mutate(visit_livestock=ifelse(is.na(visit_livestock_cattle_p),
                                visit_livestock_type_cattle,
                                visit_livestock_cattle_p))

visit_data<-visit_data%>%
  pivot_longer(cols=starts_with("visit_boots")|starts_with("visit_equipment")|starts_with("visit_wheels"),
               names_to = c("fomite_type", ".value"),
               names_pattern = "^visit_(boots|equipment|wheels)_(.*)$",
               names_repair = ~ gsub("^(cleaning|enter)", "visit_\\1", .x))%>%
  filter(!is.na(fomite_type))%>%
  #Remove equipment pathway in vehicles and wheels pathway in people
  filter(!((visit_veh&fomite_type=="equipment")|(!visit_veh&fomite_type=="wheels")))%>%
  #Filter if the fomite never enters the farm (in vehicle wheels always enter)
  filter(!(visit_enter=="never"&visit_direct)|(visit_veh))%>%
  #Filter "boots" if driver never has direct contact
  filter(!(fomite_type=="boots"&!visit_direct))%>%
  #If visit_enter is NA, it is assumed worst case scenario: "always"
  mutate(visit_enter=ifelse(is.na(visit_enter),
                            ifelse(fomite_type=="wheels", "always",
                                   as.character(visit_enter)),
                            as.character(visit_enter)),
         #Add fomite_id
         fomite_id=paste0(visit_id,"_",fomite_type),
         #Unknown pathogen plan status
         status="unk")

#PROVISIONAL for back-compatibility
if(form_version(bsg)<215){
  visit_data<-visit_data%>%
    mutate(visit_close = grepl("veh_driver_help",bsg$complete_comments)
         &(visit_type=="slaughter_veh"|
             visit_type=="calf_veh"|
             visit_type=="adult_veh"))
}

visit_own_exists<-"visit_own"%in%names(visit_data) #BACKWARDS COMPATIBILITY for #141

visit_data<-visit_data%>%
  hjoin_region(pathogen_region)%>%
  left_join(pathogen)%>%
  left_join(pathogen_animal)%>%
  left_join(animal)%>%
  mutate(surface="soil")%>%
  left_join(pathogen_surface)%>%
  #TO SOLVE PROVISIONALLY #141 visit_own = No risk
  mutate(visit_own=if(visit_own_exists) visit_own else NA,
         visit_enter=ifelse(!is.na(visit_own)&visit_own,"never", visit_enter))

# 
# ### Admin what-if
# 
## ------------------------------------------------------------------------------------------------------------
if(admin_wif){
  #Do not allow vehicles to enter in the farm perimeter
  visit_data<-wif_no_fomites(visit_data, fomites="wheels", scenario="Do not allow vehicles to enter in the farm perimeter")
  
  #Provide boots to all visitors
  visit_data<-wif_no_fomites(visit_data, fomites="boots", scenario="Provide boots to all visitors")
  
  #Do not share equipment with other farms
  visit_data<-wif_no_fomites(visit_data, fomites="equipment", scenario="Do not share equipment with other farms")
  
}

# 
# ### Evaluate I
# 
## ------------------------------------------------------------------------------------------------------------
visit_farm_expression<-c(visit_link = visit_link_expression,
                         visit_origin = origin_expression)

visit<-eval_model_expression(model_expression = visit_farm_expression,
                             data=visit_data)

visit<-get_totals(mcmodule=visit, 
                  mcnodes=c("a_inf_infectious","a_inf_pi"),
                  prefix="prev",
                  all_mcnodes=FALSE)

visit<-get_agg_totals(mcmodule=visit, 
                         mcnode=c("prev_b_inf_infectious"),
                         keys_names=c("pathogen", "farm_id", "fomite_id"),
                         agg_variates=TRUE)

visit<-get_agg_totals(mcmodule=visit, 
                         mcnode=c("prev_b_inf_pi"),
                         keys_names=c("pathogen", "farm_id", "fomite_id"),
                         agg_variates=TRUE)

# 
# ### Evaluate II
# 
## ------------------------------------------------------------------------------------------------------------
visit_indir_expression<-c(visit_indir_contact=indir_contact_expression)

visit_indir<-eval_model_expression(model_expression = visit_indir_expression,
                             prev_mcmodule = visit,
                             data=visit_data)

visit<-combine_modules(visit, visit_indir)

visit<-get_totals(mcmodule=visit,
                  mcnodes=c("a_indir_contact","a_indir_contact_pi"))

visit<-get_agg_totals(mcmodule=visit,mcnode=c("a_indir_contact"), keys_names=c("pathogen", "farm_id","scenario_id", "visit_type"), suffix="type")

visit<-get_agg_totals(mcmodule=visit,mcnode=c("a_indir_contact"), keys_names=c("pathogen", "farm_id","scenario_id", "visit_veh"), suffix="veh")

visit<-get_agg_totals(mcmodule=visit, mcnode=c("a_indir_contact_all"))

visit<-add_prefix(visit)

message("\nVISITS PATHWAY OK!")

# 
# ## Save results
# 

# 
# ## Output values
# 
# We are only going to save the probability of disease introduction at aggregated level (one year movements)
}
# ---
# title: "Outdoors introduction pathway"
# author: "Natalia Ciria"
# editor: visual
# bibliography: references.bib
# execute:
#   output: false
# ---
# 

# 
## ------------------------------------------------------------------------------------------------------------
message("\nNEIGHBOURS PATHWAY: ")

# 
# ## Neighbour farms module
# 
# ### Prepare data
# 
## ------------------------------------------------------------------------------------------------------------
neighbour_data<-get_region_code(bsg, pathway = "farm")

neighbour_data<-neighbour_data%>%
  left_join(tidy_prefix(bsg, "neighbour", rm_prefix = FALSE))%>%
  hjoin_region(pathogen_region)%>%
  left_join(pathogen)

# 
# ### Admin what-if
# 
## ------------------------------------------------------------------------------------------------------------
if(admin_wif){
  #Do not allow direct contact with neighbour farms
  neighbour_data<-wif_no_dc(neighbour_data, scenario="Do not allow direct contact with neighbour farms")
  
}

# 
# ### Evaluate
# 
## ------------------------------------------------------------------------------------------------------------
neighbour_expression<-c(neighbour_link = neighbour_link_expression,
                        dir_contact_neighbour=dir_contact_expression,
                        area_inf_neighbour=area_inf_expression)

neighbour<-eval_model_expression(model_expression = neighbour_expression,
                             data=neighbour_data)

neighbour<-at_least_one(mcmodule=neighbour, mcnodes=c("a_dir_contact","area_inf"), name="neighbour_inf")

neighbour<-get_agg_totals(mcmodule=neighbour, mcnode=c("neighbour_inf"))

neighbour<-add_prefix(neighbour)

message("\nNEIGHBOURS PATHWAY OK!")

# 
# ### Save results
# 
if(bsg$outdoors_access&&any(!bsg$outdoors_fencing_perimeter,bsg$outdoors_fencing_wildlife)){
# ---
# title: "Outdoors introduction pathway"
# author: "Natalia Ciria"
# editor: visual
# bibliography: references.bib
# execute:
#   output: false
# ---
# 

# 
## ------------------------------------------------------------------------------------------------------------
message("\nFARM WILDLIFE PATHWAY: ")

# 
# ## Wildlife area effect module
# 
# ### Prepare data
# 
## ------------------------------------------------------------------------------------------------------------
wildlife_area_data<-get_region_code(bsg, pathway = "farm")%>%
  left_join(tidy_prefix(bsg, "outdoors"))

wildlife_area_data<-wildlife_area_data%>%
  mutate(surface="soil")%>%
  left_join(pathogen_surface)%>%
  hjoin_region(wildlife_region)%>%
  hjoin_region(wildlife_pathogen_region)%>%
  left_join(wildlife_point_density)%>%
  hjoin_region(pathogen_region)%>%
  left_join(pathogen)%>%
  #Filter those that have wildlife prevalence
  filter(pathogen%in%
           unique(wildlife_pathogen_region$pathogen[wildlife_pathogen_region$wl_prev_mode>0]))

message("\nwildlife_area_data (",paste(dim(wildlife_area_data), collapse=", "),") created")

# 
# ### Prepare what-if
# 
## ------------------------------------------------------------------------------------------------------------
if(admin_wif){
  #There are no admin what-ifs in this module, but we set up wif for compatibility
  wildlife_area_data<-set_up_wif(wildlife_area_data)
}

# 
# ### Evaluate
# 
## ------------------------------------------------------------------------------------------------------------
wildlife_area_expression<-c(wildlife_area_link = wildlife_area_link_expression,                         area_inf_wildlife_area=area_inf_expression)  

wildlife<-eval_model_expression(model_expression = wildlife_area_expression,                              data=wildlife_area_data)  

wildlife<-get_agg_totals(mcmodule=wildlife, mcnode=c("area_inf"))

wildlife<-add_prefix(wildlife)

# 
# ## Wildlife contact points module
# 
# ### Prepare data
# 
## ------------------------------------------------------------------------------------------------------------
farm_census_data<-tidy_prefix(bsg, "farm_census")%>%
  pivot_longer(
      cols = starts_with("cattle"),
      names_to = "animal_category",
      values_to = "census",
      values_drop_na = TRUE
    )%>%
  left_join(animal)%>%
  mutate(livestock_units=census*livestock_units)

wildlife_water_data<-wildlife_area_data%>%
  pivot_longer(
    cols = ends_with("_p")&!mud_p&!sum_p,
    names_to = "waterpoint_type",
    values_to = "point_p",
    values_drop_na = TRUE
  )%>%
  mutate(
    animals_n=sum(farm_census_data$census),
    livestock_units=sum(farm_census_data$livestock_units),
    contact_point_type=sub("_p", "", waterpoint_type),
    mud_p=ifelse(contact_point_type%in%c("low", "high"),ifelse(mud,as.numeric(gsub("%", "", mud_p)),0),100)/100,
    contact_point_p=as.numeric(gsub("%", "", point_p))/100
  )

wildlife_water_data<-wildlife_water_data%>%
  left_join(wildlife_point)%>%
  left_join(animal_point)

message("\nwildlife_water_data (",paste(dim(wildlife_water_data), collapse=", "),") created")

# 
# ### Admin what-if
# 
## ------------------------------------------------------------------------------------------------------------
if(admin_wif){
  #Avoid mud on wateres
  wildlife_water_data<-wif_no_mud(wildlife_water_data, scenario="Avoid mud on wateres")
  wildlife_water_data<-wif_high_waterer(wildlife_water_data, scenario="Use high wateres")
}

# 
# ### Evaluate
# 
## ------------------------------------------------------------------------------------------------------------
wildlife_water_expression<-c(
  wildlife_link = wildlife_water_link_expression,
                            time_wildlife = time_expression,
                            indir_contact_wildlife=indir_contact_expression)

#debugonce(get_node_list)
wildlife_water<-eval_model_expression(
  model_expression = wildlife_water_expression,
  data=wildlife_water_data)

wildlife_water<-get_totals(mcmodule=wildlife_water, mcnode=c("a_indir_contact"))

wildlife_water<-get_agg_totals(mcmodule=wildlife_water, mcnode=c("b_indir_contact"))

wildlife_water<-add_prefix(wildlife_water)

# 
# ## Combine modules
# 
## ------------------------------------------------------------------------------------------------------------
wildlife<-combine_modules(wildlife, wildlife_water)

wildlife<-at_least_one(mcmodule=wildlife, mcnodes=c("wildlife_area_inf_agg","wildlife_water_b_indir_contact_agg"), name="wildlife_inf_agg")

message("\nFARM WILDLIFE PATHWAY OK!")

# 
# ## Save results
# 
}
# ---
# title: "Export outputs"
# author: "Natalia Ciria"
# editor: visual
# bibliography: references.bib
# execute:
#   output: false
# ---
# 

# 
# ### Combine all pathways
# 
## ------------------------------------------------------------------------------------------------------------
message("\nOUTPUTS: ")

# 
## ------------------------------------------------------------------------------------------------------------
message("\nMerging mcmodules: ")
#rm(wildlife)
#message("\nWILDLIFE MODULE REMOVED UNTIL #143 BUG IS SOLVED")
if(exists("wildlife")){
  message("\n- wildlife: YES")
  surroundings<-combine_modules(neighbour, wildlife)
  
  surroundings<-at_least_one(mcmodule=surroundings, mcnodes=c("neighbour_inf_agg","wildlife_inf_agg"), name="surroundings_inf_agg")
  
  outputs_surroundings<-c(neighbour_total="neighbour_inf_agg",
                          wildlife_total="wildlife_inf_agg")
  
  outputs_hg_surroundings<-c(farm_neigbour="neighbour_inf",
                          farm_wildlife_water="wildlife_water_b_indir_contact",
                          farm_wildlife_area= "wildlife_area_inf")

}else{
  message("\n- wildlife: NO")

  surroundings<-neighbour
  
  surroundings$node_list$surroundings_inf_agg<-
    surroundings$node_list$neighbour_inf_agg
  
  outputs_surroundings<-c(neighbour_total="neighbour_inf_agg")
  outputs_hg_surroundings<-c(farm_neigbour="neighbour_inf")

}

if(exists("visit")){
  message("\n- visit: YES")

  no_purchase<-combine_modules(surroundings, visit)
  
  no_purchase<-at_least_one(mcmodule=no_purchase, mcnodes=c("visit_a_indir_contact_all_agg","surroundings_inf_agg"), name="no_purchase_inf_agg")
  

  outputs_visits<-c(visit_total="visit_a_indir_contact_all_agg")

  outputs_hg_visits<-c(farm_visits="visit_a_indir_contact_all")
  
  outputs_agg_visits<-c(visit_type="visit_a_indir_contact_type",
                        visit_veh="visit_a_indir_contact_veh")
  
}else{
  message("\n- visit: NO")
  
  no_purchase<-surroundings
  
  no_purchase$node_list$no_purchase_inf_agg<-
    surroundings$node_list$surroundings_inf_agg
  
  outputs_visits<-c()
  outputs_hg_visits<-c()
  outputs_agg_visits<-c()

}


#If farm has purchase movements
if(exists("purchase")){
  message("\n- purchase: YES")

  no_pasture<-combine_modules(purchase, no_purchase)

  no_pasture<-at_least_one(mcmodule=no_pasture, mcnodes=c("no_purchase_inf_agg","b_entry_agg"), name="no_pasture_agg")
    
  outputs_purchase<-c(purchase_origin_prequaran="farm_b_no_detect_all_agg",
           purchase_veh_prequaran="transport_b_contact_all_agg",
           purchase_prequaran="purchase_b_inf_all_agg",
           purchase_total="b_entry_agg",
           purchase_origin="quarantine_II_farm_b_no_detect_all_weight_agg",
           purchase_veh="quarantine_II_transport_b_contact_all_weight_agg")
  
  outputs_hg_purchase<-c(
           purchase_origin="quarantine_II_farm_b_no_detect_all_weight",
           purchase_veh="quarantine_II_transport_b_contact_all_weight")
  
  outputs_agg_purchase<-c(purchase_by_mov="b_entry_mov",
                          purchase_by_mov_animal="b_entry_mov_animal")
  
}else{
  message("\n- purchase: NO")
  
  no_pasture<-no_purchase
  
  no_pasture$node_list[["no_pasture_agg"]]<-no_pasture$node_list$no_purchase_inf_agg
  
  outputs_purchase<-c()
  
  outputs_hg_purchase<-c()
  
  outputs_agg_purchase<-c()
}

#If farm has pasture movements
if(exists("pasture")){
  message("\n- pasture: YES")
  
  intro<-combine_modules(pasture, no_pasture)

  intro<-at_least_one(mcmodule=intro, mcnodes=c("no_pasture_agg","pasture_inf_agg"), name="total_agg")
  
  outputs_pasture<-c(pasture_mix="pasture_mix_b_dir_contact_all_agg",
           pasture_veh="pasture_transport_b_contact_all_agg",
           pasture_livestock="pasture_b_livestock_all_agg",
           pasture_wildlife="pasture_wildlife_inf_agg",
           pasture_total="pasture_inf_agg")
  
  outputs_hg_pasture<-c(pasture_mix="pasture_mix_b_dir_contact_all",
         pasture_veh="pasture_transport_b_contact_all")
  
  outputs_agg_pasture<-c(pasture_livestock_by_mov="pasture_b_livestocK_all_mov",
                         pasture_wildlife_by_mov="pasture_wildlife_all_mov",
                         pasture_by_mov="pasture_all_mov")

  
  if(exists("purchase")){
      intro<-at_least_one(mcmodule=intro, mcnodes=c("b_entry_agg","pasture_inf_agg"), name="mov_total_agg")
  }
  
}else{
  message("\n- pasture: NO")
  
  intro<-no_pasture
  
  intro$node_list[["total_agg"]]<-intro$node_list$no_pasture_agg
  
  outputs_pasture<-c()
  outputs_hg_pasture<-c()
  outputs_agg_pasture<-c()

}

#Reassure class
class(intro)<-"mcmodule"

# 

# 
# ### Select result mcnodes
# 
## ------------------------------------------------------------------------------------------------------------
outputs_all<-c(neighbour_total="neighbour_inf_agg",
               wildlife_total="wildlife_inf_agg",
               visit_total="visit_a_indir_contact_all_agg",
               purchase_origin="farm_b_no_detect_all_agg",
               purchase_veh="transport_b_contact_all_agg",
               purchase_prequaran="purchase_b_inf_all_agg",
               purchase_total="b_entry_agg",
               no_purchase_total="no_purchase_inf_agg",
               pasture_mix="pasture_mix_b_dir_contact_all_agg",
               pasture_veh="pasture_transport_b_contact_all_agg",
               #pasture_wildlife="purchase_b_inf_all_agg",
               pasture_total="pasture_b_livestock_all_agg",
               total="total_agg")

outputs_agg<-c(outputs_agg_visits,outputs_agg_purchase,outputs_agg_pasture)


outputs<-c(outputs_pasture,
           outputs_purchase,
           outputs_visits,
           outputs_surroundings,
           no_purchase_total="no_purchase_inf_agg",
           no_pasture_total="no_pasture_agg",
           total="total_agg")

# 

# 
# ### Output by hg (all simulations)
# 
## ------------------------------------------------------------------------------------------------------------
if(exists("hg_csv")&&hg_csv){
  outputs_hg<-c(outputs_hg_pasture,
           outputs_hg_purchase,
           outputs_hg_visits,
           outputs_hg_surroundings)
  
  hg_df<-data.frame()
for(i in 1:length(outputs_hg)){
  hg_df_i<-mc_to_long_df(intro, outputs_hg[i], hg_keys=TRUE)
  hg_df_i$mcnode<-outputs_hg[i]
  hg_df_i$pathway<-names(outputs_hg[i])
  hg_df_i$farm_id<-farm_id
  hg_df_i$timestamp<-Sys.time()
  hg_df<-bind_rows(hg_df,hg_df_i)
}
  
  write.table(hg_df, file="output_files/hg_df.csv", row.names=FALSE, col.names=!file.exists("output_files/hg_df.csv"), append=append_csv, sep=";")
  
  message(nrow(hg_df), "hg simulations of farm ", farm_id," added to output_files/hg_df.csv")
}

# 
# ### Output medians
# 
## ------------------------------------------------------------------------------------------------------------

#Create a table with the mean of the selected mcnodes
wif_median<-intro$node_list[[outputs["total"]]][["summary"]]%>%
  select(pathogen, scenario_id, median = `50%`)

names(wif_median)[names(wif_median)=="median"]<-names(outputs["total"])



for(i in 1:length(outputs[!names(outputs)=="total"])){
  wif_median_i<-intro$node_list[[outputs[i]]][["summary"]]%>%
    select(pathogen, scenario_id, median = `50%`)
  names(wif_median_i)[names(wif_median_i)=="median"]<-names(outputs[i])
  
  wif_median<-left_join(wif_median, wif_median_i, by = c("pathogen", "scenario_id"))
}

#Add outputs that have not been calculated (probability equals zero)
wif_median[names(outputs_all[!outputs_all%in%outputs])]<-0
wif_median[is.na(wif_median)]<-0



#Create a table with the mean of the selected mcnodes
wif_median<-intro$node_list[[outputs["total"]]][["summary"]]%>%
  select(pathogen, scenario_id, median = `50%`)

names(wif_median)[names(wif_median)=="median"]<-names(outputs["total"])



for(i in 1:length(outputs[!names(outputs)=="total"])){
  wif_median_i<-intro$node_list[[outputs[i]]][["summary"]]%>%
  select(pathogen, scenario_id, median = `50%`)
  names(wif_median_i)[names(wif_median_i)=="median"]<-names(outputs[i])

  wif_median<-left_join(wif_median, wif_median_i, by = c("pathogen", "scenario_id"))
}

#Add outputs that have not been calculated (probability equals zero)
wif_median[names(outputs_all[!outputs_all%in%outputs])]<-0
wif_median[is.na(wif_median)]<-0

# 
# #### Output by scenario
# 
## ------------------------------------------------------------------------------------------------------------
#Current risk table
current_median <- wif_median%>%
  filter(scenario_id=="0")%>%
  mutate(scenario_id=NULL)

#Calculate relative risk
wif_median<-wif_median%>%
  left_join(current_median[c("pathogen","total")], by="pathogen",suffix=c("_wif","_current"))%>%
  mutate(
    relative=total_wif/total_current,
    relative=ifelse(relative>=1&scenario_id!="Current",
                             1-(1*10^-6),
                             relative),
  )

#Create table to plot risk reduction
wif_median_plot<-wif_median%>%
  select(pathogen,scenario_id,total=relative)%>%
  mutate(yes=total,
         no=1-total)%>%
  pivot_longer(
    cols = c(yes, no), 
    names_to = "yesno", 
    values_to = "value",
  )


#Current risk table
current_median <- wif_median%>%
  mutate(total=total_wif)%>%
  select(pathogen, scenario_id, total, 
         pasture_mix,pasture_veh,
         purchase_origin,purchase_veh,
         visit_total, neighbour_total,wildlife_total)%>%
  pivot_longer(cols=-c(pathogen, scenario_id, total))%>%
  filter(scenario_id=="0")%>%
  transmute(pathogen,
            pathway=name,
            abs_total=total,
            rel_pathway=value/total,
            abs_pathway=value)
         

# 
## ------------------------------------------------------------------------------------------------------------
summary_median<-wif_median%>%
  arrange(relative)%>%
  transmute(Pathogen=pathogen,
            Scenario=ifelse(scenario_id=="0","Current",scenario_id),
         `Risk disease entry`=paste0(Signif(total_wif*100,2),"%"),
         `Risk animal infected from origin`=paste0(Signif(purchase_origin*100,2),"%"),
         `Risk animal infected during purchase transport`=paste0(Signif(purchase_veh*100,2),"%"),
         `Risk disease entry by visits`=paste0(Signif(visit_total*100,2),"%"),
         `Risk disease entry by neigbour farms`=paste0(Signif(neighbour_total*100,2),"%"),
          `Risk disease entry by wildlife near the farm`=paste0(Signif(wildlife_total*100,2),"%"),
         `Risk disease entry by movement to pastures`=paste0(Signif(pasture_total*100,2),"%"))

summary_median[summary_median=="NA%"]<-"-"

# 
# #### Save summary tables
# 
## ------------------------------------------------------------------------------------------------------------
save_table_summary(list(wif_median=wif_median, 
                        wif_median_plot=wif_median_plot, 
                        current_median=current_median,
                        summary_median=summary_median))

# 
# ### Print-cat summary
# 
## ------------------------------------------------------------------------------------------------------------
cat("\nResults")
print_summary<-unique(summary_median[summary_median$Scenario=="Current",
                                     !names(summary_median)%in%"Scenario"])
cat("\n\nIBR current median:\n")
print.table(print_summary[print_summary$Pathogen=="IBR",
                          !names(print_summary)%in%"Pathogen"])

cat("\n\nBVD current median:\n")
print.table(print_summary[print_summary$Pathogen=="BVD",
                          !names(print_summary)%in%"Pathogen"])

cat("\n\nTB current median:\n")
print.table(print_summary[print_summary$Pathogen=="TB",
                          !names(print_summary)%in%"Pathogen"])

assign("intro",intro,envir=parent.frame())
}
