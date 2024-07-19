################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#First use command prompt to log in to Slade using SSH

#To install package
#library(devtools)
#install_github("Exeter-Diabetes/CPRD-analysis-package")


#Load tidyverse and lubridate
library(tidyverse)

#Load aurum package
library(aurum)

################################################################################
###SETUP########################################################################

###Connecting to data and setting up analysis###################################
#Initialise connection
cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "C:/Users/rh530/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "31/10/2021")

#Connect to 'all' analysis to get all diabetes cohort table
#analysis = cprd$analysis("all")
#all_diabetes <- all_diabetes %>% analysis$cached("diabetes_cohort", unique_indexes="patid",indexes=c("gender", "dob", "dm_diag_date_all", "dm_diag_age_all", "diabetes_type"))
#Date copied: 07/06/2023

#Connect to 'Rhian_T3c' analysis
analysis = cprd$analysis("Rhian_T3c")

#Save all diabetes table
all_diabetes <- all_diabetes %>% analysis$cached("all_diabetes_cohort", unique_indexes="patid",indexes=c("gender", "dob", "dm_diag_date_all", "dm_diag_age_all", "diabetes_type"))
all_diabetes %>% count() #1138179
#Drop people without diagnosis date (missing if this date is within -30 to +90 days (inclusive) of registration start)
all_diabetes <- all_diabetes %>% filter(!is.na(dm_diag_date_all))
all_diabetes %>% count() #1077649

################################################################################
###1: GET DATA ON PANCREATIC CONDITIONS FOR WHOLE DIABETES COHORT
################################################################################

###Get raw data tables
##Primary care
acutepancreatitis_gp_raw <- cprd$tables$observation %>% inner_join(codes$acutepancreatitis) %>% analysis$cached("all_acutepancreatitis_raw", indexes= c("patid", "obsdate"))
chronicpancreatitis_gp_raw <- cprd$tables$observation %>% inner_join(codes$chronicpancreatitis) %>% analysis$cached("all_chronicpancreatitis_raw", indexes= c("patid", "obsdate"))
haemochromatosis_gp_raw <- cprd$tables$observation %>% inner_join(codes$haemochromatosis_spocc) %>% analysis$cached("all_haemochromatosis_raw", indexes= c("patid", "obsdate"))
cysticfibrosis_gp_raw <- cprd$tables$observation %>% inner_join(codes$cysticfibrosis) %>% analysis$cached("all_cysticfibrosis_raw", indexes= c("patid", "obsdate"))
pancreaticcancer_gp_raw <- cprd$tables$observation %>% inner_join(codes$pancreaticcancer_spocc) %>% analysis$cached("all_pancreaticcancer_raw", indexes= c("patid", "obsdate"))
surgicalpancreaticresection_gp_raw <- cprd$tables$observation %>% inner_join(codes$surgicalpancreaticresection) %>% analysis$cached("all_surgicalpancreaticresection_raw", indexes= c("patid", "obsdate"))

##HES
acutepancreatitis_hes_raw <- cprd$tables$hesDiagnosisEpi %>% inner_join(codes$icd10_acutepancreatitis, sql_on="LHS.ICD LIKE CONCAT(icd10,'%')") %>% analysis$cached("all_acutepancreatitis_hes", indexes= c("patid", "epistart"))
chronicpancreatitis_hes_raw <- cprd$tables$hesDiagnosisEpi %>% inner_join(codes$icd10_chronicpancreatitis, sql_on="LHS.ICD LIKE CONCAT(icd10,'%')") %>% analysis$cached("all_chronicpancreatitis_hes", indexes= c("patid", "epistart"))
haemochromatosis_hes_raw <- cprd$tables$hesDiagnosisEpi %>% inner_join(codes$icd10_haemochromatosis, sql_on="LHS.ICD LIKE CONCAT(icd10,'%')") %>% analysis$cached("all_haemochromatosis_hes", indexes= c("patid", "epistart"))
cysticfibrosis_hes_raw <- cprd$tables$hesDiagnosisEpi %>% inner_join(codes$icd10_cysticfibrosis, sql_on="LHS.ICD LIKE CONCAT(icd10,'%')") %>% analysis$cached("all_cysticfibrosis_hes", indexes= c("patid", "epistart"))
pancreaticcancer_hes_raw <- cprd$tables$hesDiagnosisEpi %>% inner_join(codes$icd10_pancreaticcancer, sql_on="LHS.ICD LIKE CONCAT(icd10,'%')") %>% analysis$cached("all_pancreaticcancer_hes", indexes= c("patid", "epistart"))
surgicalpancreaticresection_hes_proc_raw <- cprd$tables$hesProceduresEpi %>% inner_join(codes$opcs4_surgicalpancreaticresection, by=c("OPCS"="opcs4")) %>% analysis$cached("all_surgicalpancreaticresection_hes", indexes= c("patid", "evdate"))

################################################################################
#Pancreatic conditions
panc_conditions <-c("acutepancreatitis", "chronicpancreatitis", "pancreaticcancer", "haemochromatosis", "cysticfibrosis", "surgicalpancreaticresection")

#First date of condition prior to (and including) date of diabetes diagnosis
for (i in panc_conditions) {
  #Primary care
  print(paste0(i, " GP"))
  raw_gp <- paste0(i,"_gp_raw")
  pre_diag_gp_tablename <- paste0("pre_diag_",i,"_gp")
  pre_diag_gp <- get(raw_gp) %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>% inner_join(all_diabetes %>% select(patid, dm_diag_date_all) %>% distinct()) %>%
    filter(obsdate <= dm_diag_date_all) %>% group_by(patid) %>% summarise(first_date =min(obsdate,na.rm=TRUE)) %>% ungroup()
  assign(pre_diag_gp_tablename,pre_diag_gp)
  print(pre_diag_gp %>% count())
  #HES
  raw_hes <- paste0(i,"_hes_raw")
  raw_hes_op <- paste0(i,"_hes_proc_raw")
  pre_diag_hes_tablename <- paste0("pre_diag_",i,"_hes")
  if(exists(raw_hes)) {
    print(paste0(i, " HES diagnosis"))  
    pre_diag_hes <- get(raw_hes) %>% inner_join(cprd$tables$validDateLookup) %>% filter(epistart>=min_dob) %>% inner_join(all_diabetes %>% select(patid, dm_diag_date_all) %>% distinct()) %>%
      filter(epistart <= dm_diag_date_all) %>% group_by(patid) %>% summarise(first_date =min(epistart,na.rm=TRUE)) %>% ungroup()
    assign(pre_diag_hes_tablename,pre_diag_hes)
    print(pre_diag_hes %>% count())
  }
  if(exists(raw_hes_op)) {
    print(paste0(i, " HES procedure")) 
    pre_diag_hes <- get(raw_hes_op) %>% inner_join(cprd$tables$validDateLookup) %>% filter(evdate>=min_dob) %>% inner_join(all_diabetes %>% select(patid, dm_diag_date_all) %>% distinct()) %>%
      filter(evdate <= dm_diag_date_all) %>% group_by(patid) %>% summarise(first_date =min(evdate,na.rm=TRUE)) %>% ungroup()
    assign(pre_diag_hes_tablename,pre_diag_hes)
    print(pre_diag_hes %>% count())
  }
  #Combine primary care and HES, find earliest date
  pre_diag_combined_tablename <- paste0("pre_diag_",i)
  pre_diag_combined <- pre_diag_gp %>% union(pre_diag_hes) %>% group_by(patid) %>% summarise(pre_diag_date = min(first_date, na.rm = TRUE))
  assign(pre_diag_combined_tablename,pre_diag_combined)
  print(paste0("Combined GP and HES dates for ", i))
  print(pre_diag_combined %>% count())
}

pre_diag_pancreatic_disease_all <- pre_diag_acutepancreatitis %>% union(pre_diag_chronicpancreatitis) %>% union(pre_diag_pancreaticcancer) %>%
  union(pre_diag_haemochromatosis) %>% union(pre_diag_surgicalpancreaticresection) %>% union(pre_diag_cysticfibrosis) %>%
  group_by(patid) %>% summarise(pre_diag_date = min(pre_diag_date, na.rm = TRUE))

#Latest date of condition prior to (and including) date of diabetes diagnosis
for (i in panc_conditions) {
  #Primary care
  print(paste0(i, " GP"))
  raw_gp <- paste0(i,"_gp_raw")
  pre_diag_last_gp_tablename <- paste0("pre_diag_last_",i,"_gp")
  pre_diag_last_gp <- get(raw_gp) %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>% inner_join(all_diabetes %>% select(patid, dm_diag_date_all) %>% distinct()) %>%
    filter(obsdate <= dm_diag_date_all) %>% group_by(patid) %>% summarise(last_date =max(obsdate,na.rm=TRUE)) %>% ungroup()
  assign(pre_diag_last_gp_tablename,pre_diag_last_gp)
  print(pre_diag_last_gp %>% count())
  #HES
  raw_hes <- paste0(i,"_hes_raw")
  raw_hes_op <- paste0(i,"_hes_proc_raw")
  pre_diag_last_hes_tablename <- paste0("pre_diag_last_",i,"_hes")
  if(exists(raw_hes)) {
    print(paste0(i, " HES diagnosis"))  
    pre_diag_last_hes <- get(raw_hes) %>% inner_join(cprd$tables$validDateLookup) %>% filter(epistart>=min_dob) %>% inner_join(all_diabetes %>% select(patid, dm_diag_date_all) %>% distinct()) %>%
      filter(epistart <= dm_diag_date_all) %>% group_by(patid) %>% summarise(last_date =max(epistart,na.rm=TRUE)) %>% ungroup()
    assign(pre_diag_last_hes_tablename,pre_diag_last_hes)
    print(pre_diag_last_hes %>% count())
  }
  if(exists(raw_hes_op)) {
    print(paste0(i, " HES procedure")) 
    pre_diag_last_hes <- get(raw_hes_op) %>% inner_join(cprd$tables$validDateLookup) %>% filter(evdate>=min_dob) %>% inner_join(all_diabetes %>% select(patid, dm_diag_date_all) %>% distinct()) %>%
      filter(evdate <= dm_diag_date_all) %>% group_by(patid) %>% summarise(last_date =max(evdate,na.rm=TRUE)) %>% ungroup()
    assign(pre_diag_last_hes_tablename,pre_diag_last_hes)
    print(pre_diag_last_hes %>% count())
  }
  #Combine primary care and HES, find latest date
  pre_diag_last_combined_tablename <- paste0("pre_diag_last_",i)
  pre_diag_last_combined <- pre_diag_last_gp %>% union(pre_diag_last_hes) %>% group_by(patid) %>% summarise(pre_diag_last_date = max(last_date, na.rm = TRUE))
  assign(pre_diag_last_combined_tablename,pre_diag_last_combined)
  print(paste0("Combined GP and HES dates for ", i))
  print(pre_diag_last_combined %>% count())
}

pre_diag_last_pancreatic_disease_all <- pre_diag_last_acutepancreatitis %>% union(pre_diag_last_chronicpancreatitis) %>% union(pre_diag_last_pancreaticcancer) %>%
  union(pre_diag_last_haemochromatosis) %>% union(pre_diag_last_surgicalpancreaticresection) %>% union(pre_diag_last_cysticfibrosis) %>%
  group_by(patid) %>% summarise(pre_diag_last_date = max(pre_diag_last_date, na.rm = TRUE))


#First date of condition following date of diabetes diagnosis
for (i in panc_conditions) {
  #Primary care
  print(paste0(i, " GP"))
  raw_gp <- paste0(i,"_gp_raw")
  post_diag_gp_tablename <- paste0("post_diag_",i,"_gp")
  post_diag_gp <- get(raw_gp) %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>% inner_join(all_diabetes %>% select(patid, dm_diag_date_all) %>% distinct()) %>%
    filter(obsdate > dm_diag_date_all) %>% group_by(patid) %>% summarise(first_date =min(obsdate,na.rm=TRUE)) %>% ungroup()
  assign(post_diag_gp_tablename,post_diag_gp)
  print(post_diag_gp %>% count())
  #HES
  raw_hes <- paste0(i,"_hes_raw")
  raw_hes_op <- paste0(i,"_hes_proc_raw")
  post_diag_hes_tablename <- paste0("post_diag_",i,"_hes")
  if(exists(raw_hes)) {
    print(paste0(i, " HES diagnosis"))  
    post_diag_hes <- get(raw_hes) %>% inner_join(cprd$tables$validDateLookup) %>% filter(epistart>=min_dob) %>% inner_join(all_diabetes %>% select(patid, dm_diag_date_all) %>% distinct()) %>%
      filter(epistart > dm_diag_date_all) %>% group_by(patid) %>% summarise(first_date =min(epistart,na.rm=TRUE)) %>% ungroup()
    assign(post_diag_hes_tablename,post_diag_hes)
    print(post_diag_hes %>% count())
  }
  if(exists(raw_hes_op)) {
    print(paste0(i, " HES procedure")) 
    post_diag_hes <- get(raw_hes_op) %>% inner_join(cprd$tables$validDateLookup) %>% filter(evdate>=min_dob) %>% inner_join(all_diabetes %>% select(patid, dm_diag_date_all) %>% distinct()) %>%
      filter(evdate > dm_diag_date_all) %>% group_by(patid) %>% summarise(first_date =min(evdate,na.rm=TRUE)) %>% ungroup()
    assign(post_diag_hes_tablename,post_diag_hes)
    print(post_diag_hes %>% count())
  }
  #Combine primary care and HES, find ealiest date after diabetes diagnosis
  post_diag_combined_tablename <- paste0("post_diag_",i)
  post_diag_combined <- post_diag_gp %>% union(post_diag_hes) %>% group_by(patid) %>% summarise(post_diag_date = min(first_date, na.rm = TRUE))
  assign(post_diag_combined_tablename,post_diag_combined)
  print(paste0("Combined GP and HES dates for ", i))
  print(post_diag_combined %>% count())
}

post_diag_pancreatic_disease_all <- post_diag_acutepancreatitis %>% union(post_diag_chronicpancreatitis) %>% union(post_diag_pancreaticcancer) %>%
  union(post_diag_haemochromatosis) %>% union(post_diag_surgicalpancreaticresection) %>% union(post_diag_cysticfibrosis) %>%
  group_by(patid) %>% summarise(post_diag_date = min(post_diag_date, na.rm = TRUE))

################################################################################
##Joining with cohort

cohort_pancreatic_disease <- all_diabetes %>% select(patid, dm_diag_date_all) %>% distinct()
panc_conditions <- c(panc_conditions, "pancreatic_disease_all")

for(i in panc_conditions) {
  print(i)
  pre_diag_tablename <- paste0("pre_diag_",i)
  pre_diag_last_tablename <- paste0("pre_diag_last_",i)
  post_diag_tablename <- paste0("post_diag_",i)
  pre_diag_variablename <- paste0("first_",i,"_pre_diagnosis")
  pre_diag_last_variablename <- paste0("last_",i,"_pre_diagnosis")
  post_diag_variablename <- paste0("first_",i, "_post_diagnosis")
  prior_variablename <- paste0("prior_", i)
  interim_table <- paste0("all_interim_", i)
  pre_diag <- get(pre_diag_tablename)
  post_diag <- get(post_diag_tablename)
  pre_diag_last <- get(pre_diag_last_tablename)
  cohort_pancreatic_disease <- cohort_pancreatic_disease %>% left_join(pre_diag) %>% left_join(pre_diag_last) %>% left_join(post_diag) %>% 
    mutate(prior = ifelse(!is.na(pre_diag_date) | (datediff(post_diag_date, dm_diag_date_all)) <=30, 1L, 0L)) %>%
    mutate(prior = ifelse(is.na(prior), 0L, prior)) %>% rename({{pre_diag_variablename}}:= pre_diag_date) %>% rename({{pre_diag_last_variablename}}:= pre_diag_last_date) %>% 
    rename({{post_diag_variablename}}:= post_diag_date) %>% rename({{prior_variablename}}:= prior) %>% 
    analysis$cached(interim_table, unique_indexes = "patid")
}  

diabetes_pancreatic_disease  <- all_diabetes %>% left_join(cohort_pancreatic_disease, by = c("patid", "dm_diag_date_all")) %>% 
  mutate(prior_acutepancreatitis_only = ifelse(prior_acutepancreatitis ==1 & prior_chronicpancreatitis ==0 & prior_pancreaticcancer ==0 & prior_haemochromatosis ==0 & prior_surgicalpancreaticresection ==0, 1L, 0L)) %>%
  analysis$cached("all_diabetes_pancreaticdisease", unique_indexes = "patid")

diabetes_pancreatic_disease %>% count(prior_acutepancreatitis ==1) #12494
diabetes_pancreatic_disease %>% count(prior_acutepancreatitis_only ==1) #8358
diabetes_pancreatic_disease %>% count(prior_chronicpancreatitis ==1) #5661
diabetes_pancreatic_disease %>% count(prior_pancreaticcancer ==1) #825
diabetes_pancreatic_disease %>% count(prior_haemochromatosis ==1) #935
diabetes_pancreatic_disease %>% count(prior_cysticfibrosis ==1) #831 - want to exclude these people!
diabetes_pancreatic_disease %>% count(prior_surgicalpancreaticresection ==1) #1352

################################################################################
###2: REFINE COHORT
# 'type' variable:
# T3c = people with type 1 or type 2 diabetes with prior pancreatic disease OR people with secondary (non-drug induced) diabetes and with prior pancreatic disease
# T2 = people with type 2 diabetes and no prior pancreatic disease (T2 controls)
# T1 = people with type 1 diabetes and no prior pancreatic disease
################################################################################

raw_exclusion_diabetes_medcodes <- cprd$tables$observation %>% inner_join(codes$exclusion_diabetes, by="medcodeid") %>% analysis$cached("raw_exclusion_diabetes_medcodes", indexes=c("patid", "obsdate", "exclusion_diabetes_cat"))

#Find people who have code for secondary (non drug induced) diabetes
secondary_diabetes_patids <- raw_exclusion_diabetes_medcodes %>% filter(medcodeid == "356078011" | medcodeid == "189721000000113" | medcodeid == "198461000000116" | medcodeid == "8120401000006111" | medcodeid == "15518018" | medcodeid == "1230929011" | medcodeid == "967681000006119" | medcodeid == "2474730014") %>%
  distinct(patid) %>% mutate(secondary_diabetes_non_drug =1L)

#Pull together T3c/T1/T2 cohort
t3c_t2controls <- diabetes_pancreatic_disease %>% left_join(secondary_diabetes_patids, by = "patid") %>% mutate(secondary_diabetes_non_drug = ifelse(is.na(secondary_diabetes_non_drug), 0L, secondary_diabetes_non_drug)) %>%
  #Drop people with other types of diabetes and cystic fibrosis ever
  filter(diabetes_type == "type 2" | diabetes_type == "type 1" | (secondary_diabetes_non_drug ==1 & prior_pancreatic_disease_all ==1)) %>% filter(is.na(first_cysticfibrosis_pre_diagnosis) & is.na(first_cysticfibrosis_post_diagnosis)) %>% 
  #Define 'type' variable: T1, T2, or T3c
  mutate(type = ifelse(diabetes_type == "type 2" & prior_pancreatic_disease_all ==0, "T2", ifelse(diabetes_type == "type 1" & prior_pancreatic_disease_all ==0, "T1", ifelse(prior_pancreatic_disease_all ==1, "T3c", NA)))) %>%
  #Drop surgical resection only group (still want to keep those with surgery and pancreatic cancer etc)
  filter(type == "T2" | type == "T1" | prior_acutepancreatitis_only ==1 | prior_chronicpancreatitis ==1 | prior_pancreaticcancer ==1 | prior_haemochromatosis ==1) %>%
  #Define binary grouping variables for T3c subtypes using hierachy acute pancreatitis < haemochromatosis < chronic pancreatitis < pancreatic cancer
  mutate(t3c_all = prior_pancreatic_disease_all, t3c_acutepancreatitis_only = prior_acutepancreatitis_only, t3c_pancreaticcancer = prior_pancreaticcancer, t3c_chronicpancreatitis = ifelse(prior_chronicpancreatitis ==1 & prior_pancreaticcancer ==0, 1L, 0L), t3c_haemochromatosis = ifelse(prior_haemochromatosis ==1 & prior_chronicpancreatitis ==0 & prior_pancreaticcancer ==0, 1L, 0L)) %>%
  #Define a categorical variable of T3c subgroups
  mutate(t3c_subgroup = ifelse(t3c_acutepancreatitis_only ==1, "Acute pancreatitis only", ifelse(t3c_pancreaticcancer ==1, "Pancreatic cancer", ifelse(t3c_chronicpancreatitis ==1, "Chronic pancreatitis", ifelse(t3c_haemochromatosis ==1, "Haemochromatosis", ifelse(t3c_all ==0, "Type 2", NA)))))) %>%
  analysis$cached("cohort", unique_indexes = "patid")

#Count numbers
t3c_t2controls %>% count(type) #T2: 971019, T1: 73057, T3c: 15421

t3c_t2controls %>% count(prior_pancreatic_disease_all) #0 = 1044076, 1 = 15421

t3c_t2controls %>% filter(t3c_all ==1) %>% count() #15421
t3c_t2controls %>% filter(t3c_acutepancreatitis_only ==1) %>% count() #8242
t3c_t2controls %>% filter(t3c_chronicpancreatitis ==1) %>% count() #5453
t3c_t2controls %>% filter(t3c_pancreaticcancer ==1) %>% count() #819
t3c_t2controls %>% filter(t3c_haemochromatosis==1) %>% count() #907

t3c_t2controls %>% count(t3c_subgroup)

################################################################################
###3: DEFINE PANCREATIC EXOCRINE INSUFFICIENCY (PEI)
# PEI = Faecal elastase (FE1) <200 
# OR Pancreatic exocrine insufficiency diagnosis 
# OR Pancreatic enzyme replacement therapy (PERT) prescription
################################################################################

#Raw tables
PEI_gp_raw <- cprd$tables$observation %>% inner_join(codes$pancreatic_exocrine_insufficiency) %>% analysis$cached("all_PEI_raw", indexes= c("patid", "obsdate"))
FE1_gp_raw <- cprd$tables$observation %>% inner_join(codes$fe1) %>% analysis$cached("all_FE1_raw", indexes= c("patid", "obsdate"))
PERT_raw <- cprd$tables$drugIssue %>% inner_join(codes$pert) %>% analysis$cached("all_PERT_raw", indexes= c("patid", "issuedate"))

#FE1 under 200
FE1_cleaned <- FE1_gp_raw %>% filter(!is.na(obsdate)) %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>%
  filter(numunitid == 279 | numunitid == 194 | numunitid == 1760 | numunitid == 1793 | numunitid == 1812 | numunitid == 540 | numunitid == 2196 | is.na(numunitid) | numunitid == 280 | numunitid == 405 | numunitid == 4819 | numunitid == 7161) %>%
  filter(!is.na(testvalue)) %>% group_by(patid,obsdate) %>% summarise(testvalue=mean(testvalue)) %>% ungroup() %>% select(patid, obsdate, testvalue) %>% distinct()

FE1_under200 <- FE1_cleaned %>% filter(testvalue <200) %>% group_by(patid) %>% summarise(first_FE1_under200 = min(obsdate, na.rm = TRUE))

#PEI clinical codes
PEI_cleaned <- PEI_gp_raw %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date)
PEI_diag <- PEI_cleaned %>% group_by(patid) %>% summarise(first_PEI = min(obsdate, na.rm = TRUE))

#PERT prescriptions
PERT_cleaned <- PERT_raw %>% inner_join(cprd$tables$validDateLookup) %>% filter(!is.na(issuedate)) %>% filter(issuedate>=min_dob & issuedate<=gp_ons_end_date)
N_PERT <- PERT_cleaned %>% group_by(patid) %>% summarise(N_PERT_ever = n()) %>% ungroup()
PERT_prescription <- PERT_cleaned %>% group_by(patid) %>% summarise(first_PERT = min(issuedate, na.rm = TRUE)) %>% ungroup()

#Any of above
PEI_all <- FE1_under200 %>% select(patid, date = first_FE1_under200) %>% union(PEI_diag %>% select(patid, date = first_PEI)) %>% union(PERT_prescription %>% select(patid, date = first_PERT)) %>%
  group_by(patid) %>% summarise(first_PEI_all = min(date, na.rm = TRUE))

#Join with cohort
cohort_with_PEI <- t3c_t2controls %>% left_join(FE1_under200, by = "patid") %>% left_join(PEI_diag, by = "patid") %>% left_join(PERT_prescription, by = "patid") %>%
  left_join(PEI_all, by = "patid") %>% mutate(PEI_ever = ifelse(!is.na(first_PEI_all), 1L, 0L)) %>% mutate(PEI_ever = ifelse(is.na(PEI_ever), 0L, PEI_ever)) %>%
  left_join(N_PERT, by = "patid") %>% mutate(N_PERT_ever = ifelse(is.na(N_PERT_ever), 0, N_PERT_ever)) %>%
  analysis$cached("cohort_with_PEI", unique_indexes = "patid")


################################################################################
###4: JOIN COHORT WITH AT DIAGNOSIS TABLE
################################################################################
#Connect to 'at_diag' analysis to get all diabetes cohort table
#analysis = cprd$analysis("at_diag")
#all_diabetes_at_diag <- all_diabetes_at_diag %>% analysis$cached("final_merge", unique_indexes="patid")
#Date copied: 07/06/2023

#Connect back to 'Rhian_T3c' analysis
analysis = cprd$analysis("Rhian_T3c")

all_diabetes_at_diag <- all_diabetes_at_diag %>% analysis$cached("all_diabetes_at_diagnosis", unique_indexes="patid")
#everyone in this table has a diagnosis post-registration (diagnosis before reg excluded)
all_diabetes_at_diag %>% count() #771678

#Join with cohort
at_diagnosis_t3c_t2 <- cohort_with_PEI %>% select(patid, first_acutepancreatitis_pre_diagnosis:PEI_ever) %>% filter(type == "T1" | type == "T2" | type == "T3c") %>%
  inner_join(all_diabetes_at_diag, by = "patid") %>% filter(dm_diag_date >= as.Date("2004-01-01")) %>% mutate(PEI_prior = ifelse(!is.na(first_PEI_all) & first_PEI_all < dm_diag_date, 1L,0L)) %>% 
  left_join(all_diabetes_at_diag %>% select(patid, dm_diag_date) %>% inner_join(PERT_cleaned, by = "patid") %>% filter(issuedate < dm_diag_date) %>% 
              group_by(patid) %>% summarise(N_PERT_prior = n())) %>% mutate(N_PERT_prior = ifelse(is.na(N_PERT_prior), 0, N_PERT_prior)) %>%
  analysis$cached("cohort_at_diagnosis", unique_indexes = "patid")

at_diagnosis_t3c_t2 %>% count() #550079
at_diagnosis_t3c_t2 %>% count(type) # T2: 524084, T3c: 10318, T1: 15677
at_diagnosis_t3c_t2 %>% filter(type == "T3c") %>% count(diabetes_type) # 9443 T2, 616 T1, 259 other

at_diagnosis_t3c_t2 %>% filter(t3c_all ==1) %>% count() #10318
at_diagnosis_t3c_t2 %>% filter(t3c_acutepancreatitis_only ==1) %>% count() #5322
at_diagnosis_t3c_t2 %>% filter(t3c_chronicpancreatitis ==1) %>% count() #3623
at_diagnosis_t3c_t2 %>% filter(t3c_pancreaticcancer ==1) %>% count() #692
at_diagnosis_t3c_t2 %>% filter(t3c_haemochromatosis==1) %>% count() #681

################################################################################
###5: JOIN COHORT WITH MASTERMIND DATA
#Mastermind data for whole diabetes cohort from running adapted version of 10_mm_final_merge script from Github: CPRD-Cohort-Scripts
#Date ran: 08/11/2023
################################################################################

analysis = cprd$analysis("mm")
all_diabetes_1stinstance <- all_diabetes_1stinstance %>% analysis$cached("20231108_all_diabetes_1stinstance", indexes=c("patid", "dstartdate", "drugclass"))
all_diabetes_all_drug_periods <- all_diabetes_all_drug_periods %>% analysis$cached("20231108_all_diabetes_all_drug_periods", indexes=c("patid", "dstartdate", "drugclass"))

#Connect back to 'Rhian_T3c' analysis
analysis = cprd$analysis("Rhian_T3c")

#Define treatment response cohort
t3c_t2_mm <- cohort_with_PEI %>% select(patid, first_acutepancreatitis_pre_diagnosis:PEI_ever) %>% filter(type == "T1" | type == "T2" | type == "T3c") %>% inner_join(all_diabetes_1stinstance, by = "patid") %>% 
  filter(dm_diag_date_all >= as.Date("2004-01-01")) %>% mutate(PEI_prior = ifelse(!is.na(first_PEI_all) & first_PEI_all < dstartdate, 1L,0L)) %>% 
  left_join(all_diabetes_1stinstance %>% select(patid, dstartdate) %>% distinct() %>% inner_join(PERT_cleaned, by = "patid") %>% filter(issuedate < dstartdate) %>% 
              group_by(patid, dstartdate) %>% summarise(N_PERT_prior = n())) %>% mutate(N_PERT_prior = ifelse(is.na(N_PERT_prior), 0, N_PERT_prior)) %>%
  left_join(all_diabetes_1stinstance %>% select(patid, dstartdate) %>% distinct() %>% inner_join(PERT_cleaned, by = "patid") %>% filter(issuedate < dstartdate) %>% 
              group_by(patid, dstartdate) %>% summarise(latest_PERT = max(issuedate))) %>%
  analysis$cached("cohort_20231108_all_1stinstance", indexes=c("patid", "dstartdate", "drugclass"))

t3c_t2_mm %>% count() #923293
t3c_t2_mm %>% distinct(patid) %>% count() #483737

t3c_t2_mm %>% filter(type != "T1") %>% count() #901216
t3c_t2_mm %>% filter(type != "T1") %>% distinct(patid) %>% count() #466943

#Also save all drug periods table
dm_diag_dates <- t3c_t2_mm %>% select(patid, dm_diag_date_all) %>% distinct()

t3c_t2_all_drug_periods <- cohort_with_PEI %>% select(patid, first_acutepancreatitis_pre_diagnosis:PEI_ever) %>% filter(type == "T1" | type == "T2" | type == "T3c") %>%
  inner_join(all_diabetes_all_drug_periods, by = "patid") %>% left_join(dm_diag_dates, by = "patid") %>% filter(dm_diag_date_all >= as.Date("2004-01-01")) %>%
  analysis$cached("cohort_20231108_all_drug_periods", indexes=c("patid", "dstartdate", "drugclass"))

t3c_t2_all_drug_periods %>% count() #
t3c_t2_all_drug_periods %>% distinct(patid) %>% count() #

t3c_t2_all_drug_periods %>% filter(type != "T1") %>% count() #1263067
t3c_t2_all_drug_periods %>% filter(type != "T1") %>% distinct(patid) %>% count() #466943

################################################################################
###END
rm(list=ls())
