################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#First use command prompt to log in to Slade using SSH

#To install package
#library(devtools)
#install_github("Exeter-Diabetes/CPRD-analysis-package")


#Load packages
library(tidyverse)
library(MatchIt)
library(cobalt)
library(broom)
library(kableExtra)
library(patchwork)
library(ggpubr)

#Load aurum package
library(aurum)

################################################################################
###SETUP########################################################################

###Connecting to data and setting up analysis###################################
#Initialise connection
cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "C:/Users/rh530/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "31/10/2021")

#Connect to 'Rhian_T3c' analysis
analysis = cprd$analysis("Rhian_T3c")

#Load cohort
mm_all <- mm_all %>% analysis$cached("cohort_20231108_all_1stinstance", indexes=c("patid", "dstartdate", "drugclass"))

mm_all <- mm_all %>% select(patid, drugclass, type, diabetes_type, t3c_subgroup, t3c_acutepancreatitis_only, t3c_chronicpancreatitis, t3c_pancreaticcancer, t3c_haemochromatosis,
                            t3c_all, gender, ethnicity_5cat, imd2015_10, dstartdate_age, dstartdate, drugline_all, INS, numdrugs, prehba1c, hba1cresp12m, hba1cresp6m,
                            preweight, weightresp6m, weightresp12m, stopdrug_6m_3mFU, first_PEI_all, PEI_ever, PEI_prior, first_pancreatic_disease_all_post_diagnosis)

mm_all <- collect(mm_all)

#Drop T1s
mm_all <- mm_all %>% filter(type != "T1")
mm_all %>% distinct(patid, type) %>% count(type) # 9166 T3c, 457777 T2

################################################################################
##Cohort setup

mm_all <- mm_all %>% mutate(patid = as.character(patid)) %>% mutate(t3c = ifelse(type == "T3c", 1, ifelse(type == "T2", 0, NA)))

#Add year of drug initiation
mm_all$drug_year <- as.numeric(format(mm_all$dstartdate, format = "%Y"))

#Add 10 year age bands for drug initiation age
mm_all <- mm_all %>% mutate(d_age_band = ifelse(dstartdate_age < 30, "<30", ifelse(dstartdate_age <30 & dstartdate_age >= 20, "20-30", 
                                                                                   ifelse(dstartdate_age <40 & dstartdate_age >= 30, "30-40", 
                                                                                          ifelse(dstartdate_age < 50 & dstartdate_age >= 40, "40-50", 
                                                                                                 ifelse(dstartdate_age <60 & dstartdate_age >=50, "50-60", 
                                                                                                        ifelse(dstartdate_age < 70 & dstartdate_age >= 60, "60-70",
                                                                                                               ifelse(dstartdate_age <80 & dstartdate_age >=70, "70-80",
                                                                                                                      ifelse(dstartdate_age>=80, "80+", NA)))))))))

#Add 5 year bands for drug start date
mm_all <- mm_all %>% mutate(d_year_band = ifelse(dstartdate < as.Date("2021-01-01") & dstartdate >= as.Date("2016-01-01"), "2016-2020",
                                                 ifelse(dstartdate < as.Date("2016-01-01") & dstartdate >= as.Date("2011-01-01"), "2011-2015",
                                                        ifelse(dstartdate < as.Date("2011-01-01") & dstartdate >= as.Date("2006-01-01"), "2006-2010",
                                                               ifelse(dstartdate < as.Date("2006-01-01") & dstartdate >= as.Date("2001-01-01"), "2001-2005",
                                                                      ifelse(dstartdate < as.Date("2001-01-01") & dstartdate >= as.Date("1996-01-01"), "1996-2000",
                                                                             ifelse(dstartdate < as.Date("1996-01-01") & dstartdate >= as.Date("1991-01-01"), "1991-1995",
                                                                                    ifelse(dstartdate < as.Date("1991-01-01") & dstartdate >= as.Date("1986-01-01"), "1986-1990",
                                                                                           ifelse(dstartdate < as.Date("1986-01-01") & dstartdate >= as.Date("1981-01-01"), "1981-1985",
                                                                                                  ifelse(dstartdate < as.Date("1981-01-01") & dstartdate >= as.Date("1976-01-01"), "1976-1980",
                                                                                                         ifelse(dstartdate < as.Date("1976-01-01") & dstartdate >= as.Date("1971-01-01"), "1971-1975",
                                                                                                                ifelse(dstartdate < as.Date("1971-01-01") & dstartdate >= as.Date("1966-01-01"), "1966-1970", NA))))))))))))

#Add drugline category : 1,2 or 3+
mm_all <- mm_all %>% mutate(drugline_band = ifelse(drugline_all == 1, "1", ifelse(drugline_all ==2, "2", ifelse(drugline_all >= 3, "3+", NA))))

#Group IMD into quintiles
mm_all <- mm_all %>% mutate(imd_quintile = ifelse(imd2015_10 == 1 | imd2015_10 == 2, "1", ifelse(imd2015_10 == 3 | imd2015_10 == 4, "2", ifelse(imd2015_10 == 5 | imd2015_10 == 6, "3", ifelse(imd2015_10 == 7 | imd2015_10 == 8, "4", ifelse(imd2015_10 == 9 | imd2015_10 == 10, "5", "Missing")))))) %>%
  mutate(imd_quintile = ifelse(is.na(imd_quintile), "Missing", imd_quintile))

#Change ethnicity NA to 5
mm_all <- mm_all %>% mutate(ethnicity_5cat = ifelse(is.na(ethnicity_5cat), "5", ethnicity_5cat))

mm_all <- mm_all %>% mutate(gender = as.factor(gender), ethnicity_5cat = as.factor(ethnicity_5cat), imd_quintile = as.factor(imd_quintile), d_year_band = as.factor(d_year_band), d_age_band = as.factor(d_age_band), drugline_all = as.factor(drugline_all), drugline_band = as.factor(drugline_band), INS = as.factor(INS))

#Keep just drug classes interested in
mm_all <- mm_all %>% filter(drugclass == "MFN" | drugclass == "SU" | drugclass == "TZD" | drugclass == "DPP4" | drugclass == "SGLT2")
mm_all %>% distinct(patid, type) %>% count(type) # 7999 T3c, 450990 T2

#Exclude T2s who develop pancreatic disease before drug start date
mm_all <- mm_all %>% filter(type == "T3c" |(type == "T2" & is.na(first_pancreatic_disease_all_post_diagnosis)) | (type == "T2" & first_pancreatic_disease_all_post_diagnosis > dstartdate))
mm_all %>% distinct(patid, type) %>% count(type) # 7999 T3c, 449633 T2

#Drop anyone with concurrent insulin
mm_all <- mm_all %>% filter(INS != 1)
mm_all %>% distinct(patid, type) %>% count(type) # 7622 T3c, 444369 T2

#Drop anyone with T2 and PEI ever, or T3c and PEI after drug start
mm_all <- mm_all %>% filter((type == "T3c" & PEI_prior ==1) | (type == "T3c" & PEI_ever ==0 & PEI_prior ==0) | (type == "T2" & PEI_ever ==0))
mm_all %>% distinct(patid, type) %>% count(type) # 7182 T3c, 441774 T2

#Make a categorical subgroup variable: T3c PEI, T3c no PEI, or T2
mm_all <- mm_all %>% mutate(subgroup = ifelse(type == "T2", "T2", ifelse(type == "T3c" & PEI_prior ==1, "T3c PEI", ifelse(type == "T3c" & PEI_prior ==0, "T3c No PEI", NA))))

#Check counts
mm_all %>% count(subgroup) # 771795 T2, 10363 T3c No PEI, 1893 T3c PEI
mm_all %>% count(t3c_subgroup) # 771795 T2, 6932 AP, 3975 CP, 832 haemochromatosis, 517 pancreatic cancer

mm_all %>% distinct(patid, subgroup) %>% count(subgroup) # 1185 T3c PEI, 5997 T3c no PEI, 441774 T2
mm_all %>% count(type, drugclass)
mm_all %>% count(subgroup, drugclass)

################################################################################
################################################################################
###Matching

drugs <- c("SGLT2", "DPP4", "MFN", "SU", "TZD")

for (i in drugs) {
  print(i)
  exact_cohort_name <- paste0(i,"_all_exact_matches")
  matched_cohort_name <- paste0(i,"_matched_cohort")
  
  
  full_drug_cohort <- mm_all %>% filter(drugclass == i)
  full_drug_cohort %>% count(t3c)
  
  #Exact matching: sex, ethnicity, imd, drug start year (5 year bands), age (10 year bands), drugline (1,2, or 3+), on insulin
  match <- matchit(t3c ~ gender + ethnicity_5cat + imd_quintile + d_year_band + d_age_band + drugline_band, data = full_drug_cohort,
                   method = "exact")
  print(summary(match)) 
  exact_match_all <- match.data(match) %>% rename (weights_1 = weights, subclass_1 = subclass)
  assign(exact_cohort_name,exact_match_all)
  
  #Nearest neighbour match on continuous age + calendar year
  #Distance = "glm" estimates propensity scores using logistic regression
  match_ratio <- matchit(t3c ~ dstartdate_age + drug_year,
                         data = exact_match_all, exact = ~ subclass_1, method = "nearest", distance="glm", caliper=0.05, ratio = 10)
  print(summary(match_ratio))
  matched_cohort <- match.data(match_ratio)
  assign(matched_cohort_name,matched_cohort)
}

overall_matched_cohort <- SGLT2_matched_cohort %>% union(DPP4_matched_cohort) %>% union(MFN_matched_cohort) %>% union(SU_matched_cohort) %>% union(TZD_matched_cohort)

#Check counts
overall_matched_cohort %>% count(type) # 11971 T3c, 107367 T2
overall_matched_cohort %>% count(subgroup) # 107367 T2, 10118 T3c No PEI, 1853 T3c PEI
overall_matched_cohort %>% count(t3c_subgroup) # 107367 T2, 6766 AP, 3880 CP, 819 haemochromatosis, 506 pancreatic cancer

overall_matched_cohort %>% distinct(patid, subgroup) %>% count(subgroup) # 1167 T3c PEI, 5917 T3c no PEI, 97227 T2
overall_matched_cohort %>% count(subgroup, drugclass)

################################################################################
#Get cohorts with valid HbA1c responses
for (i in drugs) {
  print(i)
  matched_table_name <- paste0(i,"_matched_cohort")
  HbA1cresp_table_name <- paste0(i, "_HbA1cresp_cohort")
  
  if(exists(matched_table_name)) {
    matched_table <- get(matched_table_name)
    print(matched_table %>% count(t3c))
    
    #Drop all cases and controls with no valid HbA1c response and controls belonging to case with no valid response
    match_hba1c <- matched_table %>% 
      mutate(prehba1c = as.numeric(prehba1c), numdrugs = as.numeric(numdrugs), hba1cresp12m = as.numeric(hba1cresp12m), hba1cresp6m = as.numeric(hba1cresp6m), t3c = as.numeric(t3c)) %>%
      mutate(hba1cresp = coalesce(hba1cresp12m, hba1cresp6m)) %>% filter(!is.na(hba1cresp))
    match_hba1c_resp <- match_hba1c %>% inner_join(match_hba1c %>% filter(t3c == 1) %>% distinct(subclass)) %>% inner_join(match_hba1c %>% filter(t3c == 0) %>% distinct(subclass))
    print(match_hba1c_resp %>% count(t3c))
    
    weights_table <- match_hba1c_resp %>% group_by(drugclass, subclass, type) %>% summarise(n=n()) %>% ungroup() %>%
      pivot_wider(names_from = type, values_from = n) %>% rename(T2_total = T2, T3c_total = T3c)
    
    match_hba1c_resp <- match_hba1c_resp %>% select(-weights, -weights_1) %>% left_join(weights_table) %>%
      arrange(subclass) %>% mutate(weights = ifelse(t3c==1, T3c_total/T3c_total, ifelse(t3c==0, T3c_total/T2_total, NA)))
    
    assign(HbA1cresp_table_name, match_hba1c_resp)
    
  }
}

overall_HbA1cresp_cohort <- SGLT2_HbA1cresp_cohort %>% union(DPP4_HbA1cresp_cohort) %>% union(MFN_HbA1cresp_cohort) %>% union(SU_HbA1cresp_cohort) %>% union(TZD_HbA1cresp_cohort)

################################################################################
#Get cohorts with valid discontinuation (within 6 months with 3 months follow-up)
for (i in drugs) {
  print(i)
  matched_table_name <- paste0(i,"_matched_cohort")
  discontinuation_table_name <- paste0(i, "_discontinuation_cohort")
  
  if(exists(matched_table_name)) {
    matched_table <- get(matched_table_name)
    print(matched_table %>% count(t3c))
    
    #Drop all cases and controls with no valid discontinuation outcome and controls belonging to case with no valid outcome
    match_discont <- matched_table %>% filter(!is.na(stopdrug_6m_3mFU))
    match_discont_resp <- match_discont %>% inner_join(match_discont %>% filter(t3c == 1) %>% distinct(subclass)) %>% inner_join(match_discont %>% filter(t3c == 0) %>% distinct(subclass))
    print(match_discont_resp %>% count(t3c))
    
    weights_table <- match_discont_resp %>% group_by(drugclass, subclass, type) %>% summarise(n=n()) %>% ungroup() %>%
      pivot_wider(names_from = type, values_from = n) %>% rename(T2_total = T2, T3c_total = T3c)
    
    match_discont_resp <- match_discont_resp %>% select(-weights, -weights_1) %>% left_join(weights_table) %>%
      arrange(subclass) %>% mutate(weights = ifelse(t3c==1, T3c_total/T3c_total, ifelse(t3c==0, T3c_total/T2_total, NA)))
    
    assign(discontinuation_table_name, match_discont_resp)
    
  }
}

overall_discontinuation_cohort <- SGLT2_discontinuation_cohort %>% union(DPP4_discontinuation_cohort) %>% union(MFN_discontinuation_cohort) %>% union(SU_discontinuation_cohort) %>% union(TZD_discontinuation_cohort)

################################################################################
#Get cohorts with valid weight responses
for (i in drugs) {
  print(i)
  matched_table_name <- paste0(i,"_matched_cohort")
  weightresp_table_name <- paste0(i, "_weightresp_cohort")
  
  if(exists(matched_table_name)) {
    matched_table <- get(matched_table_name)
    print(matched_table %>% count(t3c))
    
    #Drop all cases and controls with no valid weight response and controls belonging to case with no valid response
    match_weight <- matched_table %>% 
      mutate(preweight = as.numeric(preweight), numdrugs = as.numeric(numdrugs), weightresp12m = as.numeric(weightresp12m), weightresp6m = as.numeric(weightresp6m), t3c = as.numeric(t3c)) %>%
      mutate(weightresp = coalesce(weightresp12m, weightresp6m)) %>% filter(!is.na(weightresp))
    match_weight_resp <- match_weight %>% inner_join(match_weight %>% filter(t3c == 1) %>% distinct(subclass)) %>% inner_join(match_weight %>% filter(t3c == 0) %>% distinct(subclass))
    print(match_weight_resp %>% count(t3c))
    
    weights_table <- match_weight_resp %>% group_by(drugclass, subclass, type) %>% summarise(n=n()) %>% ungroup() %>%
      pivot_wider(names_from = type, values_from = n) %>% rename(T2_total = T2, T3c_total = T3c)
    
    match_weight_resp <- match_weight_resp %>% select(-weights, -weights_1) %>% left_join(weights_table) %>%
      arrange(subclass) %>% mutate(weights = ifelse(t3c==1, T3c_total/T3c_total, ifelse(t3c==0, T3c_total/T2_total, NA)))
    
    assign(weightresp_table_name, match_weight_resp)
    
  }
}

overall_weightresp_cohort <- SGLT2_weightresp_cohort %>% union(DPP4_weightresp_cohort) %>% union(MFN_weightresp_cohort) %>% union(SU_weightresp_cohort) %>% union(TZD_weightresp_cohort)


################################################################################
##Analysis
################################################################################

###HbA1c
##For all subtypes except haemochromatosis, plus overall
drugs2 <- c(drugs, "overall")
conditions <- c("all", "acutepancreatitis_only", "chronicpancreatitis", "pancreaticcancer")

HbA1c_results_ALL <- NULL

for (i in drugs2) {
  print(i)
  HbA1cresp_table_name <- paste0(i, "_HbA1cresp_cohort")
  
  if(exists(HbA1cresp_table_name)) {
    match_hba1c_resp <- get(HbA1cresp_table_name)
    print(match_hba1c_resp %>% count(t3c))
    
    #Linear regression, adjusted for baseline HbA1c and number of current glucose-lowering therapies
    HbA1c_results <- NULL
    
    for (c in conditions) {
      print(c)
      
      prior_condition <- paste0("t3c_", c)
      
      sub_cohort <- match_hba1c_resp %>% 
        inner_join(match_hba1c_resp %>% filter(!!as.symbol(prior_condition) == 1) %>% distinct(drugclass, subclass))
      model <- lm(hba1cresp ~ subgroup + prehba1c + numdrugs, data = sub_cohort, weights = weights)
      print(summary(model))
      
      PEI <- sub_cohort %>% filter(subgroup == "T3c PEI")
      no_PEI <- sub_cohort %>% filter(subgroup == "T3c No PEI")
      t2 <- sub_cohort %>% filter(subgroup == "T2")
      
      if(i == "MFN") {
        PEI_means <- data.frame(subgroup = "T3c PEI", prehba1c = 73, numdrugs =1)
        no_PEI_means <- data.frame(subgroup = "T3c No PEI", prehba1c = 73, numdrugs =1)
        t2_means <- data.frame(subgroup = "T2", prehba1c = 73, numdrugs =1) }
      
      else {
        PEI_means <- data.frame(subgroup = "T3c PEI", prehba1c = 73, numdrugs =2)
        no_PEI_means <- data.frame(subgroup = "T3c No PEI", prehba1c = 73, numdrugs =2)
        t2_means <- data.frame(subgroup = "T2", prehba1c = 73, numdrugs =2) }
      
      # Tidy
      model_df <- tidy(model, conf.int = TRUE) %>% filter(term == "subgroupT3c PEI" | term == "subgroupT3c No PEI") %>% 
        mutate(term = ifelse(term == "subgroupT3c PEI", "PEI", ifelse(term == "subgroupT3c No PEI", "No PEI", term))) %>%
        
        # Add N's, average baseline HbA1c, average unadjusted response, and estimated average adjusted response in all groups
        mutate(av_PEI_prehba1c = mean(PEI$prehba1c, na.rm = TRUE)) %>%
        mutate(av_no_PEI_prehba1c = mean(no_PEI$prehba1c, na.rm = TRUE)) %>%
        mutate(av_t2_prehba1c = mean(t2$prehba1c, na.rm = TRUE)) %>%
        mutate(av_PEI_unadjust = mean(PEI$hba1cresp, na.rm = TRUE)) %>%
        mutate(av_no_PEI_unadjust = mean(no_PEI$hba1cresp, na.rm = TRUE)) %>%
        mutate(av_t2_unadjust = mean(t2$hba1cresp, na.rm = TRUE)) %>%
        mutate(av_PEI = predict(model, PEI_means, type = "response", interval = "confidence")) %>%
        mutate(av_no_PEI = predict(model, no_PEI_means, type = "response", interval = "confidence")) %>%
        mutate(av_t2 = predict(model, t2_means, type = "response", interval = "confidence")) %>%
        mutate(PEI = paste0(count(sub_cohort %>% filter(subgroup == "T3c PEI")))) %>%
        mutate(no_PEI = paste0(count(sub_cohort %>% filter(subgroup == "T3c No PEI")))) %>%
        mutate(t2 = paste0(count(sub_cohort %>% filter(subgroup == "T2")))) %>%
        mutate(condition = paste0(prior_condition))
      
      if (is.null(HbA1c_results))
        HbA1c_results <- model_df
      else
        HbA1c_results <- bind_rows(HbA1c_results, model_df)
      
    }
    
    HbA1c_results <- HbA1c_results %>% mutate(drug = i) %>% rename(subgroup = term)
    
    if (is.null(HbA1c_results_ALL))
      HbA1c_results_ALL <- HbA1c_results
    else
      HbA1c_results_ALL <- bind_rows(HbA1c_results_ALL, HbA1c_results)
    
  }
}


###Add haemochromatosis (no PEI group)
haemochromatosis_results_HbA1c <- NULL

for (i in drugs2) {
  print(i)
  HbA1cresp_table_name <- paste0(i, "_HbA1cresp_cohort")
  
  if(exists(HbA1cresp_table_name)) {
    match_hba1c_resp <- get(HbA1cresp_table_name)
    print(match_hba1c_resp %>% count(t3c))
    
    #Linear regression, adjusted for baseline HbA1c and number of current glucose-lowering therapies
    sub_cohort <- match_hba1c_resp %>% 
      inner_join(match_hba1c_resp %>% filter(t3c_haemochromatosis == 1 & subgroup == "T3c No PEI") %>% distinct(drugclass, subclass))
    model <- lm(hba1cresp ~ subgroup + prehba1c + numdrugs, data = sub_cohort, weights = weights)
    print(summary(model))
    
    no_PEI <- sub_cohort %>% filter(subgroup == "T3c No PEI")
    t2 <- sub_cohort %>% filter(subgroup == "T2")
    
    if(i == "MFN") {
      no_PEI_means <- data.frame(subgroup = "T3c No PEI", prehba1c = 73, numdrugs =1)
      t2_means <- data.frame(subgroup = "T2", prehba1c = 73, numdrugs =1) }
    
    else {
      no_PEI_means <- data.frame(subgroup = "T3c No PEI", prehba1c = 73, numdrugs =2)
      t2_means <- data.frame(subgroup = "T2", prehba1c = 73, numdrugs =2) }
    
    # Tidy
    model_df <- tidy(model, conf.int = TRUE) %>% filter(term == "subgroupT3c No PEI") %>% 
      mutate(term = ifelse(term == "subgroupT3c No PEI", "No PEI", term)) %>%
      
      # Add N's, average baseline HbA1c, average unadjusted response, and estimated average adjusted response in all groups
      mutate(av_no_PEI_prehba1c = mean(no_PEI$prehba1c, na.rm = TRUE)) %>%
      mutate(av_t2_prehba1c = mean(t2$prehba1c, na.rm = TRUE)) %>%
      mutate(av_no_PEI_unadjust = mean(no_PEI$hba1cresp, na.rm = TRUE)) %>%
      mutate(av_t2_unadjust = mean(t2$hba1cresp, na.rm = TRUE)) %>%
      mutate(av_no_PEI = predict(model, no_PEI_means, type = "response", interval = "confidence")) %>%
      mutate(av_t2 = predict(model, t2_means, type = "response", interval = "confidence")) %>%
      mutate(no_PEI = paste0(count(sub_cohort %>% filter(subgroup == "T3c No PEI")))) %>%
      mutate(t2 = paste0(count(sub_cohort %>% filter(subgroup == "T2")))) %>%
      mutate(condition = "t3c_haemochromatosis") %>% mutate(drug = i)
    
    if (is.null(haemochromatosis_results_HbA1c))
      haemochromatosis_results_HbA1c <- model_df
    else
      haemochromatosis_results_HbA1c <- bind_rows(haemochromatosis_results_HbA1c, model_df)
  }
}

################################################################################

###Discontinuation
##For all subtypes except haemochromatosis, plus overall
discontinuation_results_ALL <- NULL

for (i in drugs2) {
  print(i)
  discontinuation_table_name <- paste0(i, "_discontinuation_cohort")
  
  if(exists(discontinuation_table_name)) {
    match_discont_resp <- get(discontinuation_table_name)
    print(match_discont_resp %>% count(t3c))
    
    #Use logistic regression, adjusted for baseline HbA1c and number of current glucose-lowering therapies, to get odd ratio for discontinuation
    discontinuation_results <- NULL
    
    for (c in conditions) {
      print(c)
      
      prior_condition <- paste0("t3c_", c)
      
      sub_cohort <- match_discont_resp %>% 
        inner_join(match_discont_resp %>% filter(!!as.symbol(prior_condition) == 1) %>% distinct(drugclass, subclass))
      model <- glm(stopdrug_6m_3mFU ~ subgroup + prehba1c + numdrugs, data = sub_cohort, family = binomial(link = "logit"), weights = weights)
      print(summary(model))
      
      PEI <- sub_cohort %>% filter(subgroup == "T3c PEI")
      no_PEI <- sub_cohort %>% filter(subgroup == "T3c No PEI")
      t2 <- sub_cohort %>% filter(subgroup == "T2")
      
      if(i == "MFN") {
        PEI_means <- data.frame(subgroup = "T3c PEI", prehba1c = 73, numdrugs =1)
        no_PEI_means <- data.frame(subgroup = "T3c No PEI", prehba1c = 73, numdrugs =1)
        t2_means <- data.frame(subgroup = "T2", prehba1c = 73, numdrugs =1) }
      
      else {
        PEI_means <- data.frame(subgroup = "T3c PEI", prehba1c = 73, numdrugs =2)
        no_PEI_means <- data.frame(subgroup = "T3c No PEI", prehba1c = 73, numdrugs =2)
        t2_means <- data.frame(subgroup = "T2", prehba1c = 73, numdrugs =2) }
      
      # Tidy
      model_df <- tidy(model,  conf.int = TRUE, exponentiate = TRUE) %>% filter(term == "subgroupT3c PEI" | term == "subgroupT3c No PEI") %>% 
        mutate(term = ifelse(term == "subgroupT3c PEI", "PEI", ifelse(term == "subgroupT3c No PEI", "No PEI", term))) %>%
        
        # Add N's, % discontinuation unadjusted and adjusted estimate
        mutate(prob_PEI_unadjust = mean(PEI$stopdrug_6m_3mFU, na.rm = TRUE)) %>%
        mutate(prob_no_PEI_unadjust = mean(no_PEI$stopdrug_6m_3mFU, na.rm = TRUE)) %>%
        mutate(prob_t2_unadjust = mean(t2$stopdrug_6m_3mFU, na.rm = TRUE)) %>%
        mutate(PEI = paste0(count(sub_cohort %>% filter(subgroup == "T3c PEI")))) %>%
        mutate(no_PEI = paste0(count(sub_cohort %>% filter(subgroup == "T3c No PEI")))) %>%
        mutate(t2 = paste0(count(sub_cohort %>% filter(subgroup == "T2")))) %>%
        mutate(condition = paste0(prior_condition))
      
      pred_PEI <- as.data.frame(predict(model, PEI_means, type = "link", se.fit=TRUE)) %>% 
        mutate(lower = model$family$linkinv(fit - 1.96*se.fit), point.est = model$family$linkinv(fit), upper = model$family$linkinv(fit + 1.96*se.fit)) %>%
        select(lower_predict_PEI = lower, predict_est_PEI = point.est, upper_predict_PEI = upper) %>% mutate(condition = paste0(prior_condition))
      pred_no_PEI <- as.data.frame(predict(model, no_PEI_means, type = "link", se.fit=TRUE)) %>% 
        mutate(lower = model$family$linkinv(fit - 1.96*se.fit), point.est = model$family$linkinv(fit), upper = model$family$linkinv(fit + 1.96*se.fit)) %>%
        select(lower_predict_no_PEI = lower, predict_est_no_PEI = point.est, upper_predict_no_PEI = upper) %>% mutate(condition = paste0(prior_condition))
      pred_t2 <- as.data.frame(predict(model, t2_means, type = "link", se.fit=TRUE)) %>% 
        mutate(lower = model$family$linkinv(fit - 1.96*se.fit), point.est = model$family$linkinv(fit), upper = model$family$linkinv(fit + 1.96*se.fit)) %>%
        select(lower_predict_t2 = lower, predict_est_t2 = point.est, upper_predict_t2 = upper) %>% mutate(condition = paste0(prior_condition))
      
      model_df <- model_df %>% left_join(pred_PEI) %>% left_join(pred_no_PEI) %>% left_join(pred_t2)
      
      if (is.null(discontinuation_results))
        discontinuation_results <- model_df
      else
        discontinuation_results <- bind_rows(discontinuation_results, model_df)
      
    }
    
    discontinuation_results <- discontinuation_results %>% mutate(drug = i) %>% rename(subgroup = term)
    
    if (is.null(discontinuation_results_ALL))
      discontinuation_results_ALL <- discontinuation_results
    else
      discontinuation_results_ALL <- bind_rows(discontinuation_results_ALL, discontinuation_results)
    
  }
}


###Add haemochromatosis (no PEI group)
haemochromatosis_results_discont <- NULL

for (i in drugs2) {
  print(i)
  discontinuation_table_name <- paste0(i, "_discontinuation_cohort")
  
  if(exists(discontinuation_table_name)) {
    match_discont_resp <- get(discontinuation_table_name)
    print(match_discont_resp %>% count(t3c))
    
    #Use logistic regression, adjusted for baseline HbA1c and number of current glucose-lowering therapies, to get odd ratio for discontinuation
    sub_cohort <- match_discont_resp %>% 
      inner_join(match_discont_resp %>% filter(t3c_haemochromatosis == 1 & subgroup == "T3c No PEI") %>% distinct(drugclass, subclass))
    model <- glm(stopdrug_6m_3mFU ~ subgroup + prehba1c + numdrugs, data = sub_cohort, family = binomial(link = "logit"), weights = weights)
    print(summary(model))
    
    no_PEI <- sub_cohort %>% filter(subgroup == "T3c No PEI")
    t2 <- sub_cohort %>% filter(subgroup == "T2")
    
    if(i == "MFN") {
      no_PEI_means <- data.frame(subgroup = "T3c No PEI", prehba1c = 73, numdrugs =1)
      t2_means <- data.frame(subgroup = "T2", prehba1c = 73, numdrugs =1) }
    
    else {
      no_PEI_means <- data.frame(subgroup = "T3c No PEI", prehba1c = 73, numdrugs =2)
      t2_means <- data.frame(subgroup = "T2", prehba1c = 73, numdrugs =2) }
    
    # Tidy
    model_df <- tidy(model,  conf.int = TRUE, exponentiate = TRUE) %>% filter(term == "subgroupT3c No PEI") %>% 
      mutate(term = ifelse(term == "subgroupT3c No PEI", "No PEI", term)) %>%
      
      # Add N's, % discontinuation unadjusted and adjusted estimate
      mutate(prob_no_PEI_unadjust = mean(no_PEI$stopdrug_6m_3mFU, na.rm = TRUE)) %>%
      mutate(prob_t2_unadjust = mean(t2$stopdrug_6m_3mFU, na.rm = TRUE)) %>%
      mutate(no_PEI = paste0(count(sub_cohort %>% filter(subgroup == "T3c No PEI")))) %>%
      mutate(t2 = paste0(count(sub_cohort %>% filter(subgroup == "T2")))) %>%
      mutate(condition = "t3c_haemochromatosis") %>% mutate(drug = i)
    
    pred_no_PEI <- as.data.frame(predict(model, no_PEI_means, type = "link", se.fit=TRUE)) %>% 
      mutate(lower = model$family$linkinv(fit - 1.96*se.fit), point.est = model$family$linkinv(fit), upper = model$family$linkinv(fit + 1.96*se.fit)) %>%
      select(lower_predict_no_PEI = lower, predict_est_no_PEI = point.est, upper_predict_no_PEI = upper) %>% mutate(condition = "t3c_haemochromatosis")
    pred_t2 <- as.data.frame(predict(model, t2_means, type = "link", se.fit=TRUE)) %>% 
      mutate(lower = model$family$linkinv(fit - 1.96*se.fit), point.est = model$family$linkinv(fit), upper = model$family$linkinv(fit + 1.96*se.fit)) %>%
      select(lower_predict_t2 = lower, predict_est_t2 = point.est, upper_predict_t2 = upper) %>% mutate(condition = "t3c_haemochromatosis")
    
    model_df <- model_df %>% left_join(pred_no_PEI) %>% left_join(pred_t2)
    
    if (is.null(haemochromatosis_results_discont))
      haemochromatosis_results_discont <- model_df
    else
      haemochromatosis_results_discont <- bind_rows(haemochromatosis_results_discont, model_df)
  }
}

################################################################################

###Weight
##For all subtypes except haemochromatosis, plus overall
weight_results_ALL <- NULL

for (i in drugs2) {
  print(i)
  weightresp_table_name <- paste0(i, "_weightresp_cohort")
  
  if(exists(weightresp_table_name)) {
    match_weight_resp <- get(weightresp_table_name)
    print(match_weight_resp %>% count(t3c))
    
    #Linear regression, adjusted for baseline weight and number of current glucose-lowering therapies
    weight_results <- NULL
    
    for (c in conditions) {
      print(c)
      
      prior_condition <- paste0("t3c_", c)
      
      sub_cohort <- match_weight_resp %>% 
        inner_join(match_weight_resp %>% filter(!!as.symbol(prior_condition) == 1) %>% distinct(drugclass, subclass))
      model <- lm(weightresp ~ subgroup + preweight + numdrugs, data = sub_cohort, weights = weights)
      print(summary(model))
      
      PEI <- sub_cohort %>% filter(subgroup == "T3c PEI")
      no_PEI <- sub_cohort %>% filter(subgroup == "T3c No PEI")
      t2 <- sub_cohort %>% filter(subgroup == "T2")
      
      if(i == "MFN") {
        PEI_means <- data.frame(subgroup = "T3c PEI", preweight = 93, numdrugs =1)
        no_PEI_means <- data.frame(subgroup = "T3c No PEI", preweight = 93, numdrugs =1)
        t2_means <- data.frame(subgroup = "T2", preweight = 93, numdrugs =1) }
      
      else {
        PEI_means <- data.frame(subgroup = "T3c PEI", preweight = 93, numdrugs =2)
        no_PEI_means <- data.frame(subgroup = "T3c No PEI", preweight = 93, numdrugs =2)
        t2_means <- data.frame(subgroup = "T2", preweight = 93, numdrugs =2) }
      
      # Tidy
      model_df <- tidy(model, conf.int = TRUE) %>% filter(term == "subgroupT3c PEI" | term == "subgroupT3c No PEI") %>% 
        mutate(term = ifelse(term == "subgroupT3c PEI", "PEI", ifelse(term == "subgroupT3c No PEI", "No PEI", term))) %>%
        
        # Add N's, average baseline weight, average unadjusted response, and estimated average adjusted response in all groups
        mutate(av_PEI_preweight = mean(PEI$preweight, na.rm = TRUE)) %>%
        mutate(av_no_PEI_preweight = mean(no_PEI$preweight, na.rm = TRUE)) %>%
        mutate(av_t2_preweight = mean(t2$preweight, na.rm = TRUE)) %>%
        mutate(av_PEI_unadjust = mean(PEI$weightresp, na.rm = TRUE)) %>%
        mutate(av_no_PEI_unadjust = mean(no_PEI$weightresp, na.rm = TRUE)) %>%
        mutate(av_t2_unadjust = mean(t2$weightresp, na.rm = TRUE)) %>%
        mutate(av_PEI = predict(model, PEI_means, type = "response", interval = "confidence")) %>%
        mutate(av_no_PEI = predict(model, no_PEI_means, type = "response", interval = "confidence")) %>%
        mutate(av_t2 = predict(model, t2_means, type = "response", interval = "confidence")) %>%
        mutate(PEI = paste0(count(sub_cohort %>% filter(subgroup == "T3c PEI")))) %>%
        mutate(no_PEI = paste0(count(sub_cohort %>% filter(subgroup == "T3c No PEI")))) %>%
        mutate(t2 = paste0(count(sub_cohort %>% filter(subgroup == "T2")))) %>%
        mutate(condition = paste0(prior_condition))
      
      if (is.null(weight_results))
        weight_results <- model_df
      else
        weight_results <- bind_rows(weight_results, model_df)
      
    }
    
    weight_results <- weight_results %>% mutate(drug = i) %>% rename(subgroup = term)
    
    if (is.null(weight_results_ALL))
      weight_results_ALL <- weight_results
    else
      weight_results_ALL <- bind_rows(weight_results_ALL, weight_results)
    
  }
}


###Add haemochromatosis (no PEI group)
haemochromatosis_results_weight <- NULL

for (i in drugs2) {
  print(i)
  weightresp_table_name <- paste0(i, "_weightresp_cohort")
  
  if(exists(weightresp_table_name)) {
    match_weight_resp <- get(weightresp_table_name)
    print(match_weight_resp %>% count(t3c))
    
    #Linear regression, adjusted for baseline HbA1c and number of current glucose-lowering therapies
    sub_cohort <- match_weight_resp %>% 
      inner_join(match_weight_resp %>% filter(t3c_haemochromatosis == 1 & subgroup == "T3c No PEI") %>% distinct(drugclass, subclass))
    model <- lm(weightresp ~ subgroup + preweight + numdrugs, data = sub_cohort, weights = weights)
    print(summary(model))
    
    no_PEI <- sub_cohort %>% filter(subgroup == "T3c No PEI")
    t2 <- sub_cohort %>% filter(subgroup == "T2")
    
    if(i == "MFN") {
      no_PEI_means <- data.frame(subgroup = "T3c No PEI", preweight = 93, numdrugs =1)
      t2_means <- data.frame(subgroup = "T2", preweight = 93, numdrugs =1) }
    
    else {
      no_PEI_means <- data.frame(subgroup = "T3c No PEI", preweight = 93, numdrugs =2)
      t2_means <- data.frame(subgroup = "T2", preweight = 93, numdrugs =2) }
    
    # Tidy
    model_df <- tidy(model, conf.int = TRUE) %>% filter(term == "subgroupT3c No PEI") %>% 
      mutate(term = ifelse(term == "subgroupT3c No PEI", "No PEI", term)) %>%
      
      # Add N's, average baseline weight, average unadjusted response, and estimated average adjusted response in all groups
      mutate(av_no_PEI_preweight = mean(no_PEI$preweight, na.rm = TRUE)) %>%
      mutate(av_t2_preweight = mean(t2$preweight, na.rm = TRUE)) %>%
      mutate(av_no_PEI_unadjust = mean(no_PEI$weightresp, na.rm = TRUE)) %>%
      mutate(av_t2_unadjust = mean(t2$weightresp, na.rm = TRUE)) %>%
      mutate(av_no_PEI = predict(model, no_PEI_means, type = "response", interval = "confidence")) %>%
      mutate(av_t2 = predict(model, t2_means, type = "response", interval = "confidence")) %>%
      mutate(no_PEI = paste0(count(sub_cohort %>% filter(subgroup == "T3c No PEI")))) %>%
      mutate(t2 = paste0(count(sub_cohort %>% filter(subgroup == "T2")))) %>%
      mutate(condition = "t3c_haemochromatosis") %>% mutate(drug = i)
    
    if (is.null(haemochromatosis_results_weight))
      haemochromatosis_results_weight <- model_df
    else
      haemochromatosis_results_weight <- bind_rows(haemochromatosis_results_weight, model_df)
  }
}

################################################################################
##Output

###HbA1c
#Reformat HbA1c_results_ALL table so have 3 rows for each subgroup/drug: PEI, No PEI, and T2 (T2 estimate etc will be missing as is reference group)
#Also add haemochromatosis results (No PEI/ T2)
PEI_results_HbA1c <- HbA1c_results_ALL %>% filter(subgroup == "PEI") %>% mutate(av = as.numeric(av_PEI[,1]), av_lower = as.numeric(av_PEI[,2]), av_upper = as.numeric(av_PEI[,3])) %>%
  select(drug, condition, subgroup, N = PEI, av_prehba1c = av_PEI_prehba1c, av_unadjust = av_PEI_unadjust, av, av_lower, av_upper, estimate, std.error, statistic, p.value, conf.low, conf.high)
no_PEI_results_HbA1c <- HbA1c_results_ALL %>% filter(subgroup == "No PEI") %>% mutate(av = as.numeric(av_no_PEI[,1]), av_lower = as.numeric(av_no_PEI[,2]), av_upper = as.numeric(av_no_PEI[,3])) %>%
  select(drug, condition, subgroup, N = no_PEI, av_prehba1c = av_no_PEI_prehba1c, av_unadjust = av_no_PEI_unadjust, av, av_lower, av_upper, estimate, std.error, statistic, p.value, conf.low, conf.high)
T2_results_HbA1c <- HbA1c_results_ALL %>% mutate(av = as.numeric(av_t2[,1]), av_lower = as.numeric(av_t2[,2]), av_upper = as.numeric(av_t2[,3])) %>%
  select(drug, condition, N = t2, av_prehba1c = av_t2_prehba1c, av_unadjust = av_t2_unadjust, av, av_lower, av_upper) %>% mutate(subgroup = "Type 2", estimate = NA, std.error = NA, statistic = NA, p.value = NA, conf.low = NA, conf.high = NA) %>% distinct()
no_PEI_haemochromatosis_HbA1c <- haemochromatosis_results_HbA1c %>% filter(term == "No PEI") %>% mutate(av = as.numeric(av_no_PEI[,1]), av_lower = as.numeric(av_no_PEI[,2]), av_upper = as.numeric(av_no_PEI[,3])) %>%
  select(drug, condition, subgroup = term, N = no_PEI, av_prehba1c = av_no_PEI_prehba1c, av_unadjust = av_no_PEI_unadjust, av, av_lower, av_upper, estimate, std.error, statistic, p.value, conf.low, conf.high)
T2_haemochromatosis_HbA1c <- haemochromatosis_results_HbA1c %>% mutate(av = as.numeric(av_t2[,1]), av_lower = as.numeric(av_t2[,2]), av_upper = as.numeric(av_t2[,3])) %>%
  select(drug, condition, N = t2, av_prehba1c = av_t2_prehba1c, av_unadjust = av_t2_unadjust, av, av_lower, av_upper) %>% mutate(subgroup = "Type 2", estimate = NA, std.error = NA, statistic = NA, p.value = NA, conf.low = NA, conf.high = NA) %>% distinct()

HbA1c_results_ALL <- PEI_results_HbA1c %>% union(no_PEI_results_HbA1c) %>% union(T2_results_HbA1c) %>% union(no_PEI_haemochromatosis_HbA1c) %>% union(T2_haemochromatosis_HbA1c) %>%
  mutate(condition = ifelse(condition == "t3c_acutepancreatitis_only", "Acute pancreatitis only", ifelse(condition == "t3c_chronicpancreatitis", "Chronic pancreatitis", ifelse(condition == "t3c_pancreaticcancer", "Pancreatic cancer",
                                                                                                                                                                                ifelse(condition == "t3c_haemochromatosis", "Haemochromatosis", ifelse(condition == "t3c_all", "All T3c", NA))))))

HbA1c_results_tidied <- HbA1c_results_ALL %>%
  mutate(est_ci = ifelse(!is.na(estimate), paste0(round(estimate,digits =1)," (",round(conf.low,digits=1),",",round(conf.high,digits=1),")"), NA), p_value = ifelse(p.value <0.001, "<0.001", round(p.value, digits = 3))) %>%
  mutate(av = paste0(round(av,digits =1)," (",round(av_lower,digits=1),",",round(av_upper,digits=1),")")) %>%
  mutate(av_unadjust = round(av_unadjust, digits =1)) %>% mutate(av_prehba1c = round(av_prehba1c, digits =1)) %>%
  select(drug, condition, subgroup, 'N with non-missing outcome' = N, 'Average baseline HbA1c' = av_prehba1c, 'Average HbA1c response unadjusted' = av_unadjust, 'Average HbA1c response (95%CI)' = av, 'Average difference (95%CI)' = est_ci, 'p value' = p_value) %>%
  arrange(condition)

#Output HTML table
kableExtra::kable(HbA1c_results_tidied,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  column_spec(1,width="3cm") %>%
  column_spec(2,width="6cm") %>%
  column_spec(3,width="3cm") %>%
  column_spec(4,width="4cm") %>%
  column_spec(5,width="4cm") %>%
  column_spec(6,width="4cm") %>%
  column_spec(7,width="6cm") %>%
  column_spec(8,width="6cm") %>%
  column_spec(9,width="3cm") %>%
  cat(.,file="HbA1c_results.html")

#Reformat for supplementary tables
HbA1c_results_new_table <- HbA1c_results_tidied %>% filter(subgroup == "Type 2") %>% 
  select(drug, condition, 'N T2' = 'N with non-missing outcome', 'Average HbA1c response unadjusted T2' = 'Average HbA1c response unadjusted', 'Average HbA1c response (95%CI) T2' = 'Average HbA1c response (95%CI)') %>%
  left_join(HbA1c_results_tidied %>% filter(subgroup == "PEI") %>% 
              select(drug, condition, 'N PEI' = 'N with non-missing outcome', 'Average HbA1c response unadjusted PEI' = 'Average HbA1c response unadjusted', 'Average HbA1c response (95%CI) PEI' = 'Average HbA1c response (95%CI)', 'Average difference (95%CI) PEI' = 'Average difference (95%CI)')) %>%
  left_join(HbA1c_results_tidied %>% filter(subgroup == "No PEI") %>% 
              select(drug, condition, 'N No PEI' = 'N with non-missing outcome', 'Average HbA1c response unadjusted No PEI' = 'Average HbA1c response unadjusted', 'Average HbA1c response (95%CI) No PEI' = 'Average HbA1c response (95%CI)', 'Average difference (95%CI) No PEI' = 'Average difference (95%CI)'))

#Output HTML table
kableExtra::kable(HbA1c_results_new_table,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  column_spec(1,width="3cm") %>%
  column_spec(2,width="6cm") %>%
  column_spec(3,width="4cm") %>%
  column_spec(4,width="4cm") %>%
  column_spec(5,width="6cm") %>%
  column_spec(6,width="4cm") %>%
  column_spec(7,width="4cm") %>%
  column_spec(8,width="6cm") %>%
  column_spec(9,width="6cm") %>%
  column_spec(10,width="4cm") %>%
  column_spec(11,width="4cm") %>%
  column_spec(12,width="6cm") %>%
  column_spec(13,width="6cm") %>%
  cat(.,file="HbA1c_supplementary_table.html")


###Discontinuation
#Reformat discontinuation_results_ALL table so have 3 rows for each subgroup/drug: PEI, No PEI, and T2 (T2 estimate etc will be missing as is reference group)
PEI_results_discont <- discontinuation_results_ALL %>% filter(subgroup == "PEI") %>%
  select(drug, condition, subgroup, N = PEI, prob_unadjust = prob_PEI_unadjust, predict_est = predict_est_PEI, lower_predict = lower_predict_PEI, upper_predict = upper_predict_PEI, estimate, std.error, statistic, p.value, conf.low, conf.high)
no_PEI_results_discont <- discontinuation_results_ALL %>% filter(subgroup == "No PEI") %>%
  select(drug, condition, subgroup, N = no_PEI, prob_unadjust = prob_no_PEI_unadjust, predict_est = predict_est_no_PEI, lower_predict = lower_predict_no_PEI, upper_predict = upper_predict_no_PEI, estimate, std.error, statistic, p.value, conf.low, conf.high)
T2_results_discont <- discontinuation_results_ALL %>%
  select(drug, condition, N = t2, prob_unadjust = prob_t2_unadjust, predict_est = predict_est_t2, lower_predict = lower_predict_t2, upper_predict = upper_predict_t2) %>% mutate(subgroup = "Type 2", estimate = NA, std.error = NA, statistic = NA, p.value = NA, conf.low = NA, conf.high = NA) %>% distinct()
no_PEI_haemochromatosis_discont <- haemochromatosis_results_discont %>% filter(term == "No PEI") %>%
  select(drug, condition, subgroup = term, N = no_PEI, prob_unadjust = prob_no_PEI_unadjust, predict_est = predict_est_no_PEI, lower_predict = lower_predict_no_PEI, upper_predict = upper_predict_no_PEI, estimate, std.error, statistic, p.value, conf.low, conf.high)
T2_haemochromatosis_discont <- haemochromatosis_results_discont %>%
  select(drug, condition, N = t2, prob_unadjust = prob_t2_unadjust, predict_est = predict_est_t2, lower_predict = lower_predict_t2, upper_predict = upper_predict_t2) %>% mutate(subgroup = "Type 2", estimate = NA, std.error = NA, statistic = NA, p.value = NA, conf.low = NA, conf.high = NA) %>% distinct()


discontinuation_results_ALL <- PEI_results_discont %>% union(no_PEI_results_discont) %>% union(T2_results_discont) %>% union(no_PEI_haemochromatosis_discont) %>% union(T2_haemochromatosis_discont) %>%
  mutate(condition = ifelse(condition == "t3c_acutepancreatitis_only", "Acute pancreatitis only", ifelse(condition == "t3c_chronicpancreatitis", "Chronic pancreatitis", ifelse(condition == "t3c_pancreaticcancer", "Pancreatic cancer",
                                                                                                                                                                                ifelse(condition == "t3c_haemochromatosis", "Haemochromatosis", ifelse(condition == "t3c_all", "All T3c", NA))))))

discontinuation_results_tidied <- discontinuation_results_ALL %>%
  mutate(est_ci = ifelse(!is.na(estimate), paste0(round(estimate,digits =2)," (",round(conf.low,digits=2),",",round(conf.high,digits=2),")"), NA), p_value = ifelse(p.value <0.001, "<0.001", round(p.value, digits = 3))) %>%
  mutate(prob = paste0(round(predict_est*100,digits =1)," (",round(lower_predict*100,digits=1),",",round(upper_predict*100,digits=1),")")) %>%
  mutate(prob_unadjust = round(prob_unadjust *100, digits =1)) %>%
  select(drug, condition, subgroup, 'N with non-missing outcome' = N, 'Average discontinuation unadjusted, %' = prob_unadjust, 'Average discontinuation (95%CI), %' = prob, 'Odds ratio for discontinuation (95%CI)' = est_ci, 'p value' = p_value) %>%
  arrange(condition)

#Output HTML table
kableExtra::kable(discontinuation_results_tidied,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  column_spec(1,width="3cm") %>%
  column_spec(2,width="6cm") %>%
  column_spec(3,width="3cm") %>%
  column_spec(4,width="4cm") %>%
  column_spec(5,width="4cm") %>%
  column_spec(6,width="6cm") %>%
  column_spec(7,width="6cm") %>%
  column_spec(8,width="3cm") %>%
  cat(.,file="discontinuation_results.html")

#Reformat for supplementary tables
discontinuation_results_new_table <- discontinuation_results_tidied %>% filter(subgroup == "Type 2") %>% 
  select(drug, condition, 'N T2' = 'N with non-missing outcome', 'Average discontinuation unadjusted T2' = 'Average discontinuation unadjusted, %', 'Average discontinuation (95%CI) T2, %' = 'Average discontinuation (95%CI), %') %>%
  left_join(discontinuation_results_tidied %>% filter(subgroup == "PEI") %>% 
              select(drug, condition, 'N PEI' = 'N with non-missing outcome', 'Average discontinuation unadjusted PEI' = 'Average discontinuation unadjusted, %', 'Average discontinuation (95%CI) PEI, %' = 'Average discontinuation (95%CI), %', 'Odds ratio for discontinuation (95%CI) PEI' = 'Odds ratio for discontinuation (95%CI)')) %>%
  left_join(discontinuation_results_tidied %>% filter(subgroup == "No PEI") %>% 
              select(drug, condition, 'N No PEI' = 'N with non-missing outcome', 'Average discontinuation unadjusted No PEI' = 'Average discontinuation unadjusted, %', 'Average discontinuation (95%CI) No PEI, %' = 'Average discontinuation (95%CI), %', 'Odds ratio for discontinuation (95%CI) No PEI' = 'Odds ratio for discontinuation (95%CI)'))

#Output HTML table
kableExtra::kable(discontinuation_results_new_table,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  column_spec(1,width="3cm") %>%
  column_spec(2,width="6cm") %>%
  column_spec(3,width="4cm") %>%
  column_spec(4,width="4cm") %>%
  column_spec(5,width="6cm") %>%
  column_spec(6,width="4cm") %>%
  column_spec(7,width="4cm") %>%
  column_spec(8,width="6cm") %>%
  column_spec(9,width="6cm") %>%
  column_spec(10,width="4cm") %>%
  column_spec(11,width="4cm") %>%
  column_spec(12,width="6cm") %>%
  column_spec(13,width="6cm") %>%
  cat(.,file="discontinuation_supplementary_table.html")


###Weight
#Reformat weight_results_ALL table so have 3 rows for each subgroup/drug: PEI, No PEI, and T2 (T2 estimate etc will be missing as is reference group)
PEI_results_weight <- weight_results_ALL %>% filter(subgroup == "PEI") %>% mutate(av = as.numeric(av_PEI[,1]), av_lower = as.numeric(av_PEI[,2]), av_upper = as.numeric(av_PEI[,3])) %>%
  select(drug, condition, subgroup, N = PEI, av_preweight = av_PEI_preweight, av_unadjust = av_PEI_unadjust, av, av_lower, av_upper, estimate, std.error, statistic, p.value, conf.low, conf.high)
no_PEI_results_weight <- weight_results_ALL %>% filter(subgroup == "No PEI") %>% mutate(av = as.numeric(av_no_PEI[,1]), av_lower = as.numeric(av_no_PEI[,2]), av_upper = as.numeric(av_no_PEI[,3])) %>%
  select(drug, condition, subgroup, N = no_PEI, av_preweight = av_no_PEI_preweight, av_unadjust = av_no_PEI_unadjust, av, av_lower, av_upper, estimate, std.error, statistic, p.value, conf.low, conf.high)
T2_results_weight <- weight_results_ALL %>% mutate(av = as.numeric(av_t2[,1]), av_lower = as.numeric(av_t2[,2]), av_upper = as.numeric(av_t2[,3])) %>%
  select(drug, condition, N = t2, av_preweight = av_t2_preweight, av_unadjust = av_t2_unadjust, av, av_lower, av_upper) %>% mutate(subgroup = "Type 2", estimate = NA, std.error = NA, statistic = NA, p.value = NA, conf.low = NA, conf.high = NA) %>% distinct()
no_PEI_haemochromatosis_weight <- haemochromatosis_results_weight %>% filter(term == "No PEI") %>% mutate(av = as.numeric(av_no_PEI[,1]), av_lower = as.numeric(av_no_PEI[,2]), av_upper = as.numeric(av_no_PEI[,3])) %>%
  select(drug, condition, subgroup = term, N = no_PEI, av_preweight = av_no_PEI_preweight, av_unadjust = av_no_PEI_unadjust, av, av_lower, av_upper, estimate, std.error, statistic, p.value, conf.low, conf.high)
T2_haemochromatosis_weight <- haemochromatosis_results_weight %>% mutate(av = as.numeric(av_t2[,1]), av_lower = as.numeric(av_t2[,2]), av_upper = as.numeric(av_t2[,3])) %>%
  select(drug, condition, N = t2, av_preweight = av_t2_preweight, av_unadjust = av_t2_unadjust, av, av_lower, av_upper) %>% mutate(subgroup = "Type 2", estimate = NA, std.error = NA, statistic = NA, p.value = NA, conf.low = NA, conf.high = NA) %>% distinct()


weight_results_ALL <- PEI_results_weight %>% union(no_PEI_results_weight) %>% union(T2_results_weight) %>% union(no_PEI_haemochromatosis_weight) %>% union(T2_haemochromatosis_weight) %>%
  mutate(condition = ifelse(condition == "t3c_acutepancreatitis_only", "Acute pancreatitis only", ifelse(condition == "t3c_chronicpancreatitis", "Chronic pancreatitis", ifelse(condition == "t3c_pancreaticcancer", "Pancreatic cancer",
                                                                                                                                                                                ifelse(condition == "t3c_haemochromatosis", "Haemochromatosis", ifelse(condition == "t3c_all", "All T3c", NA))))))

weight_results_tidied <- weight_results_ALL %>%
  mutate(est_ci = ifelse(!is.na(estimate), paste0(round(estimate,digits =1)," (",round(conf.low,digits=1),",",round(conf.high,digits=1),")"), NA), p_value = ifelse(p.value <0.001, "<0.001", round(p.value, digits = 3))) %>%
  mutate(av = paste0(round(av,digits =1)," (",round(av_lower,digits=1),",",round(av_upper,digits=1),")")) %>%
  mutate(av_unadjust = round(av_unadjust, digits =1)) %>% mutate(av_preweight = round(av_preweight, digits =1)) %>%
  select(drug, condition, subgroup, 'N with non-missing outcome' = N, 'Average baseline weight' = av_preweight, 'Average weight response unadjusted' = av_unadjust, 'Average weight response (95%CI)' = av, 'Average difference (95%CI)' = est_ci, 'p value' = p_value) %>%
  arrange(condition)

#Output HTML table
kableExtra::kable(weight_results_tidied,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  column_spec(1,width="3cm") %>%
  column_spec(2,width="6cm") %>%
  column_spec(3,width="3cm") %>%
  column_spec(4,width="4cm") %>%
  column_spec(5,width="4cm") %>%
  column_spec(6,width="4cm") %>%
  column_spec(7,width="6cm") %>%
  column_spec(8,width="6cm") %>%
  column_spec(9,width="3cm") %>%
  cat(.,file="weight_results.html")

#Reformat for supplementary tables
weight_results_new_table <- weight_results_tidied %>% filter(subgroup == "Type 2") %>% 
  select(drug, condition, 'N T2' = 'N with non-missing outcome', 'Average weight response unadjusted T2' = 'Average weight response unadjusted', 'Average weight response (95%CI) T2' = 'Average weight response (95%CI)') %>%
  left_join(weight_results_tidied %>% filter(subgroup == "PEI") %>% 
              select(drug, condition, 'N PEI' = 'N with non-missing outcome', 'Average weight response unadjusted PEI' = 'Average weight response unadjusted', 'Average weight response (95%CI) PEI' = 'Average weight response (95%CI)', 'Average difference (95%CI) PEI' = 'Average difference (95%CI)')) %>%
  left_join(weight_results_tidied %>% filter(subgroup == "No PEI") %>% 
              select(drug, condition, 'N No PEI' = 'N with non-missing outcome', 'Average weight response unadjusted No PEI' = 'Average weight response unadjusted', 'Average weight response (95%CI) No PEI' = 'Average weight response (95%CI)', 'Average difference (95%CI) No PEI' = 'Average difference (95%CI)'))

#Output HTML table
kableExtra::kable(weight_results_new_table,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  column_spec(1,width="3cm") %>%
  column_spec(2,width="6cm") %>%
  column_spec(3,width="4cm") %>%
  column_spec(4,width="4cm") %>%
  column_spec(5,width="6cm") %>%
  column_spec(6,width="4cm") %>%
  column_spec(7,width="4cm") %>%
  column_spec(8,width="6cm") %>%
  column_spec(9,width="6cm") %>%
  column_spec(10,width="4cm") %>%
  column_spec(11,width="4cm") %>%
  column_spec(12,width="6cm") %>%
  column_spec(13,width="6cm") %>%
  cat(.,file="weight_supplementary_table.html")


################################################################################
###Plots

#Drop all other tables
rm(list=ls()[! ls() %in% c("HbA1c_results_ALL", "discontinuation_results_ALL", "weight_results_ALL")])

#Theme
plot_theme <- theme(legend.position = "bottom", legend.text = element_text(size =8), legend.title = element_blank(), legend.justification = "right") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=8), axis.text.y = element_text(size = 8),
        plot.subtitle = element_text(hjust=-0.7, lineheight = 0.25), plot.title = element_text(hjust=0.5, size =12, lineheight = 0.25), plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))

###Figure 1 - All drug all T3c, acute pancreatitis, chronic pancreatitis (HbA1c, discontinuation)

##All T3c
##HbA1c
Fig1_All_A <- ggplot((HbA1c_results_ALL %>% filter(drug == "overall" & condition == "All T3c")), aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "HbA1c change at 12 months (mmol/mol)", title = "All type 3c \n ", subtitle = "A) HbA1c response \n ") +
  ylim(-20,0) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("0.7 (0.4,1.0)", "3.5 (2.9,4.1)"), y.position = c(-15, -18),
               label.size = 2.75, tip.length = -0.02, vjust = 2.5) +
  theme(legend.position = "none") + geom_text(aes(label = paste0("N= ", N)), y =-1, size = 2.5, colour = "white")

Fig1_AP_A <- ggplot((HbA1c_results_ALL %>% filter(drug == "overall" & condition == "Acute pancreatitis only")), aes(y = av, x = subgroup)) +
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "HbA1c change at 12 months (mmol/mol)", title = "Acute pancreatitis \n ", subtitle = "") +
  ylim(-20,0) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("0.7 (0.3,1.1)", "2.0 (0.6,3.5)"), y.position = c(-15, -18),
               label.size = 2.75, tip.length = -0.02, vjust = 2.5) +
  theme(legend.position = "none") + geom_text(aes(label = paste0("N= ", N)), y =-1, size = 2.5, colour = "white")

Fig1_CP_A <-ggplot((HbA1c_results_ALL %>% filter(drug == "overall" & condition == "Chronic pancreatitis")), aes(y = av, x = subgroup)) +
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "HbA1c change at 12 months (mmol/mol)", title = "Chronic pancreatitis \n ", subtitle ="") +
  ylim(-20,0) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("1.8 (1.2,2.5)", "3.8 (3.0,4.6)"), y.position = c(-15, -18),
               label.size = 2.75, tip.length = -0.02, vjust = 2.5) +
  theme(legend.position = "none") + geom_text(aes(label = paste0("N= ", N)), y =-1, size = 2.5, colour = "white")

##Discontinuation
Fig1_All_B <- ggplot((discontinuation_results_ALL %>% filter(drug == "overall" & condition == "All T3c") %>% mutate(av = predict_est*100, av_lower = lower_predict*100, av_upper = upper_predict*100)), 
                     aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Discontinuation within 6 months (%)", title = "", subtitle = "B) Discontinuation \n ") +
  ylim(0, 45) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("OR 1.08 (0.98,1.19)", "OR 2.03 (1.73,2.36)"), y.position = c(26, 32),
               label.size = 2.75, tip.length = 0.02, vjust = -0.5) + 
  theme(legend.position = "none") + geom_text(aes(label = paste0("N= ", N)), y =2, size = 2.5, colour = "white")


Fig1_AP_B <- ggplot((discontinuation_results_ALL %>% filter(drug == "overall" & condition == "Acute pancreatitis only") %>% mutate(av = predict_est*100, av_lower = lower_predict*100, av_upper = upper_predict*100)), 
                    aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Discontinuation within 6 months (%)", title = "", subtitle = "") +
  ylim(0, 45) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("OR 1.06 (0.94,1.21)", "OR 1.15 (0.71,1.78)"), y.position = c(21, 27),
               label.size = 2.75, tip.length = 0.02, vjust = -0.5) +
  theme(legend.position = "none") + geom_text(aes(label = paste0("N= ", N)), y =2, size = 2.5, colour = "white")

Fig1_CP_B <- ggplot((discontinuation_results_ALL %>% filter(drug == "overall" & condition == "Chronic pancreatitis") %>% mutate(av = predict_est*100, av_lower = lower_predict*100, av_upper = upper_predict*100)), 
                    aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) + 
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Discontinuation within 6 months (%)", title = "", subtitle = "") +
  ylim(0, 45) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("OR 1.16 (0.98,1.39)", "OR 1.98 (1.63,2.41)"), y.position = c(26, 32),
               label.size = 2.75, tip.length = 0.02, vjust = -0.5) + geom_text(aes(label = paste0("N= ", N)), y =2, size = 2.5, colour = "white")

Fig1 <- Fig1_All_A + Fig1_AP_A + Fig1_CP_A + Fig1_All_B + Fig1_AP_B + Fig1_CP_B + plot_layout(ncol =3)
Fig1

##Save
pdf.options(reset = TRUE, onefile = TRUE)
pdf("Fig1_T3c_treatment_response_all_drugs_all_T3c_acute_chronic_pancreatitis.pdf",width=7.5,height=6)
Fig1
dev.off()


###Supplementary: All drugs, pancreatic cancer

##HbA1c
SFig_PC_A <-ggplot((HbA1c_results_ALL %>% filter(drug == "overall" & condition == "Pancreatic cancer")), aes(y = av, x = subgroup)) +
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "HbA1c change at 12 months (mmol/mol)", subtitle = "A) HbA1c response \n ") +
  ylim(-20,0) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("0.5 (-1.2,2.2)", "3.1 (1.3,4.9)"), y.position = c(-16, -19),
               label.size = 2.75, tip.length = -0.02, vjust = 2.5) +
  theme(legend.position = "none") + geom_text(aes(label = paste0("N= ", N)), y =-1, size = 2.5, colour = "white")

##Discontinuation
SFig_PC_B <- ggplot((discontinuation_results_ALL %>% filter(drug == "overall" & condition == "Pancreatic cancer") %>% mutate(av = predict_est*100, av_lower = lower_predict*100, av_upper = upper_predict*100)), 
                    aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Discontinuation within 6 months (%)", title = "", subtitle = "B) Discontinuation \n ") +
  ylim(0, 45) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("OR 1.10 (0.60,1.93)", "OR 2.60 (1.60,4.25)"), y.position = c(34, 40),
               label.size = 2.75, tip.length = 0.02, vjust = -0.5) + geom_text(aes(label = paste0("N= ", N)), y =2, size = 2.5, colour = "white")

SFig_PC <- SFig_PC_A + SFig_PC_B + plot_layout(ncol =1)
SFig_PC

##Save
pdf.options(reset = TRUE, onefile = TRUE)
pdf("SFig7_T3c_treatment_response_all_drugs_pancreatic_cancer.pdf",width=2.5,height=6)
SFig_PC
dev.off()


###Supplementary: All drugs, haemochromatosis (only no PEI vs T2)

##HbA1c
SFig_HC_A <-ggplot((HbA1c_results_ALL %>% filter(drug == "overall" & condition == "Haemochromatosis")), aes(y = av, x = subgroup)) +
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "HbA1c change at 12 months (mmol/mol)", subtitle = "A) HbA1c response \n ") +
  ylim(-20,0) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI"), xmax = c("Type 2"), label = c("-1.9 (-2.8,-0.9)"), y.position = c(-17),
               label.size = 2.75, tip.length = -0.02, vjust = 2.5) +
  geom_text(x = "PEI", y = -1, label = "*") +
  theme(legend.position = "none") + geom_text(aes(label = paste0("N= ", N)), y =-1, size = 2.5, colour = "white")

##Discontinuation
SFig_HC_B <- ggplot((discontinuation_results_ALL %>% filter(drug == "overall" & condition == "Haemochromatosis") %>% mutate(av = predict_est*100, av_lower = lower_predict*100, av_upper = upper_predict*100)), 
                    aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI"), limits = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Discontinuation within 6 months (%)", title = "", subtitle = "B) Discontinuation \n ") +
  ylim(0, 45) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI"), xmax = c("Type 2"), label = c("OR 1.04 (0.72,1.50)"), y.position = c(21),
               label.size = 2.75, tip.length = 0.02, vjust = -0.5) +
  geom_text(x = "PEI", y = 1, label = "*") + geom_text(aes(label = paste0("N= ", N)), y =2, size = 2.5, colour = "white")

SFig_HC <- SFig_HC_A + SFig_HC_B + plot_layout(ncol =1)
SFig_HC

##Save
pdf.options(reset = TRUE, onefile = TRUE)
pdf("SFig8_T3c_treatment_response_all_drugs_haemochromatosis.pdf",width=2.5,height=6)
SFig_HC
dev.off()


###Figure 2 - All T3c by drug class (HbA1c, discontinuation, weight)

plot_theme <- plot_theme + theme(legend.text = element_text(size =10), axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10))

##HbA1c
Fig2_MFN_A <- ggplot((HbA1c_results_ALL %>% filter(drug == "MFN" & condition == "All T3c")), aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "HbA1c change at 12 months (mmol/mol)", title = "Metformin \n ", subtitle = "A) HbA1c response \n ") +
  ylim(-23,0) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("0.2 (-0.1,0.5)", "4.2 (3.5,5.0)"), y.position = c(-19, -22),
               label.size = 3.25, tip.length = -0.02, vjust = 2.25) +
  theme(legend.position = "none")

Fig2_SU_A <- ggplot((HbA1c_results_ALL %>% filter(drug == "SU" & condition == "All T3c")), aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "HbA1c change at 12 months (mmol/mol)", title = "Sulphonylureas \n ", subtitle = "") +
  ylim(-23,0) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("1.6 (0.9,2.4)", "4.2 (2.8,5.5)"), y.position = c(-17, -20),
               label.size = 3.25, tip.length = -0.02, vjust = 2.25) +
  theme(legend.position = "none")

Fig2_TZD_A <- ggplot((HbA1c_results_ALL %>% filter(drug == "TZD" & condition == "All T3c")), aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "HbA1c change at 12 months (mmol/mol)", title = "TZDs \n ", subtitle = "") +
  ylim(-23,0) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("0.3 (-1.6,2.3)", "2.3 (-1.9,6.4)"), y.position = c(-17, -20),
               label.size = 3.25, tip.length = -0.02, vjust = 2.25) +
  theme(legend.position = "none")

Fig2_SGLT2_A <- ggplot((HbA1c_results_ALL %>% filter(drug == "SGLT2" & condition == "All T3c")), aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "HbA1c change at 12 months (mmol/mol)", title = "SGLT2-inhibitors \n ", subtitle = "") +
  ylim(-23,0) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("1.1 (-0.1,2.3)", "-1.4 (-3.8,1.0)"), y.position = c(-15, -18),
               label.size = 3.25, tip.length = -0.02, vjust = 2.25) +
  theme(legend.position = "none")

Fig2_DPP4_A <- ggplot((HbA1c_results_ALL %>% filter(drug == "DPP4" & condition == "All T3c")), aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "HbA1c change at 12 months (mmol/mol)", title = "DPP4-inhibitors \n ", subtitle = "") +
  ylim(-23,0) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("0.8 (0.0,1.7)", "1.3 (-0.6,3.2)"), y.position = c(-11, -14),
               label.size = 3.25, tip.length = -0.02, vjust = 2.25) +
  theme(legend.position = "none")

##Discontinuation
Fig2_MFN_B <- ggplot((discontinuation_results_ALL %>% filter(drug == "MFN" & condition == "All T3c") %>% mutate(av = predict_est*100, av_lower = lower_predict*100, av_upper = upper_predict*100)), 
                     aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Discontinuation within 6 months (%)", title = "", subtitle = "B) Discontinuation \n ") +
  ylim(0, 46) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("OR 1.26 (1.05,1.50)", "OR 2.97 (2.29,3.84)"), y.position = c(16, 22),
               label.size = 3.25, tip.length = 0.02, vjust = -0.5) +
  theme(legend.position = "none")

Fig2_SU_B <- ggplot((discontinuation_results_ALL %>% filter(drug == "SU" & condition == "All T3c") %>% mutate(av = predict_est*100, av_lower = lower_predict*100, av_upper = upper_predict*100)), 
                     aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Discontinuation within 6 months (%)", title = "", subtitle = "") +
  ylim(0, 46) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("OR 0.99 (0.84,1.17)", "OR 1.51 (1.17,1.94)"), y.position = c(28, 34),
               label.size = 3.25, tip.length = 0.02, vjust = -0.5) +
  theme(legend.position = "none")

Fig2_TZD_B <- ggplot((discontinuation_results_ALL %>% filter(drug == "TZD" & condition == "All T3c") %>% mutate(av = predict_est*100, av_lower = lower_predict*100, av_upper = upper_predict*100)), 
                     aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Discontinuation within 6 months (%)", title = "", subtitle = "") +
  ylim(0, 46) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("OR 1.04 (0.68,1.58)", "OR 1.37 (0.58,2.99)"), y.position = c(38, 44),
               label.size = 3.25, tip.length = 0.02, vjust = -0.5) +
  theme(legend.position = "none")

Fig2_SGLT2_B <- ggplot((discontinuation_results_ALL %>% filter(drug == "SGLT2" & condition == "All T3c") %>% mutate(av = predict_est*100, av_lower = lower_predict*100, av_upper = upper_predict*100)), 
                     aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Discontinuation within 6 months (%)", title = "", subtitle = "") +
  ylim(0, 46) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("OR 0.95 (0.71,1.27)", "OR 1.41 (0.81,2.36)"), y.position = c(36, 42),
               label.size = 3.25, tip.length = 0.02, vjust = -0.5) +
  theme(legend.position = "none")

Fig2_DPP4_B <- ggplot((discontinuation_results_ALL %>% filter(drug == "DPP4" & condition == "All T3c") %>% mutate(av = predict_est*100, av_lower = lower_predict*100, av_upper = upper_predict*100)), 
                     aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Discontinuation within 6 months (%)", title = "", subtitle = "") +
  ylim(0, 46) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("OR 1.14 (0.91,1.41)", "OR 1.64 (1.07,2.47)"), y.position = c(32, 38),
               label.size = 3.25, tip.length = 0.02, vjust = -0.5) +
  theme(legend.position = "none")


##Weight
Fig2_MFN_C <- ggplot((weight_results_ALL %>% filter(drug == "MFN" & condition == "All T3c")), aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Weight change at 12 months (kg)", title = "", subtitle = "C) Weight change \n ") +
  ylim(-7,5) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("-0.1 (-0.3,0.0)", "-0.4 (-0.8,-0.1)"), y.position = c(-4, -5.5),
               label.size = 3.25, tip.length = -0.02, vjust = 2.25) +
  theme(legend.position = "none")

Fig2_SU_C <- ggplot((weight_results_ALL %>% filter(drug == "SU" & condition == "All T3c")), aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Weight change at 12 months (kg)", title = "", subtitle = "") +
  ylim(-7,5) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("-0.2 (-0.4,0.0)", "-1.1 (-1.5,-0.7)"), y.position = c(2, 3.5),
               label.size = 3.25, tip.length = 0.02, vjust = -0.5) +
  theme(legend.position = "none")

Fig2_TZD_C <- ggplot((weight_results_ALL %>% filter(drug == "TZD" & condition == "All T3c")), aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Weight change at 12 months (kg)", title = "", subtitle = "") +
  ylim(-7,5) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("-0.3 (-1.0,0.4)", "-1.0 (-2.5,0.5)"), y.position = c(3, 4.5),
               label.size = 3.25, tip.length = 0.02, vjust = -0.5) +
  theme(legend.position = "none")

Fig2_SGLT2_C <- ggplot((weight_results_ALL %>% filter(drug == "SGLT2" & condition == "All T3c")), aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Weight change at 12 months (kg)", title = "", subtitle = "") +
  ylim(-7,5) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("0.0 (-0.5,0.4)", "-1.4 (-2.3,-0.6)"), y.position = c(-5, -6.5),
               label.size = 3.25, tip.length = -0.02, vjust = 2.25)

Fig2_DPP4_C <- ggplot((weight_results_ALL %>% filter(drug == "DPP4" & condition == "All T3c")), aes(y = av, x = subgroup)) + 
  geom_col(position="dodge", aes(fill = subgroup)) + scale_x_discrete(limits = c("Type 2", "No PEI", "PEI")) + scale_fill_manual(values = c("No PEI" = "#1B9E77", "PEI" = "#D95F02", "Type 2" = "#7570B3"), breaks = c("Type 2", "No PEI", "PEI")) +
  geom_errorbar(aes(ymin= av_lower, ymax=av_upper), width =.2, position=position_dodge(.9)) +
  labs(y = "Weight change at 12 months (kg)", title = "", subtitle = "") +
  ylim(-7,5) + geom_hline(aes(yintercept = 0)) + plot_theme +
  geom_bracket(xmin = c("No PEI", "PEI"), xmax = c("Type 2", "Type 2"), label = c("0.4 (0.1,0.7)", "-0.2 (-0.8,0.5)"), y.position = c(-2.5, -4),
               label.size = 3.25, tip.length = -0.02, vjust = 2.25) +
  theme(legend.position = "none")


##Combine
Fig2 <- Fig2_MFN_A + Fig2_SU_A + Fig2_TZD_A + Fig2_DPP4_A + Fig2_SGLT2_A +
  Fig2_MFN_B + Fig2_SU_B + Fig2_TZD_B + Fig2_DPP4_B + Fig2_SGLT2_B +
  Fig2_MFN_C + Fig2_SU_C + Fig2_TZD_C + Fig2_DPP4_C + Fig2_SGLT2_C +
  plot_layout(ncol =5)

##Save
pdf.options(reset = TRUE, onefile = TRUE)
pdf("Fig2_T3c_treatment_response_all_T3c_by_drug.pdf",width=12.5,height=9.5)
Fig2
dev.off()


################################################################################
###END
rm(list=ls())

