################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#First use command prompt to log in to Slade using SSH

#To install package
#library(devtools)
#install_github("Exeter-Diabetes/CPRD-analysis-package")


#Load packages
library(tidyverse)
library(MatchIt)
library(tableone)
library(kableExtra)

#Load aurum package
library(aurum)

################################################################################
###SETUP########################################################################

################################################################################
###Connecting to data and setting up analysis###################################
#Initialise connection
cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "C:/Users/rh530/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "31/10/2021")
analysis = cprd$analysis("Rhian_T3c")
################################################################################

#Load treatment response cohort
mm_all <- mm_all %>% analysis$cached("cohort_20231108_all_1stinstance", indexes=c("patid", "dstartdate", "drugclass"))

mm_all <- mm_all %>% select(patid, drugclass, type, diabetes_type, t3c_acutepancreatitis_only, t3c_chronicpancreatitis, t3c_pancreaticcancer, t3c_haemochromatosis, t3c_subgroup,
                            t3c_all, gender, ethnicity_5cat, imd2015_10, dstartdate_age, dstartdate, drugline_all, INS, numdrugs, prehba1c, prebmi, preweight, 
                            first_PEI_all, PEI_ever, PEI_prior, first_pancreatic_disease_all_post_diagnosis, alcohol_cat)

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


#Drop all other tables
rm(list=ls()[! ls() %in% c("overall_matched_cohort")])

################################################################################
#Make table 
#To include:
#Baseline characteristics: sex, age (median), ethnicity, IMD, HbA1c, BMI - at drug initiation
#T3c subtype
#Medication initiated

cohort <- overall_matched_cohort

#Setting variables to factors
#Gender
cohort$gender <- factor(cohort$gender)
levels(cohort$gender) = c("Male", "Female")
cohort$gender <- relevel(cohort$gender, ref="Female")
cohort$femalegender <- relevel(factor(cohort$gender), ref = "Male")
#Ethnicity
cohort$ethnicity_5cat <- factor(cohort$ethnicity_5cat)
levels(cohort$ethnicity_5cat) = c("White", "South Asian", "Black", "Other", "Mixed", "Unknown")
#Deprivation
cohort$imd_quintile <- factor(cohort$imd_quintile)
#T3c group
cohort$t3c_subgroup <- factor(cohort$t3c_subgroup)
#Medication initiated
cohort$drugclass <- factor(cohort$drugclass)
#Number of other current glucose lowering therapies (recode numdrugs -1)
cohort <- cohort %>% mutate(numdrugs = numdrugs -1)
cohort <- cohort %>% mutate(numdrugs_cat = ifelse(numdrugs>=2, "2+ drugs", ifelse(numdrugs ==1, "1 drug", ifelse(numdrugs == 0, "0 drugs", numdrugs))))
cohort$numdrugs_cat <- factor(cohort$numdrugs_cat)
#Alcohol consumption
cohort <- cohort %>% mutate(alcohol_cat = ifelse(is.na(alcohol_cat), "Unknown alcohol", alcohol_cat))
cohort$alcohol_cat <- factor(cohort$alcohol_cat)
cohort$alcohol_cat <- relevel(factor(cohort$alcohol_cat), ref = "Within limits")
cohort$alcohol_cat <- relevel(factor(cohort$alcohol_cat), ref = "None")


#Specify variables
all_vars <- c("gender", "femalegender", "dstartdate_age", "ethnicity_5cat", "imd_quintile", "alcohol_cat", "prehba1c", "prebmi", "preweight", "numdrugs_cat", "subgroup", "t3c_subgroup")
categorical_vars <-c("gender", "femalegender", "ethnicity_5cat", "imd_quintile", "alcohol_cat", "numdrugs_cat", "subgroup", "t3c_subgroup")


#Metformin
tableone_MFN <- CreateTableOne(vars=all_vars,data=(cohort %>% filter(drugclass == "MFN")), strata = "type", factorVars=categorical_vars, test=FALSE)
tab_MFN <-as_tibble(print(tableone_MFN, nonnormal = c("dstartdate_age", "prehba1c", "prebmi", "preweight"))) %>%
  add_column(measure=row.names(print(tableone_MFN))) %>%
  mutate(measure = trimws(measure)) %>%
  select(measure, T2, T3c) %>%
  rename('MFN T2 controls' = T2, 'MFN T3c' = T3c)
summary(tableone_MFN)

#SU
tableone_SU <- CreateTableOne(vars=all_vars,data=(cohort %>% filter(drugclass == "SU")), strata = "type", factorVars=categorical_vars, test=FALSE)
tab_SU <-as_tibble(print(tableone_SU, nonnormal = c("dstartdate_age", "prehba1c", "prebmi", "preweight"))) %>%
  add_column(measure=row.names(print(tableone_SU))) %>%
  mutate(measure = trimws(measure)) %>%
  select(measure, T2, T3c) %>%
  rename('SU T2 controls' = T2, 'SU T3c' = T3c)
summary(tableone_SU)

#TZD
tableone_TZD <- CreateTableOne(vars=all_vars,data=(cohort %>% filter(drugclass == "TZD")), strata = "type", factorVars=categorical_vars, test=FALSE)
tab_TZD <-as_tibble(print(tableone_TZD, nonnormal = c("dstartdate_age", "prehba1c", "prebmi", "preweight"))) %>%
  add_column(measure=row.names(print(tableone_TZD))) %>%
  mutate(measure = trimws(measure)) %>%
  select(measure, T2, T3c) %>%
  rename('TZD T2 controls' = T2, 'TZD T3c' = T3c)
summary(tableone_TZD)

#SGLT2
tableone_SGLT2 <- CreateTableOne(vars=all_vars,data=(cohort %>% filter(drugclass == "SGLT2")), strata = "type", factorVars=categorical_vars, test=FALSE)
tab_SGLT2 <-as_tibble(print(tableone_SGLT2, nonnormal = c("dstartdate_age", "prehba1c", "prebmi", "preweight"))) %>%
  add_column(measure=row.names(print(tableone_SGLT2))) %>%
  mutate(measure = trimws(measure)) %>%
  select(measure, T2, T3c) %>%
  rename('SGLT2 T2 controls' = T2, 'SGLT2 T3c' = T3c)
summary(tableone_SGLT2)

#DPP4
tableone_DPP4 <- CreateTableOne(vars=all_vars,data=(cohort %>% filter(drugclass == "DPP4")), strata = "type", factorVars=categorical_vars, test=FALSE)
tab_DPP4 <-as_tibble(print(tableone_DPP4, nonnormal = c("dstartdate_age", "prehba1c", "prebmi", "preweight"))) %>%
  add_column(measure=row.names(print(tableone_DPP4))) %>%
  mutate(measure = trimws(measure)) %>%
  select(measure, T2, T3c) %>%
  rename('DPP4 T2 controls' = T2, 'DPP4 T3c' = T3c)
summary(tableone_DPP4)


#Combine and tidy
tab <- tab_MFN %>% left_join(tab_SU) %>% left_join(tab_TZD) %>% left_join(tab_SGLT2) %>% left_join(tab_DPP4) %>%
  mutate(measure=ifelse(measure=="gender = Male (%)","Male", measure)) %>%
  mutate(measure=ifelse(measure=="femalegender = Female (%)","Female", measure)) %>%
  mutate(measure = ifelse(measure == "1", "1 (least deprived)", measure)) %>%
  mutate(measure = ifelse(measure == "5", "5 (most deprived)", measure)) %>%
  mutate(measure=ifelse(measure=="Missing IMD","Missing", measure)) %>%
  mutate(measure=ifelse(measure=="prehba1c (mean (SD))","Median [IQR]", measure)) %>%
  mutate(measure=ifelse(measure=="prebmi (mean (SD))","Median [IQR]", measure)) %>%
  mutate(measure=ifelse(measure=="preweight (mean (SD))","Median [IQR]", measure)) %>%
  mutate(measure=ifelse(measure=="dstartdate_age (mean (SD))","Median [IQR]", measure)) %>%
  mutate(measure = ifelse(measure == "0 drugs", "0", measure)) %>%
  mutate(measure = ifelse(measure == "1 drug", "1", measure)) %>%
  mutate(measure = ifelse(measure == "2+ drugs", "2+", measure)) %>%
  mutate(measure=ifelse(measure=="Unknown alcohol","Unknown", measure)) %>%
  filter(measure != "ethnicity_5cat (%)" & measure != "imd_quintile (%)" & measure != "t3c_subgroup (%)" & measure != "Type 2" &
           measure != "subgroup (%)" & measure != "T2" & measure != "Missing" & measure != "numdrugs_cat (%)" & measure != "alcohol_cat (%)") %>%
  mutate(measure = ifelse(measure == "T3c PEI", "PEI", ifelse(measure == "T3c No PEI", "No PEI", measure))) %>%
  rename(" " = measure)

#Outputting HTML table
kableExtra::kable(tab,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  pack_rows("Sex", 2,3, bold = TRUE) %>%
  pack_rows("Age, years",4,4, bold=TRUE) %>% 
  pack_rows("Ethnicity",5,10, bold=TRUE)  %>% 
  pack_rows("Index of multiple deprivation quintile",11,15, bold=TRUE)  %>%
  pack_rows("Alcohol consumption",16,20, bold=TRUE)  %>%
  pack_rows("HbA1c, mmol/mol",21,21, bold=TRUE)  %>%
  pack_rows("BMI, kg/m2",22,22, bold=TRUE)  %>%
  pack_rows("Weight, kg",23,23, bold=TRUE)  %>%
  pack_rows("Number of other glucose-lowering therapies being taken",24,26, bold=TRUE)  %>%
  pack_rows("PEI status", 27,28, bold = TRUE) %>%
  pack_rows("3c subtype", 29,32, bold = TRUE) %>%
  cat(.,file="STab1_matched_cohorts_characteristics_by_drug_class.html")

################################################################################
###END
rm(list=ls())
