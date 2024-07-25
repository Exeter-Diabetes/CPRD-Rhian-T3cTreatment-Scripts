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

#Load at diagnosis cohort (T3c/T2 only)
at_diagnosis <- at_diagnosis  %>% analysis$cached("cohort_at_diagnosis", unique_indexes = "patid")
cohort <- collect(at_diagnosis %>% mutate(patid=as.character(patid)) %>% filter(type != "T1"))

#Make a categorical subgroup variable: T3c PEI, T3c no PEI, or T2
cohort <- cohort %>% mutate(subgroup = ifelse(type == "T2", "T2", ifelse(type == "T3c" & PEI_prior ==1, "T3c PEI", ifelse(type == "T3c" & PEI_prior ==0, "T3c No PEI", NA))))


################################################################################
#Make table 
#To include:
#Baseline characteristics: sex, age (median), ethnicity, IMD, HbA1c, BMI - at drug initiation

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
#Group IMD into quintiles
cohort <- cohort %>% mutate(imd_quintile = ifelse(imd2015_10 == 1 | imd2015_10 == 2, "1", ifelse(imd2015_10 == 3 | imd2015_10 == 4, "2", ifelse(imd2015_10 == 5 | imd2015_10 == 6, "3", ifelse(imd2015_10 == 7 | imd2015_10 == 8, "4", ifelse(imd2015_10 == 9 | imd2015_10 == 10, "5", "Missing"))))))
cohort$imd_quintile <- factor(cohort$imd_quintile)
#T3c group
cohort$t3c_subgroup <- factor(cohort$t3c_subgroup)


#Specify variables
all_vars <- c("gender", "femalegender", "dm_diag_age", "ethnicity_5cat", "imd_quintile", "prehba1c", "prebmi", "t3c_subgroup")
categorical_vars <-c("gender", "femalegender", "ethnicity_5cat", "imd_quintile", "t3c_subgroup")

#Type 2 matched controls
tableoneT2 <- CreateTableOne(vars=all_vars,data=(cohort %>% filter(subgroup == "T2")),factorVars=categorical_vars, test=FALSE)
tab_T2 <-as_tibble(print(tableoneT2, nonnormal = c("dm_diag_age", "prehba1c", "prebmi"))) %>%
  add_column(measure=row.names(print(tableoneT2))) %>%
  mutate(measure = trimws(measure)) %>%
  rename(T2 = Overall)

#All T3c
tableoneT3c <- CreateTableOne(vars=all_vars,data=(cohort %>% filter(type == "T3c")),factorVars=categorical_vars, test=FALSE)
tab_T3c <-as_tibble(print(tableoneT3c, nonnormal = c("dm_diag_age", "prehba1c", "prebmi"))) %>%
  add_column(measure=row.names(print(tableoneT3c))) %>%
  mutate(measure = trimws(measure)) %>%
  rename(T3c = Overall)

#PEI
tableonePEI <- CreateTableOne(vars=all_vars,data=(cohort %>% filter(subgroup == "T3c PEI")),factorVars=categorical_vars, test=FALSE)
tab_PEI <-as_tibble(print(tableonePEI, nonnormal = c("dm_diag_age", "prehba1c", "prebmi"))) %>%
  add_column(measure=row.names(print(tableonePEI))) %>%
  mutate(measure = trimws(measure)) %>%
  rename(PEI = Overall)

#No PEI
tableone_noPEI <- CreateTableOne(vars=all_vars,data=(cohort %>% filter(subgroup == "T3c No PEI")),factorVars=categorical_vars, test=FALSE)
tab_no_PEI <-as_tibble(print(tableone_noPEI, nonnormal = c("dm_diag_age", "prehba1c", "prebmi"))) %>%
  add_column(measure=row.names(print(tableone_noPEI))) %>%
  mutate(measure = trimws(measure)) %>%
  rename('No PEI' = Overall)

#Combine
tab <- tab_T3c %>% select(measure, T3c) %>% left_join(tab_PEI) %>% left_join(tab_no_PEI) %>% left_join(tab_T2) %>% select(measure, T2, T3c, PEI, 'No PEI') %>%
  mutate(measure=ifelse(measure=="gender = Male (%)","Male", measure)) %>%
  mutate(measure=ifelse(measure=="femalegender = Female (%)","Female", measure)) %>%
  mutate(measure = ifelse(measure == "1", "1 (least deprived)", measure)) %>%
  mutate(measure = ifelse(measure == "5", "5 (most deprived)", measure)) %>%
  mutate(measure=ifelse(measure=="Missing IMD","Missing", measure)) %>%
  mutate(measure=ifelse(measure=="prehba1c (mean (SD))","Median [IQR]", measure)) %>%
  mutate(measure=ifelse(measure=="prebmi (mean (SD))","Median [IQR]", measure)) %>%
  mutate(measure=ifelse(measure=="dm_diag_age (mean (SD))","Median [IQR]", measure)) %>%
  rename(" " = measure) %>% filter(T3c != " ") %>% rename ('Type 2' = T2, 'Type 3c Overall' = T3c)

#Outputting HTML table
kableExtra::kable(tab,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  pack_rows("Sex", 2,3, bold = TRUE) %>%
  pack_rows("Age, years",4,4, bold=TRUE) %>% 
  pack_rows("Ethnicity",5,9, bold=TRUE)  %>% 
  pack_rows("Index of multiple deprivation quintile",10,14, bold=TRUE)  %>%
  pack_rows("HbA1c, mmol/mol",15,15, bold=TRUE)  %>%
  pack_rows("BMI, kg/m2",16,16, bold=TRUE)  %>%
  pack_rows("3c subtype", 17,20, bold = TRUE) %>%
  column_spec(1,width="12cm") %>%
  column_spec(2,width="6cm") %>%
  column_spec(3,width="6cm") %>%
  column_spec(4,width="6cm") %>%
  column_spec(5,width="6cm") %>%
  cat(.,file="Tab1_T3c_baseline_characteristics_table_at_diagnosis.html")

##Table 1 - next add treatment after diagnosis (insulin/ oral agent initiation) info

################################################################################
###END
rm(list=ls())
