################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#First use command prompt to log in to Slade using SSH

#To install package
#library(devtools)
#install_github("Exeter-Diabetes/CPRD-analysis-package")


#Load tidyverse and lubridate
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)

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

#Load at diagnosis cohort
at_diagnosis_t3c_t2 <- at_diagnosis_t3c_t2  %>% analysis$cached("cohort_at_diagnosis", unique_indexes = "patid")
cohort <- collect(at_diagnosis_t3c_t2 %>% mutate(patid=as.character(patid)) %>% filter(type != "T1"))

#Make a categorical subgroup variable: T3c PEI, T3c no PEI, or T2
cohort <- cohort %>% mutate(subgroup = ifelse(type == "T2", "T2", ifelse(type == "T3c" & PEI_prior ==1, "T3c PEI", ifelse(type == "T3c" & PEI_prior ==0, "T3c No PEI", NA))))

#Calculate survival dates/ times 
#dm_diag_ohadate is equivalent to first insulin prescription date
cohort <- cohort %>% mutate(survival_date = dm_diag_date + days(1095)) %>%
  mutate(survival_date = as.Date(ifelse(!is.na(dm_diag_ohadate) & dm_diag_ohadate < survival_date, dm_diag_ohadate, survival_date))) %>%
  mutate(survival_date = as.Date(ifelse(!is.na(gp_record_end) & gp_record_end < survival_date, gp_record_end, survival_date))) %>% 
  mutate(survival_date = as.Date(ifelse(!is.na(death_date) & death_date < survival_date, death_date, survival_date))) %>% 
  mutate(outcome = ifelse(!is.na(dm_diag_ohadate) & dm_diag_ohadate == survival_date, 1, 0)) %>% mutate(survival_time = as.numeric(survival_date - dm_diag_date)) %>% 
  mutate(survival_time_years = survival_time/365)

#Check survival time distribution for those with oral therapy outcome
cohort %>% filter(outcome ==1) %>% ggplot(aes(survival_time)) +  geom_histogram(binwidth = 5, colour= "white", fill="plum3")
cohort %>% filter(outcome ==1 & type == "T2") %>% ggplot(aes(survival_time)) +  geom_histogram(binwidth = 5, colour= "white", fill="plum3")
cohort %>% filter(outcome ==1 & subgroup == "T3c PEI") %>% ggplot(aes(survival_time)) +  geom_histogram(binwidth = 5, colour= "white", fill="plum3")
cohort %>% filter(outcome ==1 & subgroup == "T3c No PEI") %>% ggplot(aes(survival_time)) +  geom_histogram(binwidth = 5, colour= "white", fill="plum3")


#Survival analysis

#T2
t2_sfit <- survfit(Surv(survival_time_years,outcome) ~ 1, data = (cohort %>% filter(type == "T2")))
##1y
summary(t2_sfit, time =1)
100 - (summary(t2_sfit, time =1)$surv)*100
100 - (summary(t2_sfit, time =1)$lower)*100
100 - (summary(t2_sfit, time =1)$upper)*100
##3y
summary(t2_sfit, time =3)
100 - (summary(t2_sfit, time =3)$surv)*100
100 - (summary(t2_sfit, time =3)$lower)*100
100 - (summary(t2_sfit, time =3)$upper)*100

#T3c
t3c_sfit <- survfit(Surv(survival_time_years,outcome) ~ 1, data = (cohort %>% filter(type == "T3c")))
##1y
summary(t3c_sfit, time =1)
100 - (summary(t3c_sfit, time =1)$surv)*100
100 - (summary(t3c_sfit, time =1)$lower)*100
100 - (summary(t3c_sfit, time =1)$upper)*100
##3y
summary(t3c_sfit, time =3)
100 - (summary(t3c_sfit, time =3)$surv)*100
100 - (summary(t3c_sfit, time =3)$lower)*100
100 - (summary(t3c_sfit, time =3)$upper)*100

##By PEI
#PEI
t3c_PEI_sfit <- survfit(Surv(survival_time_years,outcome) ~ 1, data = (cohort %>% filter(subgroup == "T3c PEI")))
##1y
summary(t3c_PEI_sfit, time =1)
100 - (summary(t3c_PEI_sfit, time =1)$surv)*100
100 - (summary(t3c_PEI_sfit, time =1)$lower)*100
100 - (summary(t3c_PEI_sfit, time =1)$upper)*100
##3y
summary(t3c_PEI_sfit, time =3)
100 - (summary(t3c_PEI_sfit, time =3)$surv)*100
100 - (summary(t3c_PEI_sfit, time =3)$lower)*100
100 - (summary(t3c_PEI_sfit, time =3)$upper)*100

#No PEI
t3c_no_PEI_sfit <- survfit(Surv(survival_time_years,outcome) ~ 1, data = (cohort %>% filter(subgroup == "T3c No PEI")))
##1y
summary(t3c_no_PEI_sfit, time =1)
100 - (summary(t3c_no_PEI_sfit, time =1)$surv)*100
100 - (summary(t3c_no_PEI_sfit, time =1)$lower)*100
100 - (summary(t3c_no_PEI_sfit, time =1)$upper)*100
##3y
summary(t3c_no_PEI_sfit, time =3)
100 - (summary(t3c_no_PEI_sfit, time =3)$surv)*100
100 - (summary(t3c_no_PEI_sfit, time =3)$lower)*100
100 - (summary(t3c_no_PEI_sfit, time =3)$upper)*100

###KM plot
fit <- list(T2 = t2_sfit, PEI = t3c_PEI_sfit, no_PEI = t3c_no_PEI_sfit)
KM_plot_PEI <- ggsurvplot(fit = fit, 
           fun = function(x) {100 - x*100},
           palette = c("#7570B3", "#D95F02", "#1B9E77"),
           combine = TRUE,
           xlab = "Years from diagnosis",
           ylab = "Cumulative incidence of oral therapy initiation (%)",
           ylim = c(0, 70),
           xlim = c(0,3),
           break.time.by = 1,
           font.ticks = 16,
           font.x = 16,
           font.y = 16,
           font.main = 16,
           font.legend = 16,
           risk.table.font = 6,
           legend.title = "",
           legend = "bottom",
           legend.labs = c("Type 2", "PEI", "No PEI"),
           censor = F,
           conf.int = T,
           risk.table = TRUE,
           fontsize = 6,
           tables.y.text = FALSE,
           size= 1.5,
           ggtheme = theme_classic(),
           linetype = c("strata"),
           tables.theme = theme_cleantable(),
           axes.offset = FALSE)

KM_plot_PEI

#Save
pdf.options(reset = TRUE, onefile = FALSE)
pdf("SFig5_Time_to_oral_therapy_3_year_by_PEI.pdf",width=12.5,height=9)
KM_plot_PEI
dev.off()

###Cox model
cohort$subgroup <- factor(cohort$subgroup)
model_PEI <- coxph(Surv(survival_time, outcome) ~ subgroup, data = cohort)
summary(model_PEI)


##By T3c subtype
#Acute pancreatitis
AP_sfit <- survfit(Surv(survival_time_years,outcome) ~ 1, data = (cohort %>% filter(t3c_acutepancreatitis_only ==1)))
##1y
summary(AP_sfit, time =1)
100 - (summary(AP_sfit, time =1)$surv)*100
100 - (summary(AP_sfit, time =1)$lower)*100
100 - (summary(AP_sfit, time =1)$upper)*100
##3y
summary(AP_sfit, time =3)
100 - (summary(AP_sfit, time =3)$surv)*100
100 - (summary(AP_sfit, time =3)$lower)*100
100 - (summary(AP_sfit, time =3)$upper)*100

#Chronic pancreatitis
CP_sfit <- survfit(Surv(survival_time_years,outcome) ~ 1, data = (cohort %>% filter(t3c_chronicpancreatitis ==1)))
##1y
summary(CP_sfit, time =1)
100 - (summary(CP_sfit, time =1)$surv)*100
100 - (summary(CP_sfit, time =1)$lower)*100
100 - (summary(CP_sfit, time =1)$upper)*100
##3y
summary(CP_sfit, time =3)
100 - (summary(CP_sfit, time =3)$surv)*100
100 - (summary(CP_sfit, time =3)$lower)*100
100 - (summary(CP_sfit, time =3)$upper)*100

#Haemochromatosis
HC_sfit <- survfit(Surv(survival_time_years,outcome) ~ 1, data = (cohort %>% filter(t3c_haemochromatosis ==1)))
##1y
summary(HC_sfit, time =1)
100 - (summary(HC_sfit, time =1)$surv)*100
100 - (summary(HC_sfit, time =1)$lower)*100
100 - (summary(HC_sfit, time =1)$upper)*100
##3y
summary(HC_sfit, time =3)
100 - (summary(HC_sfit, time =3)$surv)*100
100 - (summary(HC_sfit, time =3)$lower)*100
100 - (summary(HC_sfit, time =3)$upper)*100

#Pancreatic cancer
PC_sfit <- survfit(Surv(survival_time_years,outcome) ~ 1, data = (cohort %>% filter(t3c_pancreaticcancer ==1)))
##1y
summary(PC_sfit, time =1)
100 - (summary(PC_sfit, time =1)$surv)*100
100 - (summary(PC_sfit, time =1)$lower)*100
100 - (summary(PC_sfit, time =1)$upper)*100
##3y
summary(PC_sfit, time =3)
100 - (summary(PC_sfit, time =3)$surv)*100
100 - (summary(PC_sfit, time =3)$lower)*100
100 - (summary(PC_sfit, time =3)$upper)*100

###KM plot
fit2 <- list(T2 = t2_sfit, AP = AP_sfit, CP = CP_sfit, HC = HC_sfit, PC = PC_sfit)
KM_plot_subtype <- ggsurvplot(fit = fit2, 
                          fun = function(x) {100 - x*100},
                          palette = c("#7570B3", "#66A61E", "#FF7F00", "#1F78B4", "#E7298A"),
                          combine = TRUE,
                          xlab = "Years from diagnosis",
                          ylab = "Cumulative incidence of oral therapy initiation (%)",
                          ylim = c(0, 70),
                          xlim = c(0,3),
                          break.time.by = 1,
                          font.ticks = 16,
                          font.x = 16,
                          font.y = 16,
                          font.main = 16,
                          font.legend = 16,
                          risk.table.font = 6,
                          legend.title = "",
                          legend = "bottom",
                          legend.labs = c("Type 2", "Acute pancreatitis", "Chronic pancreatitis", "Haemochromatosis", "Pancreatic cancer"),
                          censor = F,
                          conf.int = T,
                          risk.table = TRUE,
                          fontsize = 6,
                          tables.y.text = FALSE,
                          size= 1.5,
                          ggtheme = theme_classic(),
                          linetype = c("strata"),
                          tables.theme = theme_cleantable(),
                          axes.offset = FALSE)

KM_plot_subtype

#Save
pdf.options(reset = TRUE, onefile = FALSE)
pdf("SFig6_Time_to_oral_therapy_3_year_by_subtype.pdf",width=12.5,height=9)
KM_plot_subtype
dev.off()

###Cox model
cohort$t3c_subgroup <- factor(cohort$t3c_subgroup)
cohort$t3c_subgroup <- relevel(cohort$t3c_subgroup, ref="Type 2")
model_subtype <- coxph(Surv(survival_time, outcome) ~ t3c_subgroup, data = cohort)
summary(model_subtype)

################################################################################
###END
rm(list=ls())
