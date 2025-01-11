# CPRD-Rhian-T3c-Scripts
This repository contains the R scripts to define type 3c diabetes phenotypes and run the analysis to evaluate treatment outcomes on major oral glucose-lowering therapy classes in type 3c vs type 2 diabetes controls, using CPRD Aurum primary care data.

&nbsp;
Published Diabetes, Obesity and Metabolism 2025: https://doi.org/10.1111/dom.16163

&nbsp;

## Scripts

**Defining cohorts:**

* defining_T3c_cohort.R: Defines type 3c diabetes phenotypes based on a record of a pancreatic condition (acute pancreatitis, chronic pancreatitis, pancreatic cancer, or haemochromatosis) prior to a diagnosis of diabetes, and defines pancreatic exocrine insuffiency (PEI). 
Note: cohorts defined in this study are based on the at-diagnosis and treatment response cohorts defined in https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts

&nbsp;

**Characteristics at diabetes diagnosis/ treatment from diagnosis:**

* Baseline_table_at_diagnosis_cohort.R: generates a table of baseline clinical and sociodemographic characteritics of an incident cohort of individuals with type 3c diabetes (with PEI and without PEI) and type 2 diabetes at the diagnosis date of their diabetes. Outputs Table 1.
&nbsp;

* Time_to_insulin.R: survival analysis of time to initiation of insulin treatment from diabetes diagnosis, stratified by PEI status and underlying pancreatic condition, and outputting Kaplan Meir survival curves and unadjusted Cox proportional hazard ratios for each analysis. Outputs Supplementary Figures 3 & 4.
&nbsp;

* Time_to_oral_therapies.R: survival analysis of time to initiation of oral glucose-lowering therapy from diabetes diagnosis, stratified by PEI status and underlying pancreatic condition, and outputting Kaplan Meir survival curves and unadjusted Cox proportional hazard ratios for each analysis. Outputs Supplementary Figures 5 & 6.

&nbsp;

**Treatment response to oral therapies (metformin, sulphonylureas, thiazolidinediones [TZDs], SGLT2-inhibitors, DPP4-inhibitors):**

* Baseline_table_matched_cohort_by_drugs.R: generates a table of baseline clinical and sociodemographic characteritics of the matched treatment response cohort of individuals with type 3c diabetes and type 2 diabetes controls at the date of drug initiation, by drug class. Outputs Supplementary Table 3.
&nbsp;

* Main_treatment_response_analysis.R: matches individuals with type 3c diabetes to up to 10 type 2 diabetes controls, and runs analysis comparing treatment outcomes (HbA1c response, early treatment discontinuation, weight change) in individuals with type 3c diabetes (with PEI/ without PEI) vs type 2 diabetes controls, overall, by drug class, and by underlying pancreatic condition. Outputs Figure 1, Figure 2, Supplementary Figures 7 & 8, and Supplementary Tables 4-8.

&nbsp;

## Type 3c phenotypes

We defined type 3c diabetes based on a record of a pancreatic condition prior to or up to 30 days after a diagnosis of diabetes (any diabetes type). Pancreatic conditions were defined in CPRD primary care data and linked HES records, and included:
* Acute pancreatitis (CPRD medcodes and ICD10 codes used available in https://github.com/Exeter-Diabetes/CPRD-Codelists)
* Chronic pancreatitis (CPRD medcodes and ICD10 codes used available in https://github.com/Exeter-Diabetes/CPRD-Codelists)
* Pancreatic cancer (CPRD medcodes provided by the SPOCC Programme team https://blogs.exeter.ac.uk/spocc/, ICD10 codes used available in https://github.com/Exeter-Diabetes/CPRD-Codelists)
* Haemochromatosis (CPRD medcodes provided by the SPOCC Programme team https://blogs.exeter.ac.uk/spocc/, ICD10 codes used available in https://github.com/Exeter-Diabetes/CPRD-Codelists)

  &nbsp;

We excluded diabetes secondary to cystic fibrosis, and diabetes secondary to surgical pancreatic resection only from this study.

&nbsp;

We defined pancreatic exocrine insufficiency as any of the following in the CPRD primary care data: 
* Diagnosis code for pancreatic exocrine insufficiency (CPRD medcodes used available in https://github.com/Exeter-Diabetes/CPRD-Codelists)
* Faecal elastase-1 (FE1) test result less than 200 ug/g
* Prescription for pancreatic enzyme replacement therapy (PERT)

&nbsp;

## Glucose-lowering drug classes

The drug classes evaluated in this study were:
* Metformin
* Sulphonylureas
* Thiazolidinediones (TZDs)
* SGLT2-inhibitors
* DPP4-inhibitors

&nbsp;

## Defining outcome variables

* **HbA1c response:** Change from baseline HbA1c 12 months after drug initiation (the closest HbA1c measure to 12 months after initiation within 3-15 months) on unchanged therapy (no addition or cessation of other glucose-lowering medications, and continued prescription of the drug of interest)

* **Early treatment discontinuation:** Discontinuation of a therapy within 6 months of initiation, with the availability of at least 3 months follow-up time after discontinuation required to confirm the drug was discontinued. A gap of over 6 months in prescriptions was used to indicate a drug being stopped.

* **Weight change:** Change from baseline weight 12 months after drug initiation (the closest weight measure to 12 months after initiation within 3-15 months) on unchanged therapy (no addition or cessation of other glucose-lowering medications, and continued prescription of the drug of interest)

&nbsp;
