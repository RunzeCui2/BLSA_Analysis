library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggplot2)
library(geepack)
library(sjPlot)
library(readxl)
library(dplyr)
library(readr)



## import the data


blsa_var = read_xlsx("Data/variables_list.xlsx")
blsa_df = read_xlsx("Data/blsa_data_v1.xlsx")
blsa_cutoff = read_xlsx("Data/Biomarkers_MaxineBlech221110.xlsx")
blsa_df[blsa_df < 0] <- NA


## Abbreviate the biomarkers' name

new_names <- blsa_var$Label
abbr_new_names <- abbreviate(new_names, minlength = 10)
names(blsa_df) <- abbr_new_names
janitor::clean_names(blsa_df)

## MRmean and RMRmean

blsa_df =
  blsa_df %>%
  rename("TBM" =  "LnMs(TBdy)",
         "T1230" = "T1230ee(/)",
         "T1300" = "T1300ee(/)",
         "T1330" = "T1330ee(/)",
         "T1400" = "T1400ee(/)",
         "T1430" = "T1430ee(/)",
         "T1500" = "T1500ee(/)",
         "T1530" = "T1530ee(/)") %>%
  mutate(MRmean = rowMeans(select(., starts_with("T1")), na.rm = TRUE)) %>% 
  mutate(RMRmean = MRmean / TBM, na.rm = TRUE)
blsa_df$MRmean[is.nan(blsa_df$MRmean)] <- NA
blsa_df$RMRmean[is.nan(blsa_df$RMRmean)] <- NA

## Rename the required biomarkers

blsa_df = 
  blsa_df %>% 
  select(
    `BLSA ID`,
    AgeatVisit, 
    Gender, 
    VisitNumbr, 
    contains(c("MR", "ANTH", "bp", "TBM", "ALB", "AMY", "CREA", "VITB", "24CF", "SOD", "TSH", "TOTAL", "FREET", "SER", "CHOLE", "TRIG", "DD", "PLAS", "FIBR", "CRP", "IL", "NF", "ADIP", "LEPT", "RES", "Ghr", "BIOA", "HOMO", "ESTR", "DEHY", "INSU", "Crbdy", "bcps" )), 
    "AFGv(mg/d)", 
    "AGva120m(/", 
    "GLUCOSE(/L", 
    "HBA1C (%)", 
    "URICACID(/") %>%
  select(!contains('Hspr')) %>% 
  rename(
    'id' = `BLSA ID`,
    "AFG" = "AFGv(mg/d)",
    "AFG120" = "AGva120m(/",
    "glucose" = "GLUCOSE(/L",
    "hba1c" = "HBA1C (%)",
    "uric_acid" = "URICACID(/",
    "age" = "AgeatVisit",
    "visit" = "VisitNumbr",
    "gender" = "Gender",
    "bmi" = "ANTH:BdyMI",
    "weight" = "ANTH:W(kg)",
    "height" = "ANTH:Hghic",
    "waist_circumference" = "ANTH:Wscic",
    "systolic_bp" = "CSbp(mmHg)",
    "diastolic_bp" = "CDbp(mmHg)",
    "albumin" = "ALBUMIN(/L",
    "total_lmass" = "TBM",
    "amylase" = "AMYLASE(/L",
    "creatinine_serum" = "CREATINIS(",
    "vb12" = "VITB12(/L)",
    "ur24_cortisol" = "URINE24CF(",
    "sodium" = "SODIUM(/L)",
    "total_t3" = "TOTALT3(/L",
    "total_t4" = "TOTALT4(/L",
    "total_testosterone" = "TOTALTEST(",
    "free_t3" = "FREET3(/L)",
    "free_t4" = "FREET4(/L)",
    "serum_cortisol" = "SERUMCORT(",
    "cholesterol" = "CHOLESTER(",
    "hdl" = "HDL-CHOLE(",
    "ldl" = "LDL-CHOLE(",
    "triglycerides" = "TRIGLYCER(",
    "d_dimer" = "DDIMER(/L)",
    "pai" = "PLASMIAI1(",
    "fibrinogen" = "FIBRINOGE(",
    "crp" = "CRP(ug/mL)",
    "il6" = "IL6(pg/mL)",
    "il18" = "IL18(pg/L)",
    "il1ra" = "IL1RA(p/L)",
    "il15" = "IL15(pg/L)",
    "stnf_r1" = "STNF-RI(/L",
    "stnf_r2" = "STNF-RII(/",
    "bio_testosterone" = "BIOAVAILT(",
    "tnfa" = "TNFA(pg/L)",
    "adiponectin" = "ADIPONECT(",
    "leptin" = "LEPTIN(/L)",
    "resistin" = "RESISTIN(/",
    "ghrelin" = "Ghr(pg/ml)",
    "homocysteine" = "HOMOCYSTEu",
    "estradiol" = "ESTRADIOL(",
    "dehy_sulfate" = "DEHYDROES(",
    "insu_gf" = "INSULIGF1(",
    "corebody_temp" = "Crbdytm(F)",
    "skin_fold" = "s(bcps)(c)",
    "tsh" = "TSH")

## Functions for cutoffs
## Based on the clinical data
cutoff <- function(df, var_name, var_high = c(Inf, Inf), var_low = c(0, 0)) {
  var_name = enquo(var_name) # Makes it so I can work with the var_name more easily
  df %>% 
    mutate("{{var_name}}_cutoff" := if_else(is.na(!!var_name), 'NA', 
                                            if_else(gender == "M" & !!var_name > var_high[1] |
                                                      gender == "M" & !!var_name < var_low[1] |
                                                      gender == "F" & !!var_name > var_high[2] |
                                                      gender == "F" & !!var_name < var_low[2],
                                                    "abnormal", 
                                                    "normal")))
}
## Based on the percentage 

percent_cutoff <- function(df, var_name, perc_high = c(1.0), perc_low = c(0.0)) {
  var_name = enquo(var_name) # Makes it so I can work with the var_name more easily
  percentile_df <- df %>% 
    group_by(gender) %>% 
    summarize(quant_high = quantile(!!var_name, probs = perc_high, na.rm = TRUE),
              quant_low = quantile(!!var_name, probs = perc_low, na.rm = TRUE)) %>% 
    ungroup()
  male_cutoffs <- percentile_df %>% filter(gender == "M") %>% select(quant_high, quant_low)
  female_cutoffs <- percentile_df %>% filter(gender == "F") %>% select(quant_high, quant_low)
  
  df %>% 
    mutate("{{var_name}}_percent_cutoff" := if_else(is.na(!!var_name), 'NA', 
                                                    if_else(gender == "M" & !!var_name >    male_cutoffs$quant_high |
                                                              gender == "M" & !!var_name <  male_cutoffs$quant_low |
                                                              gender == "F" & !!var_name >  female_cutoffs$quant_high |
                                                              gender == "F" & !!var_name <  female_cutoffs$quant_low,
                                                            "abnormal", 
                                                            "normal")))
}

## customized function


adiponectin_cutoff <- function(df, var_name) {
  var_name = enquo(var_name) # Makes it so I can work with the var_name more easily
  df %>% 
    mutate("{{var_name}}_cutoff" := if_else(is.na(!!var_name), 'NA', 
                                            if_else(    gender == "M" & bmi < 25 & !!var_name > 26 |
                                                          gender == "M" & bmi < 25 & !!var_name < 4 |
                                                          gender == "F" & bmi < 25 & !!var_name > 37 |
                                                          gender == "F" & bmi < 25 & !!var_name < 5 | 
                                                          gender == "M" & bmi == (25:30) & !!var_name > 20 |
                                                          gender == "M" & bmi == (25:30) & !!var_name < 4 |
                                                          gender == "F" & bmi == (25:30) & !!var_name > 28 |
                                                          gender == "F" & bmi == (25:30) & !!var_name < 5 |
                                                          gender == "M" & bmi > 30 & !!var_name > 20 |
                                                          gender == "M" & bmi > 30 & !!var_name < 2 |
                                                          gender == "F" & bmi > 30 & !!var_name > 22 |
                                                          gender == "F" & bmi > 30 & !!var_name < 4, 
                                                        "abnormal", 
                                                        "normal")))
}

insu_cutoff <- function(df, var_name) {
  var_name = enquo(var_name) # Makes it so I can work with the var_name more easily
  df %>% 
    mutate("{{var_name}}_cutoff" := if_else(is.na(!!var_name), 'NA', 
                                            if_else(      gender == "M" & age == (18:22) & !!var_name > 442 |
                                                            gender == "M" & age == (18:22) & !!var_name < 91 |
                                                            gender == "F" & age == (18:22) & !!var_name > 370 |
                                                            gender == "F" & age == (18:22) & !!var_name < 85 | 
                                                            gender == "M" & age == (23:25) & !!var_name > 346 |
                                                            gender == "M" & age == (23:25) & !!var_name < 66 |
                                                            gender == "F" & age == (23:25) & !!var_name > 320 |
                                                            gender == "F" & age == (23:25) & !!var_name < 73 |
                                                            gender == "M" & age == (26:30) & !!var_name > 329 | 
                                                            gender == "M" & age == (26:30) & !!var_name < 60 |
                                                            gender == "F" & age == (26:30) & !!var_name > 303 |
                                                            gender == "F" & age == (26:30) & !!var_name < 66 |
                                                            gender == "M" & age == (31:35) & !!var_name > 310 | 
                                                            gender == "M" & age == (31:35) & !!var_name < 54 |
                                                            gender == "F" & age == (31:35) & !!var_name > 279 |
                                                            gender == "F" & age == (31:35) & !!var_name < 59 |
                                                            gender == "M" & age == (36:40) & !!var_name > 292 | 
                                                            gender == "M" & age == (36:40) & !!var_name < 48 |
                                                            gender == "F" & age == (36:40) & !!var_name > 258 |
                                                            gender == "F" & age == (36:40) & !!var_name < 54 |
                                                            gender == "M" & age == (41:45) & !!var_name > 275 | 
                                                            gender == "M" & age == (41:45) & !!var_name < 44 |
                                                            gender == "F" & age == (41:45) & !!var_name > 240 |
                                                            gender == "F" & age == (41:45) & !!var_name < 49 |
                                                            gender == "M" & age == (46:50) & !!var_name > 259 | 
                                                            gender == "M" & age == (46:50) & !!var_name < 40 |
                                                            gender == "F" & age == (46:50) & !!var_name > 227 |
                                                            gender == "F" & age == (46:50) & !!var_name < 44 |
                                                            gender == "M" & age == (51:55) & !!var_name > 245 | 
                                                            gender == "M" & age == (51:55) & !!var_name < 37 |
                                                            gender == "F" & age == (51:55) & !!var_name > 217 |
                                                            gender == "F" & age == (51:55) & !!var_name < 40 |
                                                            gender == "M" & age == (56:60) & !!var_name > 232 | 
                                                            gender == "M" & age == (56:60) & !!var_name < 34 |
                                                            gender == "F" & age == (56:60) & !!var_name > 208 |
                                                            gender == "F" & age == (56:60) & !!var_name < 37 |
                                                            gender == "M" & age == (61:65) & !!var_name > 220 | 
                                                            gender == "M" & age == (61:65) & !!var_name < 33 |
                                                            gender == "F" & age == (61:65) & !!var_name > 201 |
                                                            gender == "F" & age == (61:65) & !!var_name < 35 |
                                                            gender == "M" & age == (66:70) & !!var_name > 209 | 
                                                            gender == "M" & age == (66:70) & !!var_name < 32 |
                                                            gender == "F" & age == (66:70) & !!var_name > 194 |
                                                            gender == "F" & age == (66:70) & !!var_name < 34 |
                                                            gender == "M" & age == (71:75) & !!var_name > 200 | 
                                                            gender == "M" & age == (71:75) & !!var_name < 32 |
                                                            gender == "F" & age == (71:75) & !!var_name > 187 |
                                                            gender == "F" & age == (71:75) & !!var_name < 34 |
                                                            gender == "M" & age == (76:80) & !!var_name > 192 | 
                                                            gender == "M" & age == (76:80) & !!var_name < 33 |
                                                            gender == "F" & age == (76:80) & !!var_name > 182 |
                                                            gender == "F" & age == (76:80) & !!var_name < 34 |
                                                            gender == "M" & age == (81:85) & !!var_name > 185 | 
                                                            gender == "M" & age == (81:85) & !!var_name < 33 |
                                                            gender == "F" & age == (81:85) & !!var_name > 177 |
                                                            gender == "F" & age == (81:85) & !!var_name < 34 |
                                                            gender == "M" & age == (86:90) & !!var_name > 179 | 
                                                            gender == "M" & age == (86:90) & !!var_name < 33 |
                                                            gender == "F" & age == (86:90) & !!var_name > 175 |
                                                            gender == "F" & age == (86:90) & !!var_name < 33 |
                                                            gender == "M" & age >= 91 & !!var_name > 173 | 
                                                            gender == "M" & age >= 91 & !!var_name < 32 |
                                                            gender == "F" & age >= 91 & !!var_name > 179 |
                                                            gender == "F" & age >= 91 & !!var_name < 25,
                                                          "abnormal", 
                                                          "normal")))
}




leptin_cutoff <- function(df, var_name) {
  var_name = enquo(var_name) # Makes it so I can work with the var_name more easily
  df %>% 
    mutate("{{var_name}}_cutoff" := if_else(is.na(!!var_name), 'NA', 
                                            if_else(      gender == "M" & bmi == 11 & !!var_name > 0.4 |
                                                            gender == "M" & bmi == 11 & !!var_name < 0.1|
                                                            gender == "F" & bmi == 11 & !!var_name > 3.6 |
                                                            gender == "F" & bmi == 11 & !!var_name < 0.7 | 
                                                            gender == "M" & bmi == 12 & !!var_name > 0.6 |
                                                            gender == "M" & bmi == 12 & !!var_name < 0.1 |
                                                            gender == "F" & bmi == 12 & !!var_name > 4.2 |
                                                            gender == "F" & bmi == 12 & !!var_name < 0.8 |
                                                            gender == "M" & bmi == 13 & !!var_name > 0.7 | 
                                                            gender == "M" & bmi == 13 & !!var_name < 0.1 |
                                                            gender == "F" & bmi == 13 & !!var_name > 4.8 |
                                                            gender == "F" & bmi == 13 & !!var_name < 0.9 |
                                                            gender == "M" & bmi == 14 & !!var_name > 0.9 | 
                                                            gender == "M" & bmi == 14 & !!var_name < 0.1 |
                                                            gender == "F" & bmi == 14 & !!var_name > 5.6 |
                                                            gender == "F" & bmi == 14 & !!var_name < 1.0 |
                                                            gender == "M" & bmi == 15 & !!var_name > 1.1 | 
                                                            gender == "M" & bmi == 15 & !!var_name < 0.1 |
                                                            gender == "F" & bmi == 15 & !!var_name > 6.5 |
                                                            gender == "F" & bmi == 15 & !!var_name < 1.2 |
                                                            gender == "M" & bmi == 16 & !!var_name > 1.3 | 
                                                            gender == "M" & bmi == 16 & !!var_name < 0.2 |
                                                            gender == "F" & bmi == 16 & !!var_name > 7.5 |
                                                            gender == "F" & bmi == 16 & !!var_name < 1.4 |
                                                            gender == "M" & bmi == 17 & !!var_name > 1.7 | 
                                                            gender == "M" & bmi == 17 & !!var_name < 0.2 |
                                                            gender == "F" & bmi == 17 & !!var_name > 8.7 |
                                                            gender == "F" & bmi == 17 & !!var_name < 0.6 |
                                                            gender == "M" & bmi == 18 & !!var_name > 2.1 | 
                                                            gender == "M" & bmi == 18 & !!var_name < 0.2 |
                                                            gender == "F" & bmi == 18 & !!var_name > 10 |
                                                            gender == "F" & bmi == 18 & !!var_name < 1.8 |
                                                            gender == "M" & bmi == 19 & !!var_name > 2.6 | 
                                                            gender == "M" & bmi == 19 & !!var_name < 0.3 |
                                                            gender == "F" & bmi == 19 & !!var_name > 11.6 |
                                                            gender == "F" & bmi == 19 & !!var_name < 2.1 |
                                                            gender == "M" & bmi == 20 & !!var_name > 3.2 | 
                                                            gender == "M" & bmi == 20 & !!var_name < 0.4 |
                                                            gender == "F" & bmi == 20 & !!var_name > 13.4 |
                                                            gender == "F" & bmi == 20 & !!var_name < 2.4 |
                                                            gender == "M" & bmi == 21 & !!var_name > 4 | 
                                                            gender == "M" & bmi == 21 & !!var_name < 0.4 |
                                                            gender == "F" & bmi == 21 & !!var_name > 15.6 |
                                                            gender == "F" & bmi == 21 & !!var_name < 2.8 |
                                                            gender == "M" & bmi == 22 & !!var_name > 5.0 | 
                                                            gender == "M" & bmi == 22 & !!var_name < 0.5 |
                                                            gender == "F" & bmi == 22 & !!var_name > 18 |
                                                            gender == "F" & bmi == 22 & !!var_name < 3.3 |
                                                            gender == "M" & bmi == 23 & !!var_name > 6.2 | 
                                                            gender == "M" & bmi == 23 & !!var_name < 0.8 |
                                                            gender == "F" & bmi == 23 & !!var_name > 20.9 |
                                                            gender == "F" & bmi == 23 & !!var_name < 3.8 |
                                                            gender == "M" & bmi == 24 & !!var_name > 7.7 | 
                                                            gender == "M" & bmi == 24 & !!var_name < 0.9 |
                                                            gender == "F" & bmi == 24 & !!var_name > 24.2 |
                                                            gender == "F" & bmi == 24 & !!var_name < 4.4 |
                                                            gender == "M" & bmi == 25 & !!var_name > 9.6 | 
                                                            gender == "M" & bmi == 25 & !!var_name < 1.1 |
                                                            gender == "F" & bmi == 25 & !!var_name > 28 |
                                                            gender == "F" & bmi == 25 & !!var_name < 5.1 |
                                                            gender == "M" & bmi == 26 & !!var_name > 12 | 
                                                            gender == "M" & bmi == 26 & !!var_name < 1.3 |
                                                            gender == "F" & bmi == 26 & !!var_name > 32.4 |
                                                            gender == "F" & bmi == 26 & !!var_name < 5.9 |
                                                            gender == "M" & bmi == 27 & !!var_name > 14.9 | 
                                                            gender == "M" & bmi == 27 & !!var_name < 1.6 |
                                                            gender == "F" & bmi == 27 & !!var_name > 37.5 |
                                                            gender == "F" & bmi == 27 & !!var_name < 6.8 |
                                                            gender == "M" & bmi == 28 & !!var_name > 18.6 | 
                                                            gender == "M" & bmi == 28 & !!var_name < 2 |
                                                            gender == "F" & bmi == 28 & !!var_name > 43.5 |
                                                            gender == "F" & bmi == 28 & !!var_name < 7.9 |
                                                            gender == "M" & bmi == 29 & !!var_name > 23.2 | 
                                                            gender == "M" & bmi == 29 & !!var_name < 2.5 |
                                                            gender == "F" & bmi == 29 & !!var_name > 50.4 |
                                                            gender == "F" & bmi == 29 & !!var_name < 9.1 |
                                                            gender == "M" & bmi == 30 & !!var_name > 28.9 | 
                                                            gender == "M" & bmi == 30 & !!var_name < 3.2 |
                                                            gender == "F" & bmi == 30 & !!var_name > 58.3 |
                                                            gender == "F" & bmi == 30 & !!var_name < 10.6 |
                                                            gender == "M" & bmi == 31 & !!var_name > 36 | 
                                                            gender == "M" & bmi == 31 & !!var_name < 3.9 |
                                                            gender == "F" & bmi == 31 & !!var_name > 67.5 |
                                                            gender == "F" & bmi == 31 & !!var_name < 12.2 |
                                                            gender == "M" & bmi == 32 & !!var_name > 44.9 | 
                                                            gender == "M" & bmi == 32 & !!var_name < 4.9 |
                                                            gender == "F" & bmi == 32 & !!var_name > 78.2 |
                                                            gender == "F" & bmi == 32 & !!var_name < 14.1 |
                                                            gender == "M" & bmi == 33 & !!var_name > 55.8 | 
                                                            gender == "M" & bmi == 33 & !!var_name < 6.1 |
                                                            gender == "F" & bmi == 33 & !!var_name > 90.5 |
                                                            gender == "F" & bmi == 33 & !!var_name < 16.4 |
                                                            gender == "M" & bmi == 34 & !!var_name > 69.6 | 
                                                            gender == "M" & bmi == 34 & !!var_name < 7.6 |
                                                            gender == "F" & bmi == 34 & !!var_name > 105 |
                                                            gender == "F" & bmi == 34 & !!var_name < 19 |
                                                            gender == "M" & bmi == 35 & !!var_name > 86.7 | 
                                                            gender == "M" & bmi == 35 & !!var_name < 9.5 |
                                                            gender == "F" & bmi == 35 & !!var_name > 121 |
                                                            gender == "F" & bmi == 35 & !!var_name < 22 |
                                                            gender == "M" & bmi == 36 & !!var_name > 108 | 
                                                            gender == "M" & bmi == 36 & !!var_name < 11.8 |
                                                            gender == "F" & bmi == 36 & !!var_name > 141 |
                                                            gender == "F" & bmi == 36 & !!var_name < 25.4 |
                                                            gender == "M" & bmi == 37 & !!var_name > 135 | 
                                                            gender == "M" & bmi == 37 & !!var_name < 14.8 |
                                                            gender == "F" & bmi == 37 & !!var_name > Inf |
                                                            gender == "F" & bmi == 37 & !!var_name < 0 ,
                                                          "abnormal", 
                                                          "normal")))
}





dehy_cutoff <- function(df, var_name) {
  var_name = enquo(var_name) # Makes it so I can work with the var_name more easily
  df %>% 
    mutate("{{var_name}}_cutoff" := if_else(is.na(!!var_name), 'NA', 
                                            if_else(    gender == "M" & age == (18:30) & !!var_name > 728 |
                                                          gender == "M" & age == (18:30) & !!var_name < 105 |
                                                          gender == "F" & age == (18:30) & !!var_name > 377 |
                                                          gender == "F" & age == (18:30) & !!var_name < 83 | 
                                                          gender == "M" & age == (31:40) & !!var_name > 522 |
                                                          gender == "M" & age == (31:40) & !!var_name < 57 |
                                                          gender == "F" & age == (31:40) & !!var_name > 295 |
                                                          gender == "F" & age == (31:40) & !!var_name < 45 |
                                                          gender == "M" & age == (41:50) & !!var_name > 395 | 
                                                          gender == "M" & age == (41:50) & !!var_name < 34 |
                                                          gender == "F" & age == (41:50) & !!var_name > 240 |
                                                          gender == "F" & age == (41:50) & !!var_name < 27 |
                                                          gender == "M" & age == (51:60) & !!var_name > 299 | 
                                                          gender == "M" & age == (51:60) & !!var_name < 20 |
                                                          gender == "F" & age == (51:60) & !!var_name > 195 |
                                                          gender == "F" & age == (51:60) & !!var_name < 16 |
                                                          gender == "M" & age == (61:70) & !!var_name > 227 | 
                                                          gender == "M" & age == (61:70) & !!var_name < 12 |
                                                          gender == "F" & age == (61:70) & !!var_name > 159 |
                                                          gender == "F" & age == (61:70) & !!var_name < 9.7 |
                                                          gender == "M" & age >= 71 & !!var_name > 162 | 
                                                          gender == "M" & age >= 71 & !!var_name < 6.6 |
                                                          gender == "F" & age >= 71 & !!var_name > 124 |
                                                          gender == "F" & age >= 71 & !!var_name < 5.3,
                                                        "abnormal", 
                                                        "normal")))
}

## function for N_obs and N_indiv for all biomarkers

N <- function(var_name)  {
  df1 = 
    blsa_clean %>% 
    select(c(id, age, gender, visit, {{var_name}})) %>% 
    drop_na() %>% 
    nrow()
  df2 = 
    blsa_clean %>% 
    select(c(id, age, gender, visit, {{var_name}})) %>%
    mutate(id = factor(id)) %>%
    drop_na() %>% 
    group_by(id) %>%
    summarize(n = n()) %>% 
    count()
  
  tibble(df1,df2) %>% 
    rename("N_obs" = df1,
           "N_indi" = n)
}


## function for the number of abnormal

N_abn = function(df, var_name) {
  df %>% 
    select({{var_name}}) %>% 
    filter({{var_name}} == "abnormal") %>% 
    nrow()
}




## Function for percentage converter
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}

## Function for waves:

wave <- function(var_name) {
  blsa_clean %>% 
    filter(!is.na({{var_name}})) %>% 
    group_by(visit) %>% 
    count() %>% 
    nrow()
}

## log transform data frame

blsa_df$MRmean = log(blsa_df$MRmean + 1)
blsa_df$RMRmean = log(blsa_df$RMRmean + 1)
blsa_df$bmi = log(blsa_df$bmi + 1)
blsa_df$weight = log(blsa_df$weight + 1)
blsa_df$height = log(blsa_df$height + 1)
blsa_df$waist_circumference = log(blsa_df$waist_circumference + 1)
blsa_df$systolic_bp = log(blsa_df$systolic_bp + 1)
blsa_df$diastolic_bp = log(blsa_df$diastolic_bp + 1)
blsa_df$total_lmass = log(blsa_df$total_lmass + 1)
blsa_df$albumin = log(blsa_df$albumin + 1)
blsa_df$amylase = log(blsa_df$amylase + 1)
blsa_df$creatinine_serum = log(blsa_df$creatinine_serum + 1)
blsa_df$vb12 = log(blsa_df$vb12 + 1)
blsa_df$ur24_cortisol = log(blsa_df$ur24_cortisol + 1)
blsa_df$sodium = log(blsa_df$sodium + 1)
blsa_df$tsh = log(blsa_df$tsh + 1)
blsa_df$total_t3 = log(blsa_df$total_t3 + 1)
blsa_df$total_t4 = log(blsa_df$total_t4 + 1)
blsa_df$total_testosterone = log(blsa_df$total_testosterone + 1)
blsa_df$free_t3 = log(blsa_df$free_t3 + 1)
blsa_df$free_t4 = log(blsa_df$free_t4 + 1)
blsa_df$serum_cortisol = log(blsa_df$serum_cortisol + 1)
blsa_df$cholesterol = log(blsa_df$cholesterol + 1)
blsa_df$hdl = log(blsa_df$hdl + 1)
blsa_df$ldl = log(blsa_df$ldl + 1)
blsa_df$triglycerides = log(blsa_df$triglycerides + 1)
blsa_df$d_dimer = log(blsa_df$d_dimer + 1)
blsa_df$pai = log(blsa_df$pai + 1)
blsa_df$fibrinogen = log(blsa_df$fibrinogen + 1)
blsa_df$crp = log(blsa_df$crp + 1)
blsa_df$il6 = log(blsa_df$il6 + 1)
blsa_df$il18 = log(blsa_df$il18 + 1)
blsa_df$il1ra = log(blsa_df$il1ra + 1)
blsa_df$il15 = log(blsa_df$il15 + 1)
blsa_df$bio_testosterone = log(blsa_df$bio_testosterone + 1)
blsa_df$tnfa = log(blsa_df$tnfa + 1)
blsa_df$stnf_r1 = log(blsa_df$stnf_r1 + 1)
blsa_df$stnf_r2 = log(blsa_df$stnf_r2 + 1)
blsa_df$adiponectin = log(blsa_df$adiponectin + 1)
blsa_df$leptin = log(blsa_df$leptin + 1)
blsa_df$resistin = log(blsa_df$resistin + 1)
blsa_df$ghrelin = log(blsa_df$ghrelin + 1)
blsa_df$homocysteine = log(blsa_df$homocysteine + 1)
blsa_df$estradiol = log(blsa_df$estradiol + 1)
blsa_df$dehy_sulfate = log(blsa_df$dehy_sulfate + 1)
blsa_df$insu_gf = log(blsa_df$insu_gf + 1)
blsa_df$corebody_temp = log(blsa_df$corebody_temp + 1)
blsa_df$skin_fold = log(blsa_df$skin_fold + 1)
blsa_df$AFG = log(blsa_df$AFG + 1)
blsa_df$hba1c = log(blsa_df$hba1c + 1)
blsa_df$uric_acid = log(blsa_df$uric_acid + 1)

blsa_df

## Get rid of outliers

st_5_up <- function(x) {
  mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE)*5
}

st_5_low <- function(x){
  mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE)*5
}

st_5_upper <- lapply(blsa_df[, -c(1, 2, 3, 4)], st_5_up)
st_5_lower <- lapply(blsa_df[, -c(1, 2, 3, 4)], st_5_low)
outlier_name <- colnames(blsa_df[, (-c(1:4))])

for (i in outlier_name) {
  blsa_df[i][blsa_df[i] >= st_5_upper[i]] <- NA
  blsa_df[i][blsa_df[i] <= st_5_lower[i]] <- NA
}

## exp back to raw data

blsa_df$MRmean = exp(blsa_df$MRmean) - 1
blsa_df$RMRmean = exp(blsa_df$RMRmean) - 1
blsa_df$bmi = exp(blsa_df$bmi) - 1
blsa_df$weight = exp(blsa_df$weight) - 1
blsa_df$height = exp(blsa_df$height) - 1
blsa_df$waist_circumference = exp(blsa_df$waist_circumference) - 1
blsa_df$systolic_bp = exp(blsa_df$systolic_bp) - 1
blsa_df$diastolic_bp = exp(blsa_df$diastolic_bp) - 1
blsa_df$total_lmass = exp(blsa_df$total_lmass) - 1
blsa_df$albumin = exp(blsa_df$albumin) - 1
blsa_df$amylase = exp(blsa_df$amylase) - 1
blsa_df$creatinine_serum = exp(blsa_df$creatinine_serum) - 1
blsa_df$vb12 = exp(blsa_df$vb12) - 1
blsa_df$ur24_cortisol = exp(blsa_df$ur24_cortisol) - 1
blsa_df$sodium = exp(blsa_df$sodium) - 1
blsa_df$tsh = exp(blsa_df$tsh) - 1
blsa_df$total_t3 = exp(blsa_df$total_t3) - 1
blsa_df$total_t4 = exp(blsa_df$total_t4) - 1
blsa_df$total_testosterone = exp(blsa_df$total_testosterone) - 1
blsa_df$free_t3 = exp(blsa_df$free_t3) - 1
blsa_df$free_t4 = exp(blsa_df$free_t4) - 1
blsa_df$serum_cortisol = exp(blsa_df$serum_cortisol) - 1
blsa_df$cholesterol = exp(blsa_df$cholesterol) - 1
blsa_df$hdl = exp(blsa_df$hdl) - 1
blsa_df$ldl = exp(blsa_df$ldl) - 1
blsa_df$triglycerides = exp(blsa_df$triglycerides) - 1
blsa_df$d_dimer = exp(blsa_df$d_dimer) - 1
blsa_df$pai = exp(blsa_df$pai) - 1
blsa_df$fibrinogen = exp(blsa_df$fibrinogen) - 1
blsa_df$crp = exp(blsa_df$crp) - 1
blsa_df$il6 = exp(blsa_df$il6) - 1
blsa_df$il18 = exp(blsa_df$il18) - 1
blsa_df$il1ra = exp(blsa_df$il1ra) - 1
blsa_df$il15 = exp(blsa_df$il15) - 1
blsa_df$bio_testosterone = exp(blsa_df$bio_testosterone) - 1
blsa_df$tnfa = exp(blsa_df$tnfa) - 1
blsa_df$stnf_r1 = exp(blsa_df$stnf_r1) - 1
blsa_df$stnf_r2 = exp(blsa_df$stnf_r2) - 1
blsa_df$adiponectin = exp(blsa_df$adiponectin) - 1
blsa_df$leptin = exp(blsa_df$leptin) - 1
blsa_df$resistin = exp(blsa_df$resistin) - 1
blsa_df$ghrelin = exp(blsa_df$ghrelin) - 1
blsa_df$homocysteine = exp(blsa_df$homocysteine) - 1
blsa_df$estradiol = exp(blsa_df$estradiol) - 1
blsa_df$dehy_sulfate = exp(blsa_df$dehy_sulfate) - 1
blsa_df$insu_gf = exp(blsa_df$insu_gf) - 1
blsa_df$corebody_temp = exp(blsa_df$corebody_temp) - 1
blsa_df$skin_fold = exp(blsa_df$skin_fold) - 1
blsa_df$AFG = exp(blsa_df$AFG) - 1
blsa_df$hba1c = exp(blsa_df$hba1c) - 1
blsa_df$uric_acid = exp(blsa_df$uric_acid) - 1

blsa_clean = blsa_df %>% 
  select(c(id, age, gender, visit, total_lmass, MRmean, RMRmean, everything(), -c(AFG120, glucose, weight, height)))



# Biomarkers 

## Adiponectin (ug/mL):
blsa_clean = adiponectin_cutoff(blsa_clean, `adiponectin`)

## Adjusted Fasting Glucose value (mg/dl)
blsa_clean = cutoff(blsa_clean, `AFG`, var_high = c(100, 100), var_low = c(70, 70))


## ALBUMIN (g/dL)
blsa_clean = cutoff(blsa_clean, `albumin`, var_high = c(5.4, 5.4), var_low = c(3.4, 3.4))

## AMYLASE (u/L)
blsa_clean = cutoff(blsa_clean, `amylase`, var_high = c(100, 100), var_low = c(28, 28))

## BMI
blsa_clean = cutoff(blsa_clean, `bmi`, var_high = c(30, 30), var_low = c(18.5, 18.5))

## Waist circumference in cm
blsa_clean = cutoff(blsa_clean, `waist_circumference`, var_high = c(94, 80), var_low = c(0, 0))




## BIOAVAILABLE TESTOSTERONE (ng/dL)
blsa_clean = cutoff(blsa_clean, `bio_testosterone`, var_high = c(400, 19), var_low = c(110, 1))

## CHOLESTEROL (mg/dL)
blsa_clean = cutoff(blsa_clean, `cholesterol`, var_high = c(200, 200), var_low = c(0, 0))


## Composed systolic blood pressure (mmHg)
blsa_clean = cutoff(blsa_clean, `systolic_bp`, var_high = c(120, 120), var_low = c(90, 90))

## Composed Diastolic  blood pressure (mmHg)
blsa_clean = cutoff(blsa_clean, `diastolic_bp`, var_high = c(80, 80), var_low = c(60, 60))

## Core body temperature (F)
blsa_clean = cutoff(blsa_clean, `corebody_temp`, var_high = c(99, 99), var_low = c(97, 97))

## CREATININE, SERUM (mg/dL)
blsa_clean = cutoff(blsa_clean, `creatinine_serum`, var_high = c(1.35, 1.04), var_low = c(0.74, 0.59))

## CRP (ug/mL)
blsa_clean = cutoff(blsa_clean, `crp`, var_high = c(10, 10), var_low = c(0, 0))

## D DIMER (mg/mL) drop the lower bound

blsa_clean = cutoff(blsa_clean, `d_dimer`, var_high = c(500, 500), var_low = c(0, 0))

## DEHYDROEPIANDROSTERONE SULFATE (ug/dL) different by age

blsa_clean = dehy_cutoff(blsa_clean, `dehy_sulfate`)

## ESTRADIOL (ng/DL)
blsa_clean = cutoff(blsa_clean, `estradiol`, var_high = c(4, 1), var_low = c(1, 0))

## FIBRINOGEN (mg/dL)
blsa_clean = cutoff(blsa_clean, `fibrinogen`, var_high = c(393, 393), var_low = c(200, 200))

## FREE T3 (pg/mL)
blsa_clean = cutoff(blsa_clean, `free_t3`, var_high = c(4.5, 4.5), var_low = c(1.3, 1.3))

## FREE T4 (ng/dL)
blsa_clean = cutoff(blsa_clean, `free_t4`, var_high = c(2.3, 2.3), var_low = c(0.9, 0.9))

## Ghrelin (pg/ml)
blsa_clean = percent_cutoff(blsa_clean, ghrelin, perc_high = c(0.8), perc_low = c(0.0))


## HBA1C (%)
blsa_clean = cutoff(blsa_clean, `hba1c`, var_high = c(5.7, 5.7), var_low = c(4.1, 4.1))


## HDL-CHOLESTEROL (mg/dL)
blsa_clean = cutoff(blsa_clean, `hdl`, var_high = c(Inf, Inf), var_low = c(60, 60))


## HOMOCYSTEINE (umol/L) we choose age = 65
blsa_clean = cutoff(blsa_clean, `homocysteine`, var_high = c(16.3, 15.6), var_low = c(7.1, 5.6))


## IL15 (pg/ml)
blsa_clean = cutoff(blsa_clean, `il15`, var_high = c(20, 20), var_low = c(0, 0))

## IL18
blsa_clean = cutoff(blsa_clean, `il18`, var_high = c(468, 468), var_low = c(0, 0))

## IL6
blsa_clean = percent_cutoff(blsa_clean, `il6`, perc_high = c(0.8), perc_low = c(0.0))


## IL1RA (pg/mL) Google value
blsa_clean = percent_cutoff(blsa_clean, `il1ra`, perc_high = c(1.0), perc_low = c(0.2))


## INSULIN-LIKE GROWTH FACTOR 1 (ng/mL) Different based on the age

blsa_clean = insu_cutoff(blsa_clean, `insu_gf`)

## LDL-CHOLESTEROL (mg/dL)
blsa_clean = cutoff(blsa_clean, `ldl`, var_high = c(129, 129), var_low = c(100, 100))

## lean mass no cutoff

## LEPTIN (ng/mL) based on the bmi

blsa_clean = leptin_cutoff(blsa_clean, `leptin`)

## PLASMINOGEN ACTIVATOR INHIBITOR 1 (ng/mL)
blsa_clean = percent_cutoff(blsa_clean, `pai`, perc_high = c(0.8), perc_low = c(0.0))


## RESISTIN (ng/mL)
blsa_clean = percent_cutoff(blsa_clean, `resistin`, perc_high = c(0.8), perc_low = c(0.0))


## SERUM CORTISOL (ug/dL) For the morning
blsa_clean = cutoff(blsa_clean, `serum_cortisol`, var_high = c(25, 25), var_low = c(7, 7))

## skinfold (biceps) (cm)
blsa_clean = cutoff(blsa_clean, `skin_fold`, var_high = c(127, 177.8), var_low = c(63.5, 76.2))

## SODIUM (mmol/L)
blsa_clean = cutoff(blsa_clean, `sodium`, var_high = c(145, 145), var_low = c(135, 135))

## STNF-RI (pg/mL) 
blsa_clean = percent_cutoff(blsa_clean, stnf_r1, perc_high = c(0.8), perc_low = c(0.0))

## STNF-RII(pg/mL)
blsa_clean = percent_cutoff(blsa_clean, stnf_r2, perc_high = c(0.8), perc_low = c(0.0))

## TNFA (pg/mL)
blsa_clean = percent_cutoff(blsa_clean, `tnfa`, perc_high = c(0.8), perc_low = c(0.0))


## TOTAL T3 (ng/mL)
blsa_clean = cutoff(blsa_clean, `total_t3`, var_high = c(200, 200), var_low = c(80, 80))

## TOTAL T4 (ug/dL)
blsa_clean = cutoff(blsa_clean, `total_t4`, var_high = c(11.5, 11.5), var_low = c(4.7, 4.7))

## TOTAL TESTOSTERONE (ng/dL)
blsa_clean = cutoff(blsa_clean, `total_testosterone`, var_high = c(950, 60), var_low = c(240, 8))

## TRIGLYCERIDES (mg/dL)
blsa_clean = cutoff(blsa_clean, `triglycerides`, var_high = c(150, 150), var_low = c(0, 0))

## TSH
blsa_clean = cutoff(blsa_clean, `tsh`, var_high = c(4.2, 4.2), var_low = c(0.3, 0.3))

## URIC ACID (mg/dL)
blsa_clean = cutoff(blsa_clean, `uric_acid`, var_high = c(8, 6.1), var_low = c(3.7, 2.7))

## URINE24 CORTISOL, FREE (ug/24hrs)
blsa_clean = cutoff(blsa_clean, `ur24_cortisol`, var_high = c(45, 45), var_low = c(3.5, 3.5))

## VIT B12 (pg/mL)
blsa_clean = cutoff(blsa_clean, `vb12`, var_high = c(914, 914), var_low = c(180, 180))


## N_obs, N_indi and abnormal percent for all biomarkers

## MRmean
N(MRmean)

## RMRmean
N(RMRmean)

## Adiponectin (ug/mL):
N_abn(blsa_clean, adiponectin_cutoff)
N(adiponectin)
adiponectin_abn = N_abn(blsa_clean, adiponectin_cutoff) / N(adiponectin)$N_obs
adiponectin_abper = percent(adiponectin_abn)
adiponectin_abper

## Adjusted Fasting Glucose value (mg/dl)

N_abn(blsa_clean, AFG_cutoff)
N(AFG)
AFG_abn = N_abn(blsa_clean, AFG_cutoff) / N(AFG)$N_obs
AFG_abper = percent(AFG_abn)
AFG_abper


## ALBUMIN (g/dL)
N_abn(blsa_clean, albumin_cutoff)
N(albumin)
albumin_abn = N_abn(blsa_clean, albumin_cutoff) / N(albumin)$N_obs
albumin_abper = percent(albumin_abn)
albumin_abper

## AMYLASE (u/L)
N_abn(blsa_clean, amylase_cutoff)
N(amylase)
amylase_abn = N_abn(blsa_clean, amylase_cutoff) / N(amylase)$N_obs
amylase_abper = percent(amylase_abn)
amylase_abper

## BMI
N_abn(blsa_clean, bmi_cutoff)
N(bmi)
bmi_abn = N_abn(blsa_clean, bmi_cutoff) / N(bmi)$N_obs
bmi_abper = percent(bmi_abn)
bmi_abper

## Waist circumference (cm)

N_abn(blsa_clean, waist_circumference_cutoff)
N(waist_circumference)
waist_circumference_abn = N_abn(blsa_clean, waist_circumference_cutoff) / N(waist_circumference)$N_obs
waist_circumference_abper = percent(waist_circumference_abn)
waist_circumference_abper

## BIOAVAILABLE TESTOSTERONE (ng/dL)
N_abn(blsa_clean, bio_testosterone_cutoff)
N(bio_testosterone)
bio_testosterone_abn = N_abn(blsa_clean, bio_testosterone_cutoff) / N(bio_testosterone)$N_obs
bio_testosterone_abper = percent(bio_testosterone_abn)
bio_testosterone_abper

## CHOLESTEROL (mg/dL)

N_abn(blsa_clean, cholesterol_cutoff)
N(cholesterol)
cholesterol_abn = N_abn(blsa_clean, cholesterol_cutoff) / N(cholesterol)$N_obs
cholesterol_abper = percent(cholesterol_abn)
cholesterol_abper

## Composed Systolic blood pressure (mmHg)

N_abn(blsa_clean, systolic_bp_cutoff)
N(systolic_bp)
systolic_bp_abn = N_abn(blsa_clean, systolic_bp_cutoff) / N(systolic_bp)$N_obs
systolic_bp_abper = percent(systolic_bp_abn)
systolic_bp_abper

## Composed Diastolic blood pressure (mmHg)
N_abn(blsa_clean, diastolic_bp_cutoff)
N(diastolic_bp)
diastolic_bp_abn = N_abn(blsa_clean, diastolic_bp_cutoff) / N(diastolic_bp)$N_obs
diastolic_bp_abper = percent(diastolic_bp_abn)
diastolic_bp_abper

## Core body temperature (F)
N_abn(blsa_clean, corebody_temp_cutoff)
N(corebody_temp)
corebody_temp_abn = N_abn(blsa_clean, corebody_temp_cutoff) / N(corebody_temp)$N_obs
corebody_temp_abper = percent(corebody_temp_abn)
corebody_temp_abper

## CREATININE, SERUM (mg/dL)
N_abn(blsa_clean, creatinine_serum_cutoff)
N(creatinine_serum)
creatinine_serum_abn = N_abn(blsa_clean, creatinine_serum_cutoff) / N(creatinine_serum)$N_obs
creatinine_serum_abper = percent(creatinine_serum_abn)
creatinine_serum_abper

## CRP (ug/mL)
N_abn(blsa_clean, crp_cutoff)
N(crp)
crp_abn = N_abn(blsa_clean, crp_cutoff) / N(crp)$N_obs
crp_abper = percent(crp_abn)
crp_abper

## D DIMER (mg/mL)
N_abn(blsa_clean, d_dimer_cutoff)
N(d_dimer)
d_dimer_abn = N_abn(blsa_clean, d_dimer_cutoff) / N(d_dimer)$N_obs
d_dimer_abper = percent(d_dimer_abn)
d_dimer_abper

## DEHYDROEPIANDROSTERONE SULFATE (ug/dL) different by age
N_abn(blsa_clean, dehy_sulfate_cutoff)
N(dehy_sulfate)
dehy_sulfate_abn = N_abn(blsa_clean, dehy_sulfate_cutoff) / N(dehy_sulfate)$N_obs
dehy_sulfate_abper = percent(dehy_sulfate_abn)
dehy_sulfate_abper

## ESTRADIOL (ng/DL)
N_abn(blsa_clean, estradiol_cutoff)
N(estradiol)
estradiol_abn = N_abn(blsa_clean, estradiol_cutoff) / N(estradiol)$N_obs
estradiol_abper = percent(estradiol_abn)
estradiol_abper

## FIBRINOGEN (mg/dL)
N_abn(blsa_clean, fibrinogen_cutoff)
N(fibrinogen)
fibrinogen_abn = N_abn(blsa_clean, fibrinogen_cutoff) / N(fibrinogen)$N_obs
fibrinogen_abper = percent(fibrinogen_abn)
fibrinogen_abper

## FREE T3 (pg/mL)
N_abn(blsa_clean, free_t3_cutoff)
N(free_t3)
free_t3_abn = N_abn(blsa_clean, free_t3_cutoff) / N(free_t3)$N_obs
free_t3_abper = percent(free_t3_abn)
free_t3_abper

## FREE T4 (ng/dL)
N_abn(blsa_clean, free_t4_cutoff)
N(free_t4)
free_t4_abn = N_abn(blsa_clean, free_t4_cutoff) / N(free_t4)$N_obs
free_t4_abper = percent(free_t4_abn)
free_t4_abper

## Ghrelin (pg/ml)

N_abn(blsa_clean, ghrelin_percent_cutoff)
N(ghrelin)
ghrelin_abn = N_abn(blsa_clean, ghrelin_percent_cutoff) / N(ghrelin)$N_obs
ghrelin_abper = percent(ghrelin_abn)
ghrelin_abper

## HBA1C (%)

N_abn(blsa_clean, hba1c_cutoff)
N(hba1c)
hba1c_abn = N_abn(blsa_clean, hba1c_cutoff) / N(hba1c)$N_obs
hba1c_abper = percent(hba1c_abn)
hba1c_abper

## HDL-CHOLESTEROL (mg/dL)

N_abn(blsa_clean, hdl_cutoff)
N(hdl)
hdl_abn = N_abn(blsa_clean, hdl_cutoff) / N(hdl)$N_obs
hdl_abper = percent(hdl_abn)
hdl_abper

## HOMOCYSTEINE (umol/L) we choose age = 65
N_abn(blsa_clean, homocysteine_cutoff)
N(homocysteine)
homocysteine_abn = N_abn(blsa_clean, homocysteine_cutoff) / N(homocysteine)$N_obs
homocysteine_abper = percent(homocysteine_abn)
homocysteine_abper

## IL15 (pg/ml)
N_abn(blsa_clean, il15_cutoff)
N(il15)
il15_abn = N_abn(blsa_clean, il15_cutoff) / N(il15)$N_obs
il15_abper = percent(il15_abn)
il15_abper

## IL18
N_abn(blsa_clean, il18_cutoff)
N(il18)
il18_abn = N_abn(blsa_clean, il18_cutoff) / N(il18)$N_obs
il18_abper = percent(il18_abn)
il18_abper

## IL6

N_abn(blsa_clean, il6_percent_cutoff)
N(il6)
il6_abn = N_abn(blsa_clean, il6_percent_cutoff) / N(il6)$N_obs
il6_abper = percent(il6_abn)
il6_abper

## IL1RA (pg/mL) Google value

N_abn(blsa_clean, il1ra_percent_cutoff)
N(il1ra)
il1ra_abn = N_abn(blsa_clean, il1ra_percent_cutoff) / N(il1ra)$N_obs
il1ra_abper = percent(il1ra_abn)
il1ra_abper

## INSULIN-LIKE GROWTH FACTOR 1 (ng/mL) 
N_abn(blsa_clean, insu_gf_cutoff)
N(insu_gf)
insu_gf_abn = N_abn(blsa_clean, insu_gf_cutoff) / N(insu_gf)$N_obs
insu_gf_abper = percent(insu_gf_abn)
insu_gf_abper

## LDL-CHOLESTEROL (mg/dL)

N_abn(blsa_clean, ldl_cutoff)
N(ldl)
ldl_abn = N_abn(blsa_clean, ldl_cutoff) / N(ldl)$N_obs
ldl_abper = percent(ldl_abn)
ldl_abper

## lean mass no cutoff



## LEPTIN (ng/mL) based on the age
N_abn(blsa_clean, leptin_cutoff)
N(leptin)
leptin_abn = N_abn(blsa_clean, leptin_cutoff) / N(leptin)$N_obs
leptin_abper = percent(leptin_abn)
leptin_abper

## PLASMINOGEN ACTIVATOR INHIBITOR 1 (ng/mL)

N_abn(blsa_clean, pai_percent_cutoff)
N(pai)
pai_abn = N_abn(blsa_clean, pai_percent_cutoff) / N(pai)$N_obs
pai_abper = percent(pai_abn)
pai_abper

## RESISTIN (ng/mL)

N_abn(blsa_clean, resistin_percent_cutoff)
N(resistin)
resistin_abn = N_abn(blsa_clean, resistin_percent_cutoff) / N(resistin)$N_obs
resistin_abper = percent(resistin_abn)
resistin_abper

## SERUM CORTISOL (ug/dL) For the morning
N_abn(blsa_clean, serum_cortisol_cutoff)
N(serum_cortisol)
serum_cortisol_abn = N_abn(blsa_clean, serum_cortisol_cutoff) / N(serum_cortisol)$N_obs
serum_cortisol_abper = percent(serum_cortisol_abn)
serum_cortisol_abper


## skinfold (biceps) (cm), After data cleaning, all data are be considered as outliers and be removed
N_abn(blsa_clean, skin_fold_cutoff)
N(skin_fold)
skin_fold_abn = N_abn(blsa_clean, skin_fold_cutoff) / N(skin_fold)$N_obs
skin_fold_abper = percent(skin_fold_abn)
skin_fold_abper


## SODIUM (mmol/L)
N_abn(blsa_clean, sodium_cutoff)
N(sodium)
sodium_abn = N_abn(blsa_clean, sodium_cutoff) / N(sodium)$N_obs
sodium_abper = percent(sodium_abn)
sodium_abper

## STNF-RI (pg/mL) 
## percentile cutoff
N_abn(blsa_clean, stnf_r1_percent_cutoff)
N(stnf_r1)
stnf_r1_abn = N_abn(blsa_clean, stnf_r1_percent_cutoff) / N(stnf_r1)$N_obs
stnf_r1_abper = percent(stnf_r1_abn)
stnf_r1_abper

## STNF-RII (pg/mL) 
## percentile cutoff
N_abn(blsa_clean, stnf_r2_percent_cutoff)
N(stnf_r2)
stnf_r2_abn = N_abn(blsa_clean, stnf_r2_percent_cutoff) / N(stnf_r2)$N_obs
stnf_r2_abper = percent(stnf_r2_abn)
stnf_r2_abper

## TNFA (pg/mL)

N_abn(blsa_clean, tnfa_percent_cutoff)
N(tnfa)
tnfa_abn = N_abn(blsa_clean, tnfa_percent_cutoff) / N(tnfa)$N_obs
tnfa_abper = percent(tnfa_abn)
tnfa_abper

## TOTAL T3 (ng/mL)
N_abn(blsa_clean, total_t3_cutoff)
N(total_t3)
total_t3_abn = N_abn(blsa_clean, total_t3_cutoff) / N(total_t3)$N_obs
total_t3_abper = percent(total_t3_abn)
total_t3_abper

## TOTAL T4 (ng/mL)
N_abn(blsa_clean, total_t4_cutoff)
N(total_t4)
total_t4_abn = N_abn(blsa_clean, total_t4_cutoff) / N(total_t4)$N_obs
total_t4_abper = percent(total_t4_abn)
total_t4_abper

## TOTAL TESTOSTERONE (ng/dL)
N_abn(blsa_clean, total_testosterone_cutoff)
N(total_testosterone)
total_testosterone_abn = N_abn(blsa_clean, total_testosterone_cutoff) / N(total_testosterone)$N_obs
total_testosterone_abper = percent(total_testosterone_abn)
total_testosterone_abper

## TRIGLYCERIDES (mg/dL)
N_abn(blsa_clean, triglycerides_cutoff)
N(triglycerides)
triglycerides_abn = N_abn(blsa_clean, triglycerides_cutoff) / N(triglycerides)$N_obs
triglycerides_abper = percent(triglycerides_abn)
triglycerides_abper

## TSH
N_abn(blsa_clean, tsh_cutoff)
N(tsh)
tsh_abn = N_abn(blsa_clean, tsh_cutoff) / N(tsh)$N_obs
tsh_abper = percent(tsh_abn)
tsh_abper


## URIC ACID (mg/dL)
N_abn(blsa_clean, uric_acid_cutoff)
N(uric_acid)
uric_acid_abn = N_abn(blsa_clean, uric_acid_cutoff) / N(uric_acid)$N_obs
uric_acid_abper = percent(uric_acid_abn)
uric_acid_abper

## URINE24 CORTISOL, FREE (ug/24hrs)
N_abn(blsa_clean, ur24_cortisol_cutoff)
N(ur24_cortisol)
ur24_cortisol_abn = N_abn(blsa_clean, ur24_cortisol_cutoff) / N(ur24_cortisol)$N_obs
ur24_cortisol_abper = percent(ur24_cortisol_abn)
ur24_cortisol_abper

## VIT B12 (pg/mL)
N_abn(blsa_clean, vb12_cutoff)
N(vb12)
vb12_abn = N_abn(blsa_clean, vb12_cutoff) / N(vb12)$N_obs
vb12_abper = percent(vb12_abn)
vb12_abper










# Save blsa_clean

write_csv(blsa_clean, "output/data/blsa_clean.csv")







# Choose the biomarkers used in BLSA!!!!!!!!!!!!!!!


blsa_clean

range(blsa_clean$visit, na.rm = T)

new_blsa_clean = 
  blsa_clean %>% 
  select(c(id,
           visit,
           age,
           gender,
           MRmean,
           RMRmean,
           bmi, 
           systolic_bp,
           diastolic_bp,
           waist_circumference,
           hba1c,
           uric_acid,
           albumin, 
           amylase, 
           creatinine_serum,
           vb12,
           ur24_cortisol,
           sodium,
           tsh,
           total_t4,
           free_t3,
           free_t4,
           serum_cortisol,
           cholesterol,
           hdl,
           ldl,
           triglycerides,
           pai,
           crp,
           il6,
           adiponectin,
           leptin,
           resistin,
           total_testosterone,
           bio_testosterone,
           homocysteine,
           estradiol,
           dehy_sulfate,
           insu_gf,
           AFG,
           bmi_cutoff, 
           systolic_bp_cutoff,
           diastolic_bp_cutoff,
           waist_circumference_cutoff,
           hba1c_cutoff,
           uric_acid_cutoff,
           albumin_cutoff, 
           amylase_cutoff, 
           creatinine_serum_cutoff,
           vb12_cutoff,
           ur24_cortisol_cutoff,
           sodium_cutoff,
           tsh_cutoff,
           total_t4_cutoff,
           free_t3_cutoff,
           free_t4_cutoff,
           serum_cortisol_cutoff,
           cholesterol_cutoff,
           hdl_cutoff,
           ldl_cutoff,
           triglycerides_cutoff,
           pai_percent_cutoff,
           crp_cutoff,
           il6_percent_cutoff,
           adiponectin_cutoff,
           leptin_cutoff,
           resistin_percent_cutoff,
           total_testosterone_cutoff,
           bio_testosterone_cutoff,
           homocysteine_cutoff,
           estradiol_cutoff,
           dehy_sulfate_cutoff,
           insu_gf_cutoff,
           AFG_cutoff))

new_blsa_clean$nbm <- rowSums(!is.na(new_blsa_clean %>% select(bmi:AFG)))
new_blsa_clean  =
  new_blsa_clean %>% 
  mutate(aldata = '')
new_blsa_clean$aldata[new_blsa_clean$nbm >= 24] <- 1
new_blsa_clean$aldata[new_blsa_clean$nbm < 24] <- 0

new_blsa_clean = 
  new_blsa_clean %>% 
  filter(visit < 19)


## Compare and check the data by tabulated table with Dan:

sort = new_blsa_clean %>% select(c(visit, nbm, aldata))
sort = sort[order(sort$visit),]

sort0 = 
  sort %>% 
  filter(aldata == 0) %>% 
  group_by(visit) %>% 
  summarize(n0 = n())

sort1 = 
  sort %>% 
  filter(aldata == 1) %>% 
  group_by(visit) %>% 
  summarize(n1 = n())


merge(x = sort0, y = sort1, 
      by = "visit", all.x = TRUE) %>% 
  mutate(Total = n0 + n1)


## Pro-rated value (AL score):
names(new_blsa_clean)
new_blsa_clean = 
  new_blsa_clean %>% 
  filter(aldata == 1)
new_blsa_clean$nab <- rowSums(new_blsa_clean %>% select(bmi_cutoff:AFG_cutoff) == "abnormal")

new_blsa_clean = 
  new_blsa_clean %>% 
  mutate(al_score = round((nab / nbm) * 34))

# Save new_blsa_clean
write_csv(new_blsa_clean, "output/data/new_blsa_clean.csv")



## histogram of AL score:

al_score_gender_plot =
  new_blsa_clean %>% 
  ggplot(., aes(x = al_score, fill = gender)) + 
  geom_histogram(bins = 35) + 
  facet_wrap(~gender) +
  theme(legend.position = "none") +
  labs(caption = "By gender")


al_score_all_plot =  
  new_blsa_clean %>% 
  ggplot() +
  geom_histogram(aes(x = al_score), bins = 55, fill = "#008000") +       
  theme(legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) + 
  labs(caption = "All")

al_score_fig = al_score_gender_plot + al_score_all_plot
al_score_fig
ggsave("output/figures/al_score_fig.pdf")


## MR distribution
MR_gender_plot = 
  new_blsa_clean %>% 
  ggplot() +
  geom_histogram(aes(x = MRmean, fill = gender), bins = 100) +
  facet_wrap(~gender) +
  theme(legend.position = "none") +
  labs(caption = "By gender") + 
  scale_fill_manual(values = c("#EC748B", "#6EB1DE"))

MR_all_plot =  
  new_blsa_clean %>% 
  ggplot() +
  geom_histogram(aes(x = MRmean, fill = gender), bins = 100) +
  theme(legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) + 
  labs(caption = "All") + 
  scale_fill_manual(values = c("#EC748B", "#6EB1DE"))

MR_fig = MR_gender_plot + MR_all_plot
MR_fig
ggsave("output/figures/MR_fig.pdf")

## RMR distribution
RMR_gender_plot = 
  new_blsa_clean %>% 
  ggplot() +
  geom_histogram(aes(x = RMRmean, fill = gender), bins = 100) +
  facet_wrap(~gender) +
  theme(legend.position = "none") +
  labs(caption = "By gender") + 
  scale_fill_manual(values = c("#EC748B", "#6EB1DE"))

RMR_all_plot =  
  new_blsa_clean %>% 
  ggplot() +
  geom_histogram(aes(x = RMRmean, fill = gender), bins = 100) +
  theme(legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) + 
  labs(caption = "All") + 
  scale_fill_manual(values = c("#EC748B", "#6EB1DE"))

RMR_fig = RMR_gender_plot + RMR_all_plot
RMR_fig
ggsave("output/figures/RMR_fig.pdf")





## Association of AL score:
## Composite score vs Age:
prvage_gender_plot = 
  new_blsa_clean %>% 
  ggplot(., aes(x = al_score, y = age, color = gender)) + 
  geom_point(alpha = 0.6) + 
  geom_smooth(method = lm, se = F, color = "orange") +
  facet_wrap(~gender) +
  labs(caption = "By gender") +
  theme(legend.position = "none") +
  stat_cor(method = "pearson", label.x = 0, label.y = 17, 
           cor.coef.name = "rho", color = "black")

prvage_all_plot = 
  new_blsa_clean %>% 
  ggplot(., aes(x = al_score, y = age, color = gender)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = lm, se = F, color = "orange") +
  labs(caption = "All") +
  theme(legend.position = "none") + 
  stat_cor(method = "pearson", label.x = 0, label.y = 17, 
           cor.coef.name = "rho", color = "black")


prvage_fig = prvage_gender_plot / prvage_all_plot
prvage_fig
ggsave("output/figures/ALscore_age.pdf")

## Composite score vs MRmean
prvmr_gender_plot = 
  new_blsa_clean %>% 
  ggplot(., aes(x = al_score, y = MRmean, color = gender)) + 
  geom_point(alpha = 0.6) + 
  geom_smooth(method = lm, se = F, color = "orange") +
  facet_wrap(~gender) +
  labs(caption = "By gender") +
  theme(legend.position = "none") +
  stat_cor(method = "pearson", label.x = 0, label.y = 2.5, 
           cor.coef.name = "rho", color = "black")

prvmr_all_plot = 
  new_blsa_clean %>% 
  ggplot(., aes(x = al_score, y = MRmean, color = gender)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = lm, se = F, color = "orange") +
  labs(caption = "All") +
  theme(legend.position = "none") + 
  stat_cor(method = "pearson", label.x = 0, label.y = 2.5, 
           cor.coef.name = "rho", color = "black")


prvmr_fig = prvmr_gender_plot / prvmr_all_plot
prvmr_fig
ggsave("output/figures/ALscore_MRmean.pdf")



## Composite score vs RMRmean
prvrmr_gender_plot = 
  new_blsa_clean %>% 
  ggplot(., aes(x = al_score, y = RMRmean, color = gender)) + 
  geom_point(alpha = 0.6) + 
  geom_smooth(method = lm, se = F, color = "orange") +
  facet_wrap(~gender) +
  labs(caption = "By gender") +
  theme(legend.position = "none") +
  stat_cor(method = "pearson", label.x = 0, label.y = 4e-05, 
           cor.coef.name = "rho", color = "black")

prvrmr_all_plot = 
  new_blsa_clean %>% 
  ggplot(., aes(x = al_score, y = RMRmean, color = gender)) + 
  geom_point(alpha = 0.4) + 
  geom_smooth(method = lm, se = F, color = "orange") +
  labs(caption = "All") +
  theme(legend.position = "none") + 
  stat_cor(method = "pearson", label.x = 0, label.y = 4e-05, 
           cor.coef.name = "rho", color = "black")


prvrmr_fig = prvrmr_gender_plot / prvrmr_all_plot
prvrmr_fig
ggsave("output/figures/ALscore_RMRmean.pdf")




### (GEE for MR and RMR by gender and All)
new_blsa_clean.gee = 
  new_blsa_clean %>% 
  mutate(age2 = age^2) %>% 
  na.omit()

new_blsa_clean_male.gee = 
  new_blsa_clean.gee %>% 
  filter(gender == "M") %>% 
  na.omit()

new_blsa_clean_female.gee = 
  new_blsa_clean.gee %>% 
  filter(gender == "F") %>% 
  na.omit()

## AL ~ MRmean (All)
mral.all =
       geeglm(scale(MRmean) ~ scale(al_score) + age + age2 + gender, 
       data = new_blsa_clean.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable")
tab_model(mral.all, 
          show.reflvl = T, 
          show.intercept = T, 
          show.se = T,
          p.style = "numeric_stars")

## AL ~ MRmean (Male)
mral.male =
       geeglm(scale(MRmean) ~ scale(al_score) + age + age2, 
       data = new_blsa_clean_male.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable")
tab_model(mral.male, 
          show.reflvl = T, 
          show.intercept = T, 
          show.se = T,
          p.style = "numeric_stars")

## AL ~ MRmean (Female)
mral.female =
       geeglm(scale(MRmean) ~ scale(al_score) + age + age2, 
       data = new_blsa_clean_female.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable")
tab_model(mral.female, 
          show.reflvl = T, 
          show.intercept = T, 
          show.se = T,
          p.style = "numeric_stars")


## AL ~ RMRmean (All) 
rmral.all =
       geeglm(scale(RMRmean) ~ scale(al_score) + age + age2 + gender, 
       data = new_blsa_clean.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable") 
tab_model(rmral.all, 
          show.reflvl = T, 
          show.intercept = T, 
          show.se = T,
          p.style = "numeric_stars")


## AL ~ RMRmean (Male)
rmral.male = 
       geeglm(scale(RMRmean) ~ scale(al_score) + age + age2, 
       data = new_blsa_clean_male.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable")
tab_model(rmral.male, 
          show.reflvl = T, 
          show.intercept = T, 
          show.se = T,
          p.style = "numeric_stars")

## AL ~ RMRmean (Female)
rmral.female =       
       geeglm(scale(RMRmean) ~ scale(al_score) + age + age2, 
       data = new_blsa_clean_female.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable")
tab_model(rmral.female, 
          show.reflvl = T, 
          show.intercept = T, 
          show.se = T,
          p.style = "numeric_stars")

## Add BMI:
## AL ~ MRmean (all)
mral_bmi.all =
       geeglm(scale(MRmean) ~ scale(al_score) + age + age2 + gender + bmi, 
       data = new_blsa_clean.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable")
tab_model(mral_bmi.all, 
          show.reflvl = T, 
          show.intercept = T, 
          show.se = T,
          p.style = "numeric_stars")

## AL ~ MRmean (Male)
mral_bmi.male =
       geeglm(scale(MRmean) ~ scale(al_score) + age + age2 + bmi, 
       data = new_blsa_clean_male.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable")
tab_model(mral_bmi.male, 
          show.reflvl = T, 
          show.intercept = T, 
          show.se = T,
          p.style = "numeric_stars")


## AL ~ MRmean (Female)
mral_bmi.female = 
       geeglm(scale(MRmean) ~ scale(al_score) + age + age2 + bmi, 
       data = new_blsa_clean_female.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable")
tab_model(mral_bmi.female, 
          show.reflvl = T, 
          show.intercept = T, 
          show.se = T,
          p.style = "numeric_stars")


## AL ~ RMRmean (All)
rmral_bmi.all =
       geeglm(scale(RMRmean) ~ scale(al_score) + age + age2 + gender + bmi, 
       data = new_blsa_clean.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable")
tab_model(rmral_bmi.all, 
          show.reflvl = T, 
          show.intercept = T, 
          show.se = T,
          p.style = "numeric_stars")

## AL ~ RMRmean (Male)
rmral_bmi.male = 
       geeglm(scale(RMRmean) ~ scale(al_score) + age + age2 + bmi, 
       data = new_blsa_clean_male.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable")
tab_model(rmral_bmi.male, 
          show.reflvl = T, 
          show.intercept = T,
          show.se = T,
          p.style = "numeric_stars")

## AL ~ RMRmean (Female)
rmral_bmi.female = geeglm(scale(RMRmean) ~ scale(al_score) + age + age2 + bmi, 
       data = new_blsa_clean_female.gee, 
       family = gaussian(link = "identity"),
       id = id,
       corstr = "exchangeable")

tab_model(rmral_bmi.female, 
          show.reflvl = T, 
          show.intercept = T, 
          show.se = T,
          p.style = "numeric_stars")
