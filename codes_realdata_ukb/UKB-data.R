######
library(tidyverse)
library(data.table)
library(mice)
load("~/R/UKB/UKB_yogurt/data/Covariate.rdata")
load("~/R/UKB/UKB_yogurt/data/Covariate_add.rdata")
load("~/R/UKB/UKB_yogurt/data/Covariate_add2.rdata")
load("~/R/UKB/UKB_yogurt/data/Teadt2.rdata")
load("~/R/UKB/UKB_yogurt/data/Deathdt1.rdata")
load("~/R/UKB/UKB_yogurt/data/Popudata.rdata")
load("~/R/UKB/UKB_yogurt/data/Airpollution.rdata")
load("~/R/UKB/UKB_yogurt/data/physical_MET.rdata")
load("~/R/UKB/UKB_yogurt/data/Covariate_add3.rdata")



names(physical_MET)[1] = "f.eid"
(physical_MET) %>%
  merge(Covariate, by = "f.eid") %>% 
  merge(Covariate_add, by = "f.eid") %>%
  merge(Covariate_add2, by = "f.eid") %>%
  merge(Covariate_add3, by = "f.eid") %>%
  merge(Popudata, by = "f.eid")-> dtset

dtset %>% filter(!is.na(bone_density), !is.na(Type_PA_HeavyDIY), !is.na(blood_glucose)) %>% 
  dplyr::select("f.eid", "Ethnicity", "sex", "age", "edu", "Summins_activity",
         "BMI","Alcohol","Smoking_status",
         "highcholesterol","hypertension", 
         "Diabete", "fruit", "vitamin","redmeat","calcium","milk_new",
         "Type_PA_warking", "Type_PA_other", "Type_PA_strenuous",
         "Type_PA_LightDIY", "Type_PA_HeavyDIY", 
         "Type_PA_None", "blood_glucose", "bone_density", "gpc_1", "gpc_2", 
         "gpc_3", "gpc_4", "gpc_5", "gpc_6", "gpc_7", "gpc_8", "gpc_9", 
         "gpc_10", "gpc_11", "gpc_12", "gpc_13", "gpc_14", "gpc_15", "gpc_16", 
         "gpc_17", "gpc_18", "gpc_19", "gpc_20", "gpc_21", "gpc_22", "gpc_23", 
         "gpc_24", "gpc_25", "gpc_26", "gpc_27", "gpc_28", "gpc_29", "gpc_30", 
         "gpc_31", "gpc_32", "gpc_33", "gpc_34", "gpc_35", "gpc_36", "gpc_37", 
         "gpc_38", "gpc_39", "gpc_40") ->dtset1
summary(dtset1)
data = dtset1

m = 1
init = mice(data, maxit=0)
meth = init$method
predM = init$predictorMatrix
predM[, c("f.eid", "Type_PA_warking", "Type_PA_other", "Type_PA_strenuous",
          "Type_PA_LightDIY", "Type_PA_HeavyDIY", 
          "Type_PA_None", "blood_glucose", "bone_density")]=0
imp = mice(data,  predictorMatrix=predM, m=m)
complete(imp, 1) -> dt
save(dt, file = "dt.rdata")
##################


