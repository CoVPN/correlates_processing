# Sys.setenv(TRIAL = "janssen_pooled_real")
renv::activate(here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------

# load required libraries and functions
library(tidyverse)
library(here)
library(methods)
library(SuperLearner)
library(e1071)
library(glmnet)
library(kyotil)
library(argparse)
library(vimp)
library(nloptr)
library(RhpcBLASctl)
library(aucm)
library(mice)
library(conflicted)
library(gam)
library(xgboost)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
# conflict_prefer("omp_set_num_threads", "RhpcBLASctl")
load(paste0("output/", Sys.getenv("TRIAL"), "/objects_for_running_SL.rda"))
# load(paste0("output/", Sys.getenv("TRIAL"), "/plac_top2learners_SL_discreteSL.rda"))
# source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
# source(here("code", "utils.R")) # get CV-AUC for all algs

load(paste0("output/", Sys.getenv("TRIAL"), "/sl_riskscore_slfits.rda"))
# Predict on vaccine arm
if(!any(sapply(c("COVE", "ENSEMBLE"), grepl, study_name))){
  dat.ph1.vacc <- inputMod %>%
    filter(Riskscorecohortflag == 1 & Trt == 1) %>%
    # Keep only variables to be included in risk score analyses
    select(Ptid, Trt, all_of(endpoint), paste0(sub("1rscore", "", endpoint), paste0(vaccAUC_timepoint, "rauc")), 
           all_of(risk_vars), RiskscoreAUCflag) %>%
    # Drop any observation with NA values in Ptid, Trt, or endpoint!
    drop_na(Ptid, Trt, all_of(endpoint))
} else if(study_name == "COVE") {                    # BSEROPOS CHANGE MADE
  dat.ph1.vacc <- inputMod %>%
    filter(Riskscorecohortflag == 1 & Trt == 1) %>%
    # Keep only variables to be included in risk score analyses
    select(Ptid, Trt, all_of(endpoint), all_of(risk_vars)) %>%
    # Drop any observation with NA values in Ptid, Trt, or endpoint!
    drop_na(Ptid, Trt, all_of(endpoint))
  
} else {
  dat.ph1.vacc <- inputMod %>%
    filter(Riskscorecohortflag == 1 & Trt == 1) %>%
    # Keep only variables to be included in risk score analyses
    select(Ptid, Trt, all_of(endpoint), all_of(risk_vars)) %>%
    # Drop any observation with NA values in Ptid, Trt, or endpoint!
    drop_na(Ptid, Trt, all_of(endpoint))
}

X_covars2adjust_vacc <- dat.ph1.vacc %>%
  select(all_of(risk_vars))

# Impute missing values in any variable included in risk_vars using the mice package!
print("Make sure data is clean before conducting imputations!")
X_covars2adjust_vacc <- impute_missing_values(X_covars2adjust_vacc, risk_vars)

# Scale X_covars2adjust_vacc to have mean 0, sd 1 for all vars
for (a in colnames(X_covars2adjust_vacc)) {
  X_covars2adjust_vacc[[a]] <- scale(X_covars2adjust_vacc[[a]],
                                     center = mean(X_covars2adjust_vacc[[a]], na.rm = T),
                                     scale = sd(X_covars2adjust_vacc[[a]], na.rm = T)
  )
}

X_riskVars_vacc <- X_covars2adjust_vacc

pred_on_vaccine <- predict(sl_riskscore_slfits, newdata = X_riskVars_vacc, onlySL = TRUE)$pred %>%
  as.data.frame()

if(!any(sapply(c("COVE", "ENSEMBLE"), grepl, study_name))){
  vacc <- dat.ph1.vacc %>% select(Ptid, all_of(endpoint), paste0(sub("1rscore", "", endpoint), paste0(vaccAUC_timepoint, "rauc")), RiskscoreAUCflag)
} else {
  vacc <- dat.ph1.vacc %>% select(Ptid, all_of(endpoint))
}

vacc <- bind_cols(vacc, pred_on_vaccine) %>%
  rename(pred = V1) %>%
  mutate(risk_score = log(pred / (1 - pred)),
         standardized_risk_score = scale(risk_score,
                                         center = mean(risk_score, na.rm = T),
                                         scale = sd(risk_score, na.rm = T))) 

if(study_name == "COVE"){
  pred_on_plac_bseropos <- predict(sl_riskscore_slfits, 
                                   newdata = plac_bseropos %>% select(names(X_riskVars_vacc)), 
                                   onlySL = TRUE)$pred %>%
    as.data.frame()
  
  plac_bseropos <- bind_cols(plac_bseropos %>% select(Ptid, all_of(endpoint)), pred_on_plac_bseropos) %>%
    rename(pred = V1) %>%
    mutate(risk_score = log(pred / (1 - pred)),
           standardized_risk_score = scale(risk_score,
                                           center = mean(risk_score, na.rm = T),
                                           scale = sd(risk_score, na.rm = T))) 
}

if(!any(sapply(c("COVE", "ENSEMBLE"), grepl, study_name))){
  # AUC for vaccine arm computed only on cohort with RiskscoreAUCflag==1 and based off endpoint rauc!
  AUCvacc <- vacc %>% 
    filter(RiskscoreAUCflag == 1) %>%
    mutate(AUCchar = format(round(fast.auc(pred, get(paste0(sub("1rscore", "", endpoint), paste0(vaccAUC_timepoint, "rauc")))), 3), nsmall = 3)) %>%
    distinct(AUCchar)
} else {
  AUCvacc <- vacc %>% 
    mutate(AUCchar = format(round(fast.auc(pred, get(endpoint)), 3), nsmall = 3)) %>%
    distinct(AUCchar)
}


vacc <- vacc %>% mutate(AUCchar = AUCvacc$AUCchar) 

write.csv(vacc, here("output", Sys.getenv("TRIAL"), "vaccine_ptids_with_riskscores.csv"), row.names = FALSE)
write.csv(plac_bseropos, here("output", Sys.getenv("TRIAL"), "plac_bseropos_ptids_with_riskscores.csv"), row.names = FALSE)

