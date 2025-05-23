# Sys.setenv(TRIAL = "janssen_pooled_realbAb")
# Sys.setenv(TRIAL = "prevent19")
# Sys.setenv(TRIAL = "covail")
renv::activate(here::here(".."))
# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------
## load libraries and source files #############################################
library(cvAUC)
library(conflicted)
library(tidyverse)
library(vimp)
library(kyotil)
library(grid)
library(gridExtra)
library(cowplot)
library(here)
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("load", "base")
source(here("code", "study_specific_functions.R"))
source(here("code", "utils.R"))
method <- "method.CC_nloglik" # since SuperLearner relies on this to be in GlobalEnv
ggplot2::theme_set(theme_cowplot())

print("TABLES_FIGURES.R")

if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  load(file = here("output", Sys.getenv("TRIAL"), args[1], "objects_for_running_SL.rda"))
}else{
  load(file = here("output", Sys.getenv("TRIAL"), "objects_for_running_SL.rda"))
}

rm(Y, X_riskVars, weights, maxVar)

if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  load(file = here("output", Sys.getenv("TRIAL"), args[1], "cvsl_risk_placebo_cvaucs.rda"))
}else{
  load(file = here("output", Sys.getenv("TRIAL"), "cvsl_risk_placebo_cvaucs.rda"))
}


######## Table of demographic variables used to derive the risk score ##########
if(study_name == "COVAIL"){
  dat <- inputMod %>%
    filter(Riskscorecohortflag == 1) %>%
    select(all_of(risk_vars)) 
  
  risk_vars <- dat %>%
    map(~ sum(is.na(.))) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Variable Name") %>%
    mutate(V1 = paste0(V1, "/", nrow(dat), " (", format(round((V1 / nrow(dat)) * 100, 1), nsmall = 1), "%)")) %>%
    get_defs_comments_riskVars() %>%
    rename(`Total missing values` = V1) %>%
    select(`Variable Name`, Definition, `Total missing values`, Comments) 
} else {
  dat <- inputMod %>%
    filter(Riskscorecohortflag == 1 & Trt == 0) %>%
    select(all_of(risk_vars)) 
  
  risk_vars <- dat %>%
    map(~ sum(is.na(.))) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Variable Name") %>%
    mutate(V1 = paste0(V1, "/", nrow(dat), " (", format(round((V1 / nrow(dat)) * 100, 1), nsmall = 1), "%)")) %>%
    get_defs_comments_riskVars() %>%
    rename(`Total missing values` = V1) %>%
    select(`Variable Name`, Definition, `Total missing values`, Comments) 
}


if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  risk_vars %>% write.csv(here("output", Sys.getenv("TRIAL"), args[1], "risk_vars.csv"))
}else{
  risk_vars %>% write.csv(here("output", Sys.getenv("TRIAL"), "risk_vars.csv"))
}

rm(risk_vars)

######## learner-screens #######################################################
caption <- "All learner-screen combinations (28 in total) used as input to the superlearner."

if (run_prod) {
  tab <- risk_placebo_cvaucs %>%
    filter(!Learner %in% c("SL", "Discrete SL")) %>%
    select(Learner, Screen) %>%
    mutate(
      Screen = fct_relevel(Screen, c(
        "all", "glmnet", "univar_logistic_pval",
        "highcor_random"
      )),
      Learner = as.factor(Learner),
      Learner = fct_relevel(Learner, c(
        "SL.mean", "SL.glm", "SL.glm.interaction",
        "SL.glmnet", "SL.gam", 
        "SL.xgboost", "SL.ranger.imp"
      ))
    ) %>%
    arrange(Learner, Screen) %>%
    distinct(Learner, Screen) %>%
    rename("Screen*" = Screen)
} else {
  tab <- risk_placebo_cvaucs %>%
    filter(!Learner %in% c("SL", "Discrete SL")) %>%
    select(Learner, Screen) %>%
    mutate(
      Screen = fct_relevel(Screen, c(
        "all", "glmnet", "univar_logistic_pval",
        "highcor_random"
      )),
      Learner = as.factor(Learner),
      Learner = fct_relevel(Learner, c("SL.mean", "SL.glm"))
    ) %>%
    arrange(Learner, Screen) %>%
    distinct(Learner, Screen) %>%
    rename("Screen*" = Screen)
}

if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  tab %>% write.csv(here("output", Sys.getenv("TRIAL"), args[1], "learner-screens.csv"))
}else{
  tab %>% write.csv(here("output", Sys.getenv("TRIAL"), "learner-screens.csv"))
}


######## SLperformance-plac ####################################################
if(study_name=="COVE" | study_name=="MockCOVE"){
  caption <- "Performance of Superlearner and all learner-screen combinations (CV-AUCs with 95\\% CIs) for risk score analyses using placebo group and EventIndPrimaryD57 as outcome. Constraint of np/20 is applied to all learners such that no more than 6 input variables were allowed in any model."
}else if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE"){
  caption <- "Performance of Superlearner and all learner-screen combinations (CV-AUCs with 95\\% CIs) for risk score analyses using placebo group and EventIndPrimaryD29 (including those cases that are not molecularly confirmed) as outcome. Constraint of np/20 is applied to all learners such that no more than 5 input variables were allowed in any model."
}
  
sl.perf <- risk_placebo_cvaucs %>%
  mutate(AUCstr = ifelse(AUC %in% tail(sort(AUC), 1), paste0(AUCstr, "*"), AUCstr)) %>%
  select(Learner, Screen, AUCstr)

if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  sl.perf %>% write.csv(here("output", Sys.getenv("TRIAL"), args[1], "SLperformance-plac.csv"))
}else{
  sl.perf %>% write.csv(here("output", Sys.getenv("TRIAL"), "SLperformance-plac.csv")) 
}

################################################################################
# Forest plots for risk_placebo model, yd57 endpoint
options(bitmapType = "cairo")
if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  if (run_prod) {
    png(file = here("output", Sys.getenv("TRIAL"), args[1], "risk_placebo_cvaucs.png"),
        width = 2000, height = 1100)
    top_learner <- make_forest_plot_prod(risk_placebo_cvaucs)
  } else {
    png(file = here("output", Sys.getenv("TRIAL"), args[1], "risk_placebo_cvaucs.png"),
        width = 2000, height = 700)
    top_learner <- make_forest_plot_demo(risk_placebo_cvaucs)
  }
  grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot,
               ncol = 2)
  dev.off()
}else{
  if (run_prod) {
    png(file = here("output", Sys.getenv("TRIAL"), "risk_placebo_cvaucs.png"),
        width = 2000, height = 1100)
    top_learner <- make_forest_plot_prod(risk_placebo_cvaucs)
  } else {
    png(file = here("output", Sys.getenv("TRIAL"), "risk_placebo_cvaucs.png"),
        width = 2000, height = 700)
    top_learner <- make_forest_plot_demo(risk_placebo_cvaucs)
  }
  grid.arrange(top_learner$top_learner_nms_plot, top_learner$top_learner_plot,
               ncol = 2)
  dev.off()
}


################################################################################
# plot ROC curve and pred.Prob with SL, Discrete SL and top 2 best-performing individual Learners
top2_plac <- bind_rows(
  risk_placebo_cvaucs %>% arrange(-AUC) %>%
    filter(!Learner %in% c("SL", "Discrete SL")) %>%
    dplyr::slice(1:2),
  risk_placebo_cvaucs %>%
    filter(Learner == "SL"),
  risk_placebo_cvaucs %>%
    filter(Learner == "Discrete SL")
) %>%
  mutate(LearnerScreen = ifelse(Learner == "SL", "Super Learner",
                                ifelse(Learner == "Discrete SL", Learner,
                                       paste0(Learner, "_", Screen_fromRun))))

# Get cvsl fit and extract cv predictions
if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  load(file = here("output", Sys.getenv("TRIAL"), args[1], "cvsl_riskscore_cvfits.rda"))
}else{
  load(file = here("output", Sys.getenv("TRIAL"), "cvsl_riskscore_cvfits.rda"))
}

pred <- get_cv_predictions(cv_fit = cvfits[[1]], cvaucDAT = top2_plac)

# plot ROC curve
options(bitmapType = "cairo")
if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  png(file = here("output", Sys.getenv("TRIAL"), args[1], "ROCcurve_riskscore_plac.png"),
      width = 1000, height = 1000)
}else{
  png(file = here("output", Sys.getenv("TRIAL"), "ROCcurve_riskscore_plac.png"),
      width = 1000, height = 1000)
}

p1 <- plot_roc_curves(pred, cvaucDAT = top2_plac)
print(p1)
dev.off()

# plot pred prob plot
options(bitmapType = "cairo")
if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  png(file = here("output", Sys.getenv("TRIAL"), args[1], "predProb_riskscore_plac.png"),
      width = 1100, height = 1400)
}else{
  png(file = here("output", Sys.getenv("TRIAL"), "predProb_riskscore_plac.png"),
      width = 1100, height = 1400)
}

if(!any(sapply(c("COVE", "ENSEMBLE"), grepl, study_name))){
  p2 <- plot_predicted_probabilities(pred, 1)
} else {
  p2 <- plot_predicted_probabilities(pred, risk_timepoint)
}
print(p2)
dev.off()

# Use SuperLearner to generate risk scores!
if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  load(file = here("output", Sys.getenv("TRIAL"), args[1], "risk_placebo_ptids.rda"))
}else{
  load(file = here("output", Sys.getenv("TRIAL"), "risk_placebo_ptids.rda"))
}

plac <- bind_cols(
  risk_placebo_ptids,
  pred %>% filter(Learner == "SL") %>% select(pred, AUCchar)
) %>%
  mutate(risk_score = log(pred / (1 - pred)),
         standardized_risk_score = scale(risk_score,
                                         center = mean(risk_score, na.rm = T),
                                         scale = sd(risk_score, na.rm = T)))

if(study_name %in% c("VAT08m", "VAT08", "PREVENT19")){
  write.csv(plac, here("output", Sys.getenv("TRIAL"), args[1], "placebo_ptids_with_riskscores.csv"),
            row.names = FALSE)
  
  save(top2_plac, file = here("output", Sys.getenv("TRIAL"), args[1], "plac_top2learners_SL_discreteSL.rda"))
}else{
  write.csv(plac, here("output", Sys.getenv("TRIAL"), "placebo_ptids_with_riskscores.csv"),
            row.names = FALSE)
  
  save(top2_plac, file = here("output", Sys.getenv("TRIAL"), "plac_top2learners_SL_discreteSL.rda"))
}

rm(cvfits, pred, p1, p2)
