# Sys.setenv(TRIAL = "janssen_pooled_realbAb")
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
# load(paste0("output/", Sys.getenv("TRIAL"), "/objects_for_running_SL.rda"))
# load(paste0("output/", Sys.getenv("TRIAL"), "/plac_top2learners_SL_discreteSL.rda"))
# source(here("code", "sl_screens.R")) # set up the screen/algorithm combinations
# source(here("code", "utils.R")) # get CV-AUC for all algs

load(paste0("output/", Sys.getenv("TRIAL"), "/sl_riskscore_slfits.rda"))
# Predict on vaccine arm
dat.ph1.vacc <- inputMod %>%
  filter(Riskscorecohortflag == 1 & Trt == 1) %>%
  # Keep only variables to be included in risk score analyses
  select(Ptid, Trt, all_of(endpoint), all_of(risk_vars)) %>%
  # Drop any observation with NA values in Ptid, Trt, or endpoint!
  drop_na(Ptid, Trt, all_of(endpoint))

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

vacc <- bind_cols(
  dat.ph1.vacc %>% select(Ptid, all_of(endpoint)),
  pred_on_vaccine
) %>%
  rename(pred = V1) %>%
  # add AUC
  mutate(AUCchar = format(round(fast.auc(pred, get(endpoint)), 3), nsmall = 3),
         risk_score = log(pred / (1 - pred)),
         standardized_risk_score = scale(risk_score,
                            center = mean(risk_score, na.rm = T),
                            scale = sd(risk_score, na.rm = T)))

write.csv(vacc, here("output", Sys.getenv("TRIAL"), "vaccine_ptids_with_riskscores.csv"), row.names = FALSE)

# plot ROC curve on vaccinees
pred.obj <- ROCR::prediction(vacc$pred, vacc %>% pull(endpoint))
perf.obj <- ROCR::performance(pred.obj, "tpr", "fpr")

options(bitmapType = "cairo")
png(file = here("output", Sys.getenv("TRIAL"), "ROCcurve_riskscore_vacc_onlySL.png"),
    width = 1000, height = 1000)

print(data.frame(xval = perf.obj@x.values[[1]],
           yval = perf.obj@y.values[[1]],
           learner = paste0("Superlearner (", unique(vacc$AUCchar), ")")) %>% 
  ggplot(aes(x = xval, y = yval, col = learner)) +
  geom_step(lwd = 2) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.direction = "vertical",
    legend.box = "horizontal",
    legend.title=element_text(size=20),
    legend.text=element_text(size=20),
    axis.ticks.length = unit(.35, "cm"),
    axis.text = element_text(size = 23),
    axis.title = element_text(size = 30)
  ) +
  labs(x = "False Positive Rate", y = "True Positive Rate", col = "Model (AUC)") +
  geom_abline(intercept = 0, slope = 1) + 
  scale_color_manual(values = "purple"))

dev.off()

# plot pred prob plot on vaccinees
options(bitmapType = "cairo")
png(file = here("output", Sys.getenv("TRIAL"), "predProb_riskscore_vacc_onlySL.png"),
    width = 1100, height = 700)
# if(study_name == "COVE" | study_name == "MockCOVE"){
#   cases = "Post Day 57 Cases"
# }
# if(study_name == "ENSEMBLE" | study_name == "MockENSEMBLE"){
#   cases = "Post Day 29 Cases"
# }
print(vacc %>%
  mutate(Ychar = ifelse(get(endpoint) == 0, "Non-Cases", paste0("Post Day ", risk_timepoint, " Cases"))) %>%
  ggplot(aes(x = Ychar, y = pred, color = Ychar)) +
  geom_jitter(width = 0.06, size = 3, shape = 21, fill = "white") +
  geom_violin(alpha = 0.05, color = "black", lwd=1.5) +
  geom_boxplot(alpha = 0.05, width = 0.15, color = "black", outlier.size = NA, outlier.shape = NA, lwd=1.5) +
  theme_bw() +
  #scale_color_manual(values = c("#56B4E9", "#E69F00")) +
  scale_color_manual(values = c("#00468B", "#8B0000")) +
  labs(y = "Predicted probability of COVID-19 disease", x = "") +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 25),
    axis.text = element_text(size = 23),
    axis.ticks.length = unit(.35, "cm"),
    axis.title.y = element_text(size = 30)
  ))
dev.off()
