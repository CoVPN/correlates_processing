# Check if only change in input dataset is marker data. 
# If so, simply pull the earlier risk scores and add to the new processed dataset!
generate_new_riskscores <- function(){
  if(study_name == "COVAIL"){
    source(here("code", "clean_output_dir.R"))
    source(here("code", "run_cvsl_riskscore.R"))
    source(here("code", "createRDAfiles_fromSLobjects.R"))
    source(here("code", "tables_figures.R"))
    source(here("code", "constructSL_getSLweights_Modelpredictors.R"))
    source(here("code", "append_risk_score_to_data.R"))
  } else if(study_name == "ENSEMBLE" & Sys.getenv("TRIAL") == "janssen_sa_partA_3008"){
    source(here("code", "clean_output_dir.R"))
    source(here("code", "run_cvsl_riskscore.R"))
    source(here("code", "createRDAfiles_fromSLobjects.R"))
    source(here("code", "tables_figures.R"))
    source(here("code", "constructSL_getSLweights_Modelpredictors.R"))
    source(here("code", "append_risk_score_to_data.R"))
    source(here("code", "test_risk_scores.R"))
  } else {
    source(here("code", "clean_output_dir.R"))
    source(here("code", "run_cvsl_riskscore.R"))
    source(here("code", "createRDAfiles_fromSLobjects.R"))
    source(here("code", "tables_figures.R"))
    source(here("code", "constructSL_getSLweights_Modelpredictors.R"))
    source(here("code", "predict_on_vaccine.R"))
    source(here("code", "append_risk_score_to_data.R"))
    source(here("code", "performance_on_vaccine.R"))
    source(here("code", "test_risk_scores.R"))
  }
}

if(!study_name %in% c("COVE", "PROFISCOV")){
  if((file.exists(paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile_with_riskscore.RData")) |
     (startsWith(Sys.getenv("TRIAL"), "janssen") & endsWith(Sys.getenv("TRIAL"), "EUA")      & file.exists(paste0("output/janssen_pooled_EUA/inputFile_with_riskscore.RData"))) |
     (startsWith(Sys.getenv("TRIAL"), "janssen") & endsWith(Sys.getenv("TRIAL"), "partA")    & file.exists(paste0("output/janssen_pooled_partA/inputFile_with_riskscore.RData"))) | 
     (startsWith(Sys.getenv("TRIAL"), "janssen") & endsWith(Sys.getenv("TRIAL"), "partA_VL") & file.exists(paste0("output/janssen_pooled_partA/inputFile_with_riskscore.RData")))) |
     (study_name == "VAT08" & (file.exists(paste0("output/vat08_combined/inputFile_with_riskscore.RData")))) |
     (study_name == "PREVENT19" & (file.exists(here("output", "prevent19", args[1], "inputFile_with_riskscore.RData")))) 
     ){
    
      if(startsWith(Sys.getenv("TRIAL"), "janssen") & endsWith(Sys.getenv("TRIAL"), "EUA") & file.exists(paste0("output/janssen_pooled_EUA/inputFile_with_riskscore.RData"))){
        load("output/janssen_pooled_EUA/inputFile_with_riskscore.RData")
      } else if(startsWith(Sys.getenv("TRIAL"), "janssen") & endsWith(Sys.getenv("TRIAL"), "partA") & file.exists(paste0("output/janssen_pooled_partA/inputFile_with_riskscore.RData"))){
        load("output/janssen_pooled_partA/inputFile_with_riskscore.RData")
      } else if(startsWith(Sys.getenv("TRIAL"), "janssen") & endsWith(Sys.getenv("TRIAL"), "partA_VL") & file.exists(paste0("output/janssen_pooled_partA/inputFile_with_riskscore.RData"))){
        load("output/janssen_pooled_partA/inputFile_with_riskscore.RData")
      }else if(study_name %in% c("VAT08")){
        load("output/vat08_combined/inputFile_with_riskscore.RData")
      }else if(study_name == "PREVENT19"){
        load(here("output", "prevent19", args[1], "inputFile_with_riskscore.RData"))
      }else{
        load(paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile_with_riskscore.RData"))
      } 
    
    if(study_name == "COVAIL"){
      old_processed <- inputFile_with_riskscore %>%
        rename(pre.study.booster.until.studydose1.day = pre_study_booster_until_studydose1_day,
               pre.study.booster.until.studydose1.ind = pre_study_booster_until_studydose1_ind,
               primary.booster.type = primary_booster_type) %>%
        select(Ptid, Riskscorecohortflag, treatment_actual, all_of(endpoint), all_of(original_risk_vars), risk_score, standardized_risk_score)
      
      new_processed <- inputFile  %>%
        rename(pre.study.booster.until.studydose1.day = pre_study_booster_until_studydose1_day,
               pre.study.booster.until.studydose1.ind = pre_study_booster_until_studydose1_ind,
               primary.booster.type = primary_booster_type) %>%
        select(Ptid, Riskscorecohortflag, treatment_actual, all_of(endpoint), all_of(original_risk_vars))
      
      row.names(old_processed) <- 1:nrow(old_processed)
      row.names(new_processed) <- 1:nrow(new_processed)
    } else {
      old_processed <- inputFile_with_riskscore %>%
        select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), all_of(original_risk_vars), risk_score, standardized_risk_score)
      
      new_processed <- inputFile %>%
        select(Ptid, Riskscorecohortflag, Trt, all_of(endpoint), all_of(original_risk_vars))
    }
    
    
    if(all.equal(old_processed %>% select(-c(risk_score, standardized_risk_score)), new_processed) == TRUE){
      message("Variables related to risk score generation in input data have not changed. Superlearner will not be run. Risk scores from earlier run will be appended to raw data!")
      inputFile_with_riskscore <- left_join(inputFile,
                                            old_processed %>%
                                              select(Ptid, risk_score, standardized_risk_score), by = "Ptid")
      if(study_name == "PREVENT19"){
        save(inputFile_with_riskscore, file = here("output", Sys.getenv("TRIAL"), args[1], "inputFile_with_riskscore.RData"))
      } else {
        save(inputFile_with_riskscore, file = here("output", Sys.getenv("TRIAL"), "inputFile_with_riskscore.RData"))
      }
      
      # if(study_name %in% c("VAT08m", "VAT08"))
      #   args[1] = "SLnotrun"
      
    }else{
      message("Variables related to risk score generation in input data have changed! Superlearner needs to be run and new risk scores generated!")
      generate_new_riskscores()
    }
  }else{
    message(paste0("riskscore_baseline/", paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile_with_riskscore.RData"), " does not exist. Superlearner needs to be run and new risk scores generated!"))
    generate_new_riskscores()
  }
}else if(study_name == "COVE"){
  if(file.exists(paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile_with_riskscore.RData"))){
    load(paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile_with_riskscore.RData"))
    # Check 
    all.equal(names(inputFile_with_riskscore %>% select(Ptid, risk_score, standardized_risk_score)), c("Ptid", "risk_score", "standardized_risk_score"))
    print("Risk scores for Moderna real dataset were generated at Moderna's end using CoVPN Stats/SCHARP code and are already present in input file.")
    print("To suit the Stage 2 COVEBoost trial, risk scores were generated for baseline seropositives and also 15 extra subjects that were in Stage 2 but not in Stage 1 trial.")
    print("Superlearner will not be run!")
    # if(!file.exists(paste0("output/", Sys.getenv("TRIAL")))){
    #   dir.create(paste0("output/", Sys.getenv("TRIAL")))
    # }
  }else{
    message(paste0("riskscore_baseline/", paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile_with_riskscore.RData"), " does not exist. Superlearner needs to be run and new risk scores generated!"))
    generate_new_riskscores()
  }
}else if(study_name == "PROFISCOV"){
  inputFile_with_riskscore <- inputFile %>% mutate(risk_score = 1,
                                                   standardized_risk_score = NA)
  # Check 
  all.equal(names(inputFile_with_riskscore %>% select(Ptid, risk_score)), c("Ptid", "risk_score"))
  print("All subjects for Butantan dataset are assigned risk score of 1. Risk scores are not generated by design. Superlearner will not be run!")
  if(!file.exists(paste0("output/", Sys.getenv("TRIAL")))){
    dir.create(paste0("output/", Sys.getenv("TRIAL")))
  }
  save(inputFile_with_riskscore, file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile_with_riskscore.RData"))
}

# Perform sanity check if mock data!
if(study_name != "VAT08" | (study_name == "VAT08" & args[1] == "stackonly")){
  source(here("code", "unit_test.R"))
}

