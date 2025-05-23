# Sys.setenv(TRIAL = "moderna_mock")
# Sys.setenv(TRIAL = "moderna_real")
# Sys.setenv(TRIAL = "janssen_pooled_mock")
# Sys.setenv(TRIAL = "azd1222") # Astra-Zeneca
# Sys.setenv(TRIAL = "prevent19") # Novavax
# Sys.setenv(TRIAL = "vat08m") # Sanofi
# Sys.setenv(TRIAL = "vat08_combined") # Sanofi
# Sys.setenv(TRIAL = "janssen_pooled_partA") 
# Sys.setenv(TRIAL = "janssen_sa_partA_3008") 
# Sys.setenv(TRIAL = "butantan")
# Sys.setenv(TRIAL = "moderna_boost")
# Sys.setenv(TRIAL = "covail")

print("GET_RISKSCORES.R")

# # Since risk scores are generated for VAT08m and VAT08b using the combined data, 
# # if risk score code is called with study_name = VAT08b, change the study_name to VAT08m
# riskscore_called_using_STUDYNAME <- config$study_name
# if(config$study_name == "VAT08b"){
#   study_name = "VAT08m"
# }

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)

source("code/loadlibraries_readinputdata.R")

if(study_name %in% c("ENSEMBLE", "MockENSEMBLE", "PREVENT19", "AZD1222", "VAT08m", "VAT08", "PROFISCOV", "COVAIL")){
  inputFile <- inputFile %>%
    rename(Ptid = Subjectid)
}else if(study_name == "MockCOVE"){
  inputFile <- inputFile %>%
    rename(Ptid = X)
}

# Identify the risk demographic variable names that will be used to compute the risk score
# Identify the endpoint variable
if(study_name %in% c("COVE", "MockCOVE")){
  endpoint <- "EventIndPrimaryD57"
  risk_timepoint <- 57
  studyName_for_report <- "COVE"
  inputMod <- inputFile
  if(study_name %in% c("COVE")){
    risk_vars <- c(
      "MinorityInd", "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
      "Black", "Asian", "NatAmer", "PacIsl",
      "Multiracial", "Other",
      "Notreported", "Unknown",
      "HighRiskInd", "Sex", "Age", "BMI"
    )
  }

  if(study_name %in% c("MockCOVE")){ # as MinorityInd variable is absent in mock!
    risk_vars <- c(
      "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown", 
      "Black", "Asian", "NatAmer", "PacIsl",  
      "Multiracial", "Other", 
      "Notreported", "Unknown",
      "HighRiskInd", "Sex", "Age", "BMI"
    )
  }
  original_risk_vars <- risk_vars
}


if(study_name %in% c("ENSEMBLE", "MockENSEMBLE")){
  
  if(Sys.getenv("TRIAL") != "janssen_sa_partA_3008"){
    risk_vars <- c(
      "EthnicityHispanic","EthnicityNotreported", "EthnicityUnknown",
      "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
      "URMforsubcohortsampling", "HighRiskInd", "HIVinfection", 
      "Sex", "Age", "BMI",
      "Country.X1", "Country.X2", "Country.X3", "Country.X4", "Country.X5", "Country.X6", "Country.X7", 
      "Region.X1", "Region.X2", 
      "CalDtEnrollIND.X1"
    )
    
    # Store original original risk variables as well to check in check_if_SL_needs_be_run.R!
    original_risk_vars <- c(
      "EthnicityHispanic","EthnicityNotreported", "EthnicityUnknown",
      "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
      "URMforsubcohortsampling", "HighRiskInd", "HIVinfection", 
      "Sex", "Age", "BMI",
      "Country", "Region", "CalendarDateEnrollment"
    )
    
    if(run_prod){
      risk_vars <- append(risk_vars, c("CalDtEnrollIND.X2", "CalDtEnrollIND.X3"))
    }
    
    endpoint <- "EventIndPrimaryIncludeNotMolecConfirmedD29"
    risk_timepoint <- 29
    studyName_for_report <- "ENSEMBLE"
    
    # Create binary indicator variables for Country and Region
    inputMod <- inputFile %>%
      drop_na(CalendarDateEnrollment, all_of(endpoint)) %>%
      mutate(Sex.rand = sample(0:1, n(), replace = TRUE),
             Sex = ifelse(Sex %in% c(2, 3), Sex.rand, Sex), # assign Sex randomly as 0 or 1 if Sex is 2 or 3.
             Country = as.factor(Country),
             Region = as.factor(Region),
             CalDtEnrollIND = case_when(CalendarDateEnrollment < 28 ~ 0,
                                        CalendarDateEnrollment >= 28 & CalendarDateEnrollment < 56 ~ 1,
                                        CalendarDateEnrollment >= 56 & CalendarDateEnrollment < 84 ~ 2,
                                        CalendarDateEnrollment >= 84 & CalendarDateEnrollment < 112 ~ 3,
                                        CalendarDateEnrollment >= 112 & CalendarDateEnrollment < 140 ~ 4,
                                        CalendarDateEnrollment >= 140 & CalendarDateEnrollment < 168 ~ 5),
             CalDtEnrollIND = as.factor(CalDtEnrollIND)) %>%
      select(-Sex.rand)
    
    rec <- recipe(~ Country + Region + CalDtEnrollIND, data = inputMod)
    dummies <- rec %>%
      step_dummy(Country, Region, CalDtEnrollIND) %>%
      prep(training = inputMod)
    inputMod <- inputMod %>% bind_cols(bake(dummies, new_data = NULL)) 
    # %>%
    #   select(-c(Country, Region, CalDtEnrollIND))
    names(inputMod)<-gsub("\\_",".",names(inputMod))
    
    # # Create interaction variables between Region and CalDtEnrollIND
    # rec <- recipe(EventIndPrimaryIncludeNotMolecConfirmedD29 ~., data = inputMod)
    # int_mod_1 <- rec %>%
    #   step_interact(terms = ~ starts_with("Region"):starts_with("CalDtEnrollIND"))
    # int_mod_1 <- prep(int_mod_1, training = inputMod)
    # inputMod <- bake(int_mod_1, inputMod)
    # names(inputMod)<-gsub("\\_",".",names(inputMod))
    # if(run_prod){
    #   risk_vars <- append(risk_vars, c("Region.X1.x.CalDtEnrollIND.X1", "Region.X1.x.CalDtEnrollIND.X2",
    #                                    "Region.X1.x.CalDtEnrollIND.X3",
    #                                    "Region.X2.x.CalDtEnrollIND.X1", "Region.X2.x.CalDtEnrollIND.X2",
    #                                    "Region.X2.x.CalDtEnrollIND.X3"))
    # }
  } else if(Sys.getenv("TRIAL") == "janssen_sa_partA_3008"){
    
    # Restrict to South African placebo participants and use only Sex, Age and BMI as input variables for risk score!
    inputFile <- inputFile %>% filter(Country == 7, Trt == 0)
    
    risk_vars <- c("Sex", "Age", "BMI")
    
    # Store original original risk variables as well to check in check_if_SL_needs_be_run.R!
    original_risk_vars <- risk_vars
    
    endpoint <- "EventIndPrimaryIncludeNotMolecConfirmedD29"
    risk_timepoint <- 29
    studyName_for_report <- "ENSEMBLE"
    
    # Create binary indicator variables for Country and Region
    inputMod <- inputFile %>%
      drop_na(CalendarDateEnrollment, all_of(endpoint)) %>%
      mutate(Sex.rand = sample(0:1, n(), replace = TRUE),
             Sex = ifelse(Sex %in% c(2, 3), Sex.rand, Sex), # assign Sex randomly as 0 or 1 if Sex is 2 or 3.
             Country = as.factor(Country),
             Region = as.factor(Region),
             CalDtEnrollIND = case_when(CalendarDateEnrollment < 28 ~ 0,
                                        CalendarDateEnrollment >= 28 & CalendarDateEnrollment < 56 ~ 1,
                                        CalendarDateEnrollment >= 56 & CalendarDateEnrollment < 84 ~ 2,
                                        CalendarDateEnrollment >= 84 & CalendarDateEnrollment < 112 ~ 3,
                                        CalendarDateEnrollment >= 112 & CalendarDateEnrollment < 140 ~ 4,
                                        CalendarDateEnrollment >= 140 & CalendarDateEnrollment < 168 ~ 5),
             CalDtEnrollIND = as.factor(CalDtEnrollIND)) %>%
      select(-Sex.rand)
    
    rec <- recipe(~ CalDtEnrollIND, data = inputMod)
    dummies <- rec %>%
      step_dummy(CalDtEnrollIND) %>%
      prep(training = inputMod)
    inputMod <- inputMod %>% bind_cols(bake(dummies, new_data = NULL)) 
    # %>%
    #   select(-c(Country, Region, CalDtEnrollIND))
    names(inputMod)<-gsub("\\_",".",names(inputMod))
    }
}

if(study_name == "PREVENT19"){
  inputFile <- inputFile %>%
    mutate(EventIndPrimaryD1rscore = EventIndPrimaryD1,
           EventIndPrimaryD35rauc = ifelse(RiskscoreAUCflag == 1, EventIndPrimaryD35, NA)
           )
  
  risk_vars <- c(
    "Age", "Sex", "Black", "Asian", "NatAmer", "PacIsl",  
    "Multiracial", "Notreported", "Unknown",
    "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
    "Height", "Weight", "BMI", "HighRiskInd"
  )
  original_risk_vars <- risk_vars
  endpoint <- "EventIndPrimaryD1rscore"
  #endpoint <- paste0(endpoint, "rscore")
  riskscore_timepoint <- 1
  vaccAUC_timepoint <- 35
  studyName_for_report <- "PREVENT19"
  if(args[1] == "onlyUSsubjects"){
    inputMod <- inputFile %>% filter(Country == 0) # Risk score in prevent19 was derived twice. 
                                                   # First risk score (risk_score) was derived only for US subjects (SL model was trained using only US subjects).  
  } else if (args[1] == "allsubjects"){
    inputMod <- inputFile # Second risk score (risk_score2) was derived for all subjects (SL model was trained using both US and Mexican subjects). 
  }
}

if(study_name == "AZD1222"){
  inputFile <- inputFile %>%
    mutate(EventIndPrimaryD1rscore = EventIndPrimaryD1,
           EventIndPrimaryD57rauc = ifelse(RiskscoreAUCflag == 1, EventIndPrimaryD57, NA))
  risk_vars <- c(
    "Age", "Sex", "Black", "Asian", "NatAmer", "PacIsl",  
    "Multiracial", "Notreported", "Unknown",
    "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
    "BMI", "Country.X1", "Country.X2"
  )
  # Store original original risk variables as well to check in check_if_SL_needs_be_run.R!
  original_risk_vars <- c(
    "Age", "Sex", "Black", "Asian", "NatAmer", "PacIsl",  
    "Multiracial", "Notreported", "Unknown",
    "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
    "BMI", "HighRiskInd", "Country"
  )
  
  endpoint <- "EventIndPrimaryD1rscore"
  riskscore_timepoint <- 1
  vaccAUC_timepoint <- 57
  studyName_for_report <- "AZD1222"
  
  # Create binary indicator variables for Country and Region
  inputMod <- inputFile %>%
    #drop_na(all_of(endpoint)) %>%
    mutate(Country = as.factor(Country))
  
  rec <- recipe(~ Country, data = inputMod)
  dummies <- rec %>%
    step_dummy(Country) %>%
    prep(training = inputMod)
  inputMod <- inputMod %>% bind_cols(bake(dummies, new_data = NULL)) 
  # %>%
  #   select(-c(Country, Region, CalDtEnrollIND))
  names(inputMod)<-gsub("\\_",".",names(inputMod))
}

if(study_name %in% c("VAT08", "VAT08m")){
  inputFile <- inputFile %>%
    mutate(EventIndPrimaryD22rscore = EventIndPrimaryD22,
           EventIndPrimaryD43rauc = ifelse(RiskscoreAUCflag == 1, EventIndPrimaryD43, NA),
           pooled.age.grp = ifelse(Age >= 60, 1, 0))
  
  # Assign geographic region: Honduras, not Honduras for the Stage 1 trial
  # Assign geographic region: India, Mexico, Other for the Stage 2 trial
  inputFile <- inputFile %>%
    filter(Trialstage == 1) %>%
    mutate(Country.ind = case_when(Country == 3 ~ "Honduras",
                                   Country != 3 ~ "NotHonduras",
                                   TRUE ~ "Other")) %>%
    bind_rows(inputFile %>%
                filter(Trialstage == 2) %>%
                mutate(Country.ind = case_when(Country == 4 ~ "India",
                                               Country == 9 ~ "Mexico",
                                               TRUE ~ "Other")))
  
  if(study_name == "VAT08m"){
    risk_vars <- c(
      "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
      "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
      "URMforsubcohortsampling", "HighRiskInd", "HIVinfection",
      "Sex", "Age", "pooled.age.grp", "BMI", #"BMI.group", "Height", "Weight", 
      "Country.ind.NotHonduras", #"Country.ind.India", "Country.ind.Mexico", "Country.ind.Other", 
      #"USAInd",  
      "CalDtEnrollIND.X1", "CalDtEnrollIND.X2", "CalDtEnrollIND.X3", "CalDtEnrollIND.X4", "CalDtEnrollIND.X5"
    )
  } else if(study_name == "VAT08"){
    risk_vars <- c(
      "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
      "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
      "URMforsubcohortsampling", "HighRiskInd", "HIVinfection",
      "Sex", "Age", "pooled.age.grp", "BMI", #"BMI.group", "Height", "Weight", 
      "Country.ind.India", "Country.ind.Mexico", "Country.ind.NotHonduras", "Country.ind.Other", 
      #"USAInd",  
      "CalDtEnrollIND.X1", "CalDtEnrollIND.X2", "CalDtEnrollIND.X3", "CalDtEnrollIND.X4", "CalDtEnrollIND.X5"
      # ,
      # "FOI"  # added FOI and generated new risk scores upon Michal's request for MV vs BV manuscript!
    )
  }
  
  # Store original original risk variables as well to check in check_if_SL_needs_be_run.R!
  original_risk_vars <- c(
    "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
    "Black", "Asian", "NatAmer", "PacIsl", "Multiracial", "Notreported", "Unknown",
    "URMforsubcohortsampling", "HighRiskInd", "HIVinfection",
    "Sex", "Age", "pooled.age.grp", "BMI", #"BMI.group", "Height", "Weight", 
    "Country", 
    #"USAInd", 
    "CalendarDateEnrollment"#, "FOI"
  )
  
  endpoint <- "EventIndPrimaryD22rscore"
  riskscore_timepoint <- 1
  vaccAUC_timepoint <- 43
  studyName_for_report <- "VAT08m_VAT08b_combined"

  # Create binary indicator variables for Country and CalendarDateEnrollment
  inputMod <- inputFile %>%
    mutate(Country.ind = as.factor(Country.ind),
           CalDtEnrollIND = case_when(CalendarDateEnrollment < 28 ~ 0,
                                      CalendarDateEnrollment >= 28 & CalendarDateEnrollment < 56 ~ 1,
                                      CalendarDateEnrollment >= 56 & CalendarDateEnrollment < 84 ~ 2,
                                      CalendarDateEnrollment >= 84 & CalendarDateEnrollment < 112 ~ 3,
                                      CalendarDateEnrollment >= 112 & CalendarDateEnrollment < 140 ~ 4,
                                      CalendarDateEnrollment >= 140 & CalendarDateEnrollment < 168 ~ 5, 
                                      CalendarDateEnrollment >= 168 & CalendarDateEnrollment < 196 ~ 6),
           CalDtEnrollIND = as.factor(CalDtEnrollIND)) 
  
  rec <- recipe(~ Country.ind + CalDtEnrollIND, data = inputMod)
  dummies <- rec %>%
    step_dummy(Country.ind, CalDtEnrollIND) %>%
    prep(training = inputMod)
  inputMod <- inputMod %>% bind_cols(bake(dummies, new_data = NULL)) 
  names(inputMod) <- gsub("\\_", ".", names(inputMod))
}



if(study_name == "COVAIL"){
  inputFile <- inputFile %>%
    mutate(EventIndPrimaryD15rscore = EventIndPrimaryD15,
           EventIndPrimaryD15rauc = EventIndPrimaryD15,
           primary_booster_type = case_when(primary_booster_type == "J, J" ~ "J.J",
                                            primary_booster_type == "J, M" ~ "J.M",
                                            primary_booster_type == "J, P" ~ "J.P",
                                            primary_booster_type == "M, M" ~ "M.M",
                                            primary_booster_type == "M, P" ~ "M.P",
                                            primary_booster_type == "P, M" ~ "P.M",
                                            primary_booster_type == "P, P" ~ "P.P")
    ) 
  
  risk_vars <- c(
    "Age", "Age65C", "Sex", "Black", "Asian", "NatAmer", "PacIsl",  
    "Multiracial", "Unknown",
    "EthnicityHispanic", "EthnicityNotreported", 
    "pre.study.booster.until.studydose1.day", "pre.study.booster.until.studydose1.ind",
    "primary.booster.type.J.M", "primary.booster.type.J.P", "primary.booster.type.M.M",
    "primary.booster.type.M.P", "primary.booster.type.P.M", "primary.booster.type.P.P"
  )
  
  original_risk_vars <- c(
    "Age", "Age65C", "Sex", "Black", "Asian", "NatAmer", "PacIsl",  
    "Multiracial", "Unknown",
    "EthnicityHispanic", "EthnicityNotreported", 
    "pre.study.booster.until.studydose1.day", "pre.study.booster.until.studydose1.ind",
    "primary.booster.type"
  )

  endpoint <- "EventIndPrimaryD15rscore"
  #endpoint <- paste0(endpoint, "rscore")
  riskscore_timepoint <- 15
  vaccAUC_timepoint <- 15
  studyName_for_report <- "COVAIL"
  inputMod <- inputFile 
  
  rec <- recipe(~ primary_booster_type, data = inputMod)
  dummies <- rec %>%
    step_dummy(primary_booster_type) %>%
    prep(training = inputMod)
  inputMod <- inputMod %>% bind_cols(bake(dummies, new_data = NULL)) 
  names(inputMod) <- gsub("\\_", ".", names(inputMod))
}

# Check there are no NA values in Riskscorecohortflag!
if(!study_name %in% c("COVE", "PROFISCOV")){
  assertthat::assert_that(
    all(!is.na(inputMod$Riskscorecohortflag)), msg = "NA values present in Riskscorecohortflag!"
  )
  
  # Save inputFile 
  if(study_name %in% c("VAT08m", "VAT08")){
    if(!dir.exists(paste0("output/", Sys.getenv("TRIAL")))){
      dir.create(paste0("output/", Sys.getenv("TRIAL")))
      if(args[1] %in% c("bseroneg", "bseropos"))
        dir.create(paste0("output/", Sys.getenv("TRIAL"), "/", args[1]))
    }
    if((args[1] %in% c("bseroneg", "bseropos"))  &  (!dir.exists(paste0("output/", Sys.getenv("TRIAL"), "/", args[1])))){
      dir.create(paste0("output/", Sys.getenv("TRIAL"), "/", args[1]))
    }
    save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile.RData"))
    if(args[1] == "bseroneg"){
      save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/inputFile.RData"))
    }else if(args[1] == "bseropos"){
      save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/inputFile.RData"))
    }
  }else if(study_name == "PREVENT19"){
    if(!dir.exists(paste0("output/", Sys.getenv("TRIAL")))){
      dir.create(paste0("output/", Sys.getenv("TRIAL")))
      dir.create(paste0("output/", Sys.getenv("TRIAL"), "/", args[1]))
    }
    if(!dir.exists(paste0("output/", Sys.getenv("TRIAL"), "/", args[1]))){
      dir.create(paste0("output/", Sys.getenv("TRIAL"), "/", args[1]))
    }
    save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile.RData"))
    if(args[1] == "onlyUSsubjects"){
      save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/inputFile.RData"))
    }else if(args[1] == "allsubjects"){
      save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", args[1], "/inputFile.RData"))
    }
  }else{
    if(!dir.exists(paste0("output/", Sys.getenv("TRIAL")))){
      dir.create(paste0("output/", Sys.getenv("TRIAL")))
    }
    save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile.RData"))
  }
}

if(study_name == "COVE"){
  save(inputFile, file = paste0("output/", Sys.getenv("TRIAL"), "/", "inputFile.RData"))
}

# if(study_name == "VAT08" & args[1] == "stackonly"){
#   source("code/stack_bseroneg_bseropos.R")
# }

source(here("code", "check_if_SL_needs_be_run.R"))
