{
library(kyotil)
library(copcor)
library(methods)
library(dplyr)
library(digest)
set.seed(98109)
 

if(Sys.getenv("TRIAL")=="") stop("Environmental variable TRIAL not defined!!!!!!!!!!!!!!")
TRIAL=Sys.getenv("TRIAL")

config <- config::get(config = TRIAL)
for(opt in names(config))eval(parse(text = paste0(names(config[opt])," <- config[[opt]]")))
 
data_name = paste0(attr(config, "config"), "_data_processed_with_riskscore.csv")

# disabling lower level parallelization in favor of higher level of parallelization
# set parallelization in openBLAS and openMP
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1L)
stopifnot(blas_get_num_procs() == 1L)
omp_set_num_threads(1L)

verbose=Sys.getenv("VERBOSE")=="1"

# if this flag is true, then the N IgG binding antibody is reported 
# in the immuno report (but is not analyzed in the cor or cop reports).
include_bindN <- !study_name %in% c("PREVENT19","AZD1222","VAT08")
}


################################################################################
# assay limits
# hardcoded in the code for some, driven by config for others
# For bAb, IU and BAU are the same thing
# all values on BAU or IU
# LOQ can not be NA, it is needed for computing delta

if (!is.null(config$assay_metadata)) {
  # use metadata file for assay when exists
  assay_metadata = read.csv(paste0(dirname(attr(config,"file")),"/",config$assay_metadata))
  
  assays=assay_metadata$assay
  
  # created named lists for assay metadata for easier access, e.g. assay_labels_short["bindSpike"]
  assay_labels=assay_metadata$assay_label; names(assay_labels)=assays
  assay_labels_short=assay_metadata$assay_label_short; names(assay_labels_short)=assays
  uloqs=assay_metadata$uloq; names(uloqs)=assays
  lloqs=assay_metadata$lloq; names(lloqs)=assays
  llods=assay_metadata$lod; names(llods)=assays
  
} else {
  # the following part should be a copy of the same code from reporting2 repo _common.R
  names(assays)=assays # add names so that lapply results will have names
  pos.cutoffs<-llods<-lloqs<-uloqs<-c()
  lloxs=NULL
  if (study_name %in% c("COVE", "MockCOVE", "MockENSEMBLE")) {
    tmp=list(
      bindSpike=c(
        pos.cutoff=10.8424,
        LLOD = 0.3076,
        ULOD = 172226.2,
        LLOQ = 1.7968,
        ULOQ = 10155.95)
      ,
      bindRBD=c(
        pos.cutoff=14.0858,
        LLOD = 1.593648,
        ULOD = 223074,
        LLOQ = 3.4263,
        ULOQ = 16269.23)
      ,
      bindN=c( 
        pos.cutoff=23.4711,
        LLOD = 0.093744,
        ULOD = 52488,
        LLOQ = 4.4897,
        ULOQ = 574.6783)
      ,
      pseudoneutid50=c( 
        pos.cutoff=2.42,# as same lod
        LLOD = 2.42,
        ULOD = NA,
        LLOQ = 4.477,
        ULOQ = 10919)
      ,
      pseudoneutid80=c( 
        pos.cutoff=15.02,# as same lod
        LLOD = 15.02,
        ULOD = NA,
        LLOQ = 21.4786,
        ULOQ = 15368)
      ,
      liveneutmn50=c( 
        pos.cutoff=82.1*0.276,# as same lod
        LLOD = 82.11*0.276,
        ULOD = NA,
        LLOQ =  159.79*0.276,
        ULOQ = 11173.21*0.276)
    )
    
    pos.cutoffs=sapply(tmp, function(x) unname(x["pos.cutoff"]))
    llods=sapply(tmp, function(x) unname(x["LLOD"]))
    lloqs=sapply(tmp, function(x) unname(x["LLOQ"]))
    uloqs=sapply(tmp, function(x) unname(x["ULOQ"]))        
    
    
  } else if(study_name=="ENSEMBLE") {
    
    # data less than pos cutoff is set to pos.cutoff/2
    llods["bindSpike"]=NA 
    lloqs["bindSpike"]=1.7968 
    uloqs["bindSpike"]=238.1165 
    pos.cutoffs["bindSpike"]=10.8424
    
    # data less than pos cutoff is set to pos.cutoff/2
    llods["bindRBD"]=NA                 
    lloqs["bindRBD"]=3.4263                 
    uloqs["bindRBD"]=172.5755    
    pos.cutoffs["bindRBD"]=14.0858
    
    # data less than lod is set to lod/2
    llods["ADCP"]=11.57
    lloqs["ADCP"]=8.87
    uloqs["ADCP"]=211.56
    pos.cutoffs["ADCP"]=11.57# as same lod
    
    llods["bindN"]=0.093744
    lloqs["bindN"]=4.4897
    uloqs["bindN"]=574.6783
    pos.cutoffs["bindN"]=23.4711
    
    # the limits below are different for EUA and Part A datasets
    if (contain(attr(config, "config"), "EUA")) {
      # EUA data
      
      # data less than lloq is set to lloq/2
      llods["pseudoneutid50"]=NA  
      lloqs["pseudoneutid50"]=42*0.0653  #2.7426
      uloqs["pseudoneutid50"]=9484*0.0653 # 619.3052
      pos.cutoffs["pseudoneutid50"]=lloqs["pseudoneutid50"]
      
      # repeat for two synthetic markers that are adapted to SA and LA
      llods["pseudoneutid50sa"]=NA  
      lloqs["pseudoneutid50sa"]=42*0.0653  #2.7426
      uloqs["pseudoneutid50sa"]=9484*0.0653 # 619.3052
      pos.cutoffs["pseudoneutid50sa"]=lloqs["pseudoneutid50sa"]
      
      llods["pseudoneutid50la"]=NA  
      lloqs["pseudoneutid50la"]=42*0.0653  #2.7426
      uloqs["pseudoneutid50la"]=9484*0.0653 # 619.3052
      pos.cutoffs["pseudoneutid50la"]=lloqs["pseudoneutid50la"]
      
      
    } else if (contain(attr(config, "config"), "partA")) {
      # complete part A data
      
      # data less than lloq is set to lloq/2
      llods["pseudoneutid50"]=NA  
      lloqs["pseudoneutid50"]=75*0.0653  #4.8975
      uloqs["pseudoneutid50"]=12936*0.0653 # 844.7208
      pos.cutoffs["pseudoneutid50"]=lloqs["pseudoneutid50"]
    }
    
  } else if(study_name=="PREVENT19") {
    # Novavax
    
    # data less than lloq is set to lloq/2 in the raw data
    llods["bindSpike"]=NA 
    lloqs["bindSpike"]=150.4*0.0090 # 1.3536
    uloqs["bindSpike"]=770464.6*0.0090 # 6934.181
    pos.cutoffs["bindSpike"]=10.8424 # use same as COVE
    
    # data less than lloq is set to lloq/2
    llods["bindRBD"]=NA  
    lloqs["bindRBD"]=1126.7*0.0272  #30.6
    uloqs["bindRBD"]=360348.7*0.0272 # 9801
    pos.cutoffs["bindRBD"]=lloqs["bindRBD"]
    
    # data less than lod is set to lod/2 in the raw data
    llods["pseudoneutid50"]=2.612 # 40 * 0.0653
    lloqs["pseudoneutid50"]=51*0.0653 # 3.3303
    uloqs["pseudoneutid50"]=127411*0.0653 # 8319.938
    pos.cutoffs["pseudoneutid50"]=llods["pseudoneutid50"]
    
    llods["bindN"]=0.093744
    lloqs["bindN"]=4.4897
    uloqs["bindN"]=574.6783
    pos.cutoffs["bindN"]=23.4711
    
    # data less than pos is set to pos/2 by Yiwen
    llods["bindNVXIgG"]=200
    lloqs["bindNVXIgG"]=200
    uloqs["bindNVXIgG"]=2904275
    pos.cutoffs["bindN"]=500
    
    llods["bindNVXIgGIU"]=200/22
    lloqs["bindNVXIgGIU"]=200/22
    uloqs["bindNVXIgGIU"]=2904275/22
    pos.cutoffs["bindNIU"]=500/22
    
    llods["ACE2"]=10
    lloqs["ACE2"]=10
    uloqs["ACE2"]=Inf
    pos.cutoffs["ACE2"]=10
    
    
  } else if(study_name=="AZD1222") {
    
    # data less than lloq is set to lloq/2 in the raw data, Nexelis
    llods["bindSpike"]=NA 
    lloqs["bindSpike"]=62.8*0.0090 # 0.5652
    uloqs["bindSpike"]=238528.4*0.0090 # 2146.756
    pos.cutoffs["bindSpike"]=10.8424 # use same as COVE
    
    # data less than lod is set to lod/2
    llods["pseudoneutid50"]=2.612  
    lloqs["pseudoneutid50"]=56*0.0653 # 3.6568
    uloqs["pseudoneutid50"]=47806*0.0653 # 3121.732
    pos.cutoffs["pseudoneutid50"]=llods["pseudoneutid50"]
    
    # bindN info missing in SAP
    
  } else if(study_name=="HVTN705") {
    
    # get uloqs and lloqs from config
    # config$uloqs is a list before this processing
    if (!is.null(config$uloqs)) uloqs=sapply(config$uloqs, function(x) ifelse(is.numeric(x), x, Inf))  else uloqs=sapply(assays, function(a) Inf)
    if (!is.null(config$lloxs)) lloxs=sapply(config$lloxs, function(x) ifelse(is.numeric(x), x, NA))   else lloxs=sapply(assays, function(a) NA)
    names(uloqs)=assays # this is necessary because config$uloqs does not have names
    names(lloxs)=assays
    
  } else if(study_name=="PROFISCOV") { # Butantan
    
    # lod and lloq are the same
    # data less than lod is set to lloq/2
    
    #SARS-CoV-2 Spike           49 70,000 696 49
    #SARS-CoV-2 Spike (P.1)     32 36,000 463 32
    #SARS-CoV-2 Spike (B.1.351) 72 21,000 333 72
    #SARS-CoV-2 Spike (B.1.1.7) 70 47,000 712 70
    
    lloqs["bindSpike"] <- llods["bindSpike"] <- 49*0.0090 # 0.441
    uloqs["bindSpike"]=70000*0.0090 # 630
    pos.cutoffs["bindSpike"]=696*0.0090 # 15.0
    
    lloqs["bindSpike_P.1"] <- llods["bindSpike_P.1"] <- 32*0.0090 
    uloqs["bindSpike_P.1"]=36000*0.0090 
    pos.cutoffs["bindSpike_P.1"]=463*0.0090 
    
    lloqs["bindSpike_B.1.351"] <- llods["bindSpike_B.1.351"] <- 72*0.0090 
    uloqs["bindSpike_B.1.351"]=21000*0.0090 
    pos.cutoffs["bindSpike_B.1.351"]=333*0.0090 
    
    lloqs["bindSpike_B.1.1.7"] <- llods["bindSpike_B.1.1.7"] <- 70*0.0090 
    uloqs["bindSpike_B.1.1.7"]=47000*0.0090 
    pos.cutoffs["bindSpike_B.1.1.7"]=712*0.0090 
    
    #SARS-CoV-2 S1 RBD           35  30,000 1264 35
    #SARS-CoV-2 S1 RBD (P.1)     91  10,000 572  91
    #SARS-CoV-2 S1 RBD (B.1.351) 53  6,300  368  53
    #SARS-CoV-2 S1 RBD (B.1.1.7) 224 20,000 1111 224
    
    lloqs["bindRBD"] <- llods["bindRBD"] <- 35*0.0272 
    uloqs["bindRBD"]=30000*0.0272 # 630
    pos.cutoffs["bindRBD"]=1264*0.0272 # 15.0
    
    lloqs["bindRBD_P.1"] <- llods["bindRBD_P.1"] <- 91*0.0272 
    uloqs["bindRBD_P.1"]=10000*0.0272 
    pos.cutoffs["bindRBD_P.1"]=572*0.0272 
    
    lloqs["bindRBD_B.1.351"] <- llods["bindRBD_B.1.351"] <- 53*0.0272 
    uloqs["bindRBD_B.1.351"]=6300*0.0272 
    pos.cutoffs["bindRBD_B.1.351"]=368*0.0272 
    
    lloqs["bindRBD_B.1.1.7"] <- llods["bindRBD_B.1.1.7"] <- 224*0.0272 
    uloqs["bindRBD_B.1.1.7"]=20000*0.0272 
    pos.cutoffs["bindRBD_B.1.1.7"]=1111*0.0272 
    
    #SARS-CoV-2 Nucleocapsid 46 80,000 7015 46
    
    lloqs["bindN"] <- llods["bindN"] <- 46*0.00236 
    uloqs["bindN"]=80000*0.00236 
    pos.cutoffs["bindN"]=7015*0.00236 
    
    #LVMN
    llods["liveneutmn50"]=27.56 
    lloqs["liveneutmn50"]=27.84
    uloqs["liveneutmn50"]=20157.44 
    pos.cutoffs["liveneutmn50"]=llods["liveneutmn50"] 
    
  } else if(study_name=="COVEBoost") { 
    # nothing to do, but this is needed so that _common.R can be called for making risk score
    
  } else stop("unknown study_name 1") # not necessary if there is assay metadata file
  
}


###############################################################################
# figure labels and titles for markers
###############################################################################

markers <- c(outer(times[which(times %in% c("B", paste0("Day", config$timepoints)))], assays, "%.%"))

# race labeling
labels.race <- c(
  "White", 
  "Black or African American",
  "Asian", 
  if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & startsWith(attr(config, "config"),"janssen_la")) "Indigenous South American" else "American Indian or Alaska Native",
  "Native Hawaiian or Other Pacific Islander", 
  "Multiracial",
  if ((study_name=="COVE" | study_name=="MockCOVE")) "Other", 
  "Not reported and unknown"
)

# ethnicity labeling
labels.ethnicity <- c(
  "Hispanic or Latino", "Not Hispanic or Latino",
  "Not reported and unknown"
)

labels.assays.short <- c("Anti N IgG (BAU/ml)", 
                         "Anti Spike IgG (BAU/ml)", 
                         "Anti RBD IgG (BAU/ml)", 
                         "Pseudovirus-nAb cID50", 
                         "Pseudovirus-nAb cID80", 
                         "Live virus-nAb cMN50")
names(labels.assays.short) <- c("bindN",
  "bindSpike",
  "bindRBD",
  "pseudoneutid50",
  "pseudoneutid80",
  "liveneutmn50")

# hacky fix for tabular, since unclear who else is using
# the truncated labels.assays.short later
labels.assays.short.tabular <- labels.assays.short

labels.time <- c("Day 1", paste0("Day ", config$timepoints), paste0("D", config$timepoints, " fold-rise over D1"), "D57 fold-rise over D29")

names(labels.time) <- c("B", paste0("Day", config$timepoints), paste0("Delta", config$timepoints, "overB"), "Delta57over29")

# axis labeling
labels.axis <- outer(
  rep("", length(times)),
  labels.assays.short[assays],
  "%.%"
)
labels.axis <- as.data.frame(labels.axis)
rownames(labels.axis) <- times

labels.assays <- c("Binding Antibody to Spike", 
                   "Binding Antibody to RBD",
                   "PsV Neutralization 50% Titer",
                   "PsV Neutralization 80% Titer",
                   "WT LV Neutralization 50% Titer")

names(labels.assays) <- c("bindSpike", 
                          "bindRBD", 
                          "pseudoneutid50",
                          "pseudoneutid80",
                          "liveneutmn50")

# title labeling
labels.title <- outer(
  labels.assays[assays],
  ": " %.% c("Day 1", paste0("Day ", config$timepoints), paste0("D", config$timepoints, " fold-rise over D1"), "D57 fold-rise over D29"),
  paste0
)
labels.title <- as.data.frame(labels.title)
colnames(labels.title) <- times
# NOTE: hacky solution to deal with changes in the number of markers
rownames(labels.title)[seq_along(assays)] <- assays
labels.title <- as.data.frame(t(labels.title))

# creating short and long labels
labels.assays.short <- labels.axis[1, ]
labels.assays.long <- labels.title

# baseline stratum labeling
if (study_name=="COVE" | study_name=="MockCOVE") {
    Bstratum.labels <- c(
      "Age >= 65",
      "Age < 65, At risk",
      "Age < 65, Not at risk"
    )
    
} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    Bstratum.labels <- c(
      "Age < 60, Not at risk",
      "Age < 60, At risk",
      "Age >= 60, Not at risk",
      "Age >= 60, At risk"
    )
    
} else if (study_name %in% c("PREVENT19","AZD1222","NVX_UK302")) {
    Bstratum.labels <- c(
      "Age >= 65",
      "Age < 65"
    )

} else if (study_name %in% c("VAT08")) {
    Bstratum.labels <- c(
      "Age >= 60",
      "Age < 60"
    )

} else if (study_name %in% c("PROFISCOV", "COVAIL", "ILiAD_IB202P", "ILiAD_IB201P")) {
  Bstratum.labels <- c(
    "All"
  )
  
} else if(study_name %in% c("COVEBoost", "HVTN705", "NextGen_Mock")) { 
    # nothing to do, but this is needed so that _common.R can be called for making risk score
  
} else stop("unknown study_name 2")



# baseline stratum labeling
if (study_name=="COVE" | study_name=="MockCOVE") {
    demo.stratum.labels <- c(
      "Age >= 65, URM",
      "Age < 65, At risk, URM",
      "Age < 65, Not at risk, URM",
      "Age >= 65, White non-Hisp",
      "Age < 65, At risk, White non-Hisp",
      "Age < 65, Not at risk, White non-Hisp"
    )

} else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
    demo.stratum.labels <- c(
      "US URM, Age 18-59, Not at risk",
      "US URM, Age 18-59, At risk",
      "US URM, Age >= 60, Not at risk",
      "US URM, Age >= 60, At risk",
      "US White non-Hisp, Age 18-59, Not at risk",
      "US White non-Hisp, Age 18-59, At risk",
      "US White non-Hisp, Age >= 60, Not at risk",
      "US White non-Hisp, Age >= 60, At risk",
      "Latin America, Age 18-59, Not at risk",
      "Latin America, Age 18-59, At risk",
      "Latin America, Age >= 60, Not at risk",
      "Latin America, Age >= 60, At risk",
      "South Africa, Age 18-59, Not at risk",
      "South Africa, Age 18-59, At risk",
      "South Africa, Age >= 60, Not at risk",
      "South Africa, Age >= 60, At risk"
    )

} else if (study_name=="PREVENT19") {
    demo.stratum.labels <- c(
      "US White non-Hisp, Age 18-64, Not at risk",
      "US White non-Hisp, Age 18-64, At risk",
      "US White non-Hisp, Age >= 65, Not at risk",
      "US White non-Hisp, Age >= 65, At risk",
      "US URM, Age 18-64, Not at risk",
      "US URM, Age 18-64, At risk",
      "US URM, Age >= 65, Not at risk",
      "US URM, Age >= 65, At risk",
      "Mexico, Age 18-64",
      "Mexico, Age >= 65"
    )

} else if (study_name=="AZD1222") {
    demo.stratum.labels <- c(
      "US White non-Hisp, Age 18-64",
      "US White non-Hisp, Age >= 65",
      "US URM, Age 18-64",
      "US URM, Age >= 65",
      "Non-US, Age 18-64",
      "Non-US, Age >= 65"
    )

} else if (study_name %in% c("VAT08")) {
#    demo.stratum.labels <- c(
#      "Not HND, Age 18-59",
#      "Not HND, Age >= 60",
#      "HND, Age 18-59",
#      "HND, Age >= 60",
#      "USA, Age 18-59",
#      "USA, Age >= 60",
#      "JPN, Age 18-59",
#      "JPN, Age >= 60"
#    )

    # in this partial dataset, we need to collapse "Not HND, US or JPN, senior" and "HND, senior" due to sparsity
    demo.stratum.labels <- c(
      "Not HND, Age 18-59",
      "Not USA or JPN, Age >= 60",
      "HND, Age 18-59",
      "USA, Age 18-59",
      "USA, Age >= 60",
      "JPN, Age 18-59",
      "JPN, Age >= 60"
    )

} else if (study_name %in% c("PROFISCOV")) {
    demo.stratum.labels <- c("All")

} else if (study_name=="HVTN705") {
    # do nothing

} else if(study_name %in% c("COVEBoost", "HVTN705", "NextGen_Mock", "ILiAD_IB202P", "ILiAD_IB201P","COVAIL")) {
  # nothing to do, but this is needed so that _common.R can be called for making risk score

} else stop("unknown study_name 3")


labels.regions.ENSEMBLE =c("0"="Northern America", "1"="Latin America", "2"="Southern Africa")
regions.ENSEMBLE=0:2
names(regions.ENSEMBLE)=labels.regions.ENSEMBLE

labels.countries.ENSEMBLE=c("0"="United States", "1"="Argentina", "2"="Brazil", "3"="Chile", "4"="Columbia", "5"="Mexico", "6"="Peru", "7"="South Africa")
countries.ENSEMBLE=0:7
names(countries.ENSEMBLE)=labels.countries.ENSEMBLE


###############################################################################
# reproduciblity options
###############################################################################

# NOTE: used in appendix.Rmd to store digests of input raw/processed data files
# hash algorithm picked based on https://csrc.nist.gov/projects/hash-functions
hash_algorithm <- "sha256"


###############################################################################
# theme options
###############################################################################

# fixed knitr chunk options
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  out.width = "80%",
  out.extra = "",
  fig.pos = "H",
  fig.show = "hold",
  fig.align = "center",
  fig.width = 6,
  fig.asp = 0.618,
  fig.retina = 0.8,
  dpi = 600,
  echo = FALSE,
  message = FALSE,
  warning = FALSE
)

# global options
options(
  digits = 6,
  #scipen = 999,
  dplyr.print_min = 6,
  dplyr.print_max = 6,
  crayon.enabled = FALSE,
  bookdown.clean_book = TRUE,
  knitr.kable.NA = "NA",
  repos = structure(c(CRAN = "https://cran.rstudio.com/"))
)

# no complaints from installation warnings
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

# overwrite options by output type
if (knitr:::is_html_output()) {
  #options(width = 80)

  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}
if (knitr:::is_latex_output()) {
  #knitr::opts_chunk$set(width = 67)
  #options(width = 67)
  options(cli.unicode = TRUE)

  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}

# create and set global ggplot theme
# borrowed from https://github.com/tidymodels/TMwR/blob/master/_common.R
theme_transparent <- function(...) {
  # use black-white theme as base
  ret <- ggplot2::theme_bw(...)

  # modify with transparencies
  trans_rect <- ggplot2::element_rect(fill = "transparent", colour = NA)
  ret$panel.background  <- trans_rect
  ret$plot.background   <- trans_rect
  ret$legend.background <- trans_rect
  ret$legend.key        <- trans_rect

  # always have legend below
  ret$legend.position <- "bottom"
  return(ret)
}

library(ggplot2)
theme_set(theme_transparent())
theme_update(
  text = element_text(size = 25),
  axis.text.x = element_text(colour = "black", size = 30),
  axis.text.y = element_text(colour = "black", size = 30)
)

# custom ggsave function with updated defaults
ggsave_custom <- function(filename = default_name(plot),
                          height= 15, width = 21, ...) {
  ggsave(filename = filename, height = height, width = width, ...)
}



preprocess=function(dat_raw, study_name) {
    dat_proc=dat_raw
    
    if(is_ows_trial){
        dat_proc=subset(dat_proc, !is.na(Bserostatus))
    }
    
    
    if(study_name=="ENSEMBLE") {
        # EventTimePrimaryIncludeNotMolecConfirmedD29 are the endpoint of interest and should be used to compute weights
        dat_proc$EventTimePrimaryD29=dat_proc$EventTimePrimaryIncludeNotMolecConfirmedD29
        dat_proc$EventIndPrimaryD29 =dat_proc$EventIndPrimaryIncludeNotMolecConfirmedD29
        dat_proc$EventTimePrimaryD1 =dat_proc$EventTimePrimaryIncludeNotMolecConfirmedD1
        dat_proc$EventIndPrimaryD1  =dat_proc$EventIndPrimaryIncludeNotMolecConfirmedD1
        
    } else if (startsWith(study_name,"VAT08")) {
      # needed for risk score as mapped data does not have these variables defined
      dat_proc$EventTimePrimaryD43=dat_proc$EventTimeFirstInfectionD43
      dat_proc$EventIndPrimaryD43 =dat_proc$EventIndFirstInfectionD43
      dat_proc$EventTimePrimaryD22=dat_proc$EventTimeFirstInfectionD22
      dat_proc$EventIndPrimaryD22 =dat_proc$EventIndFirstInfectionD22
      dat_proc$EventTimePrimaryD1 =dat_proc$EventTimeFirstInfectionD1
      dat_proc$EventIndPrimaryD1  =dat_proc$EventIndFirstInfectionD1
      
    } else if (study_name=='COVAIL') {
      dat_proc$EventTimePrimaryD15 =dat_proc$COVIDtimeD22toD181
      dat_proc$EventIndPrimaryD15  =dat_proc$COVIDIndD22toD181
    }
    
    
    # Mar 29, 24 skip this check because things are now more complicated, the dataset does not always have EventTimePrimaryDxx
    # for(tp in timepoints) {
    #   empty = nrow(dat_proc[is.na(dat_proc[["EventTimePrimaryD"%.%tp]]), ])
    #   if (empty>0) {
    #     warning("there are "%.%empty%.%" ptids without event time. subset to ptids without missing EventTimePrimaryD"%.%tp%.%"\n")
    #   }
    #   dat_proc=dat_proc[!is.na(dat_proc[["EventTimePrimaryD"%.%tp]]), ]
    # }
    
    
    # define earlyendpoint 
    if (study_name=='COVAIL') {
      dat_proc$EarlyinfectionD15=dat_proc$EarlyendpointD15
      
    } else if (TRIAL %in% c("prevent19", "azd1222_bAb", "azd1222", "janssen_partA_VL")) {
      for(tp in timepoints) {
        dat_proc[["EarlyendpointD"%.%tp]] <- with(dat_proc, 
          ifelse(get("EarlyinfectionD"%.%tp)==1 | (EventIndPrimaryD1==1 & EventTimePrimaryD1 < get("NumberdaysD1toD"%.%tp) + 7),1,0))
      }
    } else {
      # do nothing
      # EarlyendpointDxx is not necessary anymore. For example, in prevent19_stage2, we use EarlyinfectionDxx to define ph1.
    }
    

    # ENSEMBLE only, since we are not using this variable to define Riskscorecohortflag and we are not doing D29start1 analyses for other trials
    if (study_name %in% c("MockENSEMBLE", "ENSEMBLE")) {
        dat_proc[["EarlyendpointD29start1"]]<- with(dat_proc, ifelse(get("EarlyinfectionD29start1")==1| (EventIndPrimaryD1==1 & EventTimePrimaryD1 < get("NumberdaysD1toD29") + 1),1,0))
    } else {
#        # commented out on Aug 25, 2022 because it is dangerous
#        # this is not necessary, but it is kept here to make the hash checks for mock datasets happy
#        dat_proc$EarlyinfectionD29start1=dat_proc$EarlyinfectionD29 
    }
   
    dat_proc    
}
