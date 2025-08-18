Sys.setenv(TRIAL = "covail")
source(here::here("_common.R"))
# no need to run renv::activate(here::here()) b/c .Rprofile exists

{
library(tidyverse)
library(glue)
library(Hmisc) # wtd.quantile, cut2
library(mice)
library(dplyr)
library(here)
library(mdw)

begin=Sys.time()
}


########################################################################################################
# 1. read mapped data with risk score added

{
  # # load risk score from running risk analysis
  # load(file = paste0('riskscore_baseline/output/',TRIAL,'/inputFile_with_riskscore.RData'))
  # dat_proc <- inputFile_with_riskscore    
  
  # load risk score from a file
  dat.risk = read.csv("/trials/covpn/COVAILcorrelates/analysis/correlates/adata/risk_score.csv")
  # read mapped data
  dat_raw = read.csv(mapped_data)
  dat_proc = preprocess(dat_raw, study_name)   
  names(dat_proc)[[1]]="Ptid"
  dat_proc$risk_score = dat.risk$risk_score[match(dat_proc$Pti, dat.risk$Ptid)]
  dat_proc$standardized_risk_score = dat.risk$standardized_risk_score[match(dat_proc$Ptid, dat.risk$Ptid)]
  
  # bring in imputed variant column
  dat.lineage = read.csv('/trials/covpn/COVAILcorrelates/analysis/correlates/adata/lineages/covail_lineages_export_v1.csv')
  dat_proc$COVIDlineage = dat.lineage$inf1.lineage[match(dat_proc$Ptid, dat.lineage$ptid)]
  dat_proc$COVIDlineageObserved = !dat.lineage$inf1.imputed[match(dat_proc$Ptid, dat.lineage$ptid)]
  # check NA
  stopifnot(!any(is.na(dat_proc$COVIDlineage[dat_proc$ph1.D15==1 & dat_proc$COVIDIndD22toD181==1])))
  stopifnot(!any(is.na(dat_proc$COVIDlineage[dat_proc$ph1.D29==1 & dat_proc$COVIDIndD36toD181==1])))
  # this is not true: !any(is.na(dat_proc$COVIDlineage[dat_proc$ph1.D15==1 & dat_proc$AsympInfectIndD15to181==1]))
  
  # bring in FOI
  dat.foi = read.csv('/trials/covpn/COVAILcorrelates/analysis/correlates/adata/covail_foi_v2.csv')
  dat_proc$FOIoriginal = dat.foi$foi[match(dat_proc$Ptid, dat.foi$ptid)]
  dat_proc$FOIstandardized = scale(dat_proc$FOIoriginal)
  # check NA
  stopifnot(!any(is.na(dat_proc$FOI[dat_proc$ph1.D15==1])))
  stopifnot(!any(is.na(dat_proc$FOI[dat_proc$ph1.D29==1])))
  
  assay_metadata=read.csv(config$assay_metadata)
  assays=assay_metadata$assay
  tmp=assay_metadata$assay[8:nrow(assay_metadata)]
  tcellsubsets = sub("_Wuhan.N", "", tmp[endsWith(tmp, ".N")]); tcellsubsets
  S=tmp[endsWith(tmp, ".S") & startsWith(tmp, "c")]
  S1=c(tcellsubsets%.%"_COV2.CON.S1", tcellsubsets%.%"_BA.4.5.S1")
  S2=c(tcellsubsets%.%"_COV2.CON.S2", tcellsubsets%.%"_BA.4.5.S2")
  N=tcellsubsets%.%"_Wuhan.N"
  tcellvv=c(S1, S2, N)
  
  nAb = setdiff(assays[startsWith(assays, "pseudoneutid50_")], c('pseudoneutid50_MDW'))
  frnt=assays[startsWith(assays, "frnt")]
  
  # special treatment for seromyx
  # there are 198 assays, which are not tracked in assay metadata
  tmp=read.csv("/trials/covpn/COVAILcorrelates/download_data/covail/SeromYxDataJuly23_2025.csv")
  tmp$assay= paste0(tmp$tgttype, "_", tmp$vstrain)
  seromyx=sort(unique(tmp$assay))
  # len(unique(dat$subjid))
  
}


########################################################################################################
# define Senior and race/ethnicity
{
  colnames(dat_proc)[colnames(dat_proc)=="Subjectid"] <- "Ptid" 
  dat_proc <- dat_proc %>% mutate(age.geq.65 = as.integer(Age >= 65))
  dat_proc$Senior = as.integer(dat_proc$Age>=switch(study_name, COVE=65, MockCOVE=65, ENSEMBLE=60, MockENSEMBLE=60, PREVENT19=65, AZD1222=65, VAT08=60, PROFISCOV=NA, COVAIL=65, NVX_UK302=65, stop("unknown study_name 1")))
  
  # ethnicity labeling
  dat_proc$ethnicity <- ifelse(dat_proc$EthnicityHispanic == 1, labels.ethnicity[1], labels.ethnicity[2])
  dat_proc$ethnicity[dat_proc$EthnicityNotreported == 1 | dat_proc$EthnicityUnknown == 1] <- labels.ethnicity[3]
  dat_proc$ethnicity <- factor(dat_proc$ethnicity, levels = labels.ethnicity)
  
  dat_proc$race = 0 # not applicable, but has to define a value so that the next chunk of code can run

  dat_proc$WhiteNonHispanic <- NA
  # WhiteNonHispanic=1 IF race is White AND ethnicity is not Hispanic
  dat_proc$WhiteNonHispanic <-
    ifelse(dat_proc$race == "White" &
             dat_proc$ethnicity == "Not Hispanic or Latino", 1,
           dat_proc$WhiteNonHispanic
    )
  # WhiteNonHispanic=0 IF race is not "white or unknown" OR ethnicity is Hispanic
  dat_proc$WhiteNonHispanic <-
    ifelse(!dat_proc$race %in% c(labels.race[1], last(labels.race)) |
             dat_proc$ethnicity == "Hispanic or Latino", 0,
           dat_proc$WhiteNonHispanic
    )
  dat_proc$MinorityInd = 1-dat_proc$WhiteNonHispanic
  # set NA to 0 in both WhiteNonHispanic and MinorityInd. This means the opposite things for how NA's are interpreted and that gave us the option to use one or the other
  dat_proc$WhiteNonHispanic[is.na(dat_proc$WhiteNonHispanic)] = 0
  dat_proc$MinorityInd[is.na(dat_proc$MinorityInd)] = 0

}




###############################################################################
# 3. stratum variables
{
dat_proc$Bstratum = 1 # there are no demographics stratum for subcohort sampling
names(Bstratum.labels) <- Bstratum.labels
  
dat_proc$demo.stratum = 1 # # there are no demographics stratum for subcohort sampling

# Stratum 1-17 nnaive, 51-67 naive
dat_proc <- dat_proc %>% mutate(tps.stratum = arm)
dat_proc <- dat_proc %>% mutate(tps.stratum.tcell = arm + strtoi(paste0(naive), base = 2) * 50)
if (!is.null(dat_proc$tps.stratum)) table(dat_proc$tps.stratum)


# Wstratum, 1 ~ max(tps.stratum), max(tps.stratum)+1, ..., max(tps.stratum)+4. 
# Used to compute sampling weights. 
# Differs from tps stratum in that case is a separate stratum within each of the four groups defined by Trt and Bserostatus
# A case will have a Wstratum even if its tps.stratum is NA
# The case is defined using EventIndPrimaryD29
dat_proc$Wstratum = dat_proc$tps.stratum

# cases are stratified by stage and vaccine-proximal/vaccine-distal. cases post D181 are sampled as controls
# 91-98 case strata
dat_proc$Wstratum.tcell = dat_proc$tps.stratum.tcell
dat_proc$Wstratum.tcell[with(dat_proc, COVIDIndD22toD91 ==1 & stage==1)]=91 
dat_proc$Wstratum.tcell[with(dat_proc, COVIDIndD92toD181==1 & stage==1)]=92
dat_proc$Wstratum.tcell[with(dat_proc, COVIDIndD22toD91 ==1 & stage==2)]=93 
dat_proc$Wstratum.tcell[with(dat_proc, COVIDIndD92toD181==1 & stage==2)]=94
dat_proc$Wstratum.tcell[with(dat_proc, COVIDIndD22toD91 ==1 & stage==3)]=95 
dat_proc$Wstratum.tcell[with(dat_proc, COVIDIndD92toD181==1 & stage==3)]=96
dat_proc$Wstratum.tcell[with(dat_proc, COVIDIndD22toD91 ==1 & stage==4)]=97 
dat_proc$Wstratum.tcell[with(dat_proc, COVIDIndD92toD181==1 & stage==4)]=98

dat_proc$Wstratum.frnt = dat_proc$Wstratum.tcell
dat_proc$Wstratum.seromyx  = dat_proc$Wstratum.tcell
dat_proc$Wstratum.xassays  = dat_proc$Wstratum.tcell

if (!is.null(dat_proc$Wstratum.tcell)) table(dat_proc$Wstratum.tcell) # variables may be named other than Wstratum
}


################################################################################
# 4. Define ph1, ph2, and weights
# Note that Wstratum may have NA if any variables to form strata has NA
{
# the whole cohort is treated as ph1 and ph2
dat_proc$TwophasesampIndD15 = dat_proc$ph1.D15 
dat_proc$TwophasesampIndD29 = dat_proc$ph1.D29
  
# PP = no violation + marker available at d1 and d15
# Immunemarkerset = PP & no infection between enrollment and D15+6
# ph1.D15 = Immunemarkerset & arm!=3
dat_proc[["ph2.D15"]]=dat_proc$ph1.D15
dat_proc[["wt.D15"]] = 1
dat_proc[["ph2.D92"]]=dat_proc$ph1.D92
dat_proc[["wt.D92"]] = 1
dat_proc[["ph2.D29"]]=dat_proc$ph1.D29
dat_proc[["wt.D29"]] = 1

# filter out arm 16 and 17
dat_proc$ph1.D15.tcell = ifelse (dat_proc$ph1.D15==1 & dat_proc$arm<=15, 1, 0)
dat_proc$ph1.D92.tcell = ifelse (dat_proc$ph1.D92==1 & dat_proc$arm<=15, 1, 0)
dat_proc$ph1.D15.frnt = dat_proc$ph1.D15.tcell
dat_proc$ph1.D92.frnt = dat_proc$ph1.D92.tcell 
dat_proc$ph1.D15.seromyx = dat_proc$ph1.D15.tcell
dat_proc$ph1.D92.seromyx = dat_proc$ph1.D92.tcell 
dat_proc$ph1.D15.xassays = dat_proc$ph1.D15.tcell
dat_proc$ph1.D92.xassays = dat_proc$ph1.D92.tcell 


# N is not included in the def of TwophasesampInd b/c ptids with S also have N, at almost all time points
dat_proc$TwophasesampIndD15.tcell  = apply(dat_proc, 1, function (x) any(!is.na(x[glue("Day15{c(S1,S2)}")])))
dat_proc$TwophasesampIndB.tcell  = apply(dat_proc, 1, function (x) any(!is.na(x[glue("B{c(S1,S2)}")])))
mytable(dat_proc$TwophasesampIndB.tcell, dat_proc$TwophasesampIndD15.tcell, dat_proc$COVIDIndD22toend)

# make two phase indicator 1 if either B or D15 are present
dat_proc$TwophasesampIndD15.tcell = dat_proc$TwophasesampIndB.tcell | dat_proc$TwophasesampIndD15.tcell
# remove TwophasesampIndB.tcell
dat_proc$TwophasesampIndB.tcell = NULL

dat_proc$TwophasesampIndD15.frnt  = apply(dat_proc, 1, function (x) any(!is.na(x[c(glue("B{frnt}"), glue("Day15{frnt}"))])))
dat_proc$TwophasesampIndD15.seromyx  = apply(dat_proc, 1, function (x) any(!is.na(x[c(glue("B{seromyx}"), glue("Day15{seromyx}"))])))


# this shows that cases from D182 on are sampled like controls 
dat_proc$case.period=NA
dat_proc$case.period[dat_proc$COVIDIndD22toend==1]=3
dat_proc$case.period[dat_proc$COVIDIndD92toD181==1]=2
dat_proc$case.period[dat_proc$COVIDIndD22toD91 ==1]=1

with(dat_proc[dat_proc$ph1.D15.tcell==1,], mytable(case.period, TwophasesampIndD15.tcell))

# use TwophasesampIndD15.tcell for both time points
for (tp in c(15,92)) {
  dat_proc[["ph2.D"%.%tp%.%".tcell"]] = dat_proc[["ph1.D"%.%tp%.%".tcell"]] & dat_proc[["TwophasesampIndD15.tcell"]]
  # 6/12/25 bug fix: Wstratum => Wstratum.tcell
  dat_proc = add.wt(dat_proc, ph1="ph1.D"%.%tp%.%".tcell", ph2="ph2.D"%.%tp%.%".tcell", Wstratum="Wstratum.tcell", 
                    wt="wt.D"%.%tp%.%".tcell", verbose=T) 
  
  dat_proc[["ph2.D"%.%tp%.%".frnt"]] = dat_proc[["ph1.D"%.%tp%.%".frnt"]] & dat_proc[["TwophasesampIndD15.frnt"]]
  dat_proc = add.wt(dat_proc, ph1="ph1.D"%.%tp%.%".frnt",  ph2="ph2.D"%.%tp%.%".frnt", Wstratum="Wstratum.frnt", 
                    wt="wt.D"%.%tp%.%".frnt", verbose=T) 
  
  dat_proc[["ph2.D"%.%tp%.%".seromyx"]] = dat_proc[["ph1.D"%.%tp%.%".seromyx"]] & dat_proc[["TwophasesampIndD15.seromyx"]]
  dat_proc = add.wt(dat_proc, ph1="ph1.D"%.%tp%.%".seromyx",  ph2="ph2.D"%.%tp%.%".seromyx", Wstratum="Wstratum.seromyx", 
                    wt="wt.D"%.%tp%.%".seromyx", verbose=T) 
  
  dat_proc[["ph2.D"%.%tp%.%".xassays"]] = dat_proc[["ph2.D"%.%tp%.%".seromyx"]] & dat_proc[["ph2.D"%.%tp%.%".tcell"]]
  dat_proc = add.wt(dat_proc, ph1="ph1.D"%.%tp%.%".xassays",  ph2="ph2.D"%.%tp%.%".xassays", Wstratum="Wstratum.xassays", 
                    wt="wt.D"%.%tp%.%".xassays", verbose=T) 
  
}

}


###############################################################################
# 5. impute missing biomarkers in ph2 (assay imputation)

# NO MORE IMPUTATIONS can be done in this section
# any additional imputation should be done in section 11 so that casedeletion results are reproducible

#     impute vaccine and placebo, baseline pos and neg, separately
#     use all assays (not bindN)
#     use baseline, each time point, but not Delta


{
# nAb data
  
# loop through the time points
# first impute (B, D29, D57) among TwophasesampIndD57==1
# next impute (B, D29) among TwophasesampIndD29==1

#### impute D15 BA.4/BA.5 nAb among ph1.D15

dat.tmp.impute <- subset(dat_proc, ph1.D15 == 1)

# first, do it based on a linear regression with D29 BA.4/BA.5 nAb and number of days between the D15 and D29 visit
fit = lm(Day15pseudoneutid50_BA.4.BA.5~Day29pseudoneutid50_BA.4.BA.5, dat.tmp.impute)
predicted = predict(fit, subset(dat.tmp.impute, select=Day29pseudoneutid50_BA.4.BA.5))
dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5 = ifelse(is.na(dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5), predicted, dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5)

# there are still some ptids with missing D15 BA4.5 b/c not all have D29 BA4.5. impute with BA1: mypairs(dat.tmp.impute["Day15"%.%assays[1:5]])
fit = lm(Day15pseudoneutid50_BA.4.BA.5~Day15pseudoneutid50_BA.1, dat.tmp.impute)
predicted = predict(fit, subset(dat.tmp.impute, select=Day15pseudoneutid50_BA.1))
dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5 = ifelse(is.na(dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5), predicted, dat.tmp.impute$Day15pseudoneutid50_BA.4.BA.5)

# save a copy for Lauren before assigning imputed values
dat_proc$Day15pseudoneutid50_BA.4.BA.5_unimputed = dat_proc$Day15pseudoneutid50_BA.4.BA.5
# populate dat_proc with the imputed values
imp.markers='Day15pseudoneutid50_BA.4.BA.5'
dat_proc[dat_proc[["ph1.D15"]]==1, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["ph1.D15"]]==1, "Ptid"], dat.tmp.impute$Ptid), ]

assertthat::assert_that(
  all(complete.cases(dat_proc[dat_proc[["ph1.D15"]] == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)


#### none missing at D29 among those at risk at D29
dat.tmp.impute <- subset(dat_proc, ph1.D15==1 & COVIDtimeD22toD181>NumberdaysD15toD29 & AsympInfectIndD15to29==0)
with(dat.tmp.impute, print(table(!is.na(get("Day29"%.%assays[1])), !is.na(get("Day15"%.%assays[1])))))
# thus, no missingness actually

#### none missing at D91 among those at risk at D91
dat.tmp.impute <- subset(dat_proc, ph1.D15==1 & COVIDtimeD22toD181>NumberdaysD15toD91 & AsympInfectIndD15to91==0)
with(dat.tmp.impute, print(table(!is.na(get("Day91"%.%assays[1])), !is.na(get("Day15"%.%assays[1])))))
# thus, no missingness actually

#### none missing at D181 among those at risk at D181
dat.tmp.impute <- subset(dat_proc, ph1.D15==1 & COVIDtimeD22toD181>NumberdaysD15toD181 & AsympInfectIndD15to181==0)
with(dat.tmp.impute, print(table(!is.na(get("Day181"%.%assays[1])), !is.na(get("Day15"%.%assays[1])))))
# thus, no missingness actually
}


{
# T cell data

# before imputation, log transform markers because distributions are skewed. we want models for imputation work on transformed variables
  # for FS, imputation is done on log scale. later, in step 6, they will be transformed back to the linear scale
tmp=c(outer(c("B", "Day15", "Day91", "Day181"), 
            tcellvv,
            "%.%"))
dat_proc[tmp] = log10 (dat_proc[tmp])
  
dat_proc$arm.factor = as.factor(dat_proc$arm)

# impute S1, S2, and N-stim T cell markers at B and D15 together. not impute markers at Day91 or D181
# impute different arms and naive/nnaive together, but pass arm and naive as covariates

n.imp=1
tp=15

# if (ii==1) {
#   imp.markers=c(outer(c("B", "Day"%.%tp), 
#                       setdiff(tcellvv, c("cd8_FS_Wuhan.N", "cd8_FS_COV2.CON.S1","cd8_FS_BA.4.5.S1", "cd8_FS_COV2.CON.S2", "cd8_FS_BA.4.5.S2")), 
#                       "%.%"))
# } else {
  imp.markers=c(outer(c("B", "Day"%.%tp), tcellvv, "%.%"))
# }

# add arm and naive to the imputation dataset
imp.markers =  c(imp.markers, "arm.factor", "naive")

dat.tmp.impute <- subset(dat_proc, get("TwophasesampIndD15.tcell") == 1)

imp <- dat.tmp.impute %>% select(all_of(imp.markers))     

# there is one -Inf in a FS
imp[imp==-Inf]=NA

if(any(is.na(imp))) {
  # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
  dat.tmp.impute[, imp.markers] <- mice::complete(imp, action = 1)
}                

# missing markers imputed properly?
assertthat::assert_that(
  all(complete.cases(dat.tmp.impute[, imp.markers])),
  msg = "missing markers imputed properly?"
)    

# populate dat_proc imp.markers with the imputed values
dat_proc[dat_proc[["TwophasesampIndD15.tcell"]]==1, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["TwophasesampIndD15.tcell"]]==1, "Ptid"], dat.tmp.impute$Ptid), ]

assertthat::assert_that(
  all(complete.cases(dat_proc[dat_proc[["TwophasesampIndD15.tcell"]] == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)




# imputing 0/1 variables with mice is extremely slow, so we do it with glm

# there are no resp variables for FS
imp.markers=c(outer(c("B", "Day"%.%tp), tcellvv[!grepl("FS", tcellvv)], "%.%")) %.%"_resp"

dat.tmp.impute <- subset(dat_proc, get("TwophasesampIndD15.tcell") == 1)

for (a in imp.markers) {
  kp = is.na(dat.tmp.impute[[a]])
  if (any(kp)) {
    f=as.formula(glue("{a} ~ {sub('_resp','',a)} + arm.factor + naive"))
    fit = glm(f, dat.tmp.impute, family=binomial)
    dat.tmp.impute[[a]][kp] = rbinom(sum(kp), 1, prob=predict(fit, dat.tmp.impute[kp,], type="response"))
  } else next
}  

# missing markers imputed properly?
assertthat::assert_that(
  all(complete.cases(dat.tmp.impute[, imp.markers])),
  msg = "missing markers imputed properly?"
)    

# populate dat_proc imp.markers with the imputed values
dat_proc[dat_proc[["TwophasesampIndD15.tcell"]]==1, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["TwophasesampIndD15.tcell"]]==1, "Ptid"], dat.tmp.impute$Ptid), ]

assertthat::assert_that(
  all(complete.cases(dat_proc[dat_proc[["TwophasesampIndD15.tcell"]] == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)

}



{
# FRNT B and D15

# impute FRNT50, FRNT80 at B and D15 together.
# impute different arms and naive/nnaive together, but pass arm and naive as covariates

n.imp=1
tp=15

imp.markers=c(outer(c("B", "Day"%.%tp), frnt, "%.%"))

# add arm and naive to the imputation dataset
imp.markers =  c(imp.markers, "arm.factor", "naive")

dat.tmp.impute <- subset(dat_proc, get("TwophasesampIndD15.frnt") == 1)

imp <- dat.tmp.impute %>% select(all_of(imp.markers))     

# there is one -Inf in a FS
imp[imp==-Inf]=NA

if(any(is.na(imp))) {
  # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
  dat.tmp.impute[, imp.markers] <- mice::complete(imp, action = 1)
}                

# missing markers imputed properly?
assertthat::assert_that(
  all(complete.cases(dat.tmp.impute[, imp.markers])),
  msg = "missing markers imputed properly?"
)    

# populate dat_proc imp.markers with the imputed values
dat_proc[dat_proc[["TwophasesampIndD15.frnt"]]==1, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["TwophasesampIndD15.frnt"]]==1, "Ptid"], dat.tmp.impute$Ptid), ]

assertthat::assert_that(
  all(complete.cases(dat_proc[dat_proc[["TwophasesampIndD15.frnt"]] == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)

}



{
  # FRNT D91
  
  n.imp=1
  
  imp.markers=c(outer(c("Day15", "Day91"), frnt, "%.%"))
  
  # add arm and naive to the imputation dataset
  imp.markers =  c(imp.markers, "arm.factor", "naive", "COVIDIndD22toD181")
  
  dat.tmp.impute <- subset(dat_proc, get("TwophasesampIndD15.frnt") == 1)
  
  imp <- dat.tmp.impute %>% select(all_of(imp.markers))     
  
  # there is one -Inf in a FS
  imp[imp==-Inf]=NA
  
  if(any(is.na(imp))) {
    # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
    imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
    dat.tmp.impute[, imp.markers] <- mice::complete(imp, action = 1)
  }                
  
  # missing markers imputed properly?
  assertthat::assert_that(
    all(complete.cases(dat.tmp.impute[, imp.markers])),
    msg = "missing markers imputed properly?"
  )    
  
  # populate dat_proc imp.markers with the imputed values
  dat_proc[dat_proc[["TwophasesampIndD15.frnt"]]==1, imp.markers] <-
    dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["TwophasesampIndD15.frnt"]]==1, "Ptid"], dat.tmp.impute$Ptid), ]
  
  assertthat::assert_that(
    all(complete.cases(dat_proc[dat_proc[["TwophasesampIndD15.frnt"]] == 1, imp.markers])),
    msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
  )
  
}



{
  # FRNT D181
  
  n.imp=1
  
  imp.markers=c(outer(c("Day15", "Day91", "Day181"), frnt, "%.%"))
  
  # add arm and naive to the imputation dataset
  imp.markers =  c(imp.markers, "arm.factor", "naive", "COVIDIndD22toD181")
  
  # define a new flag to exclude the 24 ptids without the D181 markers
  dat_proc$Day181frnt_missing_cnt = with(dat_proc, is.na(Day181frnt50_D614G) + is.na(Day181frnt50_BA.1) + is.na(Day181frnt50_BA.2) + is.na(Day181frnt50_BA.4.BA.5) + is.na(Day181frnt80_D614G) + is.na(Day181frnt80_BA.1) + is.na(Day181frnt80_BA.2) + is.na(Day181frnt80_BA.4.BA.5) )
  dat_proc$TwophasesampIndD181.frnt = dat_proc$TwophasesampIndD15.frnt & !(dat_proc$Day181frnt_missing_cnt==8 & (dat_proc$COVIDIndD92toD181 == 0 & !is.na(dat_proc$COVIDIndD92toD181)))
  dat_proc$ph2.D181.frnt = dat_proc$ph2.D15.frnt & dat_proc$TwophasesampIndD181.frnt
  
  dat.tmp.impute <- subset(dat_proc, get("TwophasesampIndD181.frnt") == 1)
  
  imp <- dat.tmp.impute %>% select(all_of(imp.markers))     
  
  # there is one -Inf in a FS
  imp[imp==-Inf]=NA
  
  if(any(is.na(imp))) {
    # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
    imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
    dat.tmp.impute[, imp.markers] <- mice::complete(imp, action = 1)
  }                
  
  # missing markers imputed properly?
  assertthat::assert_that(
    all(complete.cases(dat.tmp.impute[, imp.markers])),
    msg = "missing markers imputed properly?"
  )    
  
  # populate dat_proc imp.markers with the imputed values
  dat_proc[dat_proc[["TwophasesampIndD181.frnt"]]==1, imp.markers] <-
    dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["TwophasesampIndD181.frnt"]]==1, "Ptid"], dat.tmp.impute$Ptid), ]
  
  assertthat::assert_that(
    all(complete.cases(dat_proc[dat_proc[["TwophasesampIndD181.frnt"]] == 1, imp.markers])),
    msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
  )
  
}



###############################################################################
# 6. transformation of the markers

# create S-stimulated markers by summing up S1 and S2 on the anti log scale
for (a in c(tcellsubsets%.%"_COV2.CON.S", tcellsubsets%.%"_BA.4.5.S")) {
  for (tp in c("B","Day15","Day91","Day181")) {
    dat_proc[[tp%.%a]] = log10 (10^dat_proc[[tp%.%a%.%"1"]] + 10^dat_proc[[tp%.%a%.%"2"]])
  }
}

# transform FS back to the linear scale
fs=tcellsubsets[grepl("FS", tcellsubsets)]
for (a in c(outer(fs, c("_COV2.CON.S", "_BA.4.5.S", "_COV2.CON.S1", "_BA.4.5.S1", "_COV2.CON.S2", "_BA.4.5.S2", "_Wuhan.N"), paste0))) {
  for (tp in c("B","Day15","Day91","Day181")) {
    dat_proc[[tp%.%a]] = 10^dat_proc[[tp%.%a]]
  }
}

# add pos call columns for S markers as the OR of S1 and S2
tmp=c(outer(c("B", "Day15", "Day91", "Day181"), S[!grepl("FS", S)], "%.%"))
for (a in tmp) dat_proc[[glue("{a}_resp")]] = pmax(dat_proc[[glue("{a}1_resp")]], dat_proc[[glue("{a}2_resp")]])


# # at 20% threshold, pos calls for BA.4.5.S and COV2.CON.S markers are identical
# COV2.CON.S = S[endsWith(S, "COV2.CON.S")]
# BA.4.5.S = S[endsWith(S, "BA.4.5.S")]
# stopifnot(all(sub("COV2.CON.S","",COV2.CON.S) == sub("BA.4.5.S","",BA.4.5.S))) # make sure ordering is the same
# 
# tmp=c(outer(c("B", "Day15", "Day91", "Day181"), COV2.CON.S, "%.%"))
# pos1 = sapply(tmp%.%"_resp", function(x) mean(dat_proc[[x]], na.rm=T))
# tmp=c(outer(c("B", "Day15", "Day91", "Day181"), BA.4.5.S, "%.%"))
# pos2 = sapply(tmp%.%"_resp", function(x) mean(dat_proc[[x]], na.rm=T))
# all((pos1>=0.2) == (pos2>=0.2))




###############################################################################
# 7. add mdw scores for nAb
{
kp = dat_proc$ph1.D15==1

# myboxplot(dat_proc[kp, c("B"%.%assays[1:5], "Day15"%.%assays[1:5])], names=sub("pseudoneutid50_", "", rep(assays[1:5],2)))
# mypairs(dat_proc["Day15"%.%assays[1:5]])
# mypairs(dat_proc["B"%.%assays[1:5]])
# corplot(Day15pseudoneutid50_D614G~Bpseudoneutid50_D614G, dat_proc)
# sapply(dat_proc[kp, c("B"%.%assays[1:5], "Day15"%.%assays[1:5])], sd)


# use Day15 to derive weights
mdw.wt.nAb=tryCatch({
  tree.weight(cor(dat_proc[kp, "Day15"%.%nAb], use='complete.obs'))
}, error = function(err) {
  print(err$message)
  rep(1/length(nAb), length(nAb))
})
# using Day15 and using B lead to similar weights: BA1 and BA4 BA5 together takes about half of the weight
print(mdw.wt.nAb)
write.csv(mdw.wt.nAb, file = here("data_clean", "csv", TRIAL%.%"_nAb_mdw_weights.csv"))

# apply to all time points
for (t in c("B", "Day15", "Day29", "Day91", "Day181")) {
  dat_proc[, t%.%'pseudoneutid50_MDW'] = as.matrix(dat_proc[, t%.%nAb]) %*% mdw.wt.nAb
}

# pc1
pc1.wt.nAb = princomp(dat_proc[kp, "Day15"%.%nAb])$loadings[,1] # alt for when there are no weights
for (t in c("B", "Day15", "Day29", "Day91", "Day181")) {
  dat_proc[, t%.%'pseudoneutid50_PC1'] = as.matrix(dat_proc[, t%.%nAb]) %*% mdw.wt.nAb
}


}


{# add pc1 and mdw scores for SL. compute score for N and NN separately
# compute scores for all ptids, but weights are based on for Stage 1 + 2 + 3 one-dose mRNA arms.


get.pc1.wt = function(markers) {
  # the first column of markers is weight
  stopifnot(complete.cases(markers))
  cw <- cov.wt(markers[,-1], wt = markers[,1], cor = FALSE)
  eig <- eigen(cw$cov)
  eig$vectors[, 1]
  # princomp(markers[-1])$loadings[,1] # alt for when there are no weights
}

get.mdw.wt = function(markers) {
  # the first column of markers is weight
  stopifnot(complete.cases(markers))
  cw <- cov.wt(markers[,-1], wt = markers[,1], cor = TRUE)
  tree.weight(cw$cor, plot=F)
}


marker_sets = c("fc_spike", "fc_N", "bAb_spike", "bAb_N", "cd4_spike", "cd8_spike")

fc_assays = subset(assay_metadata, panel=="Fc", assay, drop=T)
fc_spike = fc_assays[endsWith(fc_assays, "S")]
fc_N = fc_assays[endsWith(fc_assays, "N")]

bAb_assays = subset(assay_metadata, panel=="bAb", assay, drop=T)
bAb_spike = bAb_assays[endsWith(bAb_assays, "S")]
bAb_N = bAb_assays[endsWith(bAb_assays, "N")]

cd4_spike = c(S,N)[startsWith(S,"cd4")]
cd8_spike = c(S,N)[startsWith(S,"cd8")]
dat.tmp = subset(dat_proc, ph1.D15 & TrtonedosemRNA==1 & !arm %in% c(16,17) & naive==0)
# filter by pos rate among NN
pos1 = sapply("Day15"%.%cd4_spike%.%"_resp", function(x) sum(dat.tmp[[x]] * dat.tmp$ph2.D15.tcell * dat.tmp$wt.D15.tcell, na.rm=T)/sum(dat.tmp$ph1.D15.tcell))
cd4_spike = cd4_spike[pos1>=0.1]
# filter by pos rate among NN
pos1 = sapply("Day15"%.%cd8_spike%.%"_resp", function(x) sum(dat.tmp[[x]] * dat.tmp$ph2.D15.tcell * dat.tmp$wt.D15.tcell, na.rm=T)/sum(dat.tmp$ph1.D15.tcell))
cd8_spike = cd8_spike[pos1>=0.1]

# 8: at B, seven_PC1 plus (i) CD4+ T cells IFNg and/or IL-2 Index N peptides and (ii) CD8+ T cells IFNg and/or IL-2 Index N peptides

seven_PC1 = c('pseudoneutid50',marker_sets) %.% "_PC1"
seven_PC1 = seven_PC1[!contain(seven_PC1, "_N_")]
seven_MDW = c('pseudoneutid50',marker_sets) %.% "_MDW"
seven_MDW = seven_MDW[!contain(seven_MDW, "_N_")]

seven_PC1_N = c(c('pseudoneutid50',marker_sets) %.% "_PC1", "cd4_IFNg.IL2_Wuhan.N")
seven_MDW_N = c(c('pseudoneutid50',marker_sets) %.% "_MDW", "cd4_IFNg.IL2_Wuhan.N")


for (assay_group in c(marker_sets,"seven_PC1","seven_MDW","seven_PC1_N","seven_MDW_N")) {
  # define weights using one dose mRNA arms in stage 1 and 2. define weights separately for naive and non naive
  
  # use B and D15 to compute weight for N and spike markers, respectively
  tp.wt = ifelse(endsWith(assay_group,"_N"), "B", "Day15")
  
  # pc1
  pc1.wts = sapply (0:1, function(.naive) {
    get.pc1.wt(subset(dat_proc, naive==.naive & ph2.D15.xassays==1)[, c("wt.D15.xassays", glue("{tp.wt}{get(assay_group)}"))])
  })
  for (tp in c("B", if (!endsWith(assay_group,"_N")) "Day15")) {
    dat_proc[[glue("{tp}{assay_group}_PC1")]] = ifelse(dat_proc$naive==0,
                                                       c(as.matrix(dat_proc[,glue("{tp}{get(assay_group)}")]) %*% pc1.wts[,1]),
                                                       c(as.matrix(dat_proc[,glue("{tp}{get(assay_group)}")]) %*% pc1.wts[,2]))
  }
  
  # mdw 
  mdw.wts = sapply (0:1, function(.naive) {
    get.mdw.wt(subset(dat_proc, naive==.naive & ph2.D15.xassays==1)[, c("wt.D15.xassays", glue("{tp.wt}{get(assay_group)}"))])
  })
  for (tp in c("B", if (!endsWith(assay_group,"_N")) "Day15")) {
    dat_proc[[glue("{tp}{assay_group}_MDW")]] = ifelse(dat_proc$naive==0,
                                                       c(as.matrix(dat_proc[,glue("{tp}{get(assay_group)}")]) %*% mdw.wts[,1]),
                                                       c(as.matrix(dat_proc[,glue("{tp}{get(assay_group)}")]) %*% mdw.wts[,2]))
  }
  
}
# remove 12 unwanted from 16 that has seven_ and only keep 4: Bseven_PC1_N_PC1, Bseven_MDW_N_MDW, Day15seven_PC1_PC1, Day15seven_MDW_MDW
for (tp in c("B","Day15")) {
  dat_proc[[glue("{tp}seven_PC1_MDW")]] = NULL
  dat_proc[[glue("{tp}seven_MDW_PC1")]] = NULL
  dat_proc[[glue("{tp}seven_PC1_N_MDW")]] = NULL
  dat_proc[[glue("{tp}seven_MDW_N_PC1")]] = NULL
}
dat_proc[[glue("Bseven_PC1_PC1")]] = NULL
dat_proc[[glue("Bseven_MDW_MDW")]] = NULL
dat_proc[[glue("Day15seven_PC1_N_PC1")]] = NULL
dat_proc[[glue("Day15seven_MDW_N_MDW")]] = NULL

}


###############################################################################
# 8. add fold change markers
# assuming data has been censored at the lower limit
# thus no need to do, say, lloq censoring
# but there is a need to do uloq censoring before computing delta
{
  

# mdw scores delta are computed as weighted average of delta, not as delta of mdw
# include S1, S2, S, N tcell markers
assays1=union(assays, tcellvv) 
tmp=list()
for (a in assays1) {
  for (t in c("B", paste0("Day", config$timepoints)) ) {
    tmp[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] > log10(uloqs[marker.name.to.assay(t %.% a)]), log10(uloqs[marker.name.to.assay(t %.% a)]), dat_proc[[t %.% a]])
  }
}
tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame
for (tp in rev(timepoints)) {
  dat_proc["Delta"%.%tp%.%"overB" %.% assays1] <- tmp["Day"%.%tp %.% assays1] - tmp["B" %.% assays1]
}   
if(two_marker_timepoints) {
  dat_proc["Delta"%.%timepoints[2]%.%"over"%.%timepoints[1] %.% assays1] <- tmp["Day"%.% timepoints[2]%.% assays1] - tmp["Day"%.%timepoints[1] %.% assays1]
}


# also need D29 delta for sanofi arms
assays1=c(setdiff(nAb, "pseudoneutid50Duke_BA.2.12.1"), "pseudoneutid50_MDW")
tmp=list()
for (a in assays1) {
  for (t in c("B", "Day29") ) {
    tmp[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat_proc[[t %.% a]])
  }
}
tmp=as.data.frame(tmp) # cannot subtract list from list, but can subtract data frame from data frame
dat_proc["Delta29overB" %.% assays1] <- tmp["Day29" %.% assays1] - tmp["B" %.% assays1]
}


###############################################################################
# 9. add discrete/trichotomized markers
{

## ID50
  
# mRNA arms, D15, B, fold change
dat_proc$tmp = with(dat_proc, ph1.D15 & TrtonedosemRNA==1) 
# neut
assays1 = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
all.markers = c("B"%.%assays1, "Day15"%.%assays1, "Delta15overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D15", verbose=T)

# Sanofi arms, Day29, fold change from B
dat_proc$tmp = with(dat_proc, ph1.D29 & TrtSanofi==1)
# neut
assays1 = c(nAb, "pseudoneutid50_MDW")
all.markers = c("Day29"%.%assays1, "Delta29overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D29", verbose=F)


## frnt, both mRNA and Sanofi, D15, B, fold change
dat_proc$tmp = with(dat_proc, ph2.D15.frnt) 
assays1 = frnt
all.markers = c("B"%.%assays1, "Day15"%.%assays1, "Delta15overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D15.frnt", verbose=T)

# no day29 frnt data yet
# # frnt, Sanofi, Day29, fold change from B
# dat_proc$tmp = with(dat_proc, ph2.D15.frnt & TrtSanofi==1)
# assays1 = frnt
# all.markers = c("Day29"%.%assays1, "Delta29overB"%.%assays1)
# dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D29", verbose=F)


## T cell markers for both mRNA and Sanofi, D15, B, fold change
dat_proc$tmp = with(dat_proc, ph1.D15) 
assays1=c(S, S1, S2, N)
all.markers = c("B"%.%assays1, "Day15"%.%assays1, "Delta15overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D15.tcell", verbose=T)


## Seromyx markers for both mRNA and Sanofi, D15, B, fold change
dat_proc$tmp = with(dat_proc, ph2.D15.seromyx) 
assays1=seromyx
all.markers = c("B"%.%assays1, "Day15"%.%assays1, "Delta15overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D15.seromyx", verbose=T)


# remove the temp ph2 column
dat_proc$tmp = NULL


# subset on subset_variable
if(!is.null(config$subset_variable) & !is.null(config$subset_value)){
  if(subset_value != "All") {
    include_in_subset <- dat_proc[[subset_variable]] == subset_value
    dat_proc <- dat_proc[include_in_subset, , drop = FALSE]
  }
}
}

###############################################################################
# 10. impute covariates if necessary



###############################################################################
# 11. special handling 


#### create an indicator SubcohortInd.casedeletion for immunogenicity studies

# loop through naive/nnaive, arm
todelete=c()
for (i in 1:2) {
  if (i==1) dat=subset(dat_proc, naive==1 & ph1.D15.tcell) # Naive
  if (i==2) dat=subset(dat_proc, naive==0 & ph1.D15.tcell) # NN
  
  dat$casep = dat$case.period
  dat$casep[is.na(dat$case.period)] = 3
  dat$casep = factor(dat$casep, c(1,2,3))
  
  arms=unique(dat$arm)
  for (a in arms) {
    dat.1=subset(dat, arm==a)
    
    tab = mytable(dat.1$ph2.D15.tcell, dat.1$casep); tab
    # sampling ratios among booster proximal cases, booster distal cases, controls; could be NaN, if so, set to 0
    p.case.prox = tab[2,"1"]/sum(tab[,"1"]); if (is.nan(p.case.prox)) p.case.prox=0
    p.case.dist = tab[2,"2"]/sum(tab[,"2"]); if (is.nan(p.case.dist)) p.case.dist=0
    p.ctrl      = tab[2,"3"]/sum(tab[,"3"])
    c(p.case.prox, p.case.dist, p.ctrl)
    
    # delete proximal cases 
    if (p.case.prox>p.ctrl) {
      n=round(sum(tab[,"1"]) * p.ctrl)
      Num=tab[2,"1"]
      if(n<Num) {
        ptids = subset(dat.1, casep==1 & ph2.D15.tcell, Ptid, drop=T)
        stopifnot(len(ptids)==Num)
        todelete = c(todelete, sample(ptids, Num-n))
      }
    }
    
    # delete distal cases 
    if (p.case.dist>p.ctrl) {
      n=round(sum(tab[,"2"]) * p.ctrl)
      Num=tab[2,"2"]
      if(n<Num) {
        ptids = subset(dat.1, casep==2 & ph2.D15.tcell, Ptid, drop=T)
        stopifnot(len(ptids)==Num)
        todelete = c(todelete, sample(ptids, Num-n))
      }
    }
  }
  print(len(todelete))
}

dat_proc$SubcohortInd.casedeletion = dat_proc$ph2.D15.tcell # T/F
dat_proc$SubcohortInd.casedeletion[dat_proc$Ptid %in% todelete] = F

# validation
dat.tmp=dat_proc[dat_proc$ph1.D15.tcell==1,]
mytable(dat.tmp$ph2.D15.tcell, dat.tmp$COVIDIndD22toD181, dat.tmp$naive)
mytable(dat.tmp$SubcohortInd.casedeletion, dat.tmp$COVIDIndD22toD181, dat.tmp$naive)



#### impute FRNT50 and FRNT80 at B and D15 for the 3 ptids in ph2.D15.xassays but not in ph2.D15.frnt

# impute FRNT50, FRNT80 at B and D15 together.
# impute different arms and naive/nnaive together, but pass arm and naive as covariates

n.imp=1
tp=15

imp.markers=c(outer(c("B", "Day"%.%tp), c(frnt, nAb), "%.%"))

# add arm and naive to the imputation dataset
imp.markers =  c(imp.markers, "arm.factor", "naive")

dat.tmp.impute <- subset(dat_proc, get("ph2.D15.xassays") == 1)

imp <- dat.tmp.impute %>% select(all_of(imp.markers))     

# there is one -Inf in a FS
imp[imp==-Inf]=NA

if(any(is.na(imp))) {
  # diagnostics = FALSE , remove_collinear=F are needed to avoid errors due to collinearity
  imp <- imp %>% mice(m = n.imp, printFlag = FALSE, seed=1, diagnostics = FALSE , remove_collinear = FALSE)            
  dat.tmp.impute[, imp.markers] <- mice::complete(imp, action = 1)
}                

# missing markers imputed properly?
assertthat::assert_that(
  all(complete.cases(dat.tmp.impute[, imp.markers])),
  msg = "missing markers imputed properly?"
)    

# populate dat_proc imp.markers with the imputed values
dat_proc[dat_proc[["ph2.D15.xassays"]]==1, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["ph2.D15.xassays"]]==1, "Ptid"], dat.tmp.impute$Ptid), ]

assertthat::assert_that(
  all(complete.cases(dat_proc[dat_proc[["ph2.D15.xassays"]] == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)
  
# add the corresponding fold change markers
mytable(dat_proc$Delta15overBfrnt50_BA.1 == dat_proc$Day15frnt50_BA.1 - dat_proc$Bfrnt50_BA.1)
dat_proc$Delta15overBfrnt50_BA.1 = dat_proc$Day15frnt50_BA.1 - dat_proc$Bfrnt50_BA.1
mytable(dat_proc$Delta15overBfrnt50_BA.1 == dat_proc$Day15frnt50_BA.1 - dat_proc$Bfrnt50_BA.1)
# the table after the operation sees 3 more values of delta


#### add NLPC1 made from fsdam

tmp=read.csv("../covail_xassays_NLPC1.csv")
dat_proc$Day15seven_PC1_NLPC1 = tmp$Day15seven_PC1_NLPC1 [match(dat_proc$Ptid, tmp$Ptid)]
dat_proc$Bseven_PC1_N_NLPC1   = tmp$Bseven_PC1_N_NLPC1   [match(dat_proc$Ptid, tmp$Ptid)]



#### impute T cell _vacresp markers

n.imp=1
tp=15

imp.markers=names(dat_proc)
imp.markers = imp.markers[contain(imp.markers,"_vacresp")]
print(imp.markers)

dat.tmp.impute <- subset(dat_proc, get("TwophasesampIndD15.tcell") == 1)

for (a in imp.markers) {
  kp = is.na(dat.tmp.impute[[a]])
  if (any(kp)) {
    amarker = sub('_vacresp','',a)
    amarker = sub("Day15","",amarker)
    f=as.formula(glue("{a} ~ B{amarker} + Day15{amarker} + arm.factor + naive"))
    fit = glm(f, dat.tmp.impute, family=binomial)
    dat.tmp.impute[[a]][kp] = rbinom(sum(kp), 1, prob=predict(fit, dat.tmp.impute[kp,], type="response"))
  } else next
}  

# missing markers imputed properly?
assertthat::assert_that(
  all(complete.cases(dat.tmp.impute[, imp.markers])),
  msg = "missing markers imputed properly?"
)    

# populate dat_proc imp.markers with the imputed values
dat_proc[dat_proc[["TwophasesampIndD15.tcell"]]==1, imp.markers] <-
  dat.tmp.impute[imp.markers][match(dat_proc[dat_proc[["TwophasesampIndD15.tcell"]]==1, "Ptid"], dat.tmp.impute$Ptid), ]

assertthat::assert_that(
  all(complete.cases(dat_proc[dat_proc[["TwophasesampIndD15.tcell"]] == 1, imp.markers])),
  msg = "imputed values of missing markers merged properly for all individuals in the two phase sample?"
)


# correct FRNT values???

# vv=names(dat_proc)
# vv=vv[contain(vv, "frnt") & contain(vv, "_")]
# vv=vv[!endsWith(vv,"cat") & !startsWith(vv,"Delta")]
# vv=setdiff(vv, "Day181frnt_missing_cnt")
# 
# sapply (vv, function(a) sum(dat_proc[[a]]<log10(20) & dat_proc[[a]]>1, na.rm=T))
# sapply (vv, function(a) min(dat_proc[[a]], na.rm=T))


###############################################################################
# digest check

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
    tmp = switch(TRIAL,
         covail = "1b78c30d75ae006b327846ef7f588a71",
         NA)    
    if (!is.na(tmp)) assertthat::validate_that(digest(dat_proc[order(names(dat_proc))])==tmp, 
      msg = "--------------- WARNING: failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))])%.%' ----------------')    
}

data_name = paste0(TRIAL, "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")

if (!dir.exists("data_clean/csv")) dir.create("data_clean/csv")

write_csv(dat_proc, file = here("data_clean", "csv", data_name))



print("run time: "%.%format(Sys.time()-begin, digits=1))
