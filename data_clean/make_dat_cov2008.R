Sys.setenv(TRIAL = "cov2008")
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
  # read mapped data
  dat_raw = read.csv(mapped_data)
  dat_proc = preprocess(dat_raw, study_name)   
  names(dat_proc)[[1]]="Ptid"
  
  assay_metadata=read.csv(config$assay_metadata)
  assays=assay_metadata$assay
  
}


########################################################################################################
# define Senior and race/ethnicity
{
  colnames(dat_proc)[colnames(dat_proc)=="Subjectid"] <- "Ptid" 
  dat_proc <- dat_proc %>% mutate(age.geq.65 = as.integer(Age >= 65))
  dat_proc$Senior = as.integer(dat_proc$Age>=switch(study_name, COVE=65, MockCOVE=65, ENSEMBLE=60, MockENSEMBLE=60, PREVENT19=65, AZD1222=65, VAT08=60, PROFISCOV=NA, COV2008=65, NVX_UK302=65, stop("unknown study_name 1")))
  
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

dat_proc <- dat_proc %>% mutate(tps.stratum = Trt)
if (!is.null(dat_proc$tps.stratum)) table(dat_proc$tps.stratum)

dat_proc$Wstratum.tcell = dat_proc$tps.stratum
if (!is.null(dat_proc$Wstratum.tcell)) table(dat_proc$Wstratum.tcell) # variables may be named other than Wstratum
}


################################################################################
# 4. Define ph1, ph2, and weights
# Note that Wstratum may have NA if any variables to form strata has NA
{
dat_proc$ph1.D15 = ifelse(dat_proc$PPEFL=='Y',1,0)
dat_proc$ph1.D15.tcell = dat_proc$ph1.D15

dat_proc$ph2.D15.tcell=ifelse(!is.na(dat_proc$BTerminally_Diff_CD4_Any_Cov2_S_IFNg_OR_IL2) & !is.na(dat_proc$Day15Terminally_Diff_CD4_Any_Cov2_S_IFNg_OR_IL2), 1, 0)

tp=15
dat_proc = add.wt(dat_proc, ph1="ph1.D"%.%tp%.%".tcell", ph2="ph2.D"%.%tp%.%".tcell", Wstratum="Wstratum.tcell", 
                  wt="wt.D"%.%tp%.%".tcell", verbose=T) 


}


###############################################################################
# 5. impute missing biomarkers in ph2 (assay imputation)





###############################################################################
# 6. transformation of the markers





###############################################################################
# 7. add mdw scores for nAb






###############################################################################
# 8. add fold change markers
# assuming data has been censored at the lower limit
# thus no need to do, say, lloq censoring
# but there is a need to do uloq censoring before computing delta

{
# mdw scores delta are computed as weighted average of delta, not as delta of mdw
# include S1, S2, S, N tcell markers
assays1=assays 
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


}


###############################################################################
# 9. add discrete/trichotomized markers
{
dat_proc$tmp = with(dat_proc, ph2.D15.tcell) 
assays1 = assays
all.markers = c("B"%.%assays1, "Day15"%.%assays1, "Delta15overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D15.tcell", verbose=T)


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




###############################################################################
# digest check

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
    tmp = switch(TRIAL,
         cov2008 = "441aba2840ca1415583d29bae9718649",
         NA)    
    if (!is.na(tmp)) assertthat::validate_that(digest(dat_proc[order(names(dat_proc))])==tmp, 
      msg = "--------------- WARNING: failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))])%.%' ----------------')    
}

data_name = paste0(TRIAL, "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")

if (!dir.exists("data_clean/csv")) dir.create("data_clean/csv")

write_csv(dat_proc, file = here("data_clean", "csv", data_name))



print("run time: "%.%format(Sys.time()-begin, digits=1))
