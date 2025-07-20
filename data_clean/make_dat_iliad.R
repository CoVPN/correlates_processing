Sys.setenv(TRIAL = "iliad_ib202p")
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
# 1. read mapped data 

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
{}




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

dat_proc$Wstratum = dat_proc$tps.stratum
}


################################################################################
# 4. Define ph1, ph2, and weights
# Note that Wstratum may have NA if any variables to form strata has NA
{
# the whole cohort is treated as ph1 and ph2
dat_proc$TwophasesampIndC0 <- dat_proc$ph1.C0 <- 1

dat_proc[["ph2.C0"]]=dat_proc$ph1.C0
dat_proc[["wt.C0"]] = 1

}


###############################################################################
# 5. impute missing biomarkers in ph2 (assay imputation)




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
  
# mRNA arms, C0, B, fold change
dat_proc$tmp = with(dat_proc, ph1.C0 & TrtonedosemRNA==1) 
# neut
assays1 = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "pseudoneutid50_Beta", "pseudoneutid50_BA.1", "pseudoneutid50_BA.4.BA.5", "pseudoneutid50_MDW")
all.markers = c("B"%.%assays1, "Day15"%.%assays1, "Delta15overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.C0", verbose=T)

# Sanofi arms, Day29, fold change from B
dat_proc$tmp = with(dat_proc, ph1.D29 & TrtSanofi==1)
# neut
assays1 = c(nAb, "pseudoneutid50_MDW")
all.markers = c("Day29"%.%assays1, "Delta29overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D29", verbose=F)


## frnt, both mRNA and Sanofi, C0, B, fold change
dat_proc$tmp = with(dat_proc, ph2.C0.frnt) 
assays1 = frnt
all.markers = c("B"%.%assays1, "Day15"%.%assays1, "Delta15overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.C0.frnt", verbose=T)

# no day29 frnt data yet
# # frnt, Sanofi, Day29, fold change from B
# dat_proc$tmp = with(dat_proc, ph2.C0.frnt & TrtSanofi==1)
# assays1 = frnt
# all.markers = c("Day29"%.%assays1, "Delta29overB"%.%assays1)
# dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.D29", verbose=F)


## T cell markers for both mRNA and Sanofi, C0, B, fold change
dat_proc$tmp = with(dat_proc, ph1.C0) 
assays1=c(S, S1, S2, N)
all.markers = c("B"%.%assays1, "Day15"%.%assays1, "Delta15overB"%.%assays1)
dat_proc = add.trichotomized.markers (dat_proc, all.markers, ph2.col.name="tmp", wt.col.name="wt.C0.tcell", verbose=T)


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
# special handling 

# create an indicator SubcohortInd.casedeletion for immunogenicity studies

# loop through naive/nnaive, arm
todelete=c()
for (i in 1:2) {
  if (i==1) dat=subset(dat_proc, naive==1 & ph1.C0.tcell) # Naive
  if (i==2) dat=subset(dat_proc, naive==0 & ph1.C0.tcell) # NN
  
  dat$casep = dat$case.period
  dat$casep[is.na(dat$case.period)] = 3
  dat$casep = factor(dat$casep, c(1,2,3))
  
  arms=unique(dat$arm)
  for (a in arms) {
    dat.1=subset(dat, arm==a)
    
    tab = mytable(dat.1$ph2.C0.tcell, dat.1$casep); tab
    # sampling ratios among booster proximal cases, booster distal cases, controls; could be NaN, if so, set to 0
    p.case.prox = tab[2,"1"]/sum(tab[,"1"]); if (is.nan(p.case.prox)) p.case.prox=0
    p.case.dist = tab[2,"2"]/sum(tab[,"2"]); if (is.nan(p.case.dist)) p.case.dist=0
    p.ctrl      = tab[2,"3"]/sum(tab[,"3"])
    c(p.case.prox, p.case.dist, p.ctrl)
    
    # delete proximal cases 
    if (p.case.prox>p.ctrl) {
      n=round(sum(tab[,"1"]) * p.ctrl)
      N=tab[2,"1"]
      if(n<N) {
        ptids = subset(dat.1, casep==1 & ph2.C0.tcell, Ptid, drop=T)
        stopifnot(len(ptids)==N)
        todelete = c(todelete, sample(ptids, N-n))
      }
    }
    
    # delete distal cases 
    if (p.case.dist>p.ctrl) {
      n=round(sum(tab[,"2"]) * p.ctrl)
      N=tab[2,"2"]
      if(n<N) {
        ptids = subset(dat.1, casep==2 & ph2.C0.tcell, Ptid, drop=T)
        stopifnot(len(ptids)==N)
        todelete = c(todelete, sample(ptids, N-n))
      }
    }
  }
  print(len(todelete))
}

dat_proc$SubcohortInd.casedeletion = dat_proc$ph2.C0.tcell # T/F
dat_proc$SubcohortInd.casedeletion[dat_proc$Ptid %in% todelete] = F

# validation
dat.tmp=dat_proc[dat_proc$ph1.C0.tcell==1,]
mytable(dat.tmp$ph2.C0.tcell, dat.tmp$COVIDIndD22toD181, dat.tmp$naive)
mytable(dat.tmp$SubcohortInd.casedeletion, dat.tmp$COVIDIndD22toD181, dat.tmp$naive)


###############################################################################
# digest check

library(digest)
if(Sys.getenv ("NOCHECK")=="") {    
    tmp = switch(TRIAL,
         covail = "cc1fc791c8e384ba9b1c96790da26930",
         NA)    
    if (!is.na(tmp)) assertthat::validate_that(digest(dat_proc[order(names(dat_proc))])==tmp, 
      msg = "--------------- WARNING: failed make_dat_proc digest check. new digest "%.%digest(dat_proc[order(names(dat_proc))])%.%' ----------------')    
}

data_name = paste0(TRIAL, "_data_processed_", format(Sys.Date(), "%Y%m%d"), ".csv")

if (!dir.exists("data_clean/csv")) dir.create("data_clean/csv")

write_csv(dat_proc, file = here("data_clean", "csv", data_name))



print("run time: "%.%format(Sys.time()-begin, digits=1))
