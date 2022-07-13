###### quantifying uncertainty in WDNR walleye PEs (don't really have this for our bass PEs) & Fishscapes lmb efCPUE (don't really have this for WDNR walleye)
###### 2022-07-12
###### cd and sej

##########################
#### WDNR walleye PEs ####
##########################
# added by stuart to load file from Colin
dat=read.csv("~/Documents/Research/People/Students/current/Mosley_Camille/hyperstabilitySurveyManuscript/FisheriesResearch_revisions3/PEuncertaintyFromColin/WAEpopest2021.csv",header=TRUE) 

library(tidyverse)
#dat=WAEpopest2021[-1,]

uncert=dat%>%
  select(MWBCODE,REFDATE,NUMBER,NUMHIGH,NUMLOW,NUMACRE,NUMAHIGH,NUMALOW,CVAPE,L95APE,CVTPE)

##### code added by Stuart
# assuming PE's are lognormally distributed - calculate SD from CV then back out sigma and mu
uncert=uncert[-1,]
uncert$SD=uncert$NUMBER*uncert$CVAPE
uncert$VAR=uncert$SD^2
sigma=sqrt(log((uncert$VAR/uncert$NUMBER^2)+1))
mu=log(uncert$NUMBER)-(sigma^2)/2

plot(mu,sigma)
summary(lm(sigma~mu))  # sigma weakly decreases with mean, so assuming this doesn't matter

mean(sigma,na.rm=TRUE) #0.1451

###############################
#### fishscapes lmb efCPUE ####
###############################
library(MFEUtilities)
dbdir="~/Documents/Research/MFE/database"
db="MFEdb.db"


# load fish samples and fish info
fs=dbTable("FISH_SAMPLES")
fi=dbTable("FISH_INFO")

# lakes we have PEs for
lakes=c("HT","BY","DY","SE","SV","UG","WB","AR","BA","FD","JS","LC","LH","LR","TO","WC","WS")

# subset fish samples to fishscapes 1/2 mile electrofishing
fsfs=fs[fs$projectID==37,]
fsfsdnr=fsfs[fsfs$metadataID%in%c("FishscapesSurvey.0.5mile.20180606","FishscapesSurvey.1.5mile.20180607"),]
fsfsdnrpe=fsfsdnr[fsfsdnr$lakeID%in%lakes,]

# create vector to store mean and sd of log efCPE for each lake
efCPE=numeric(length(lakes))
SDefCPE=numeric(length(lakes))
for(i in 1:length(lakes)){
  cur=fsfsdnrpe[fsfsdnrpe$lakeID==lakes[i],] # pull only samples from the ith lake
  cur=cur[cur$useCPUE=="yes",]
  
  cur_efCPE=numeric(nrow(cur))     # create vector to store individual trip (sample) efCPE for the  ith lake
  for(j in 1:nrow(cur)){
    curFI=fi[fi$sampleID==cur$sampleID[j],]
    cur_efCPE[j]=sum(curFI$otu=="largemouth_bass")/as.numeric(cur$distanceShocked[j])
  }
  
  SDefCPE[i]=sd(log(cur_efCPE))
  efCPE[i]=mean(log(cur_efCPE))
}

names(efCPE)=lakes
names(SDefCPE)=lakes
# remove SE and BA - no DNR-style  electrofishing; WS only  caught bass on one  ef trip; WB only had  one EF sample
SDefCPE=SDefCPE[is.finite(efCPE)]
efCPE=efCPE[is.finite(efCPE)]

efCPE=efCPE[is.finite(SDefCPE)]
SDefCPE=SDefCPE[is.finite(SDefCPE)]

plot(efCPE,SDefCPE) # no scaling with mean
summary(lm(SDefCPE~efCPE))

mean(SDefCPE)  # 0.57