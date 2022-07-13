#Cleaned up hyperstability survey code 
#looking at hyperstability of angling CPUE as a function of electrofishing CPUE
#9-24-2021
# CLM, CJD, SEJ, CS
rm(list=ls())

# load function to load data from google drive
source("gdriveURL.R")
library(dplyr)

### Data Wrangle ####

######## angling CPUE
# load creel data from google drive
creel1=gdriveURL("https://drive.google.com/open?id=1lxUd742QZMXDQunyFBnENKMYZ1XNM_Pc")
creel2=gdriveURL("https://drive.google.com/open?id=1UYhbGH28WXjmi-4BzhfwO4KYwrBCNO2Q")
creel=rbind(creel1,creel2)

# reduce to columns we care about
creel=creel[,c(1,3,6,12,18,25:26,30,36,37,38)]

# calculate effort
# add zeroes to times with only 2 or 3 digits
creel$timeStart[nchar(creel$timeStart)==3]=paste("0",creel$timeStart[nchar(creel$timeStart)==3],sep="")
creel$timeStart[nchar(creel$timeStart)==2]=paste("00",creel$timeStart[nchar(creel$timeStart)==2],sep="")
creel$timeStart[creel$timeStart=="0"]="0000"
creel=creel[creel$timeStart!="1",]  # 4 entries with "1", so we don't know start time

creel$timeEnd[nchar(creel$timeEnd)==3]=paste("0",creel$timeEnd[nchar(creel$timeEnd)==3],sep="")
creel$timeEnd[nchar(creel$timeStart)==2]=paste("00",creel$timeEnd[nchar(creel$timeEnd)==2],sep="")
creel$timeEnd[creel$timeEnd=="0"]="0000"

creel$boatHrs=0
# remove rows when end time is less than start time (assumes the boat was out over midnight)
creel=creel[strptime(creel$timeEnd,format="%H%M")>=strptime(creel$timeStart,format="%H%M"),]
# calculate difference of time in hours for rows where end time is greater than start time (fishing occurred in one day only)
creel$boatHrs[strptime(creel$timeEnd,format="%H%M")>=strptime(creel$timeStart,format="%H%M")]=as.numeric(difftime(strptime(creel$timeEnd[strptime(creel$timeEnd,format="%H%M")>=strptime(creel$timeStart,format="%H%M")],format="%H%M"),strptime(creel$timeStart[strptime(creel$timeEnd,format="%H%M")>=strptime(creel$timeStart,format="%H%M")],format="%H%M"),units="hours"))

# removing rows with a non-zero notFishingAmt because we don't know what it means to be non-zero...
creel=creel[creel$notFishingAmt==0,]

# remove rows with non-integer anglersAmt 
creel=creel[!grepl(".",creel$anglersAmt,fixed=TRUE),]

# remove rows with anglersAmt above 10? (arbitrary choice for now)
creel=creel[creel$anglersAmt<=10,]

# get angler hours of effort from party size and boat hours
creel$anglerHrs=creel$boatHrs*creel$anglersAmt

# remove rows with no species code
creel=creel[!is.na(creel$fishSpeciesCode),]

# remove rows with NA for caughtAmt
creel=creel[!is.na(creel$caughtAmt),]

# remove no effort (anglerHrs==0) rows
creel=creel[creel$anglerHrs>0,]

# calculate angling CPUE
creel$anglingCPUE=creel$caughtAmt/creel$anglerHrs

# removing instances of CPUE >=30 (this seems impossible...)
creel=creel[creel$anglingCPUE<30,]

# remove ice fishing... -> results are robust to inclusion or removal of ice fishing
#creel$DOY=as.numeric(strftime(strptime(creel$dateSample,format="%d-%b-%y"),format="%j"))
#creel=creel[creel$DOY<320,]
#creel=creel[creel$DOY>120,]

write.table(creel,"cleanedWDNRcreel.csv",sep=",",row.names=FALSE)

####### electrofishing abundance
bassEF=gdriveURL("https://drive.google.com/open?id=11v8FbT2wnKx_CqUfxu_V9r_8fyCfcdD2")
bassEF=bassEF[,c(1,3,5,13,19,27:29)]
bassEF$CPEkm=bassEF$CPEmile/1.60934   # convert fish per mile to fish per km
bassEF$distanceShockedKm=bassEF$distanceShockedMiles*0.621371 # convert miles to km
lake_yearBASSef= bassEF %>%
  group_by(WBIC,species,surveyYear,county) %>%
  summarize(meanEF_CPEkm=mean(CPEkm),
            totalFishCaught=sum(totalNumberCaughtFish),
            totalDistShockedKm=sum(distanceShockedKm),
            totalHoursSampled=sum(numberHoursSampled),
            std=sd(CPEkm),
            N=n())
lake_yearBASSef=as.data.frame(lake_yearBASSef)


panEF=gdriveURL("https://drive.google.com/open?id=1QIqCBQ9gbOgRFUJQbnokwwTZJi5VZZIR")
panEF=panEF[,c(1,3,5,13,19,25:27)]
panEF$CPEkm=panEF$CPEmile/1.60934   # convert fish per mile to fish per km
panEF$distanceShockedKm=panEF$distanceShockedMiles*0.621371 # convert miles to km
lake_yearPANef= panEF %>%
  group_by(WBIC,species,surveyYear,county) %>%
  summarize(meanEF_CPEkm=mean(CPEkm),
            totalFishCaught=sum(totalNumberCaughtFish),
            totalDistShockedKm=sum(distanceShockedKm),
            totalHoursSampled=sum(numberHoursSampled),
            std=sd(CPEkm),
            N=n())
lake_yearPANef=as.data.frame(lake_yearPANef)

walleyeEF=gdriveURL("https://drive.google.com/open?id=1DPRROWv6Cf_fP6Z-kE9ZgUfdf_F_jSNT")
walleyeEF=walleyeEF[,c(1,3,5,13,19,23:24,27)]
walleyeEF$CPEkm=walleyeEF$CPEmile/1.60934   # convert fish per mile to fish per km
walleyeEF$distanceShockedKm=walleyeEF$distanceShockedMiles*0.621371 # convert miles to km
#remove commas from total fish caught
walleyeEF$totalNumberCaughtFish=as.numeric(gsub(",","",walleyeEF$totalNumberCaughtFish))
lake_yearWALLef= walleyeEF %>%
  group_by(WBIC,species,surveyYear,county) %>%
  summarize(meanEF_CPEkm=mean(CPEkm),
            totalFishCaught=sum(totalNumberCaughtFish),
            totalDistShockedKm=sum(distanceShockedKm),
            totalHoursSampled=sum(numberHoursSampled),
            std=sd(CPEkm),
            N=n())
lake_yearWALLef=as.data.frame(lake_yearWALLef)

write.table(lake_yearBASSef,"cleanedWDNRbassEF.csv",sep=",",row.names=FALSE)
write.table(lake_yearPANef,"cleanedWDNRpanEF.csv",sep=",",row.names=FALSE)
write.table(lake_yearWALLef,"cleanedWDNRwalleyeEF.csv",sep=",",row.names=FALSE)

##********* lake characteristics using Hansen et al. dataset
# compile lake information from Gretchen Hansen
hansenFiles=list.files("~/Documents/Research/People/Students/current/Mosley_Camille/codeHelp/hyperstability_analyses/HansenData/")   # the 1st and last files have year-specific data; rest have a single observtion for 14,364 lakes
i=4
cur=read.csv(paste("~/Documents/Research/People/Students/current/Mosley_Camille/codeHelp/hyperstability_analyses/HansenData/",hansenFiles[i],sep=""),header=TRUE)
name=gsub("Hansen_WI_","",hansenFiles[i])
name=gsub(".csv","",name)
colnames(cur)[2:ncol(cur)]=paste(name,colnames(cur)[2:ncol(cur)],sep="_")
hansen=cur
hansen=hansen[order(hansen$wbic),]

for(i in c(2:3,5,7:9)){
  cur=read.csv(paste("~/Documents/Research/People/Students/current/Mosley_Camille/codeHelp/hyperstability_analyses/HansenData/",hansenFiles[i],sep=""),header=TRUE)
  name=gsub("Hansen_WI_","",hansenFiles[i])
  name=gsub(".csv","",name)
  cur=cur[order(cur$wbic),]
  print(sum(cur$wbic==hansen$wbic))
  colnames(cur)=paste(name,colnames(cur),sep="_")
  hansen=cbind(hansen,cur[,2:ncol(cur)])
}

dim(hansen)
colnames(hansen)[22]="lakeOrder"
colnames(hansen)[1]="WBIC"
colnames(hansen)

write.table(hansen,"cleanedHansenetalLakeinfo.csv",row.names=FALSE,sep=",")

####  efCPUE vs PEs

##### walleye
#### looking at WDNR walleye PE vs efCPUE ####

#join by WBIC, surveyYear
library(stringr)
setwd("~/Documents/Research/People/Students/current/Mosley_Camille/hyperstabilitySurveyManuscript/2021-08-23")
WallPE=read.csv("WalleyeData_Age0toAge1MortalitywithPE_Zebro_9_18_20.csv",skip=2)
WallPE=WallPE[,c(2,6,10:11)]
colnames(WallPE)=c('WBIC','surveyYear','PE','Density')

WallPE_2018=read.csv("WalleyeWDNR_18-20_PEs.csv") #REFDATE = Survey year, NumberAcre is the density of adult walleye
WallPE_2018=WallPE_2018[,c(1,6,15,8)]
colnames(WallPE_2018)=c('WBIC','surveyYear','PE','Density')

WallPE=full_join(WallPE,WallPE_2018)

write.table(WallPE,"cleanedWDNRwalleyePE.csv",sep=",",row.names=FALSE)

##### bass
library(MFEUtilities)
#change this path to where you stored the current database on your own computer
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

# create vector to store average efCPE for each lake
efCPE=numeric(length(lakes))
for(i in 1:length(lakes)){
  cur=fsfsdnrpe[fsfsdnrpe$lakeID==lakes[i],] # pull only samples from the ith lake
  cur=cur[cur$useCPUE=="yes",]
  
  cur_efCPE=numeric(nrow(cur))     # create vector to store individual trip (sample) efCPE for the  ith lake
  for(j in 1:nrow(cur)){
    curFI=fi[fi$sampleID==cur$sampleID[j],]
    cur_efCPE[j]=sum(curFI$otu=="largemouth_bass")/as.numeric(cur$distanceShocked[j])
  }
  
  efCPE[i]=mean(cur_efCPE)
}

names(efCPE)=lakes
names(efCPEsd)=lakes
names(efCPElog_sd)=lakes

# remove SE and BA - no DNR-style  electrofishing
efCPE=efCPE[!is.na(efCPE)]
efCPEsd=efCPEsd[!is.na(efCPE)]
efCPElog_sd=efCPElog_sd[!is.na(efCPE)]

bassPE2018=read.csv("~/Documents/Research/Fishscapes/hyperstability/hsSurvey/fishscapes2018_peSum_20180914.csv",row.names=1)
bassPE2019=read.csv("~/Documents/Research/Fishscapes/hyperstability/hsSurvey/2019PEs.csv")
bassPE2019=bassPE2019[1:7,]

li=dbTable("LAKES")

li2019=li[li$lakeID%in%bassPE2019$lakeID,]

bassPE2019$WBIC=li2019$WBIC[match(bassPE2019$lakeID,li2019$lakeID)]

bassPE2019$shoreLength=hansen$lakeChar_Length[match(bassPE2019$WBIC,hansen$WBIC)]/1000

bassPE2019$fishPerKmShoreline=bassPE2019$nHat/bassPE2019$shoreLength

bassPE=data.frame(lakeID=c(bassPE2018$lakeID,bassPE2019$lakeID),
                  PE=c(bassPE2018$nHat,bassPE2019$nHat),
                  densityShoreline=c(bassPE2018$fishPerKmShoreline,bassPE2019$fishPerKmShoreline))
   
bassPE=bassPE[bassPE$lakeID%in%names(efCPE),]

bassPE=bassPE[order(bassPE$lakeID),]
efCPE=efCPE[order(names(efCPE))]
sum(bassPE$lakeID==names(efCPE))

bassPE$efCPE=efCPE


plot(bassPE$densityShoreline,bassPE$efCPE)

bassEFvPEfit=lm(bassPE$efCPE~bassPE$densityShoreline)
summary(bassEFvPEfit)
text(bassPE$densityShoreline,bassPE$efCPE,bassPE$lakeID)


# Day Lake and  White Sand Lake have super  low efCPE for their PE
# I think this  is because of super  low  conductivity and poor  performance
# of  electrofishing in those  two lakes

##****** look at conductivity!!!
biocomplex=read.csv("biocomplexity_spC.csv")
spc=tapply(biocomplex$SpC,biocomplex$lakeName,FUN=mean,na.rm=TRUE)

lakes=dbTable("LAKES")
PElakes=lakes[lakes$lakeID%in%bassPE$lakeID,]

PElakes[,1:2]

# conductivity  from driver logs
bassPE$spc=c(115,191,9.4,53.5,24,181,119,89,149.5,236,49.7,115,80,170,7.4)


# SE and BA removed already due to no DNR-style EF to generate efCPUE estimate
# removing WS (0.06), DY (0.1) for SPC  < 10
bassPEfinal=bassPE[-c(3,15),]

li_final=li[li$lakeID%in%bassPEfinal$lakeID,]
li_final=li_final[order(li_final$lakeID),]
li_final$surfaceArea[3]=136 # from DNR website
bassPEfinal=bassPEfinal[order(bassPEfinal$lakeID),]

sum(li_final$lakeID==bassPEfinal$lakeID)

bassPEfinal$densityAreal=bassPEfinal$PE/(li_final$surfaceArea*10000/1e6) # km^2


#### get SD for beta correction
efCPE=efCPE[names(efCPE)%in%bassPEfinal$lakeID]

write.table(bassPEfinal,"cleanedMFEbassPE.csv",row.names=FALSE,sep=",")
