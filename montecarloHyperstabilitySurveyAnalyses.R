#### getting estimates of measurement error for efCPUE
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

# create vector to store average efCPE for each lake
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
  
  SDefCPE[i]=sd(cur_efCPE)  
  efCPE[i]=mean(cur_efCPE)
}

names(efCPE)=lakes
names(SDefCPE)=lakes
# remove SE and BA - no DNR-style  electrofishing
SDefCPE=SDefCPE[!is.na(efCPE)]
efCPE=efCPE[!is.na(efCPE)]

mean(SDefCPE/efCPE,na.rm=TRUE)  #0.63
range(SDefCPE/efCPE,na.rm=TRUE)


# doing corrections from Hansen et al. 2005

cur=BWPJoin[BWPJoin$species=="LARGEMOUTH BASS",]
cur=BWPJoin[BWPJoin$species=="SMALLMOUTH BASS",]
cur=BWPJoin[BWPJoin$species=="WALLEYE",]
cur=BWPJoin[BWPJoin$species=="BLACK CRAPPIE",]
cur=BWPJoin[BWPJoin$species=="YELLOW PERCH",]
cur=BWPJoin[BWPJoin$species=="BLUEGILL",]

simDens=matrix(NA,nrow(cur),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,cur$meanEF_CPEkm[i],cur$meanEF_CPEkm[i]*0.63)
  #simDens[i,]=rnorm(1000,cur$meanEF_CPEkm[i],cur$meanEF_CPEkm[i]*0.75)
  #simDens[i,]=rnorm(1000,cur$meanEF_CPEkm[i],cur$meanEF_CPEkm[i]*1)
}

bmc=numeric(1000)
for(i in 1:length(bmc)){
  logX=log(simDens[simDens[,i]>0,i])
  logY=log(cur$meanCPUE[simDens[,i]>0])
  fit=lm(logY~logX)
  bmc[i]=coef(fit)[2]
}

bols=coef(lm(cur$logCPUE~cur$logAbun))[2]

bbc=bols+(bols-bmc)

bols
mean(bbc)
median(bbc)
sort(bbc)[25]
sort(bbc)[975]
hist(bbc)



#### run monte carlo bias correction approach with mixed effect species model...
library(lme4)

simDens=matrix(NA,nrow(BWPJoin),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,BWPJoin$meanEF_CPEkm[i],BWPJoin$meanEF_CPEkm[i]*0.63)
  #simDens[i,]=rnorm(1000,BWPJoin$meanEF_CPEkm[i],BWPJoin$meanEF_CPEkm[i]*0.75)
  #simDens[i,]=rnorm(1000,BWPJoin$meanEF_CPEkm[i],BWPJoin$meanEF_CPEkm[i]*1)
}

BWP_species<-lmer(logCPUE~logAbun+species+logAbun:species+(1|WBICfactor),data=BWPJoin)
params=unlist(as.data.frame(coef(BWP_species)$WBICfactor)[1,])
bols=c(params[2],params[2]+params[8:12])

store=matrix(NA,1000,7)
for(i in 1:nrow(store)){
  logX=log(simDens[simDens[,i]>0,i])
  logY=log(BWPJoin$meanCPUE[simDens[,i]>0])
  toFit=data.frame(logX=logX,logY=logY,species=BWPJoin$species[simDens[,i]>0],WBICfactor=BWPJoin$WBICfactor[simDens[,i]>0])
  curBWP_null<-lmer(logY~logX+species+(1|WBICfactor),data=toFit)
  curBWP_species<-lmer(logY~logX+species+logX:species+(1|WBICfactor),data=toFit)
  curLRT=anova(curBWP_null,curBWP_species)
  
  params=unlist(as.data.frame(coef(curBWP_species)$WBICfactor)[1,])
  bmc=c(params[2],params[2]+params[8:12])
  
  store[i,1]=curLRT$`Pr(>Chisq)`[2]
  store[i,2:7]=bols+(bols-bmc)
}
colnames(store)=c("LRTp",paste("bbc",c("BC","BG","LMB","SMB","WALL","YP"),sep="_"))

mean(store[,1])
median(store[,1])
sum(store[,1]>0.05)


colMeans(store[,2:7])

colSums(store[,2:7]>1)



#### habitat tests
## BLACK CRAPPIE - a few warnings about convergence code 3 from bobyqa
cur=BWPlc[BWPlc$species=="BLACK CRAPPIE",] # N = 66

simDens=matrix(NA,nrow(cur),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,cur$meanEF_CPEkm[i],cur$meanEF_CPEkm[i]*0.63)
}

store=matrix(NA,1000,7)
for(i in 1:nrow(store)){
  logX=log(simDens[simDens[,i]>0,i])
  logY=log(cur$meanCPUE[simDens[,i]>0])
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity[simDens[,i]>0],area=cur$logLake[simDens[,i]>0],ripdev=cur$riparian_Developed[simDens[,i]>0],complex2=cur$ShoreComp2[simDens[,i]>0],area2=cur$LkArea2[simDens[,i]>0],ripdev2=cur$RipDev2[simDens[,i]>0],WBICfactor=cur$WBICfactor[simDens[,i]>0])
  
  ### fit six models and their nulls
  
  # riparian development
  dev1_null<-lmer(logY~logX+ripdev+(1|WBICfactor),data=toFit)
  dev1<-lmer(logY~logX+ripdev+logX:ripdev+(1|WBICfactor),data=toFit)
    
  dev2_null<-lmer(logY~logX+ripdev+ripdev2+(1|WBICfactor),data=toFit)
  dev2<-lmer(logY~logX+ripdev+ripdev2+ripdev:logX+ripdev2:logX+(1|WBICfactor),data=toFit)
    
  #shoreline complexity
  shor1_null<-lmer(logY~logX+complex+(1|WBICfactor),data=toFit)
  shor1<-lmer(logY~logX+complex+logX:complex+(1|WBICfactor),data=toFit)
  
  
  shor2_null<-lmer(logY~logX+complex+complex2+(1|WBICfactor),data=toFit)
  shor2<-lmer(logY~logX+complex+complex2+complex:logX+complex2:logX+(1|WBICfactor),data=toFit)
    
  #lake area
  area1_null<-lmer(logY~logX+area+(1|WBICfactor),data=toFit)
  area1<-lmer(logY~logX+area+logX:area+(1|WBICfactor),data=toFit)
  
  
  area2_null<-lmer(logY~logX+area+area2+(1|WBICfactor),data=toFit)
  area2<-lmer(logY~logX+area+area2+area:logX+area2:logX+(1|WBICfactor),data=toFit)
  

  # run six LRT and store p-values
  ripdev1lrt=anova(dev1_null,dev1)
  ripdev2lrt=anova(dev2_null,dev2)
  
  complex1lrt=anova(shor1_null,shor1) 
  complex2lrt=anova(shor2_null,shor2) 
  
  area1lrt=anova(area1_null,area1)
  area2lrt=anova(area2_null,area2)
  
  # both of these had median p-values below alpha=0.05/36
  rd1Vrd2=anova(dev1,dev2)
  
  store[i,1]=ripdev1lrt$`Pr(>Chisq)`[2]
  store[i,2]=ripdev2lrt$`Pr(>Chisq)`[2]
  store[i,3]=complex1lrt$`Pr(>Chisq)`[2]
  store[i,4]=complex2lrt$`Pr(>Chisq)`[2]
  store[i,5]=area1lrt$`Pr(>Chisq)`[2]
  store[i,6]=area2lrt$`Pr(>Chisq)`[2]
  store[i,7]=rd1Vrd2$`Pr(>Chisq)`[2]
}
BCstore=store

apply(BCstore,2,median)
colSums(BCstore>(0.05))/1000

## BLUEGILL - couple of warnings about convergence code
cur=BWPlc[BWPlc$species=="BLUEGILL",] # N = 89

simDens=matrix(NA,nrow(cur),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,cur$meanEF_CPEkm[i],cur$meanEF_CPEkm[i]*0.63)
}

store=matrix(NA,1000,6)
for(i in 1:nrow(store)){
  logX=log(simDens[simDens[,i]>0,i])
  logY=log(cur$meanCPUE[simDens[,i]>0])
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity[simDens[,i]>0],area=cur$logLake[simDens[,i]>0],ripdev=cur$riparian_Developed[simDens[,i]>0],complex2=cur$ShoreComp2[simDens[,i]>0],area2=cur$LkArea2[simDens[,i]>0],ripdev2=cur$RipDev2[simDens[,i]>0],WBICfactor=cur$WBICfactor[simDens[,i]>0])
  
  ### fit six models and their nulls
  
  # riparian development
  dev1_null<-lmer(logY~logX+ripdev+(1|WBICfactor),data=toFit)
  dev1<-lmer(logY~logX+ripdev+logX:ripdev+(1|WBICfactor),data=toFit)
  
  dev2_null<-lmer(logY~logX+ripdev+ripdev2+(1|WBICfactor),data=toFit)
  dev2<-lmer(logY~logX+ripdev+ripdev2+ripdev:logX+ripdev2:logX+(1|WBICfactor),data=toFit)
  
  #shoreline complexity
  shor1_null<-lmer(logY~logX+complex+(1|WBICfactor),data=toFit)
  shor1<-lmer(logY~logX+complex+logX:complex+(1|WBICfactor),data=toFit)
  
  
  shor2_null<-lmer(logY~logX+complex+complex2+(1|WBICfactor),data=toFit)
  shor2<-lmer(logY~logX+complex+complex2+complex:logX+complex2:logX+(1|WBICfactor),data=toFit)
  
  #lake area
  area1_null<-lmer(logY~logX+area+(1|WBICfactor),data=toFit)
  area1<-lmer(logY~logX+area+logX:area+(1|WBICfactor),data=toFit)
  
  
  area2_null<-lmer(logY~logX+area+area2+(1|WBICfactor),data=toFit)
  area2<-lmer(logY~logX+area+area2+area:logX+area2:logX+(1|WBICfactor),data=toFit)
  
  
  # run six LRT and store p-values
  ripdev1lrt=anova(dev1_null,dev1)
  ripdev2lrt=anova(dev2_null,dev2)#no
  
  complex1lrt=anova(shor1_null,shor1) #no
  complex2lrt=anova(shor2_null,shor2) #no
  
  area1lrt=anova(area1_null,area1)#no
  area2lrt=anova(area2_null,area2)#no
  
  store[i,1]=ripdev1lrt$`Pr(>Chisq)`[2]
  store[i,2]=ripdev2lrt$`Pr(>Chisq)`[2]
  store[i,3]=complex1lrt$`Pr(>Chisq)`[2]
  store[i,4]=complex2lrt$`Pr(>Chisq)`[2]
  store[i,5]=area1lrt$`Pr(>Chisq)`[2]
  store[i,6]=area2lrt$`Pr(>Chisq)`[2]
}
BGstore=store

apply(BGstore,2,median)
colSums(BGstore>(0.05))/1000

## LARGEMOUTH BASS -  singular warnings...
cur=BWPlc[BWPlc$species=="LARGEMOUTH BASS",] # N = 172

simDens=matrix(NA,nrow(cur),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,cur$meanEF_CPEkm[i],cur$meanEF_CPEkm[i]*0.63)
}

store=matrix(NA,1000,6)
for(i in 1:nrow(store)){
  logX=log(simDens[simDens[,i]>0,i])
  logY=log(cur$meanCPUE[simDens[,i]>0])
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity[simDens[,i]>0],area=cur$logLake[simDens[,i]>0],ripdev=cur$riparian_Developed[simDens[,i]>0],complex2=cur$ShoreComp2[simDens[,i]>0],area2=cur$LkArea2[simDens[,i]>0],ripdev2=cur$RipDev2[simDens[,i]>0],WBICfactor=cur$WBICfactor[simDens[,i]>0])
  
  ### fit six models and their nulls
  
  # riparian development
  dev1_null<-lmer(logY~logX+ripdev+(1|WBICfactor),data=toFit)
  dev1<-lmer(logY~logX+ripdev+logX:ripdev+(1|WBICfactor),data=toFit)
  
  dev2_null<-lmer(logY~logX+ripdev+ripdev2+(1|WBICfactor),data=toFit)
  dev2<-lmer(logY~logX+ripdev+ripdev2+ripdev:logX+ripdev2:logX+(1|WBICfactor),data=toFit)
  
  #shoreline complexity
  shor1_null<-lmer(logY~logX+complex+(1|WBICfactor),data=toFit)
  shor1<-lmer(logY~logX+complex+logX:complex+(1|WBICfactor),data=toFit)
  
  
  shor2_null<-lmer(logY~logX+complex+complex2+(1|WBICfactor),data=toFit)
  shor2<-lmer(logY~logX+complex+complex2+complex:logX+complex2:logX+(1|WBICfactor),data=toFit)
  
  #lake area
  area1_null<-lmer(logY~logX+area+(1|WBICfactor),data=toFit)
  area1<-lmer(logY~logX+area+logX:area+(1|WBICfactor),data=toFit)
  
  
  area2_null<-lmer(logY~logX+area+area2+(1|WBICfactor),data=toFit)
  area2<-lmer(logY~logX+area+area2+area:logX+area2:logX+(1|WBICfactor),data=toFit)
  
  
  # run six LRT and store p-values
  ripdev1lrt=anova(dev1_null,dev1)
  ripdev2lrt=anova(dev2_null,dev2)#no
  
  complex1lrt=anova(shor1_null,shor1) #no
  complex2lrt=anova(shor2_null,shor2) #no
  
  area1lrt=anova(area1_null,area1)#no
  area2lrt=anova(area2_null,area2)#no
  
  store[i,1]=ripdev1lrt$`Pr(>Chisq)`[2]
  store[i,2]=ripdev2lrt$`Pr(>Chisq)`[2]
  store[i,3]=complex1lrt$`Pr(>Chisq)`[2]
  store[i,4]=complex2lrt$`Pr(>Chisq)`[2]
  store[i,5]=area1lrt$`Pr(>Chisq)`[2]
  store[i,6]=area2lrt$`Pr(>Chisq)`[2]
}
LMBstore=store

apply(LMBstore,2,median)
colSums(LMBstore>(0.05))/1000

## SMALLMOUTH BASS - no warnings, no errors
cur=BWPlc[BWPlc$species=="SMALLMOUTH BASS",] # N = 155

simDens=matrix(NA,nrow(cur),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,cur$meanEF_CPEkm[i],cur$meanEF_CPEkm[i]*0.63)
}

store=matrix(NA,1000,6)
for(i in 1:nrow(store)){
  logX=log(simDens[simDens[,i]>0,i])
  logY=log(cur$meanCPUE[simDens[,i]>0])
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity[simDens[,i]>0],area=cur$logLake[simDens[,i]>0],ripdev=cur$riparian_Developed[simDens[,i]>0],complex2=cur$ShoreComp2[simDens[,i]>0],area2=cur$LkArea2[simDens[,i]>0],ripdev2=cur$RipDev2[simDens[,i]>0],WBICfactor=cur$WBICfactor[simDens[,i]>0])
  
  ### fit six models and their nulls
  
  # riparian development
  dev1_null<-lmer(logY~logX+ripdev+(1|WBICfactor),data=toFit)
  dev1<-lmer(logY~logX+ripdev+logX:ripdev+(1|WBICfactor),data=toFit)
  
  dev2_null<-lmer(logY~logX+ripdev+ripdev2+(1|WBICfactor),data=toFit)
  dev2<-lmer(logY~logX+ripdev+ripdev2+ripdev:logX+ripdev2:logX+(1|WBICfactor),data=toFit)
  
  #shoreline complexity
  shor1_null<-lmer(logY~logX+complex+(1|WBICfactor),data=toFit)
  shor1<-lmer(logY~logX+complex+logX:complex+(1|WBICfactor),data=toFit)
  
  
  shor2_null<-lmer(logY~logX+complex+complex2+(1|WBICfactor),data=toFit)
  shor2<-lmer(logY~logX+complex+complex2+complex:logX+complex2:logX+(1|WBICfactor),data=toFit)
  
  #lake area
  area1_null<-lmer(logY~logX+area+(1|WBICfactor),data=toFit)
  area1<-lmer(logY~logX+area+logX:area+(1|WBICfactor),data=toFit)
  
  
  area2_null<-lmer(logY~logX+area+area2+(1|WBICfactor),data=toFit)
  area2<-lmer(logY~logX+area+area2+area:logX+area2:logX+(1|WBICfactor),data=toFit)
  
  
  # run six LRT and store p-values
  ripdev1lrt=anova(dev1_null,dev1)
  ripdev2lrt=anova(dev2_null,dev2)#no
  
  complex1lrt=anova(shor1_null,shor1) #no
  complex2lrt=anova(shor2_null,shor2) #no
  
  area1lrt=anova(area1_null,area1)#no
  area2lrt=anova(area2_null,area2)#no
  
  store[i,1]=ripdev1lrt$`Pr(>Chisq)`[2]
  store[i,2]=ripdev2lrt$`Pr(>Chisq)`[2]
  store[i,3]=complex1lrt$`Pr(>Chisq)`[2]
  store[i,4]=complex2lrt$`Pr(>Chisq)`[2]
  store[i,5]=area1lrt$`Pr(>Chisq)`[2]
  store[i,6]=area2lrt$`Pr(>Chisq)`[2]
}
SMBstore=store

apply(SMBstore,2,median)
colSums(SMBstore>(0.05))/1000

## WALLEYE - no warnings, no errors
cur=BWPlc[BWPlc$species=="WALLEYE",] # N = 302

simDens=matrix(NA,nrow(cur),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,cur$meanEF_CPEkm[i],cur$meanEF_CPEkm[i]*0.63)
}

store=matrix(NA,1000,6)
for(i in 1:nrow(store)){
  logX=log(simDens[simDens[,i]>0,i])
  logY=log(cur$meanCPUE[simDens[,i]>0])
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity[simDens[,i]>0],area=cur$logLake[simDens[,i]>0],ripdev=cur$riparian_Developed[simDens[,i]>0],complex2=cur$ShoreComp2[simDens[,i]>0],area2=cur$LkArea2[simDens[,i]>0],ripdev2=cur$RipDev2[simDens[,i]>0],WBICfactor=cur$WBICfactor[simDens[,i]>0])
  
  ### fit six models and their nulls
  
  # riparian development
  dev1_null<-lmer(logY~logX+ripdev+(1|WBICfactor),data=toFit)
  dev1<-lmer(logY~logX+ripdev+logX:ripdev+(1|WBICfactor),data=toFit)
  
  dev2_null<-lmer(logY~logX+ripdev+ripdev2+(1|WBICfactor),data=toFit)
  dev2<-lmer(logY~logX+ripdev+ripdev2+ripdev:logX+ripdev2:logX+(1|WBICfactor),data=toFit)
  
  #shoreline complexity
  shor1_null<-lmer(logY~logX+complex+(1|WBICfactor),data=toFit)
  shor1<-lmer(logY~logX+complex+logX:complex+(1|WBICfactor),data=toFit)
  
  
  shor2_null<-lmer(logY~logX+complex+complex2+(1|WBICfactor),data=toFit)
  shor2<-lmer(logY~logX+complex+complex2+complex:logX+complex2:logX+(1|WBICfactor),data=toFit)
  
  #lake area
  area1_null<-lmer(logY~logX+area+(1|WBICfactor),data=toFit)
  area1<-lmer(logY~logX+area+logX:area+(1|WBICfactor),data=toFit)
  
  
  area2_null<-lmer(logY~logX+area+area2+(1|WBICfactor),data=toFit)
  area2<-lmer(logY~logX+area+area2+area:logX+area2:logX+(1|WBICfactor),data=toFit)
  
  
  # run six LRT and store p-values
  ripdev1lrt=anova(dev1_null,dev1)
  ripdev2lrt=anova(dev2_null,dev2)#no
  
  complex1lrt=anova(shor1_null,shor1) #no
  complex2lrt=anova(shor2_null,shor2) #no
  
  area1lrt=anova(area1_null,area1)#no
  area2lrt=anova(area2_null,area2)#no
  
  store[i,1]=ripdev1lrt$`Pr(>Chisq)`[2]
  store[i,2]=ripdev2lrt$`Pr(>Chisq)`[2]
  store[i,3]=complex1lrt$`Pr(>Chisq)`[2]
  store[i,4]=complex2lrt$`Pr(>Chisq)`[2]
  store[i,5]=area1lrt$`Pr(>Chisq)`[2]
  store[i,6]=area2lrt$`Pr(>Chisq)`[2]
}
WALLstore=store

apply(WALLstore,2,median)
colSums(WALLstore>(0.05))/1000


## YELLOW PERCH - no warnings, no errors
cur=BWPlc[BWPlc$species=="YELLOW PERCH",] # N = 88

simDens=matrix(NA,nrow(cur),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,cur$meanEF_CPEkm[i],cur$meanEF_CPEkm[i]*0.63)
}

store=matrix(NA,1000,6)
for(i in 1:nrow(store)){
  logX=log(simDens[simDens[,i]>0,i])
  logY=log(cur$meanCPUE[simDens[,i]>0])
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity[simDens[,i]>0],area=cur$logLake[simDens[,i]>0],ripdev=cur$riparian_Developed[simDens[,i]>0],complex2=cur$ShoreComp2[simDens[,i]>0],area2=cur$LkArea2[simDens[,i]>0],ripdev2=cur$RipDev2[simDens[,i]>0],WBICfactor=cur$WBICfactor[simDens[,i]>0])
  
  ### fit six models and their nulls
  
  # riparian development
  dev1_null<-lmer(logY~logX+ripdev+(1|WBICfactor),data=toFit)
  dev1<-lmer(logY~logX+ripdev+logX:ripdev+(1|WBICfactor),data=toFit)
  
  dev2_null<-lmer(logY~logX+ripdev+ripdev2+(1|WBICfactor),data=toFit)
  dev2<-lmer(logY~logX+ripdev+ripdev2+ripdev:logX+ripdev2:logX+(1|WBICfactor),data=toFit)
  
  #shoreline complexity
  shor1_null<-lmer(logY~logX+complex+(1|WBICfactor),data=toFit)
  shor1<-lmer(logY~logX+complex+logX:complex+(1|WBICfactor),data=toFit)
  
  
  shor2_null<-lmer(logY~logX+complex+complex2+(1|WBICfactor),data=toFit)
  shor2<-lmer(logY~logX+complex+complex2+complex:logX+complex2:logX+(1|WBICfactor),data=toFit)
  
  #lake area
  area1_null<-lmer(logY~logX+area+(1|WBICfactor),data=toFit)
  area1<-lmer(logY~logX+area+logX:area+(1|WBICfactor),data=toFit)
  
  
  area2_null<-lmer(logY~logX+area+area2+(1|WBICfactor),data=toFit)
  area2<-lmer(logY~logX+area+area2+area:logX+area2:logX+(1|WBICfactor),data=toFit)
  
  
  # run six LRT and store p-values
  ripdev1lrt=anova(dev1_null,dev1)
  ripdev2lrt=anova(dev2_null,dev2)#no
  
  complex1lrt=anova(shor1_null,shor1) #no
  complex2lrt=anova(shor2_null,shor2) #no
  
  area1lrt=anova(area1_null,area1)#no
  area2lrt=anova(area2_null,area2)#no
  
  store[i,1]=ripdev1lrt$`Pr(>Chisq)`[2]
  store[i,2]=ripdev2lrt$`Pr(>Chisq)`[2]
  store[i,3]=complex1lrt$`Pr(>Chisq)`[2]
  store[i,4]=complex2lrt$`Pr(>Chisq)`[2]
  store[i,5]=area1lrt$`Pr(>Chisq)`[2]
  store[i,6]=area2lrt$`Pr(>Chisq)`[2]
}
YWPstore=store

apply(YWPstore,2,median)
colSums(YWPstore>(0.05))/1000


#make a box and whisker or violin plot of p-values for those worth talking about
# BC ripdev1, BC ripdev2, LMB area2

boxplot(log10(c(BCstore[,1],BCstore[,2],BCstore[,7],LMBstore[,6]))~rep(c('BC - linear','BC - quadratic','BC lin. vs. quad.','LMB quadratic'),each=1000),xlab="",ylab="LRT p-value")






###############  remaking figures 2 and 3
rm(list=ls())

setwd("~/Documents/Research/People/Students/current/Mosley_Camille/hyperstabilitySurveyManuscript/hyperstabilitySurvey/")

# load function to load data from google drive
library(dplyr)

# load cleaned creel data
creel=read.csv("cleanedWDNRcreel.csv")

Bcreel=creel[creel$fishSpeciesCode%in%c("W11","W12"),]
Wcreel=creel[creel$fishSpeciesCode=="X22",]
Pcreel=creel[creel$fishSpeciesCode%in%c("X15","W14","W09"),]

# calculate average angling CPUE and sample size for each lake-year-species combination
lake_yearCPUE=creel %>%
  group_by(WBIC,fishSpeciesCode,surveyYear,county) %>%
  summarize(meanCPUE=mean(anglingCPUE),
            N=n())
lake_yearCPUE=as.data.frame(lake_yearCPUE)

####### load cleaned electrofishing data
lake_yearBASSef=read.csv('cleanedWDNRbassEF.csv')
lake_yearPANef=read.csv('cleanedWDNRpanEF.csv')
lake_yearWALLef=read.csv('cleanedWDNRwalleyeEF.csv')

lake_yearCPUE$species=""
lake_yearCPUE$species[lake_yearCPUE$fishSpeciesCode=="X22"]="WALLEYE"
lake_yearCPUE$species[lake_yearCPUE$fishSpeciesCode=="W11"]="SMALLMOUTH BASS"
lake_yearCPUE$species[lake_yearCPUE$fishSpeciesCode=="W12"]="LARGEMOUTH BASS"
lake_yearCPUE$species[lake_yearCPUE$fishSpeciesCode=="X15"]="YELLOW PERCH"
lake_yearCPUE$species[lake_yearCPUE$fishSpeciesCode=="W14"]="BLACK CRAPPIE"
lake_yearCPUE$species[lake_yearCPUE$fishSpeciesCode=="W09"]="BLUEGILL"

# trim species without EF data (can we get other species EF data?)
lake_yearCPUE=lake_yearCPUE[lake_yearCPUE$species!="",]

bassJoin=left_join(lake_yearBASSef,lake_yearCPUE,by=c("WBIC"="WBIC","species"="species","surveyYear"="surveyYear", "county"="county"))
bassJoin=bassJoin[!is.na(bassJoin$meanCPUE),]

panJoin=left_join(lake_yearPANef,lake_yearCPUE,by=c("WBIC"="WBIC","species"="species","surveyYear"="surveyYear", "county"="county"))
panJoin=panJoin[!is.na(panJoin$meanCPUE),]


wallJoin=left_join(lake_yearWALLef,lake_yearCPUE,by=c("WBIC"="WBIC","species"="species","surveyYear"="surveyYear", "county"="county"))
wallJoin=wallJoin[!is.na(wallJoin$meanCPUE),]

bassJoin$logCPUE=log(bassJoin$meanCPUE)
bassJoin$logAbun=log(bassJoin$meanEF_CPEkm)
bassJoin<- bassJoin[!is.na(bassJoin$logCPUE),]
bassJoin<- bassJoin[!is.na(bassJoin$logAbun),]
bassJoin<- bassJoin[is.finite(bassJoin$logCPUE),]

wallJoin$logCPUE=log(wallJoin$meanCPUE)
wallJoin$logAbun=log(wallJoin$meanEF_CPEkm)
wallJoin<- wallJoin[is.finite(wallJoin$logCPUE),]

panJoin$logCPUE=log(panJoin$meanCPUE)
panJoin$logAbun=log(panJoin$meanEF_CPEkm)
panJoin<- panJoin[is.finite(panJoin$logCPUE),]

BWJoin=full_join(bassJoin,wallJoin)
BWPJoin=full_join(BWJoin,panJoin)
BWPJoin$WBICfactor=as.factor(BWPJoin$WBIC)


#2)
library(lme4)

simDens=matrix(NA,nrow(BWPJoin),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,BWPJoin$meanEF_CPEkm[i],BWPJoin$meanEF_CPEkm[i]*0.63)
  #simDens[i,]=rnorm(1000,BWPJoin$meanEF_CPEkm[i],BWPJoin$meanEF_CPEkm[i]*0.75)
  #simDens[i,]=rnorm(1000,BWPJoin$meanEF_CPEkm[i],BWPJoin$meanEF_CPEkm[i]*1)
}

BWP_species<-lmer(logCPUE~logAbun+species+logAbun:species+(1|WBICfactor),data=BWPJoin)
params=unlist(as.data.frame(coef(BWP_species)$WBICfactor)[1,])
ols=params

store=matrix(NA,1000,13)
for(i in 1:nrow(store)){
  logX=log(simDens[simDens[,i]>0,i])
  logY=log(BWPJoin$meanCPUE[simDens[,i]>0])
  toFit=data.frame(logX=logX,logY=logY,species=BWPJoin$species[simDens[,i]>0],WBICfactor=BWPJoin$WBICfactor[simDens[,i]>0])
  curBWP_null<-lmer(logY~logX+species+(1|WBICfactor),data=toFit)
  curBWP_species<-lmer(logY~logX+species+logX:species+(1|WBICfactor),data=toFit)
  curLRT=anova(curBWP_null,curBWP_species)
  
  params=unlist(as.data.frame(coef(curBWP_species)$WBICfactor)[1,])
  
  store[i,1]=curLRT$`Pr(>Chisq)`[2]
  store[i,2:13]=params
}

colnames(store)=c("LRTp",names(params))

bmc=cbind(store[,3],store[,3]+store[,9:13])

qmc=cbind(store[,2],store[,2]+store[,4:8])


bc_bmc=bmc*0
bc_qmc=qmc*0

bols=c(ols[2],ols[2]+ols[8:12])
qols=c(ols[1],ols[1]+ols[3:7])

for(i in 1:ncol(bc_bmc)){
  bc_bmc[,i]=bols[i]+(bols[i]-bmc[,i])
  bc_qmc[,i]=qols[i]+(qols[i]-qmc[,i])
  #bols+(bols-bmc)
}

colnames(bc_bmc)=c("BC","BG","LMB","SMB","WAL","YP")
colnames(bc_qmc)=c("BC","BG","LMB","SMB","WAL","YP")


hist(bc_qmc[,1])
abline(v=qols[1])
abline(v=mean(bc_qmc[,1]),col='red')


median_bcBetas=apply(bc_bmc,2,median)
median_bcQs=apply(bc_qmc,2,median)

bcBetas_025=apply(bc_bmc,2,quantile,probs=0.025)
bcBetas_975=apply(bc_bmc,2,quantile,probs=0.975)

plot(1:6,median_bcBetas,col='black',pch=16,cex=1,ylim=c(0,1),xaxt='n',xlab="",ylab=expression(beta))
arrows(1:6,bcBetas_025,1:6,bcBetas_975,length=0.025,angle=90,code=3)
axis(side=1,at=1:6,labels=c("black\ncrappie","bluegill","largemouth\nbass","smallmouth\nbass","walleye","yellow\nperch"))

abundPred=seq(0.001,150,0.01)
logAP=log(abundPred)


plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="BLACK CRAPPIE"],BWPJoin$meanCPUE[BWPJoin$species=="BLACK CRAPPIE"],xlab="relative abundance",ylab="angler CPUE",main="black crappie",col='darkgrey',pch=16)
lines(abundPred,exp(median_bcQs[1]+median_bcBetas[1]*logAP),lwd=2)
#lines(abundPred,exp(qols[1]+bols[1]*logAP),lwd=2,col='red')
text(40,2.5,expression(paste(beta,"=0.228")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="BLACK CRAPPIE"],BWPJoin$meanCPUE[BWPJoin$species=="BLACK CRAPPIE"],xlab="relative abundance",ylab="angler CPUE",main="black crappie",col='darkgrey',pch=16,log="xy")
lines(abundPred,exp(median_bcQs[1]+median_bcBetas[1]*logAP),lwd=2)


plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="BLUEGILL"],BWPJoin$meanCPUE[BWPJoin$species=="BLUEGILL"],xlab="relative abundance",ylab="angler CPUE",main="bluegill",col='darkgrey',pch=16)
lines(abundPred,exp(median_bcQs[2]+median_bcBetas[2]*logAP),lwd=2)
#lines(abundPred,exp(qols[2]+bols[2]*logAP),lwd=2,col='red')
text(120,5,expression(paste(beta,"=0.268")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="BLUEGILL"],BWPJoin$meanCPUE[BWPJoin$species=="BLUEGILL"],xlab="relative abundance",ylab="angler CPUE",main="bluegill",col='darkgrey',pch=16,log="xy")
lines(abundPred,exp(median_bcQs[2]+median_bcBetas[2]*logAP),lwd=2)


plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="LARGEMOUTH BASS"],BWPJoin$meanCPUE[BWPJoin$species=="LARGEMOUTH BASS"],xlab="relative abundance",ylab="angler CPUE",main="largemouth bass",col='darkgrey',pch=16)
lines(abundPred,exp(median_bcQs[3]+median_bcBetas[3]*logAP),lwd=2)
#lines(abundPred,exp(qols[3]+bols[3]*logAP),lwd=2,col='red')
text(50,2.25,expression(paste(beta,"=0.465")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="LARGEMOUTH BASS"],BWPJoin$meanCPUE[BWPJoin$species=="LARGEMOUTH BASS"],xlab="relative abundance",ylab="angler CPUE",main="largemouth bass",col='darkgrey',pch=16,log="xy")
lines(abundPred,exp(median_bcQs[3]+median_bcBetas[3]*logAP),lwd=2)


plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="SMALLMOUTH BASS"],BWPJoin$meanCPUE[BWPJoin$species=="SMALLMOUTH BASS"],xlab="relative abundance",ylab="angler CPUE",main="smallmouth bass",col='darkgrey',pch=16)
lines(abundPred,exp(median_bcQs[4]+median_bcBetas[4]*logAP),lwd=2)
#lines(abundPred,exp(qols[4]+bols[4]*logAP),lwd=2,col='red')
text(15,0.2,expression(paste(beta,"=0.301")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="SMALLMOUTH BASS"],BWPJoin$meanCPUE[BWPJoin$species=="SMALLMOUTH BASS"],xlab="relative abundance",ylab="angler CPUE",main="smallmouth bass",col='darkgrey',pch=16,log="xy")
lines(abundPred,exp(median_bcQs[4]+median_bcBetas[4]*logAP),lwd=2)


plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="WALLEYE"],BWPJoin$meanCPUE[BWPJoin$species=="WALLEYE"],xlab="relative abundance",ylab="angler CPUE",main="walleye",col='darkgrey',pch=16)
lines(abundPred,exp(median_bcQs[5]+median_bcBetas[5]*logAP),lwd=2)
#lines(abundPred,exp(qols[5]+bols[5]*logAP),lwd=2,col='red')
text(125,1.5,expression(paste(beta,"=0.836")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="WALLEYE"],BWPJoin$meanCPUE[BWPJoin$species=="WALLEYE"],xlab="relative abundance",ylab="angler CPUE",main="walleye",col='darkgrey',pch=16,log="xy")
lines(abundPred,exp(median_bcQs[5]+median_bcBetas[5]*logAP),lwd=2)



plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="YELLOW PERCH"],BWPJoin$meanCPUE[BWPJoin$species=="YELLOW PERCH"],xlab="relative abundance",ylab="angler CPUE",main="yellow perch",col='darkgrey',pch=16)
lines(abundPred,exp(median_bcQs[6]+median_bcBetas[6]*logAP),lwd=2)
#lines(abundPred,exp(qols[6]+bols[6]*logAP),lwd=2,col='red')
text(27,2.5,expression(paste(beta,"=0.068")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="YELLOW PERCH"],BWPJoin$meanCPUE[BWPJoin$species=="YELLOW PERCH"],xlab="relative abundance",ylab="angler CPUE",main="yellow perch",col='darkgrey',pch=16,log="xy")
lines(abundPred,exp(median_bcQs[6]+median_bcBetas[6]*logAP),lwd=2)


### just predicted curves
plot(abundPred,exp(median_bcQs[1]+median_bcBetas[1]*logAP),type='l',lwd=2,xlab="electrofishing CPUE",ylab="predicted angler CPUE",xlim=c(0,150),ylim=c(0,4),col='purple')
lines(abundPred,exp(median_bcQs[2]+median_bcBetas[2]*logAP),lwd=2,col='blue')
lines(abundPred,exp(median_bcQs[3]+median_bcBetas[3]*logAP),lwd=2,col='darkgreen')
lines(abundPred,exp(median_bcQs[4]+median_bcBetas[4]*logAP),lwd=2,col='brown')
lines(abundPred,exp(median_bcQs[5]+median_bcBetas[5]*logAP),lwd=2,col='red')
lines(abundPred,exp(median_bcQs[6]+median_bcBetas[6]*logAP),lwd=2,col='goldenrod')
legend('topleft',c('walleye','largemouth bass','smallmouth bass','black crappie','bluegill','yellow perch'),lty=1,col=c('red','darkgreen','brown','purple','blue','goldenrod'))
