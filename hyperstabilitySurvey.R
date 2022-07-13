#Cleaned up hyperstability survey code 
#looking at hyperstability of angling CPUE as a function of electrofishing CPUE
#7-12-2022
# CLM, CJD, SEJ, CS
rm(list=ls())

# load packages
library(dplyr)
library(lme4)

#############################################################
#### loading data and quantifying sample sizes (Table 1) ####
#############################################################

# load cleaned creel data - produced by dataCollationCleaning.R from WDNR files
creel=read.csv("cleanedWDNRcreel.csv")

Bcreel=creel[creel$fishSpeciesCode%in%c("W11","W12"),]
Wcreel=creel[creel$fishSpeciesCode=="X22",]
Pcreel=creel[creel$fishSpeciesCode%in%c("X15","W14","W09"),]

#### CREEL sample sizes (Table 1)
LMBcreel=Bcreel[Bcreel$fishSpeciesCode=="W11",]
length(unique(LMBcreel$surveyYear))
length(unique(LMBcreel$WBIC))
length(unique(paste(LMBcreel$surveyYear,LMBcreel$WBIC,sep="_")))
SMBcreel=Bcreel[Bcreel$fishSpeciesCode=="W12",]
length(unique(SMBcreel$surveyYear))
length(unique(SMBcreel$WBIC))
length(unique(paste(SMBcreel$surveyYear,SMBcreel$WBIC,sep="_")))
BGcreel=Pcreel[Pcreel$fishSpeciesCode=="W09",]
length(unique(BGcreel$surveyYear))
length(unique(BGcreel$WBIC))
length(unique(paste(BGcreel$surveyYear,BGcreel$WBIC,sep="_")))
BCcreel=Pcreel[Pcreel$fishSpeciesCode=="W14",]
length(unique(BCcreel$surveyYear))
length(unique(BCcreel$WBIC))
length(unique(paste(BCcreel$surveyYear,BCcreel$WBIC,sep="_")))
YWPcreel=Pcreel[Pcreel$fishSpeciesCode=="X15",]
length(unique(YWPcreel$surveyYear))
length(unique(YWPcreel$WBIC))
length(unique(paste(YWPcreel$surveyYear,YWPcreel$WBIC,sep="_")))

length(unique(Wcreel$surveyYear))
length(unique(Wcreel$WBIC))
length(unique(paste(Wcreel$surveyYear,Wcreel$WBIC,sep="_")))


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

#### efCPUE sample sizes (Table 1)
LMBef=lake_yearBASSef[lake_yearBASSef$species=="LARGEMOUTH BASS",]
length(unique(LMBef$surveyYear))
length(unique(LMBef$WBIC))
length(unique(paste(LMBef$surveyYear,LMBef$WBIC,sep="_")))
SMBef=lake_yearBASSef[lake_yearBASSef$species=="SMALLMOUTH BASS",]
length(unique(SMBef$surveyYear))
length(unique(SMBef$WBIC))
length(unique(paste(SMBef$surveyYear,SMBef$WBIC,sep="_")))
BGef=lake_yearPANef[lake_yearPANef$species=="BLUEGILL",]
length(unique(BGef$surveyYear))
length(unique(BGef$WBIC))
length(unique(paste(BGef$surveyYear,BGef$WBIC,sep="_")))
BCef=lake_yearPANef[lake_yearPANef$species=="BLACK CRAPPIE",]
length(unique(BCef$surveyYear))
length(unique(BCef$WBIC))
length(unique(paste(BCef$surveyYear,BCef$WBIC,sep="_")))
YWPef=lake_yearPANef[lake_yearPANef$species=="YELLOW PERCH",]
length(unique(YWPef$surveyYear))
length(unique(YWPef$WBIC))
length(unique(paste(YWPef$surveyYear,YWPef$WBIC,sep="_")))

length(unique(lake_yearWALLef$surveyYear))
length(unique(lake_yearWALLef$WBIC))
length(unique(paste(lake_yearWALLef$surveyYear,lake_yearWALLef$WBIC,sep="_")))


##### merge data sets from angling CPUE and electrofishing CPUE to get exact lake-year matches
# convert fishSpeciesCode in lake_yearCPUE to species (name from ef stuff)
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

table(lake_yearCPUE$species) 
table(lake_yearCPUE$WBIC) # lakes with only one lake survey year (1610300,)
nrow(lake_yearBASSef)
nrow(bassJoin)
nrow(lake_yearPANef)
nrow(panJoin)
nrow(lake_yearWALLef)
nrow(wallJoin)

#creating log CPUE and log N columns, removing na's and infinite values so glm can run and we can make fit
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


#### both CREEL and efCPUE sample sizes (Table 1)
LMBboth=bassJoin[bassJoin$species=="LARGEMOUTH BASS",]
length(unique(LMBboth$surveyYear))
length(unique(LMBboth$WBIC))
length(unique(paste(LMBboth$surveyYear,LMBboth$WBIC,sep="_")))
SMBboth=bassJoin[bassJoin$species=="SMALLMOUTH BASS",]
length(unique(SMBboth$surveyYear))
length(unique(SMBboth$WBIC))
length(unique(paste(SMBboth$surveyYear,SMBboth$WBIC,sep="_")))
BGboth=panJoin[panJoin$species=="BLUEGILL",]
length(unique(BGboth$surveyYear))
length(unique(BGboth$WBIC))
length(unique(paste(BGboth$surveyYear,BGboth$WBIC,sep="_")))
BCboth=panJoin[panJoin$species=="BLACK CRAPPIE",]
length(unique(BCboth$surveyYear))
length(unique(BCboth$WBIC))
length(unique(paste(BCboth$surveyYear,BCboth$WBIC,sep="_")))
YWPboth=panJoin[panJoin$species=="YELLOW PERCH",]
length(unique(YWPboth$surveyYear))
length(unique(YWPboth$WBIC))
length(unique(paste(YWPboth$surveyYear,YWPboth$WBIC,sep="_")))

length(unique(wallJoin$surveyYear))
length(unique(wallJoin$WBIC))
length(unique(paste(wallJoin$surveyYear,wallJoin$WBIC,sep="_")))

#combining all fish species
BWJoin=full_join(bassJoin,wallJoin)
BWPJoin=full_join(BWJoin,panJoin)
BWPJoin$WBICfactor=as.factor(BWPJoin$WBIC)

# load lake variables
hansen=read.csv("cleanedHansenetalLakeinfo.csv")

# join lake data with hyperstability data
BWPlc=left_join(BWPJoin, hansen, by="WBIC")

# calculate a couple other lake shape metrics
BWPlc$shorelineComplexity=BWPlc$lakeChar_Length/(2*sqrt(pi*BWPlc$lakeChar_Area)) #ratio of the shoreline length (i.e. perimeter) to the perimeter of an equally sized circle
BWPlc$logLake=log10(BWPlc$lakeChar_Area)

BWPlc$ShoreComp2=BWPlc$shorelineComplexity^2
BWPlc$LkArea2=BWPlc$logLake^2
BWPlc$RipDev2=BWPlc$riparian_Developed^2


#### all of CREEL, efCPUE, and lake characteristics sample sizes (Table 1)
LMBall=BWPlc[BWPlc$species=="LARGEMOUTH BASS",]
length(unique(LMBall$surveyYear))
length(unique(LMBall$WBIC))
length(unique(paste(LMBall$surveyYear,LMBall$WBIC,sep="_")))
SMBall=BWPlc[BWPlc$species=="SMALLMOUTH BASS",]
length(unique(SMBall$surveyYear))
length(unique(SMBall$WBIC))
length(unique(paste(SMBall$surveyYear,SMBall$WBIC,sep="_")))
BGall=BWPlc[BWPlc$species=="BLUEGILL",]
length(unique(BGall$surveyYear))
length(unique(BGall$WBIC))
length(unique(paste(BGall$surveyYear,BGall$WBIC,sep="_")))
BCall=BWPlc[BWPlc$species=="BLACK CRAPPIE",]
length(unique(BCall$surveyYear))
length(unique(BCall$WBIC))
length(unique(paste(BCall$surveyYear,BCall$WBIC,sep="_")))
YWPall=BWPlc[BWPlc$species=="YELLOW PERCH",]
length(unique(YWPall$surveyYear))
length(unique(YWPall$WBIC))
length(unique(paste(YWPall$surveyYear,YWPall$WBIC,sep="_")))
Wall=BWPlc[BWPlc$species=="WALLEYE",]
length(unique(Wall$surveyYear))
length(unique(Wall$WBIC))
length(unique(paste(Wall$surveyYear,Wall$WBIC,sep="_")))

nrow(BWPlc)
length(unique(BWPlc$WBIC))
length(unique(BWPlc$county))
range(BWPlc$surveyYear)


###################################
#### PEs vs. efCPUE (Figure 1) ####
###################################

##### walleye
#### looking at WDNR walleye PE vs efCPUE ####

WallPE=read.csv("cleanedWDNRwalleyePE.csv")


# join lake data with PE data to calculate fish per shoreline  distance
WallPElc=left_join(WallPE, hansen[,c(1:3)], by="WBIC")
WallPElc=WallPElc[!is.na(WallPElc$PE),]
WallPElc=WallPElc[!is.na(WallPElc$lakeChar_Length),]
WallPElc$lakeChar_Length=WallPElc$lakeChar_Length/1000  # m to km
WallPElc$densityShoreline=WallPElc$PE/WallPElc$lakeChar_Length  # fish per km  shoreline

cpueCheck_WDNR=full_join(lake_yearWALLef,WallPElc,by=c("WBIC","surveyYear"))
cpueCheck_WDNR=cpueCheck_WDNR[!is.na(cpueCheck_WDNR$Density),]
cpueCheck_WDNR=cpueCheck_WDNR[!is.na(cpueCheck_WDNR$meanEF_CPEkm),] #159 obs
cpueCheck_WDNR$wbicFactor=as.factor(cpueCheck_WDNR$WBIC)

# convert to fish /km2
cpueCheck_WDNR$Density=cpueCheck_WDNR$Density/0.00404686
densityFit<-lmer(meanEF_CPEkm~Density+(1|wbicFactor),data=cpueCheck_WDNR)
shoreDensityFit<-lmer(meanEF_CPEkm~densityShoreline+(1|wbicFactor),data=cpueCheck_WDNR)

summary(densityFit)
summary(shoreDensityFit)

plot(cpueCheck_WDNR$Density,cpueCheck_WDNR$meanEF_CPEkm,col='darkgrey',pch=16)
abline(a=mean(coefficients(densityFit)$wbicFactor[,1]),b=mean(coefficients(densityFit)$wbicFactor[,2]),lwd=2)
plot(cpueCheck_WDNR$densityShoreline,cpueCheck_WDNR$meanEF_CPEkm,col='darkgrey',pch=16)
abline(a=mean(coefficients(shoreDensityFit)$wbicFactor[,1]),b=mean(coefficients(shoreDensityFit)$wbicFactor[,2]),lwd=2)

nullWalleye=lmer(meanEF_CPEkm~(1|wbicFactor),data=cpueCheck_WDNR)

anova(densityFit,nullWalleye)
anova(shoreDensityFit,nullWalleye)

cor(cpueCheck_WDNR$meanEF_CPEkm,cpueCheck_WDNR$Density)
cor(cpueCheck_WDNR$meanEF_CPEkm,cpueCheck_WDNR$densityShoreline)

##### bass
bassPEfinal=read.csv("cleanedMFEbassPE.csv")

plot(bassPEfinal$densityShoreline,bassPEfinal$efCPE,col='darkgrey',pch=16)

bassEFvPEfitShoreline=lm(bassPEfinal$efCPE~bassPEfinal$densityShoreline)
summary(bassEFvPEfitShoreline)
abline(bassEFvPEfitShoreline,lwd=2)

plot(bassPEfinal$densityAreal,bassPEfinal$efCPE,col='darkgrey',pch=16)

bassEFvPEfitAreal=lm(bassPEfinal$efCPE~bassPEfinal$densityAreal)
summary(bassEFvPEfitAreal)
abline(bassEFvPEfitAreal,lwd=2)

cor(bassPEfinal$efCPE,bassPEfinal$densityShoreline)
cor(bassPEfinal$efCPE,bassPEfinal$densityAreal)

## estimating beta  for EF vs PE relationships

# Walleye
simDens=matrix(NA,nrow(cpueCheck_WDNR),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,log(cpueCheck_WDNR$densityShoreline[i]),0.1451)  # sd based WDNR walleye PEs
  #simDens[i,]=rnorm(1000,log(cpueCheck_WDNR$Density[i]),0.1451)  # sd based WDNR walleye PEs
}

bmc=numeric(1000)
for(i in 1:length(bmc)){
  logX=simDens[,i]
  logY=log(cpueCheck_WDNR$meanEF_CPEkm)
  fit=lmer(logY~logX+(1|cpueCheck_WDNR$wbicFactor))
  bmc[i]=fit@pp$delb[2]#coef(fit)[2]
}

logX=log(cpueCheck_WDNR$densityShoreline)
logY=log(cpueCheck_WDNR$meanEF_CPEkm)
bols=(lmer(logY~logX+(1|cpueCheck_WDNR$wbicFactor)))@pp$delb[2]

bbc=bols+(bols-bmc)

bols
mean(bbc) #0.78
median(bbc)
sort(bbc)[25]   #0.74
sort(bbc)[975]  #0.82
hist(bbc)


# bass
simDens=matrix(NA,nrow(bassPEfinal),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,log(bassPEfinal$densityShoreline[i]),0.1451)  # sd based on WDNR walleye PEs
  #simDens[i,]=rnorm(1000,log(bassPEfinal$densityAreal[i]),0.1451)  # sd based on WDNR walleye PEs
}

bmc=numeric(1000)
for(i in 1:length(bmc)){
  logX=simDens[,i]
  logY=log(bassPEfinal$efCPE)
  fit=lm(logY~logX)
  bmc[i]=coef(fit)[2]
}

logX=log(bassPEfinal$densityShoreline)
logY=log(bassPEfinal$efCPE)
bols=coef(lm(logY~logX))[2]

bbc=bols+(bols-bmc)

bols
mean(bbc)
median(bbc)
sort(bbc)[25]
sort(bbc)[975]
hist(bbc)


#############################################################
#### Species effects on beta (Figures 2 & 3 and Table 3) ####
#############################################################

# using bias correction from Hansen et al. 2005, but  using lognormal distribution (normal distribution of logged efCPUE)
simDens=matrix(NA,nrow(BWPJoin),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,log(BWPJoin$meanEF_CPEkm[i]),0.57)
}

BWP_species<-lmer(logCPUE~logAbun+species+logAbun:species+(1|WBICfactor),data=BWPJoin)
params=unlist(as.data.frame(coef(BWP_species)$WBICfactor)[1,])
ols=params

store=matrix(NA,1000,13)
for(i in 1:nrow(store)){
  logX=simDens[,i]
  logY=log(BWPJoin$meanCPUE)
  toFit=data.frame(logX=logX,logY=logY,species=BWPJoin$species,WBICfactor=BWPJoin$WBICfactor)
  curBWP_null<-lmer(logY~logX+species+(1|WBICfactor),data=toFit)
  curBWP_species<-lmer(logY~logX+species+logX:species+(1|WBICfactor),data=toFit)
  curLRT=anova(curBWP_null,curBWP_species)
  
  params=unlist(as.data.frame(coef(curBWP_species)$WBICfactor)[1,])
  
  store[i,1]=curLRT$`Pr(>Chisq)`[2]
  store[i,2:13]=params
}

colnames(store)=c("LRTp",names(params))

mean(store[,1])
median(store[,1])
sum(store[,1]>0.05)


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
text(40,2.5,expression(paste(beta,"=0.216")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="BLACK CRAPPIE"],BWPJoin$meanCPUE[BWPJoin$species=="BLACK CRAPPIE"],xlab="relative abundance",ylab="angler CPUE",main="black crappie",col='darkgrey',pch=16,log="xy")
lines(abundPred,exp(median_bcQs[1]+median_bcBetas[1]*logAP),lwd=2)


plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="BLUEGILL"],BWPJoin$meanCPUE[BWPJoin$species=="BLUEGILL"],xlab="relative abundance",ylab="angler CPUE",main="bluegill",col='darkgrey',pch=16)
lines(abundPred,exp(median_bcQs[2]+median_bcBetas[2]*logAP),lwd=2)
#lines(abundPred,exp(qols[2]+bols[2]*logAP),lwd=2,col='red')
text(120,5,expression(paste(beta,"=0.25")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="BLUEGILL"],BWPJoin$meanCPUE[BWPJoin$species=="BLUEGILL"],xlab="relative abundance",ylab="angler CPUE",main="bluegill",col='darkgrey',pch=16,log="xy")
lines(abundPred,exp(median_bcQs[2]+median_bcBetas[2]*logAP),lwd=2)


plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="LARGEMOUTH BASS"],BWPJoin$meanCPUE[BWPJoin$species=="LARGEMOUTH BASS"],xlab="relative abundance",ylab="angler CPUE",main="largemouth bass",col='darkgrey',pch=16)
lines(abundPred,exp(median_bcQs[3]+median_bcBetas[3]*logAP),lwd=2)
#lines(abundPred,exp(qols[3]+bols[3]*logAP),lwd=2,col='red')
text(50,2.25,expression(paste(beta,"=0.442")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="LARGEMOUTH BASS"],BWPJoin$meanCPUE[BWPJoin$species=="LARGEMOUTH BASS"],xlab="relative abundance",ylab="angler CPUE",main="largemouth bass",col='darkgrey',pch=16,log="xy")
lines(abundPred,exp(median_bcQs[3]+median_bcBetas[3]*logAP),lwd=2)


plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="SMALLMOUTH BASS"],BWPJoin$meanCPUE[BWPJoin$species=="SMALLMOUTH BASS"],xlab="relative abundance",ylab="angler CPUE",main="smallmouth bass",col='darkgrey',pch=16)
lines(abundPred,exp(median_bcQs[4]+median_bcBetas[4]*logAP),lwd=2)
#lines(abundPred,exp(qols[4]+bols[4]*logAP),lwd=2,col='red')
text(15,0.2,expression(paste(beta,"=0.278")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="SMALLMOUTH BASS"],BWPJoin$meanCPUE[BWPJoin$species=="SMALLMOUTH BASS"],xlab="relative abundance",ylab="angler CPUE",main="smallmouth bass",col='darkgrey',pch=16,log="xy")
lines(abundPred,exp(median_bcQs[4]+median_bcBetas[4]*logAP),lwd=2)


plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="WALLEYE"],BWPJoin$meanCPUE[BWPJoin$species=="WALLEYE"],xlab="relative abundance",ylab="angler CPUE",main="walleye",col='darkgrey',pch=16)
lines(abundPred,exp(median_bcQs[5]+median_bcBetas[5]*logAP),lwd=2)
#lines(abundPred,exp(qols[5]+bols[5]*logAP),lwd=2,col='red')
text(125,1.5,expression(paste(beta,"=0.762")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="WALLEYE"],BWPJoin$meanCPUE[BWPJoin$species=="WALLEYE"],xlab="relative abundance",ylab="angler CPUE",main="walleye",col='darkgrey',pch=16,log="xy")
lines(abundPred,exp(median_bcQs[5]+median_bcBetas[5]*logAP),lwd=2)



plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="YELLOW PERCH"],BWPJoin$meanCPUE[BWPJoin$species=="YELLOW PERCH"],xlab="relative abundance",ylab="angler CPUE",main="yellow perch",col='darkgrey',pch=16)
lines(abundPred,exp(median_bcQs[6]+median_bcBetas[6]*logAP),lwd=2)
#lines(abundPred,exp(qols[6]+bols[6]*logAP),lwd=2,col='red')
text(27,2.5,expression(paste(beta,"=0.063")))

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


######################################################
#### habitat effects (Table 4 & S2 and Figure S4) ####
######################################################

## BLACK CRAPPIE
cur=BWPlc[BWPlc$species=="BLACK CRAPPIE",] # N = 66

simDens=matrix(NA,nrow(cur),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,log(cur$meanEF_CPEkm[i]),0.57)
}

store=matrix(NA,1000,7)
for(i in 1:nrow(store)){
  logX=simDens[,i]
  logY=log(cur$meanCPUE)
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity,area=cur$logLake,ripdev=cur$riparian_Developed,complex2=cur$ShoreComp2,area2=cur$LkArea2,ripdev2=cur$RipDev2,WBICfactor=cur$WBICfactor)
  
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


## BLUEGILL
cur=BWPlc[BWPlc$species=="BLUEGILL",] # N = 89

simDens=matrix(NA,nrow(cur),1000)
for(i in 1:nrow(simDens)){
  simDens[i,]=rnorm(1000,log(cur$meanEF_CPEkm[i]),0.57)
}

store=matrix(NA,1000,6)
for(i in 1:nrow(store)){
  logX=simDens[,i]
  logY=log(cur$meanCPUE)
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity,area=cur$logLake,ripdev=cur$riparian_Developed,complex2=cur$ShoreComp2,area2=cur$LkArea2,ripdev2=cur$RipDev2,WBICfactor=cur$WBICfactor)
  
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
  simDens[i,]=rnorm(1000,log(cur$meanEF_CPEkm[i]),0.57)
}

store=matrix(NA,1000,6)
for(i in 1:nrow(store)){
  logX=simDens[,i]
  logY=log(cur$meanCPUE)
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity,area=cur$logLake,ripdev=cur$riparian_Developed,complex2=cur$ShoreComp2,area2=cur$LkArea2,ripdev2=cur$RipDev2,WBICfactor=cur$WBICfactor)
  
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
  simDens[i,]=rnorm(1000,log(cur$meanEF_CPEkm[i]),0.57)
}

store=matrix(NA,1000,6)
for(i in 1:nrow(store)){
  logX=simDens[,i]
  logY=log(cur$meanCPUE)
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity,area=cur$logLake,ripdev=cur$riparian_Developed,complex2=cur$ShoreComp2,area2=cur$LkArea2,ripdev2=cur$RipDev2,WBICfactor=cur$WBICfactor)
  
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
  simDens[i,]=rnorm(1000,log(cur$meanEF_CPEkm[i]),0.57)
}

store=matrix(NA,1000,6)
for(i in 1:nrow(store)){
  logX=simDens[,i]
  logY=log(cur$meanCPUE)
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity,area=cur$logLake,ripdev=cur$riparian_Developed,complex2=cur$ShoreComp2,area2=cur$LkArea2,ripdev2=cur$RipDev2,WBICfactor=cur$WBICfactor)
  
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
  simDens[i,]=rnorm(1000,log(cur$meanEF_CPEkm[i]),0.57)
}

store=matrix(NA,1000,6)
for(i in 1:nrow(store)){
  logX=simDens[,i]
  logY=log(cur$meanCPUE)
  toFit=data.frame(logX=logX,logY=logY,complex=cur$shorelineComplexity,area=cur$logLake,ripdev=cur$riparian_Developed,complex2=cur$ShoreComp2,area2=cur$LkArea2,ripdev2=cur$RipDev2,WBICfactor=cur$WBICfactor)
  
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
#BC ripdev2, LMB area2

boxplot(log10(c(BCstore[,2],LMBstore[,6]))~rep(c('BC - quadratic for riparian development','LMB quadratic for lake area'),each=1000),xlab="",ylab="LRT p-value")


# looking at models that consistently outperformed the null model
cur=BWPlc[BWPlc$species=="BLACK CRAPPIE",]
toFit=data.frame(logX=log(cur$meanEF_CPEkm),logY=log(cur$meanCPUE),complex=cur$shorelineComplexity,area=cur$logLake,ripdev=cur$riparian_Developed,complex2=cur$ShoreComp2,area2=cur$LkArea2,ripdev2=cur$RipDev2,WBICfactor=cur$WBICfactor)
BCripdev2<-lmer(logY~logX+ripdev+ripdev2+ripdev:logX+ripdev2:logX+(1|WBICfactor),data=toFit)
summary(BCripdev2)

a=20.0861 # quadratic interaction with abund
b=-3.7379 # linear interaction with abun
c=0.2368 #just abun effect
x=seq(0,0.4,0.01)

plot(x,a*x*x+b*x+c,type="l",xlab="riparian developed",ylab="beta", main="crappie riparian development effect")
-b/(2*a) # x coordinate of vertex



cur=BWPlc[BWPlc$species=="LARGEMOUTH BASS",]
toFit=data.frame(logX=log(cur$meanEF_CPEkm),logY=log(cur$meanCPUE),complex=cur$shorelineComplexity,area=cur$logLake,ripdev=cur$riparian_Developed,complex2=cur$ShoreComp2,area2=cur$LkArea2,ripdev2=cur$RipDev2,WBICfactor=cur$WBICfactor)
LMBarea2<-lmer(logY~logX+area+area2+area:logX+area2:logX+(1|WBICfactor),data=toFit)
summary(LMBarea2)

a=-0.29336  #quadratic interaction with abund
b=3.91775 # linear interaction with abun
c=-12.60165 #just abun effect
x=seq(2,8,0.1)

plot(x,a*x*x+b*x+c,type="l",xlab="log lake area",ylab="beta", main="LMB - lake area affect")
     
10^(-b/(2*a))/1e6 # x coordinate of vertex in km^2


######################################################################
#### generating lake-year table for data availability (Figure S2) ####
######################################################################
WBICs=sort(unique(BWPJoin$WBIC))
years=sort(unique(BWPJoin$surveyYear))

whatWeHave=paste(BWPJoin$WBIC,BWPJoin$surveyYear,sep="_")

sampTable=matrix(0,length(WBICs),length(years))
for(i in 1:length(WBICs)){
  for(j in 1:length(years)){
    if(paste(WBICs[i],years[j],sep="_")%in%whatWeHave){
      sampTable[i,j]=1
    }
  }
}

dotchart(rep(years[1],length(WBICs)),pt.cex=sampTable[,1]*0.5,labels=WBICs,xlab="year",ylab="WBIC",pch=16,xlim=c(min(years),max(years)),cex=0.5,lcolor='white',main="all")
for(i in 2:length(years)){
  points(rep(years[i],length(WBICs)),1:length(WBICs),cex=sampTable[,i]*0.5,pch=16)
  
}


#walleye
whatWeHaveW=paste(BWPJoin$WBIC[BWPJoin$species=="WALLEYE"],BWPJoin$surveyYear[BWPJoin$species=="WALLEYE"],sep="_")

sampTableW=matrix(0,length(WBICs),length(years))
for(i in 1:length(WBICs)){
  for(j in 1:length(years)){
    if(paste(WBICs[i],years[j],sep="_")%in%whatWeHaveW){
      sampTableW[i,j]=1
    }
  }
}

dotchart(rep(years[1],length(WBICs)),pt.cex=sampTableW[,1]*0.5,labels=WBICs,xlab="year",ylab="WBIC",pch=16,xlim=c(min(years),max(years)),cex=0.5,lcolor='white',main="walleye")
for(i in 2:length(years)){
  points(rep(years[i],length(WBICs)),1:length(WBICs),cex=sampTableW[,i]*0.5,pch=16)
  
}

#LMB
whatWeHaveLMB=paste(BWPJoin$WBIC[BWPJoin$species=="LARGEMOUTH BASS"],BWPJoin$surveyYear[BWPJoin$species=="LARGEMOUTH BASS"],sep="_")

sampTableLMB=matrix(0,length(WBICs),length(years))
for(i in 1:length(WBICs)){
  for(j in 1:length(years)){
    if(paste(WBICs[i],years[j],sep="_")%in%whatWeHaveLMB){
      sampTableLMB[i,j]=1
    }
  }
}

dotchart(rep(years[1],length(WBICs)),pt.cex=sampTableLMB[,1]*0.5,labels=WBICs,xlab="year",ylab="WBIC",pch=16,xlim=c(min(years),max(years)),cex=0.5,lcolor='white',main="largemouth bass")
for(i in 2:length(years)){
  points(rep(years[i],length(WBICs)),1:length(WBICs),cex=sampTableLMB[,i]*0.5,pch=16)
  
}


#SMB
whatWeHaveSMB=paste(BWPJoin$WBIC[BWPJoin$species=="SMALLMOUTH BASS"],BWPJoin$surveyYear[BWPJoin$species=="SMALLMOUTH BASS"],sep="_")

sampTableSMB=matrix(0,length(WBICs),length(years))
for(i in 1:length(WBICs)){
  for(j in 1:length(years)){
    if(paste(WBICs[i],years[j],sep="_")%in%whatWeHaveSMB){
      sampTableSMB[i,j]=1
    }
  }
}

dotchart(rep(years[1],length(WBICs)),pt.cex=sampTableSMB[,1]*0.5,labels=WBICs,xlab="year",ylab="WBIC",pch=16,xlim=c(min(years),max(years)),cex=0.5,lcolor='white',main="smallmouth bass")
for(i in 2:length(years)){
  points(rep(years[i],length(WBICs)),1:length(WBICs),cex=sampTableSMB[,i]*0.5,pch=16)
  
}

#BC
whatWeHaveBC=paste(BWPJoin$WBIC[BWPJoin$species=="BLACK CRAPPIE"],BWPJoin$surveyYear[BWPJoin$species=="BLACK CRAPPIE"],sep="_")

sampTableBC=matrix(0,length(WBICs),length(years))
for(i in 1:length(WBICs)){
  for(j in 1:length(years)){
    if(paste(WBICs[i],years[j],sep="_")%in%whatWeHaveBC){
      sampTableBC[i,j]=1
    }
  }
}

dotchart(rep(years[1],length(WBICs)),pt.cex=sampTableBC[,1]*0.5,labels=WBICs,xlab="year",ylab="WBIC",pch=16,xlim=c(min(years),max(years)),cex=0.5,lcolor='white',main="black crappie")
for(i in 2:length(years)){
  points(rep(years[i],length(WBICs)),1:length(WBICs),cex=sampTableBC[,i]*0.5,pch=16)
  
}


#BG
whatWeHaveBG=paste(BWPJoin$WBIC[BWPJoin$species=="BLUEGILL"],BWPJoin$surveyYear[BWPJoin$species=="BLUEGILL"],sep="_")

sampTableBG=matrix(0,length(WBICs),length(years))
for(i in 1:length(WBICs)){
  for(j in 1:length(years)){
    if(paste(WBICs[i],years[j],sep="_")%in%whatWeHaveBG){
      sampTableBG[i,j]=1
    }
  }
}

dotchart(rep(years[1],length(WBICs)),pt.cex=sampTableBG[,1]*0.5,labels=WBICs,xlab="year",ylab="WBIC",pch=16,xlim=c(min(years),max(years)),cex=0.5,lcolor='white',main="bluegill")
for(i in 2:length(years)){
  points(rep(years[i],length(WBICs)),1:length(WBICs),cex=sampTableBG[,i]*0.5,pch=16)
  
}


#YP
whatWeHaveYP=paste(BWPJoin$WBIC[BWPJoin$species=="YELLOW PERCH"],BWPJoin$surveyYear[BWPJoin$species=="YELLOW PERCH"],sep="_")

sampTableYP=matrix(0,length(WBICs),length(years))
for(i in 1:length(WBICs)){
  for(j in 1:length(years)){
    if(paste(WBICs[i],years[j],sep="_")%in%whatWeHaveYP){
      sampTableYP[i,j]=1
    }
  }
}

dotchart(rep(years[1],length(WBICs)),pt.cex=sampTableYP[,1]*0.5,labels=WBICs,xlab="year",ylab="WBIC",pch=16,xlim=c(min(years),max(years)),cex=0.5,lcolor='white',main="yellow perch")
for(i in 2:length(years)){
  points(rep(years[i],length(WBICs)),1:length(WBICs),cex=sampTableYP[,i]*0.5,pch=16)
  
}


#############################################################################################
#### how does hyperstable EF impact ability to detect angler hyperstability? (Figure S3) ####
#############################################################################################
x=runif(100,50,1000)
efCPE=0.15*x^0.78   # beta from WDNR walleye data
anglerBetas=seq(0.5,1.5,0.05)
emergentBeta=numeric(length(anglerBetas))
for(i in 1:length(emergentBeta)){
  CPUE=0.02*x^anglerBetas[i]
  fit=lm(log(CPUE)~log(efCPE))
  emergentBeta[i]=coef(fit)[2]
}

plot(anglerBetas,emergentBeta,xlim=c(0.4,2),ylim=c(0.4,2),xlab="known beta angling", ylab="observed beta from angler-electrofishing CPUE")
abline(a=0,b=1)









#########################
####  power analyses ####
#########################

#creelN varies from 391-421
#efN varies from  554-948
#combinedN varies from 66-302

# betas vary from 0.056 - 0.63


# walleye angCPUE: 0.0046-2.11
# walleye efCPUE: 0.41-160
# sd of residuals: 0.29

# LMB angCPUE: 0.0059-2.51
# LMB efCPUE: 0.06-60.77
# sd of residals: 0.32

# SMB angCPUE: 0.0198-1.12
# SMB efCPUE: 0.062-18.89
# sd of residuals: 0.178

# BLG range angCPUE: 0.23-5.67
# BLG efCPUE: 0.15-143.54
# sd of residuals: 1.19

# YWP angCPUE: 0.061-2.98
# YWP efCPUE: 0.006-34.64
# sd of residuals: 0.679

# BC angCPUE: 0.006-2.83
# BC efCPUE: 0.025-44.32
# BC sd of residuals: 0.568

meanResidSD=mean(c(0.29,0.32,0.178,1.19,0.678,0.568))
minEFobs=min(c(0.41,0.06,0.062,0.15,0.006,0.025))
maxEFobs=max(c(160,60.77,18.89,143.54,34.64,44.32))

## simulate efCPUE v. angCPUE data
#     -assumes range of efCPUE and angCPUE that cover observed range
#     -q is approximately representative
#     -SD is average of observed
simCPUE<-function(N,beta,SD=meanResidSD,minEF=minEFobs,maxEF=maxEFobs,q=0.5){
  EF=runif(N,minEF,maxEF)
  ang=rnorm(N,mean=q*EF^beta,sd=SD)
  ang[ang<0]=q*EF[ang<0]^beta
  
  fit=glm(log(ang)~log(EF))
  
  point=summary(fit)$coefficients[2,1]
  lower=summary(fit)$coefficients[2,1]-1.96*summary(fit)$coefficients[2,2]
  upper=summary(fit)$coefficients[2,1]+1.96*summary(fit)$coefficients[2,2]
  
  return(c(point,lower,upper))
}

####  overall effect of N on ability to detect beta different from one (Question #1)
Ntotal=seq(25,350,25)
betas=seq(0.5,1,0.05)
Nreps=25

# object for storing output
outputB=array(NA,c(length(Ntotal),length(betas),Nreps))
outputUpper=array(NA,c(length(Ntotal),length(betas),Nreps))
outputLower=array(NA,c(length(Ntotal),length(betas),Nreps))
outputDiff1=array(NA,c(length(Ntotal),length(betas),Nreps))

for(i in 1:length(Ntotal)){
  for(j in 1:length(betas)){
    for(k in 1:Nreps){
      print(c(i,j,k))
      curSim=simCPUE(N=Ntotal[i],beta=betas[j])
      
      outputB[i,j,k]=curSim[1]
      outputLower[i,j,k]=curSim[2]
      outputUpper[i,j,k]=curSim[3]
      outputDiff1[i,j,k]=ifelse(((curSim[2]< 1  & curSim[3] <1) | (curSim[2]>1 & curSim[3]>1)),0,1)
    }
  }
}

# proportion of beta  not different from zero
sum(outputDiff1)/(length(Ntotal)*length(betas)*Nreps)

# proportion of true beta in CIs
sum(outputInCI)/(length(Ntotal)*length(betas)*Nreps)

# ability to detect hyperstability, sum not different from 1
notDiff=apply(outputDiff1,MARGIN=c(1,2),FUN=sum)/Nreps
filled.contour(Ntotal,betas,notDiff,main="prop. not different from 1",xlab="N",ylab="known beta")


# ability to accurately detect beta
library(wesanderson)
paletteN=wes_palette("Zissou1",length(Ntotal),type="continuous")
colN=rep(rep(paletteN,length(betas)),Nreps)

plot(rep(rep(betas,each=length(Ntotal),Nreps)),outputB,ylim=c(0,1.2),xlab="known beta",ylab="Beta",col=colN)
abline(a=0,b=1,lwd=2,col='black')
legend('bottomright',paste("N=",c(25,125,250,350),sep=""),pch=21,col=paletteN[c(1,5,10,14)],box.lty=0)

outputDiff1_N1perlake=data.frame(propDiff=1-rowSums(outputDiff1[1,,])/25,N=Ntotal[1],beta=betas)
for(i in 2:length(Ntotal)){
  outputDiff1_N1perlake=rbind(outputDiff1_N1perlake,
                              data.frame(propDiff=1-rowSums(outputDiff1[i,,])/25,N=Ntotal[i],beta=betas))
}

summary(lm(propDiff~N*beta,data=outputDiff1_N1perlake))

summary(lm(propDiff~N*beta,data=outputDiff1_N1perlake[outputDiff1_N1perlake$beta<1,]))





#### effect of Nlake, given N on ability to detect beta different from 1 (Question 2)
## simulate efCPUE v. angCPUE data w/ repeated lakes
#     -assumes range of efCPUE and angCPUE that cover observed range
#     -q is approximately representative
#     -SD is average of observed
simCPUE_Nlakes<-function(N,beta,Nlakes=N,SD=meanResidSD,minEF=minEFobs,maxEF=maxEFobs,q=0.5){
  EF=runif(N,minEF,maxEF)
  ang=rnorm(N,mean=q*EF^beta,sd=SD)
  ang[ang<0]=q*EF[ang<0]^beta
  
  lakes=paste("lake",rep(1:Nlakes,ceiling(N/Nlakes))[1:N],sep="_")
  
  fit=lmer(log(ang)~log(EF)+(1|lakes))
  
  point=fit@beta[2]
  lower=summary(fit)$coef[2,1]-1.96*summary(fit)$coef[2,2]
  upper=summary(fit)$coef[2,1]+1.96*summary(fit)$coef[2,2]
  
  return(c(point,lower,upper))
}

#  overall effect of N on ability to detect beta different from one
Ntotal=c(25,125,250,350)
betas=seq(0.5,1,0.05)
lakeProp=seq(0.1,0.8,length.out=10)
Nreps=25

# object for storing output
out=list()

out[[1]]=list(outputB=array(NA,c(length(lakeProp),length(betas),Nreps)),
              outputUpper=array(NA,c(length(lakeProp),length(betas),Nreps)),
              outputLower=array(NA,c(length(lakeProp),length(betas),Nreps)),
              outputDiff1=array(NA,c(length(lakeProp),length(betas),Nreps)))
out[[2]]=out[[1]]
out[[3]]=out[[1]]
out[[4]]=out[[1]]


for(i in 1:length(lakeProp)){
  for(j in 1:length(betas)){
    for(k in 1:Nreps){
      print(c(i,j,k))
      curSim1=simCPUE_Nlakes(N=Ntotal[1],Nlakes=ceiling(Ntotal[1]*lakeProp[i]),beta=betas[j])
      out[[1]]$outputB[i,j,k]=curSim1[1]
      out[[1]]$outputLower[i,j,k]=curSim1[2]
      out[[1]]$outputUpper[i,j,k]=curSim1[3]
      out[[1]]$outputDiff1[i,j,k]=ifelse(((curSim1[2]< 1  & curSim1[3] <1) | (curSim1[2]>1 & curSim1[3]>1)),0,1)
      
      curSim2=simCPUE_Nlakes(N=Ntotal[2],Nlakes=ceiling(Ntotal[2]*lakeProp[i]),beta=betas[j])
      out[[2]]$outputB[i,j,k]=curSim2[1]
      out[[2]]$outputLower[i,j,k]=curSim2[2]
      out[[2]]$outputUpper[i,j,k]=curSim2[3]
      out[[2]]$outputDiff1[i,j,k]=ifelse(((curSim2[2]< 1  & curSim2[3] <1) | (curSim2[2]>1 & curSim2[3]>1)),0,1)
      
      curSim3=simCPUE_Nlakes(N=Ntotal[3],Nlakes=ceiling(Ntotal[3]*lakeProp[i]),beta=betas[j])
      out[[3]]$outputB[i,j,k]=curSim3[1]
      out[[3]]$outputLower[i,j,k]=curSim3[2]
      out[[3]]$outputUpper[i,j,k]=curSim3[3]
      out[[3]]$outputDiff1[i,j,k]=ifelse(((curSim3[2]< 1  & curSim3[3] <1) | (curSim3[2]>1 & curSim3[3]>1)),0,1)
      
      curSim4=simCPUE_Nlakes(N=Ntotal[4],Nlakes=ceiling(Ntotal[4]*lakeProp[i]),beta=betas[j])
      out[[4]]$outputB[i,j,k]=curSim4[1]
      out[[4]]$outputLower[i,j,k]=curSim4[2]
      out[[4]]$outputUpper[i,j,k]=curSim4[3]
      out[[4]]$outputDiff1[i,j,k]=ifelse(((curSim4[2]< 1  & curSim4[3] <1) | (curSim4[2]>1 & curSim4[3]>1)),0,1)
    }
  }
}

# lots of singular fits....
#calculate proportion we detect difference from 1
out1_propDiff1summ=matrix(NA,10,11)
out2_propDiff1summ=out1_propDiff1summ
out3_propDiff1summ=out1_propDiff1summ
out4_propDiff1summ=out1_propDiff1summ

for(i in 1:10){
  for(j in 1:11){
    out1_propDiff1summ[i,j]=1-sum(out[[1]]$outputDiff1[i,j,])/25 # proportion different from 1
    out2_propDiff1summ[i,j]=1-sum(out[[2]]$outputDiff1[i,j,])/25
    out3_propDiff1summ[i,j]=1-sum(out[[3]]$outputDiff1[i,j,])/25
    out4_propDiff1summ[i,j]=1-sum(out[[4]]$outputDiff1[i,j,])/25
  }
}

plot(jitter(ceiling(lakeProp*25),1),out1_propDiff1summ[,1],ylim=c(0,1),col=palette[1],xlab="Number of lakes",ylab="Proportion of sims that were hyperstable",main="N=25")
for(i in 2:length(betas)){
  points(jitter(ceiling(lakeProp*25),1),out1_propDiff1summ[,i],col=palette[i])
}
legend('bottomleft',paste("beta=",c(0.5,0.7,0.9,1),sep=""),pch=21,col=palette[c(1,5,9,11)],box.lty=0)

plot(jitter(ceiling(lakeProp*125),1),out2_propDiff1summ[,1],ylim=c(0,1),col=palette[1],xlab="Number of lakes",ylab="Proportion of sims that were hyperstable",main="N=125")
for(i in 2:length(betas)){
  points(jitter(ceiling(lakeProp*125),1),out2_propDiff1summ[,i],col=palette[i])
}
legend('bottomleft',paste("beta=",c(0.5,0.7,0.9,1),sep=""),pch=21,col=palette[c(1,5,9,11)],box.lty=0)

plot(jitter(ceiling(lakeProp*250),1),out3_propDiff1summ[,1],ylim=c(0,1),col=palette[1],xlab="Number of lakes",ylab="Proportion of sims that were hyperstable",main="N=250")
for(i in 2:length(betas)){
  points(jitter(ceiling(lakeProp*250),1),out3_propDiff1summ[,i],col=palette[i])
}
legend('bottomleft',paste("beta=",c(0.5,0.7,0.9,1),sep=""),pch=21,col=palette[c(1,5,9,11)],box.lty=0)

plot(jitter(ceiling(lakeProp*350),1),out4_propDiff1summ[,1],ylim=c(0,1),col=palette[1],xlab="Number of lakes",ylab="Proportion of sims that were hyperstable",main="N=350")
for(i in 2:length(betas)){
  points(jitter(ceiling(lakeProp*350),1),out4_propDiff1summ[,i],col=palette[i])
}
legend('bottomleft',paste("beta=",c(0.5,0.7,0.9,1),sep=""),pch=21,col=palette[c(1,5,9,11)],box.lty=0)

propLakeDIFF=data.frame(propDiff=out1_propDiff1summ[1,],prop=rep(lakeProp[1],11),beta=betas,N=25)
propLakeDIFF=rbind(propLakeDIFF,data.frame(propDiff=out2_propDiff1summ[1,],prop=lakeProp[1],beta=betas,N=125))
propLakeDIFF=rbind(propLakeDIFF,data.frame(propDiff=out3_propDiff1summ[1,],prop=lakeProp[1],beta=betas,N=250))
propLakeDIFF=rbind(propLakeDIFF,data.frame(propDiff=out4_propDiff1summ[1,],prop=lakeProp[1],beta=betas,N=350))
for(i in 2:10){
  propLakeDIFF=rbind(propLakeDIFF,data.frame(propDiff=out1_propDiff1summ[i,],prop=lakeProp[i],beta=betas,N=25))
  propLakeDIFF=rbind(propLakeDIFF,data.frame(propDiff=out2_propDiff1summ[i,],prop=lakeProp[i],beta=betas,N=125))
  propLakeDIFF=rbind(propLakeDIFF,data.frame(propDiff=out3_propDiff1summ[i,],prop=lakeProp[i],beta=betas,N=250))
  propLakeDIFF=rbind(propLakeDIFF,data.frame(propDiff=out4_propDiff1summ[i,],prop=lakeProp[i],beta=betas,N=350))
  
}

summary(lm(propDiff~prop*beta*N,data=propLakeDIFF)) # only N matters; next most important is beta*N interaction...
summary(lm(propDiff~prop*beta*N,data=propLakeDIFF[propLakeDIFF$beta<1,])) # only N matters; next most important is beta*N interaction...


####  Effect of N and difference in betas to detect difference between two species (Question 3)

simCPUE_2species<-function(N,beta1,beta2,SD=meanResidSD,minEF=minEFobs,maxEF=maxEFobs,q=0.5){
  sp=rep(c("sp1","sp2"),each=N)
  
  EF1=runif(N,minEF,maxEF)
  EF2=runif(N,minEF,maxEF)
  
  ang1=rnorm(N,mean=q*EF1^beta1,sd=SD)
  ang1[ang1<0]=q*EF1[ang1<0]^beta1
  
  ang2=rnorm(N,mean=q*EF2^beta2,sd=SD)
  ang2[ang2<0]=q*EF2[ang2<0]^beta2
  
  ang=c(ang1,ang2)
  EF=c(EF1,EF2)
  
  fit=glm(log(ang)~log(EF)*sp)
  
  if(summary(fit)$coefficients[4,4]<0.05){
    return(1)    
  }else{
    return(0)
  }
}

#  overall effect of N on ability to detect beta different from one
Ntotal=seq(12,180,12)
betaDiff=seq(0,0.5,0.05)
Nreps=25

# object for storing output
outputDiff=array(NA,c(length(Ntotal),length(betaDiff),Nreps))

for(i in 1:length(Ntotal)){
  for(j in 1:length(betaDiff)){
    for(k in 1:Nreps){
      print(c(i,j,k))
      curSim=simCPUE_2species(N=Ntotal[i],beta1=0.5,beta2=0.5+betaDiff[j])
      
      outputDiff[i,j,k]=curSim
    }
  }
}

notDiff=matrix(NA,15,11)
for(i in 1:15){
  for(j in 1:11){
    notDiff[i,j]=1-sum(outputDiff[i,j,])/25   #prop not detect species diff
  }
}

filled.contour(Ntotal*2,betaDiff,notDiff,main="prop. didn't detect species difference",xlab="N",ylab="difference in beta")


#### detecting habitat effect... (Question 4)
simCPUEhabitat<-function(N,betaMin=0.1,betaMax,SD=meanResidSD,minEF=minEFobs,maxEF=maxEFobs,q=0.5){
  EF=runif(N,minEF,maxEF)
  
  habitat=runif(N,0,1)
  betaSlope=betaMax-betaMin
  beta=habitat*betaSlope+betaMin
  
  ang=rnorm(N,mean=q*EF^beta,sd=SD)
  ang[ang<0]=q*EF[ang<0]^beta[ang<0]
  
  fit=glm(log(ang)~log(EF)*habitat)
  
  if(summary(fit)$coefficients[4,4]<0.05){
    return(1)
  }else{
    return(0)
  }
}

Ntotal=seq(25,775,50)
betaMaxs=seq(0.1,0.9,0.1)
Nreps=25

# object for storing output
outputDiff=array(NA,c(length(Ntotal),length(betaMaxs),Nreps))

for(i in 1:length(Ntotal)){
  for(j in 1:length(betaMaxs)){
    for(k in 1:Nreps){
      print(c(i,j,k))
      curSim=simCPUEhabitat(N=Ntotal[i],betaMax=betaMaxs[j])
      
      outputDiff[i,j,k]=curSim
    }
  }
}

notDiff=matrix(NA,length(Ntotal),length(betaMaxs))
for(i in 1:length(Ntotal)){
  for(j in 1:length(betaMaxs)){
    notDiff[i,j]=1-sum(outputDiff[i,j,])/25   #prop not detect species diff
  }
}

filled.contour(Ntotal,betaMaxs,notDiff,main="prop. didn't detect linear habitat effect",xlab="N",ylab="difference in beta")


## detecting QUADRATIC habitat effect...
simCPUEhabitat<-function(N,betaMin=0.1,betaMax,SD=meanResidSD,minEF=minEFobs,maxEF=maxEFobs,q=0.5){
  EF=runif(N,minEF,maxEF)
  
  habitat=runif(N,0,1)
  habitat2=habitat*habitat
  
  h=0.5
  k=betaMax
  a=(betaMin-betaMax)/(-h)^2
  
  beta=a*(habitat-h)^2+k
  
  ang=rnorm(N,mean=q*EF^beta,sd=SD)
  ang[ang<0]=q*EF[ang<0]^beta[ang<0]
  
  fit=glm(log(ang)~log(EF)*habitat+log(EF)*habitat2)
  
  if(summary(fit)$coefficients[6,4]<0.05){
    return(1)
  }else{
    return(0)
  }
}

Ntotal=seq(25,775,50)
betaMaxs=seq(0.1,0.9,0.1)
Nreps=25

# object for storing output
outputDiff=array(NA,c(length(Ntotal),length(betaMaxs),Nreps))

for(i in 1:length(Ntotal)){
  for(j in 1:length(betaMaxs)){
    for(k in 1:Nreps){
      print(c(i,j,k))
      curSim=simCPUEhabitat(N=Ntotal[i],betaMax=betaMaxs[j])
      
      outputDiff[i,j,k]=curSim
    }
  }
}

notDiff=matrix(NA,length(Ntotal),length(betaMaxs))
for(i in 1:length(Ntotal)){
  for(j in 1:length(betaMaxs)){
    notDiff[i,j]=1-sum(outputDiff[i,j,])/25   #prop not detect species diff
  }
}

filled.contour(Ntotal,betaMaxs,notDiff,main="prop. didn't detect quadratic habitat effect",xlab="N",ylab="difference in beta")


