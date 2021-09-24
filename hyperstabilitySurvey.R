#Cleaned up hyperstability survey code 
#looking at hyperstability of angling CPUE as a function of electrofishing CPUE
#8-19-2021
# CLM, CJD, SEJ, CS
rm(list=ls())

# load function to load data from google drive
library(dplyr)

# load cleaned creel data
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

### Quantifying hyperstability ####

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
library(lme4)

BWJoin=full_join(bassJoin,wallJoin)
BWPJoin=full_join(BWJoin,panJoin)
BWPJoin$WBICfactor=as.factor(BWPJoin$WBIC)
BWP_null<-lmer(logCPUE~logAbun+species+(1|WBICfactor),data=BWPJoin)
BWP_species<-lmer(logCPUE~logAbun+species+logAbun:species+(1|WBICfactor),data=BWPJoin)
anova(BWP_null,BWP_species)

# get standard errors of species-specific betas...
vc=vcov(BWP_species)

var_betas=c(diag(vc)[2],diag(vc)[2]+diag(vc)[8:12]+2*vc[8:12,2])

N_betas=c(length(unique(paste(BCboth$surveyYear,BCboth$WBIC,sep="_"))),
          length(unique(paste(BGboth$surveyYear,BGboth$WBIC,sep="_"))),
          length(unique(paste(LMBboth$surveyYear,LMBboth$WBIC,sep="_"))),
          length(unique(paste(SMBboth$surveyYear,SMBboth$WBIC,sep="_"))),
          length(unique(paste(wallJoin$surveyYear,wallJoin$WBIC,sep="_"))),
          length(unique(paste(YWPboth$surveyYear,YWPboth$WBIC,sep="_"))))

se_betas=sqrt(var_betas)  

params=unlist(as.data.frame(coef(BWP_species)$WBICfactor)[1,])
betas=c(params[2],params[2]+params[8:12])

plot(1:6,betas,col='blue',pch=16,cex=2,ylim=c(0,0.7),xaxt='n',xlab="",ylab=expression(beta))
arrows(1:6,betas-1.96*se_betas,1:6,betas+1.96*se_betas,length=0.05,angle=90,code=3)
axis(side=1,at=1:6,labels=c("black\ncrappie","bluegill","largemouth\nbass","smallmouth\nbass","walleye","yellow\nperch"))

abundPred=seq(0.001,150,0.01)
logAP=log(abundPred)

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="BLACK CRAPPIE"],BWPJoin$meanCPUE[BWPJoin$species=="BLACK CRAPPIE"],xlab="relative abundance",ylab="angler CPUE",main="black crappie",col='darkgrey')
lines(abundPred,exp(params[1]+params[2]*logAP),lwd=2)
text(40,2.5,expression(paste(beta,"=0.196")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="BLUEGILL"],BWPJoin$meanCPUE[BWPJoin$species=="BLUEGILL"],xlab="relative abundance",ylab="angler CPUE",main="bluegill",col='darkgrey')
lines(abundPred,exp(params[1]+params[3]+(params[2]+params[8])*logAP),lwd=2)
text(120,5,expression(paste(beta,"=0.221")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="LARGEMOUTH BASS"],BWPJoin$meanCPUE[BWPJoin$species=="LARGEMOUTH BASS"],xlab="relative abundance",ylab="angler CPUE",main="largemouth bass",col='darkgrey')
lines(abundPred,exp(params[1]+params[4]+(params[2]+params[9])*logAP),lwd=2)
text(50,2.25,expression(paste(beta,"=0.411")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="SMALLMOUTH BASS"],BWPJoin$meanCPUE[BWPJoin$species=="SMALLMOUTH BASS"],xlab="relative abundance",ylab="angler CPUE",main="smallmouth bass",col='darkgrey')
lines(abundPred,exp(params[1]+params[5]+(params[2]+params[10])*logAP),lwd=2)
text(15,0.2,expression(paste(beta,"=0.242")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="WALLEYE"],BWPJoin$meanCPUE[BWPJoin$species=="WALLEYE"],xlab="relative abundance",ylab="angler CPUE",main="walleye",col='darkgrey')
lines(abundPred,exp(params[1]+params[6]+(params[2]+params[11])*logAP),lwd=2)
text(125,1.5,expression(paste(beta,"=0.630")))

plot(BWPJoin$meanEF_CPEkm[BWPJoin$species=="YELLOW PERCH"],BWPJoin$meanCPUE[BWPJoin$species=="YELLOW PERCH"],xlab="relative abundance",ylab="angler CPUE",main="yellow perch",col='darkgrey')
lines(abundPred,exp(params[1]+params[7]+(params[2]+params[12])*logAP),lwd=2)
text(27,2.5,expression(paste(beta,"=0.056")))

### just predicted curves
plot(abundPred,exp(params[1]+params[2]*logAP),type='l',lwd=2,xlab="electrofishing CPUE",ylab="predicted angler CPUE",xlim=c(0,150),ylim=c(0,4),col='purple')
lines(abundPred,exp(params[1]+params[3]+(params[2]+params[8])*logAP),lwd=2,col='blue')
lines(abundPred,exp(params[1]+params[4]+(params[2]+params[9])*logAP),lwd=2,col='darkgreen')
lines(abundPred,exp(params[1]+params[5]+(params[2]+params[10])*logAP),lwd=2,col='brown')
lines(abundPred,exp(params[1]+params[6]+(params[2]+params[11])*logAP),lwd=2,col='red')
lines(abundPred,exp(params[1]+params[7]+(params[2]+params[12])*logAP),lwd=2,col='goldenrod')
legend('topleft',c('walleye','largemouth bass','smallmouth bass','black crappie','bluegill','yellow perch'),lty=1,col=c('red','darkgreen','brown','purple','blue','goldenrod'))

# load lake variables
hansen=read.csv("cleanedHansenetalLakeinfo.csv")

# join lake data with hyperstability data
BWPlc=left_join(BWPJoin, hansen, by="WBIC")

# calculate a couple other lake shape metrics
BWPlc$shorelineComplexity=BWPlc$lakeChar_Length/(2*sqrt(pi*BWPlc$lakeChar_Area)) #ratio of the shoreline length (i.e. perimeter) to the perimeter of an equally sized circle
BWPlc$logLake=log10(BWPlc$lakeChar_Area)

#lake characteristic fits for linear and quadratic model comparisons
BWPlc$RipDev2=BWPlc$riparian_Developed^2
BWPlc$LkArea2=BWPlc$logLake^2
BWPlc$ShoreComp2=BWPlc$shorelineComplexity^2


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

###### 6 species x 3 habitat covariates x 2 linear vs. quadratic = 36 tests
###### Bonferonni means that alpha = 0.05/36 = 0.0014


##### walleye fits
cur=BWPlc[BWPlc$species=="WALLEYE",] #N = 302, overall no sig results 
#riparian developed
dev1_null<-lmer(logCPUE~logAbun+riparian_Developed+(1|WBICfactor),data=cur)
dev1<-lmer(logCPUE~logAbun+riparian_Developed+logAbun:riparian_Developed+(1|WBICfactor),data=cur)
anova(dev1_null,dev1) #no

dev2_null<-lmer(logCPUE~logAbun+riparian_Developed+RipDev2+(1|WBICfactor),data=cur)
dev2<-lmer(logCPUE~logAbun+riparian_Developed+RipDev2+riparian_Developed:logAbun+RipDev2:logAbun+(1|WBICfactor),data=cur)
anova(dev2_null,dev2) #no

#shoreline complexity
shor1_null<-lmer(logCPUE~logAbun+shorelineComplexity+(1|WBICfactor),data=cur)
shor1<-lmer(logCPUE~logAbun+shorelineComplexity+logAbun:shorelineComplexity+(1|WBICfactor),data=cur)
anova(shor1_null,shor1) #no

shor2_null<-lmer(logCPUE~logAbun+shorelineComplexity+ShoreComp2+(1|WBICfactor),data=cur)
shor2<-lmer(logCPUE~logAbun+shorelineComplexity+ShoreComp2+shorelineComplexity:logAbun+ShoreComp2:logAbun+(1|WBICfactor),data=cur)
anova(shor2_null,shor2)#no

#lake area
area1_null<-lmer(logCPUE~logAbun+logLake+(1|WBICfactor),data=cur)
area1<-lmer(logCPUE~logAbun+logLake+logAbun:logLake+(1|WBICfactor),data=cur) 
anova(area1_null,area1)# no

area2_null<-lmer(logCPUE~logAbun+logLake+LkArea2+(1|WBICfactor),data=cur)
area2<-lmer(logCPUE~logAbun+logLake+LkArea2+logLake:logAbun+LkArea2:logAbun+(1|WBICfactor),data=cur)
anova(area2_null,area2) #no

##### yellow perch
cur=BWPlc[BWPlc$species=="YELLOW PERCH",] # N = 88
#riparian developed
dev1_null<-lmer(logCPUE~logAbun+riparian_Developed+(1|WBICfactor),data=cur)
dev1<-lmer(logCPUE~logAbun+riparian_Developed+logAbun:riparian_Developed+(1|WBICfactor),data=cur)
anova(dev1_null,dev1) #no

dev2_null<-lmer(logCPUE~logAbun+riparian_Developed+RipDev2+(1|WBICfactor),data=cur)
dev2<-lmer(logCPUE~logAbun+riparian_Developed+RipDev2+riparian_Developed:logAbun+RipDev2:logAbun+(1|WBICfactor),data=cur)
anova(dev2_null,dev2)#no

#shoreline complexity
shor1_null<-lmer(logCPUE~logAbun+shorelineComplexity+(1|WBICfactor),data=cur)
shor1<-lmer(logCPUE~logAbun+shorelineComplexity+logAbun:shorelineComplexity+(1|WBICfactor),data=cur)
anova(shor1_null,shor1) #no

shor2_null<-lmer(logCPUE~logAbun+shorelineComplexity+ShoreComp2+(1|WBICfactor),data=cur)
shor2<-lmer(logCPUE~logAbun+shorelineComplexity+ShoreComp2+shorelineComplexity:logAbun+ShoreComp2:logAbun+(1|WBICfactor),data=cur)
anova(shor2_null,shor2) #no

#lake area
area1_null<-lmer(logCPUE~logAbun+logLake+(1|WBICfactor),data=cur)
area1<-lmer(logCPUE~logAbun+logLake+logAbun:logLake+(1|WBICfactor),data=cur)
anova(area1_null,area1)#no

area2_null<-lmer(logCPUE~logAbun+logLake+LkArea2+(1|WBICfactor),data=cur)
area2<-lmer(logCPUE~logAbun+logLake+LkArea2+logLake:logAbun+LkArea2:logAbun+(1|WBICfactor),data=cur)
anova(area2_null,area2)#no

##### largemouth bass
cur=BWPlc[BWPlc$species=="LARGEMOUTH BASS",]# N = 172 try lm bc of errors
#riparian developed
dev1_null<-lm(logCPUE~logAbun+riparian_Developed+(1:WBICfactor),data=cur)
dev1<-lm(logCPUE~logAbun+riparian_Developed+logAbun:riparian_Developed+(1:WBICfactor),data=cur)
anova(dev1_null,dev1)#no

dev2_null<-lm(logCPUE~logAbun+riparian_Developed+RipDev2+(1:WBICfactor),data=cur)
dev2<-lm(logCPUE~logAbun+riparian_Developed+RipDev2+riparian_Developed:logAbun+RipDev2:logAbun+(1:WBICfactor),data=cur)
anova(dev2_null,dev2)#no

#shoreline complexity
shor1_null<-lm(logCPUE~logAbun+shorelineComplexity+(1:WBICfactor),data=cur)
shor1<-lm(logCPUE~logAbun+shorelineComplexity+logAbun:shorelineComplexity+(1:WBICfactor),data=cur)
anova(shor1_null,shor1)#no

shor2_null<-lm(logCPUE~logAbun+shorelineComplexity+ShoreComp2+(1:WBICfactor),data=cur)
shor2<-lm(logCPUE~logAbun+shorelineComplexity+ShoreComp2+shorelineComplexity:logAbun+ShoreComp2:logAbun+(1:WBICfactor),data=cur)
anova(shor2_null,shor2)#no

#lake area
area1_null<-lm(logCPUE~logAbun+logLake+(1:WBICfactor),data=cur)
area1<-lm(logCPUE~logAbun+logLake+logAbun:logLake+(1:WBICfactor),data=cur)
anova(area1_null,area1)#no

area2_null<-lm(logCPUE~logAbun+logLake+LkArea2+(1:WBICfactor),data=cur)
area2<-lm(logCPUE~logAbun+logLake+LkArea2+logLake:logAbun+LkArea2:logAbun+(1:WBICfactor),data=cur)
anova(area2_null,area2)#yes 0.0001842

##### smallmouth bass
cur=BWPlc[BWPlc$species=="SMALLMOUTH BASS",]# N = 155
#riparian developed
dev1_null<-lmer(logCPUE~logAbun+riparian_Developed+(1|WBICfactor),data=cur)
dev1<-lmer(logCPUE~logAbun+riparian_Developed+logAbun:riparian_Developed+(1|WBICfactor),data=cur)
anova(dev1_null,dev1)#no

dev2_null<-lmer(logCPUE~logAbun+riparian_Developed+LkArea2+(1|WBICfactor),data=cur)
dev2<-lmer(logCPUE~logAbun+riparian_Developed+LkArea2+riparian_Developed:logAbun+LkArea2:logAbun+(1|WBICfactor),data=cur)
anova(dev2_null,dev2)#no

#shoreline complexity
shor1_null<-lmer(logCPUE~logAbun+shorelineComplexity+(1|WBICfactor),data=cur)
shor1<-lmer(logCPUE~logAbun+shorelineComplexity+logAbun:shorelineComplexity+(1|WBICfactor),data=cur)
anova(shor1_null,shor1)#no p = 0.02598

shor2_null<-lmer(logCPUE~logAbun+shorelineComplexity+ShoreComp2+(1|WBICfactor),data=cur)
shor2<-lmer(logCPUE~logAbun+shorelineComplexity+ShoreComp2+shorelineComplexity:logAbun+ShoreComp2:logAbun+(1|WBICfactor),data=cur)
anova(shor2_null,shor2)# no p = 0.03467
anova(shor1,shor2)#no

#lake area
area1_null<-lmer(logCPUE~logAbun+logLake+(1|WBICfactor),data=cur)
area1<-lmer(logCPUE~logAbun+logLake+logAbun:logLake+(1|WBICfactor),data=cur)
anova(area1_null,area1)#no

area2_null<-lmer(logCPUE~logAbun+logLake+LkArea2+(1|WBICfactor),data=cur)
area2<-lmer(logCPUE~logAbun+logLake+LkArea2+logLake:logAbun+LkArea2:logAbun+(1|WBICfactor),data=cur)
anova(area2_null,area2)#no

##### bluegill
cur=BWPlc[BWPlc$species=="BLUEGILL",]# N = 89
#riparian developed
dev1_null<-lmer(logCPUE~logAbun+riparian_Developed+(1|WBICfactor),data=cur)
dev1<-lmer(logCPUE~logAbun+riparian_Developed+logAbun:riparian_Developed+(1|WBICfactor),data=cur)
anova(dev1_null,dev1)#no

dev2_null<-lmer(logCPUE~logAbun+riparian_Developed+RipDev2+(1|WBICfactor),data=cur)
dev2<-lmer(logCPUE~logAbun+riparian_Developed+RipDev2+riparian_Developed:logAbun+RipDev2:logAbun+(1|WBICfactor),data=cur)
anova(dev2_null,dev2)#no p = 0.03

#shoreline complexity
shor1_null<-lmer(logCPUE~logAbun+shorelineComplexity+(1|WBICfactor),data=cur)
shor1<-lmer(logCPUE~logAbun+shorelineComplexity+logAbun:shorelineComplexity+(1|WBICfactor),data=cur)
anova(shor1_null,shor1)# no

shor2_null<-lmer(logCPUE~logAbun+shorelineComplexity+ShoreComp2+(1|WBICfactor),data=cur)
shor2<-lmer(logCPUE~logAbun+shorelineComplexity+ShoreComp2+shorelineComplexity:logAbun+ShoreComp2:logAbun+(1|WBICfactor),data=cur)
anova(shor2_null,shor2)#no

#lake area
area1_null<-lmer(logCPUE~logAbun+logLake+(1|WBICfactor),data=cur)
area1<-lmer(logCPUE~logAbun+logLake+logAbun:logLake+(1|WBICfactor),data=cur)
anova(area1_null,area1)#no

area2_null<-lmer(logCPUE~logAbun+logLake+LkArea2+(1|WBICfactor),data=cur)
area2<-lmer(logCPUE~logAbun+logLake+LkArea2+logLake:logAbun+LkArea2:logAbun+(1|WBICfactor),data=cur)
anova(area2_null,area2)#no

##### black crappie
cur=BWPlc[BWPlc$species=="BLACK CRAPPIE",] # N = 66
#riparian developed
dev1_null<-lmer(logCPUE~logAbun+riparian_Developed+(1|WBICfactor),data=cur)
dev1<-lmer(logCPUE~logAbun+riparian_Developed+logAbun:riparian_Developed+(1|WBICfactor),data=cur)
anova(dev1_null,dev1)# no p = 0.01

dev2_null<-lmer(logCPUE~logAbun+riparian_Developed+RipDev2+(1|WBICfactor),data=cur)
dev2<-lmer(logCPUE~logAbun+riparian_Developed+RipDev2+riparian_Developed:logAbun+RipDev2:logAbun+(1|WBICfactor),data=cur)
anova(dev2_null,dev2)# very sig. p = 3.902e-06 

#shoreline complexity
shor1_null<-lmer(logCPUE~logAbun+shorelineComplexity+(1|WBICfactor),data=cur)
shor1<-lmer(logCPUE~logAbun+shorelineComplexity+logAbun:shorelineComplexity+(1|WBICfactor),data=cur)
anova(shor1_null,shor1)#no

shor2_null<-lmer(logCPUE~logAbun+shorelineComplexity+ShoreComp2+(1|WBICfactor),data=cur)
shor2<-lmer(logCPUE~logAbun+shorelineComplexity+ShoreComp2+shorelineComplexity:logAbun+ShoreComp2:logAbun+(1|WBICfactor),data=cur)
anova(shor2_null,shor2)#no

#lake area
area1_null<-lmer(logCPUE~logAbun+logLake+(1|WBICfactor),data=cur)
area1<-lmer(logCPUE~logAbun+logLake+logAbun:logLake+(1|WBICfactor),data=cur)
anova(area1_null,area1)# no

area2_null<-lmer(logCPUE~logAbun+logLake+LkArea2+(1|WBICfactor),data=cur)
area2<-lmer(logCPUE~logAbun+logLake+LkArea2+logLake:logAbun+LkArea2:logAbun+(1|WBICfactor),data=cur)
anova(area2_null,area2)# no 


#ploting results for significant habitat relationships with beta

#Black crappie
a=20.0861 #interaction with abund
b=-3.7379 # linear interaction with abun
c=0.2368 #just abun effect
x=seq(0,0.4,0.01)

plot(x,a*x*x+b*x+c,type="l",xlab="riparian developed",ylab="beta", main="crappie - riparian developed") # x axis is variable, y axis is 

-b/(2*a) # x coordinate of vertex

#Largemouth Bass
a=-0.29336  #interaction with abund
b=3.91775 # linear interaction with abun
c=-12.60165 #just abun effect
x=seq(4.5,8,0.5)

plot(x,a*x*x+b*x+c,type="l",xlab="log lake area",ylab="beta", main="LMB - lake area") # x axis is variable, y axis is 

10^(-b/(2*a))/1e6 # x coordinate of vertex in km^2


####  efCPUE vs PEs

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

Wdensity=tapply(cpueCheck_WDNR$Density,cpueCheck_WDNR$WBIC,FUN=mean)
Wshoredensity=tapply(cpueCheck_WDNR$densityShoreline,cpueCheck_WDNR$WBIC,FUN=mean)
WefCPUE=tapply(cpueCheck_WDNR$meanEF_CPEkm,cpueCheck_WDNR$WBIC,FUN=mean)

densityFit=lm(WefCPUE~Wdensity)
shoreDensityFit=lm(WefCPUE~Wshoredensity)

summary(densityFit)
summary(shoreDensityFit)

plot(Wdensity,WefCPUE)
abline(densityFit,lwd=2)
plot(Wshoredensity,WefCPUE)
abline(shoreDensityFit,lwd=2)


##### bass
bassPEfinal=read.csv("cleanedMFEbassPE.csv")

plot(bassPEfinal$densityShoreline~bassPEfinal$efCPE)

bassEFvPEfitShoreline=lm(bassPEfinal$densityShoreline~bassPEfinal$efCPE)
summary(bassEFvPEfitShoreline)
abline(bassEFvPEfitShoreline,lwd=2)


plot(bassPEfinal$densityAreal~bassPEfinal$efCPE)

bassEFvPEfitAreal=lm(bassPEfinal$densityAreal~bassPEfinal$efCPE)
summary(bassEFvPEfitAreal)
abline(bassEFvPEfitAreal,lwd=2)
