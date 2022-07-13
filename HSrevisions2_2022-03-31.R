##### Analyses for Hyperstability survey revision #2  at Fisheries Research
##### 2022-03-31
##### SEJ


####Simulation study for effect of noise in true abundance (PE) vs. efCPUE

## create a function to simulate a efCPUE vs angCPUE dataset
simHS<-function(Nlakes=100,minAbund=15,maxAbund=800,EFq=0.1,ANGq=0.1,efSDperc=0.2,trueBeta=0.4,angEffort=100,plot=FALSE,returnALL=FALSE){
  # randomly generate  true abundances
  trueAbund=runif(Nlakes,min=minAbund,max=maxAbund)

  # generate angler CPUE  based on true abundances; assume poisson, but could do NB...
  angCPUE=numeric(length(trueAbund))
  for(i in 1:length(angCPUE)){
    angCPUE[i]=mean(rpois(angEffort,lambda=ANGq*trueAbund[i]^trueBeta))
  }

  # "EF sample" - assumed proportional, but with noise
  efCPUE=trueAbund*EFq
  efCPUE=rnorm(length(efCPUE),mean=efCPUE,sd=efCPUE*efSDperc)
  efCPUE[efCPUE<1]=NA
  #efCPUE[efCPUE<0]=minAbund*EFq  ### not sure what to do beyond this...


  # plots
  if(plot){
    par(mfrow=c(2,2))
    plot(trueAbund,efCPUE,xlab="true fish abundance",ylab="efCPUE")
    plot(trueAbund,angCPUE,xlab="true fish abundance",ylab="angCPUE")

    plot(efCPUE,angCPUE,xlab="efCPUE",ylab="angCPUE")
  }

  # estimate beta from simulated data
  fit=glm(log(angCPUE)~log(efCPUE))

  pars=coef(fit)

  if(returnALL){
    return(list(trueAbund=trueAbund,angCPUE=angCPUE,efCPUE=efCPUE,pars=pars))  
  }else{
    return(pars)
  }
}

## vary SD of efCPUE systematically  and look at how well beta is estimated at 3 betas
efSDpercents=seq(0.05,1,length.out=20)
reps=50
storeBs_stable=matrix(NA,length(efSDpercents),reps)
storeBs_prop=storeBs_stable
storeBs_depleted=storeBs_stable

for(i in 1:length(efSDpercents)){
  for(j in 1:reps){
    storeBs_stable[i,j]=simHS(efSDperc=efSDpercents[i],trueBeta=0.4)[2]
    storeBs_prop[i,j]=simHS(efSDperc=efSDpercents[i],trueBeta=1)[2]
    storeBs_depleted[i,j]=simHS(efSDperc=efSDpercents[i],trueBeta=1.4)[2]
  }  
}


plot(rep(efSDpercents,reps),storeBs_stable,xlab="electrofishing SD percent",ylab="beta",ylim=c(0,0.5))
points(efSDpercents,rep(0.4,length(efSDpercents)),col='red',pch=15)  
abline(h=1,lty=3)

plot(rep(efSDpercents,reps),storeBs_prop,xlab="electrofishing SD percent",ylab="beta",ylim=c(0,1.2))
points(efSDpercents,rep(1,length(efSDpercents)),col='red',pch=15)  
abline(h=1,lty=3)

plot(rep(efSDpercents,reps),storeBs_depleted,xlab="electrofishing SD percent",ylab="beta",ylim=c(0,1.5))
points(efSDpercents,rep(1.4,length(efSDpercents)),col='red',pch=15)  
abline(h=1,lty=3)  

#### we see biased low (to very low) beta estimates as efSDpercent increases...
lowSD=simHS(efSDperc=0.05,trueBeta=0.4,returnALL=TRUE)
summary(lm(lowSD$trueAbund~lowSD$efCPUE))
lowSD$pars

highSD=simHS(efSDperc=0.75,trueBeta=0.4,returnALL=TRUE)
summary(lm(highSD$trueAbund~highSD$efCPUE))
highSD$pars

plot(lowSD$trueAbund,lowSD$efCPUE,xlim=c(0,800),ylim=c(0,250),xlab="abundance",ylab="efCPUE")
points(highSD$trueAbund,highSD$efCPUE,col='red')

plot(lowSD$efCPUE,lowSD$angCPUE,xlim=c(0,250),ylim=c(0,2),xlab="efCPUE",ylab="angler CPUE")
lines(0:250,exp(lowSD$pars[1])*(0:250)^lowSD$pars[2],lwd=2)
points(highSD$efCPUE,highSD$angCPUE,col='red')
lines(0:250,exp(highSD$pars[1])*(0:250)^highSD$pars[2],lwd=2,col='red')



#########  How to get efSDpercent from our data?
#########  this seems to work... and sd of scaled residuals is the efSDpercent
X=runif(10000,15,800)
f=0.5
q=0.1

efCPUE=rnorm(length(X),q*X,f*q*X)

plot(X,efCPUE)

fit=lm(efCPUE~X)

resids=efCPUE-predict(fit)

residsScaled=resids/predict(fit)

hist(residsScaled)

sd(residsScaled)


### Where are we on sdEFpercent axis?
setwd("~/Documents/Research/People/Students/current/Mosley_Camille/hyperstabilitySurveyManuscript/hyperstabilitySurvey/")
WallPE=read.csv("cleanedWDNRwalleyePE.csv")
hansen=read.csv("cleanedHansenetalLakeinfo.csv")

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

## calculate efSDpercent for walleye areal and shoreline density
sd(residuals(densityFit)/predict(densityFit)) # 0.72

sd(residuals(shoreDensityFit)/predict(shoreDensityFit)) # 0.78

## bass
bassPEfinal=read.csv("cleanedMFEbassPE.csv")

bassEFvPEfitShoreline=lm(bassPEfinal$efCPE~bassPEfinal$densityShoreline)
summary(bassEFvPEfitShoreline)

bassEFvPEfitAreal=lm(bassPEfinal$efCPE~bassPEfinal$densityAreal)
summary(bassEFvPEfitAreal)

sd(residuals(bassEFvPEfitShoreline)/predict(bassEFvPEfitShoreline))  # 0.5
sd(residuals(bassEFvPEfitAreal)/predict(bassEFvPEfitAreal))  # 0.45



######### ********** UPDATE 05/10/2022: this bias problem is actually  well-known in statistics and hyperstability literature...
#### Hanson et al. 2005 adjusts for this bias based on Shardlow et al. 1985

### unbaised beta = observed beta*(sdE^2/(sdX^2-sdE^2)), where sdE is measurement error and sdX is standard deviation of observed values of X

#sdX^2=3.318618 for log transformed across all species...
#sdX^2=293.2603 for untransformed across all species

# I don't really know how to estimate sdE without looking at a bunch of efCPUE from teh same lake at the same abundance...
# try using variance of efCPUE~abundance regression residuals???

#sdE^2=326.6635 for walleye with areal density; lambda = 326.6635/(293.2603-326.6635) = -9.779407...
#sdE^2=416.1411 for walleye with shoreline density

#sdE^2=30.87245 for bass with areal density
#sdE^2=31.97029 for bass with shoreline density



######################
######  power analysis
######################

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

###  overall effect of N on ability to detect beta different from one
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

# effect of N and beta (N effect is subtle)
library(wesanderson)
palette=wes_palette("Zissou1",length(betas),type="continuous")
Bcols=rep(rep(palette,each=length(Ntotal)),Nreps)

plot(rep(Ntotal,length(betas)*Nreps),outputB,ylim=c(0,1.2),xlab="N",ylab="Beta",col=Bcols)
abline(h=1,lty=3,lwd=2)

# effect of beta
paletteN=wes_palette("Zissou1",length(Ntotal),type="continuous")
colN=rep(rep(paletteN,length(betas)),Nreps)

plot(rep(rep(betas,each=length(Ntotal),Nreps)),outputB,ylim=c(0,1.2),xlab="known beta",ylab="Beta",col=colN)
abline(a=0,b=1,lwd=2,col='black')
legend('bottomright',paste("N=",c(25,125,250,350),sep=""),pch=21,col=paletteN[c(1,5,10,14)],box.lty=0)


# sum not different from 1
notDiff=apply(outputDiff1,MARGIN=c(1,2),FUN=sum)/Nreps
filled.contour(Ntotal,betas,notDiff,main="prop. not different from 1",xlab="N",ylab="known beta")

# effect of N
par(mfrow=c(2,2))
toPlot=floor(seq(1,nrow(outputB),length.out=4))

plot(rep(betas,Nreps),outputB[toPlot[1],,],ylim=c(0,1.2),xlab="known beta",ylab="Beta",main=paste("N=",Ntotal[toPlot[1]],sep=""))
points(betas,betas,col='red',pch=15)
abline(h=1,lwd=2,lty=3,col='red')

plot(rep(betas,Nreps),outputB[toPlot[2],,],ylim=c(0,1.2),xlab="known beta",ylab="Beta",main=paste("N=",Ntotal[toPlot[2]],sep=""))
points(betas,betas,col='red',pch=15)
abline(h=1,lwd=2,lty=3,col='red')

plot(rep(betas,Nreps),outputB[toPlot[3],,],ylim=c(0,1.2),xlab="known beta",ylab="Beta",main=paste("N=",Ntotal[toPlot[3]],sep=""))
points(betas,betas,col='red',pch=15)
abline(h=1,lwd=2,lty=3,col='red')

plot(rep(betas,Nreps),outputB[toPlot[4],,],ylim=c(0,1.2),xlab="known beta",ylab="Beta",main=paste("N=",Ntotal[toPlot[4]],sep=""))
points(betas,betas,col='red',pch=15)
abline(h=1,lwd=2,lty=3,col='red')


outputDiff1_N1perlake=data.frame(propDiff=1-rowSums(outputDiff1[1,,])/25,N=Ntotal[1],beta=betas)
for(i in 2:length(Ntotal)){
  outputDiff1_N1perlake=rbind(outputDiff1_N1perlake,
                              data.frame(propDiff=1-rowSums(outputDiff1[i,,])/25,N=Ntotal[i],beta=betas))
}

summary(lm(propDiff~N*beta,data=outputDiff1_N1perlake))

summary(lm(propDiff~N*beta,data=outputDiff1_N1perlake[outputDiff1_N1perlake$beta<1,]))

# do we detect  the correct beta?
plot(Ntotal,outputB)




# effect of Nlake, given N on ability to detect beta different from 1
library(lme4)

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

###  overall effect of N on ability to detect beta different from one
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
plot(rep(lakeProp,length(betas)*Nreps),out[[1]]$outputB)


plot(rep(lakeProp,length(betas)*Nreps),out[[2]]$outputB)

plot(rep(lakeProp,length(betas)*Nreps),out[[3]]$outputB)

plot(rep(lakeProp,length(betas)*Nreps),out[[4]]$outputB)


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


#####  Effect of N and difference in betas to detect difference between two species

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

###  overall effect of N on ability to detect beta different from one
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


###### detecting habitat effect...
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


###### detecting QUADRATIC habitat effect...
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

