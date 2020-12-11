#======================================================================================
#Hierarchical occupancy analysis of EUFL from Avon Park Air Force Range 2019/2020
#Acoustic detectors deployed from 3 to 25 days - "primary" occasions

#            Rob Schorr, Colorado Natural Heritage Program
#             Colorado State University, robert.schorr@colostate.edu
#=======================================================================================
#
###set working directory
getwd()
setwd("C:/1Rob_Documents/EUFL bonneted bats APAFR/EUFL_occupancy/data/")
#setwd("~/GitHub/bonneted_bat_occupancy")

###packages needed
library("R2jags")
library(dplyr)
library(mcmcplots)

### EUFL acoustic detector data
eufl <- read.csv("EUFL_hier_detections_somesitesremoved.csv", header = TRUE )

head(eufl)                                    #look at data
nSites <- length(unique(eufl$site))           #nSites = number of locations
nOcc <- max(eufl$Primary)                     #nOcc = number or primary occasions
max(nOcc)                                     #number of occasions
max(nSites)                                   #number of sites

#====================Temporal covariate data and normalizing=========================
#site-specific nightly covariates

#insect noise levels
noise<- read.csv("EUFL_noise.csv", header = TRUE)
noisenorm <- noise %>% 
  mutate_if (is.numeric, scale)               #scaling (mean = 0) even if there are character fields
summary(noisenorm)
dplyr::setdiff(eufl$site,noisenorm$Site)      #are all the sites in both data sets?
head(noise)
noisenorm <- as.matrix(noisenorm)             #Jags needs it as a matrix

#Brazilian free-tailed bat acoustic recordings
TABR<- read.csv("EUFL_TABRnoise.csv", header = TRUE)
TABRnorm <- TABR %>% 
  mutate_if (is.numeric, scale)
summary(TABRnorm)
dplyr::setdiff(eufl$site,TABRnorm$Site) 
TABRnorm <- as.matrix(TABRnorm)

#minimum nightly temperature
minT<- read.csv("EUFL_Tmin.csv", header = TRUE)
minTnorm <- minT %>% 
  mutate_if (is.numeric, scale)
summary(minTnorm)
dplyr::setdiff(eufl$site,minTnorm$Site) #are all the sites in both data sets?
minTnorm <- as.matrix(minTnorm)

#maximum nightly temperature
maxT<- read.csv("EUFL_Tmax.csv", header = TRUE)
maxTnorm <- maxT %>% 
  mutate_if (is.numeric, scale)
summary(maxTnorm)
dplyr::setdiff(eufl$site,maxTnorm$Site) #are all the sites in both data sets?
maxTnorm <- as.matrix(maxTnorm)

#mean nightly temperature
meanT<- read.csv("EUFL_Tmean.csv", header = TRUE)
meanTnorm <- meanT %>% 
  mutate_if (is.numeric, scale)
summary(meanTnorm)
dplyr::setdiff(eufl$site,meanTnorm$Site) #are all the sites in both data sets?
meanTnorm <- as.matrix(meanTnorm)

#moonphase
scoto<- read.csv("EUFL_scoto.csv", header = TRUE)
scotonorm <- scoto %>% 
  mutate_if (is.numeric, scale)
summary(scotonorm)
dplyr::setdiff(eufl$site,scotonorm$Site) #are all the sites in both data sets?
scotonorm <- as.matrix(scotonorm)

#Julian date - based on Bailey et al. 2017, JMamm 98:1586-1593
julian <- read.csv("Julian_dates.csv", header = TRUE)
juliannorm <- julian %>% 
  mutate_if (is.numeric, scale)
dplyr::setdiff(eufl$site, juliannorm$Site)
juliannorm <- as.matrix(juliannorm)

#Peak activity months based on Braun de Torrez et al. 2020
highmonths <- read.csv("High_activity_months.csv", header = TRUE)
dplyr::setdiff(eufl$site, highmonths$Site)
highmonths <- as.numeric(highmonths)
highmonths <- as.matrix(highmonths)

#==========Site Covariates===========================================================
#data
site_cov <- read.csv("EUFL_site_covariates2.csv", na.strings = "", header = TRUE)
head(site_cov)                                #look at data
dplyr::setdiff(eufl$site,site_cov$Point_ID)  #what sites are found in only 1 dataset?

#dichotomy of which kind of waterbody was closer: marsh or swamp
marsh <- site_cov$Marsh    
swamp <- site_cov$Swamp

#was there a colony nearby Yes/No
RCW <- site_cov$RCWyes     

#habitat categories from 1999 veg mapping effort because this is most commonly used
#habitat data for APAFR
GrassPrairie <- site_cov$GrassPrairie    #grasslands or prairies
SwampMarsh <- site_cov$SwampMarsh        #swamp or marsh
Oak <- site_cov$Oak                      #oak woodlands
Flatwoods<- site_cov$Flatwoods           #pine flatwoods
Plantations <- site_cov$Plantations      #pine plantations
Scrub <- site_cov$Scrub                  #scrublands

#measurements of landscape features from the site of detector deployment
normDSF <- scale(site_cov$DSF)            #days since fire
normDSF <- normDSF[-c(500,501,502,503),] #had to remove 4 extra rows

normNOG <- scale(site_cov$NOG)            #distance to orange grove
normNOG <- normNOG[-c(500,501,502,503),]

normTHM <- scale(site_cov$THmin)          #minimum tree height
normTHM <- normTHM[-c(500,501,502,503),]

normCAM <- scale(site_cov$CAmin)          #minimum crown area
normCAM <- normCAM[-c(500,501,502,503),]

normCAMax <- scale(site_cov$CAmax)        #max crown area
normCAMax <- normCAMax[-c(500,501,502,503),]

normPCC <- scale(site_cov$PCC)            #percent canopy cover 1km
normPCC <- normPCC[-c(500,501,502,503),]

normTCC <- scale(site_cov$TCC)            #total canopy cover 100m
normTCC <- normTCC[-c(500,501,502,503),]

normNWET <- scale(site_cov$NWET)            #distance to nearest wetland
normNWET <- normNWET[-c(500,501,502,503),]

normNRCW <- scale(site_cov$NRCW)            #dist to nearest RCW cluster
normNRCW <- normNRCW[-c(500,501,502,503),]

normANWET <- scale(site_cov$ANWET)            #area of nearest wetland
normANWET <- normANWET[-c(500,501,502,503),]

normNTREE <- scale(site_cov$NTREE)            #dist to nearest tree
normNTREE <- normNTREE[-c(500,501,502,503),]

normHNTREE <- scale(site_cov$HNTREE)            #ht of nearest tree
normHNTREE <- normHNTREE[-c(500,501,502,503),]

normCRMIN <- scale(site_cov$CRmin)            #percent min canopy radius 1km
normCRMIN <- normCRMIN[-c(500,501,502,503),]

normCAMax <- scale(site_cov$CRmax)            #percent min canopy radius 1km
normCAMax <- normCAMax[-c(500,501,502,503),]

normCAMin <- scale(site_cov$CAmin)            #percent min canopy radius 1km
normCAMin <- normCAMin[-c(500,501,502,503),]

#================Encounter history construction=========================================
###encounter history for occupancy
EH <- matrix( 0, nSites, nOcc )       #create encounter history matrix of 0s

eufl$SiteNum <- as.numeric(as.factor(eufl$site)) #making a numeric site num for sites

for( i in 1:nrow( eufl ) ){
  if( !is.na(eufl$Detection[i]) ){
    if( eufl$Detection[i] == 1 ){
      EH[ eufl$SiteNum[i], eufl$Primary[i] ] <- 1
    }
  }
 }


# data - populating EH and providing JAGs key for covariates
nSite <- nrow(EH)
occ.data <- list( y=EH, nSite=nSites, nocc=nOcc, swamp=swamp, 
                  marsh=marsh, GrassPrairie=GrassPrairie, SwampMarsh=SwampMarsh,
                  Oak=Oak, Flatwoods=Flatwoods, Plantations=Plantations,
                  Scrub=Scrub, RCW=RCW, noise=noisenorm, TABR=TABRnorm, 
                  scoto=scotonorm,meanT=meanTnorm, minT=minTnorm, maxT=maxTnorm, 
                  DSF=normDSF, NOG=normNOG, THM=normTHM, CAM=normCAM, 
                  CAMAX=normCAMax, PCC=normPCC, TCC=normTCC,NWET=normNWET, 
                  NRCW=normNRCW, ANWET=normANWET, NTREE=normNTREE, HNTREE=normHNTREE,
                  CRMIN=normCRMIN, CAMIN=normCAMin, julian=juliannorm, highmonths=highmonths)
 
#======================Simple models of detectability holding Psi(.)========================

#######################################
#psi(.)p(.)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    b0.p = runif( 1, -3, 3 ),
    z=z
  )
}

# parameters to track
occ.parm <- c( "b0.psi", "b0.p", "N", "mean.psi", "mean.p" )

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_pdot.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#########################################
#psi(.)p(t9)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    b1.p = runif( 1, -3, 3 ),
    b2.p = runif( 1, -3, 3 ),
    b3.p = runif( 1, -3, 3 ),
    b4.p = runif( 1, -3, 3 ),
    b5.p = runif( 1, -3, 3 ),
    b6.p = runif( 1, -3, 3 ),
    b7.p = runif( 1, -3, 3 ),
    b8.p = runif( 1, -3, 3 ),
    b9.p = runif( 1, -3, 3 ),
   
    z=z
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "b02.p", "b03.p", 
    "b04.p", "b05.p", "b06.p", "b07.p", "b08.p", 
    "N", "mean.psi", 
    "mean.p.t01", "mean.p.t02", "mean.p.t03", "mean.p.t04", 
    "mean.p.t05", "mean.p.t06", "mean.p.t07", "mean.p.t08",
    "mean.p.t09")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_ptime_9.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#########################################
#psi(.)p(swamp)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    z=z
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "N", "mean.psi", "mean.p", "mean.p.swamp" )

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_pswamp.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#########################################
#psi(.)p(marsh)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    z=z
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "N", "mean.psi", "mean.p", "mean.p.marsh" )

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_pmarsh.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#########################################
#psi(.)p(RCW)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    z=z
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "N", "mean.psi", "mean.p", "mean.p.RCW" )

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_pRCW.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

###################################################
#psi(.)p(habitats)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    b05.p = runif( 1, -3, 3 ),
    
    z=z
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "b02.p", "b03.p",
"b04.p","b05.p", "N", "mean.psi", "mean.p",
"mean.p.GrassPrairie", "mean.p.Flatwoods", 
"mean.p.Oak", "mean.p.Plantations", "mean.p.Scrub", "mean.p.SwampMarsh")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_phabitat.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#================site- and time-specific covariates for detectability===========
#########################################
#psi(.)p(insect noise)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),

    z=z,
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "N", 
               "mean.psi", "mean.p", "mu.noise", "sd.noise")

ni <- 500
nb <- 100
nt <- 1
nc <- 3 

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_ptnoise.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#########################################
#psi(.)p(TABR noise) - Brazilian free-tailed bat calls

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "N", 
               "mean.psi", "mean.p", "mu.TABR", "sd.TABR")

ni <- 500
nb <- 100
nt <- 1
nc <- 3 

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_pTABRnoise.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#########################################
#psi(.)p(lowtemp)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.minT = runif(1, -3, 3),
    sd.minT = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "N", 
               "mean.psi", "mean.p", "mu.minT", "sd.minT")

ni <- 500
nb <- 100
nt <- 1
nc <- 3 

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_pTmin.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#########################################
#psi(.)p(hightemp)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "N", 
               "mean.psi", "mean.p", "mu.maxT", "sd.maxT")

ni <- 500
nb <- 100
nt <- 1
nc <- 3 

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_pTmax.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#########################################
#psi(.)p(meantemp)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.meanT = runif(1, -3, 3),
    sd.meanT = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "N", 
               "mean.psi", "mean.p", "mu.meanT", "sd.meanT")

ni <- 500
nb <- 100
nt <- 1
nc <- 3 

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_pTmean.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#########################################
#psi(.)p(scotophase) - moonphase

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.meanT = runif(1, -3, 3),
    sd.meanT = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "N", 
               "mean.psi", "mean.p", "mu.scoto", "sd.scoto")

ni <- 500
nb <- 100
nt <- 1
nc <- 3 

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_pscoto.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#########################################
#psi(.)p(Julian) - based on Bailey et al. 2017

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.julian = runif(1, -3, 3),
    sd.julian = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "N", 
               "mean.psi", "mean.p", "mu.julian", "sd.julian")

ni <- 500
nb <- 100
nt <- 1
nc <- 3 

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_pJulian.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#########################################
#psi(.)p(high activity months = Feb, July, Aug) - based on Braun de Torrez et al. 2020

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b0.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.highactivity = runif(1, -3, 3),
    sd.highactivity = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b0.psi", "b00.p", "b01.p", "N", 
               "mean.psi", "mean.p", "mu.highactivity", "sd.highactivity")

ni <- 500
nb <- 100
nt <- 1
nc <- 3 

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_pHighActivity.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#===========combining detectability covariates with merit===============================
###############################################################
#psi(.)p(habitats+enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
 
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    b05.p = runif( 1, -3, 3 ),
    b06.p = runif( 1, -3, 3 ),
    b07.p = runif( 1, -3, 3 ),
    b08.p = runif( 1, -3, 3 ),
   
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "b05.p","b06.p", "b07.p", "b08.p",
               "N",   
               "mean.p.TABR", "mean.p.Flatwoods", "mean.p.GrassPrairie",  
               "mean.p.maxT", "mean.p.noise", "mean.p.Oak", "mean.p.Plantations", 
                "mean.p.Scrub", "mean.p.SwampMarsh",
               "mean.psi")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psidot_phabitat+enviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

#====================Modeling occupancy===============================================
###############################################################
#psi(RCW)p(habitats+enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    b05.p = runif( 1, -3, 3 ),
    b06.p = runif( 1, -3, 3 ),
    b07.p = runif( 1, -3, 3 ),
    b08.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "b05.p","b06.p", "b07.p", "b08.p",
               "N",   
               "mean.p.TABR", "mean.p.Flatwoods", "mean.p.GrassPrairie",  
               "mean.p.maxT", "mean.p.noise", "mean.p.Oak", "mean.p.Plantations", 
               "mean.p.Scrub", "mean.p.SwampMarsh",
               "mean.psi", "mean.psi.RCW")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psiRCW_phabitat+enviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)


###############################################################
#psi(days since fire)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.DSF")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_days_since_fire_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

###############################################################
#psi(dist to orange grove)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.NOG")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_dist_orangegrove_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

###############################################################
#psi(min tree height)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.THM")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_min_tree_height_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

################################################################
#psi(RCW)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.RCW")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_RCW_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

################################################################
#psi(min crown area)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.CAM")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_min_canopy_area_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

################################################################
#psi(max crown area)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.CAMax")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_max_crown_area_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

###############################################################
#psi(habitats)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    b02.psi = runif( 1, -3, 3 ),
    b03.psi = runif( 1, -3, 3 ),
    b04.psi = runif( 1, -3, 3 ),
    b05.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b02.psi", "b03.psi", "b04.psi", "b05.psi", 
               "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", "mean.p.noise", 
               "mean.psi.Flatwoods", "mean.psi.GrassPrairie", "mean.psi.Oak", 
               "mean.psi.Plantations", "mean.psi.Scrub",  "mean.psi.SwampMarsh")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_habitats_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

################################################################
#psi(percent canopy cover over 1km)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5),
    mu.PCC = runif( 1, -3, 3 ),
    sd.PCC = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.PCC")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_percent_canopy_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

################################################################
#psi(canopy cover over 100m radius)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.TCC")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_total_canopy_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

################################################################
#psi(dist to nearest wetland)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.NWET")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_dist_wetland_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)


################################################################
#psi(dist to nearest RCW cluster)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.NRCW")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_dist_RCW_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

################################################################
#psi(area of nearest wetland)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.ANWET")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_area_wetland_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

###############################################################
#psi(dist to nearest tree)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.NTREE")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_dist_tree_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

###############################################################
#psi(height of nearest tree)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.HNTREE")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_ht_tree_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

###############################################################
#psi(minimum canopy radius 100km)p(enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    
    b00.p = runif( 1, -3, 3 ),
    b01.p = runif( 1, -3, 3 ),
    b02.p = runif( 1, -3, 3 ),
    b03.p = runif( 1, -3, 3 ),
    b04.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.scoto = runif(1, -3, 3),
    sd.scoto = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5),
    mu.CRMIN = runif(1, -3, 3),
    sd.CRMIN = runif(1, 0.5, 1.5)
    
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b00.p", "b01.p", "b02.p", "b03.p", "b04.p",
               "N",   
               "mean.p.TABR", "mean.p.scoto", "mean.p.maxT", 
               "mean.p.noise", 
               "mean.psi", "mean.psi.CRMIN")

ni <- 500
nb <- 100
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_min_radius_penviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

###############################################################
#psi(full)p(habitats+enviro)

# initial values
z <- apply( EH, 1, max, na.rm=T )
occ.init <- function(){
  list(
    b00.psi = runif( 1, -3, 3 ),
    b01.psi = runif( 1, -3, 3 ),
    b02.psi = runif( 1, -3, 3 ),
    b03.psi = runif( 1, -3, 3 ),
    b04.psi = runif( 1, -3, 3 ),
    b05.psi = runif( 1, -3, 3 ),
    b06.psi = runif( 1, -3, 3 ),
    b07.psi = runif( 1, -3, 3 ),
    b08.psi = runif( 1, -3, 3 ),
    b09.psi = runif( 1, -3, 3 ),
    b10.psi = runif( 1, -3, 3 ),
    b11.psi = runif( 1, -3, 3 ),
    b12.psi = runif( 1, -3, 3 ),
    b13.psi = runif( 1, -3, 3 ),
    b14.psi = runif( 1, -3, 3 ),
    b15.psi = runif( 1, -3, 3 ),
    b16.psi = runif( 1, -3, 3 ),
    b17.psi = runif( 1, -3, 3 ),
    
    b0.p = runif( 1, -3, 3 ),
    b1.p = runif( 1, -3, 3 ),
    b2.p = runif( 1, -3, 3 ),
    b3.p = runif( 1, -3, 3 ),
    b4.p = runif( 1, -3, 3 ),
    b5.p = runif( 1, -3, 3 ),
    b6.p = runif( 1, -3, 3 ),
    b7.p = runif( 1, -3, 3 ),
    b8.p = runif( 1, -3, 3 ),
    b9.p = runif( 1, -3, 3 ),
    
    z=z,
    mu.TABR = runif(1, -3, 3),
    sd.TABR = runif(1, 0.5, 1.5),
    mu.maxT = runif(1, -3, 3),
    sd.maxT = runif(1, 0.5, 1.5),
    mu.noise = runif(1, -3, 3),
    sd.noise = runif(1, 0.5, 1.5),
    mu.julian = runif(1, -3, 3),
    sd.julian = runif(1, 0.5, 1.5)
  )
}

# parameters to track

occ.parm <- c( "b00.psi", "b01.psi", "b02.psi", "b03.psi", "b04.psi", "b05.psi",
               "b06.psi", "b07.psi", "b08.psi", "b09.psi", "b10.psi", "b11.psi",
               "b12.psi", "b13.psi", "b14.psi", "b15.psi", "b16.psi", "b17.psi",
               "b0.p", "b1.p", "b2.p", "b3.p", "b4.p",
               "b5.p","b6.p", "b7.p", "b8.p", "b9.p",
               "N",   
               "mean.psi.ANWET", "mean.psi.CAMAX", "mean.psi.CRMIN", "mean.psi.DSF", 
               "mean.psi.Flatwoods", "mean.psi.GrassPrairie",  "mean.psi.HNTREE", 
               "mean.psi.NOG", "mean.psi.NRCW", "mean.psi.NTREE","mean.psi.NWET",  
               "mean.psi.Oak", "mean.psi.PCC", "mean.psi.Plantations", 
               "mean.psi.RCW", "mean.psi.Scrub", "mean.psi.SwampMarsh", "mean.psi.TCC", 
               
               "mean.p.TABR", "mean.p.Flatwoods", "mean.p.GrassPrairie", "mean.p.julian", 
               "mean.p.maxT", "mean.p.noise", "mean.p.Oak", "mean.p.Plantations", 
               "mean.p.Scrub", "mean.p.SwampMarsh")

ni <- 10000
nb <- 2000
nt <- 1
nc <- 3

# run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here( "jags", "occ_psi_full_phabitat+enviro.txt" ),
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

