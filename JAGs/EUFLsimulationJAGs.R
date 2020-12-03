###EUFL occupancy model run
##

#data

occ.data <- list (y=EH, nSite=nsites, nocc=nocc, x=x1, x2=x2, x3=x3, x4=x4) #occupancy data

#######################################################################
#psi_dot_p_dot
#initial values

z <- apply (EH, 1, max, na.rm=T ) #initialize true state *z*
              #for EH, by row, look at the max value, and ignore NA
occ.init <- function(){
  list(
    b0.psi = runif(1, -3, 3),
    b0.p = runif(1, -3, 3),
    z=z
  )
}

#parameters to track

occ.parm <- c("b0.psi", "b0.p", "N", "mean.psi", "mean.p")

ni <- 8000
nb <- 1000
nt <- 1
nc <- 3

#run the model
occ.result <- R2jags::jags(
              occ.data,
              occ.init,
              occ.parm,
              here::here("jags", "occ_psidot_pdot.txt"),
              n.chains = nc,
              n.iter = ni,
              n.burnin = nb,
              n.thin = nt
              )

##################################################################################
#####psi_x1x2x3x4_p_dot
#initial values

z <- apply (EH, 1, max, na.rm=T ) #initialize true state *z*
#for EH, by row, look at the max value, and ignore NA
occ.init <- function(){
  list(
    b0.psi = runif(1, -3, 3),
    b1.psi = runif(1, -3, 3),
    b2.psi = runif(1, -3, 3),
    b3.psi = runif(1, -3, 3),
    b4.psi = runif(1, -3, 3),
    b0.p = runif(1, -3, 3),
    z=z,
    x1 = runif(1, -3, 3),
    x2 = runif(1, -3, 3),
    x3 = runif(1, -3, 3),
    x4 = runif(1, -3, 3)
    
  )
}

#parameters to track

occ.parm <- c("b0.psi", "b1.psi","b2.psi","b3.psi","b4.psi","b0.p", "N", "mean.psi", 
              "mean.X1.psi", "mean.X2.psi", "mean.X3.psi", "mean.X4.psi", "mean.p")

ni <- 8000
nb <- 1000
nt <- 1
nc <- 3

#run the model
occ.result <- R2jags::jags(
  occ.data,
  occ.init,
  occ.parm,
  here::here("jags", "occ_psi_X1X2X3X4_pdot.txt"),
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  n.thin = nt
)
