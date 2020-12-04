#################
#EUFL simulation

nsites <- 500
nocc <- 8

#covariates
x1 <- rnorm(nsites, 0, 1)
x2 <- rnorm(nsites, 0, 1)
x3 <- rnorm(nsites, 0, 1)
x4 <- rnorm(nsites, 0, 1)

#true occupancy model
Ltrue.psi <- -0.8 - 1*x1 + 0.5*x2 + 1.2*x3 - 0.5*x4

true.psi <-plogis(Ltrue.psi)

true.p <- 0.08

EH <- matrix( 0, nsites, nocc)

for( i in 1:nsites){                  #run through all sites 1-nsites
  if (runif(1) < true.psi){        #if random is less than true.psi then occupied
    for (t in 1:nocc){                #run through all occasions 1-nocc
      EH[i,t] <- rbinom(1, 1, true.p) #for each encounter, binomial 0/1 for detection
    }
  }
}