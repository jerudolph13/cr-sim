
###################################################################################################
#
# Purpose: Generate time-varying data with competing events
#          Estimate counterfactual contrasts of interest in competing events analysis
#
# Author: Jacqueline Rudolph (Credit to Young and Moodie for DGM)
#
# Last Update: 01 Mar 2021
#
##################################################################################################

# Read in packages
packages <- c("survival", "nnet", "cmprsk", "tidyverse", "data.table", "parallel")
for (package in packages) {
  library(package, character.only=T)
}

# Define parameters and functions
n <- 1e6           # Number of subjects
N <- 6             # Number of time points per subject
lambda <- 0.1      # Baseline rate of the outcome

expit <- function(x) {1/(1+exp(-x))}

# Prepare data set to hold simulation results
truth <- data.frame(
  outcome=c("Y1", "Y2", "Y1.noY2"),
  r1=rep(NA, 3), 
  r0=rep(NA, 3),
  rd=rep(NA, 3),
  stringsAsFactors = FALSE
)


##################################################################################################
## Data generation

## This code generates data from a structural nested model 
## compatible with a marginal structural model.
## It is based on Jessica Young's algorithm, published here:
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3635680/
##
## This code, which extends Young's algorithm to multiple outcomes
## was written by Erica Moodie, published here:
## https://www.ncbi.nlm.nih.gov/pubmed/24272681

# Begin the data-generation loop
simulation <- function (exposure, K) {

  # This is the matrix of parameters of interest, possibly different at each interval
  psi.mat <- matrix(0, nrow=K, ncol=N+1)

  # Here are the effect sizes for the K=2 causes
  psi.mat[1, ] <- log(2)    # X results in more Y1
  if (K==2) {
    psi.mat[2, ] <- log(3)  # X results in more Y2 
  }

  # Here the (untreated) all-cause rate is set to lambda, with lambda/K per cause.
  gamma.vec <- rep(log(lambda/K))
  muK <- sum(exp(gamma.vec)) 
  X<-L<-ID<-Y<-Z<-Tv<-Int<-XLast<-LLast<-LFirst <- numeric()
  T0.vec<-T.vec<-Y.vec<-Z.vec <- rep(0, n)

  # Here are the coefficients determining the mediation and treatment assignment mechanisms.
  bevec <- c(log(3/7), 2, log(0.5), log(1.5)) # Used to generate time-varying confounder
  alvec <- c(log(2/7), 0.5, 0.5, log(4))      # Used to generate exposure

  # cval is used to introduce the confounding
  cval <- 30

    for (i in 1:n) {
      # Generate the counterfactual (untreated) survival time
      T0 <- rexp(1, lambda) # Generate T0 from an exponential dist with constant rate=lambda
      Ival <- as.numeric(T0 < cval)
      
      # Begin the interval-by-interval simulation
      m <- 0
      mu.tot <- 0
      X.vec<-L.vec<-XLast.vec<-LLast.vec<-LFirst.vec <- rep(0, N+1)
      
      # Implement Young's algorithm with multiple causes
      # Generate the survival time, then the cause
      while (muK*T0 > mu.tot & m <= N) {
        if (m == 0) {
          # First interval
          eta <- bevec[1] + bevec[2]*Ival + bevec[3]*0 + bevec[4]*0
          pval <- expit(eta)
          L.vec[m+1] <- rbinom(1, 1, pval)
          eta <- alvec[1] + alvec[2]*L.vec[m+1] + alvec[3]*0 + alvec[4]*0
          pval <- expit(eta)
          if (is.null(exposure)) {X.vec[m+1] <- rbinom(1, 1, pval)}
          else {X.vec[m+1] <- exposure}
    
          XLast.vec[m+1] <- 0; LLast.vec[m+1] <- 0 
          LFirst.vec <- rep(L.vec[m+1], N + 1) 
        } else {
          # Subsequent intervals
          eta <- bevec[1] + bevec[2]*Ival + bevec[3]*X.vec[m] + bevec[4]*L.vec[m]
          pval <- expit(eta)
          L.vec[m+1] <- rbinom(1, 1, pval)
          eta <- alvec[1] + alvec[2]*L.vec[m + 1]+alvec[3]*L.vec[m] + alvec[4]*X.vec[m]
          pval <- expit(eta)
          if (is.null(exposure)) {X.vec[m+1] <- rbinom(1, 1, pval)}
          else {X.vec[m+1] <- exposure}
          XLast.vec[m+1] <- X.vec[m]; LLast.vec[m+1] <- L.vec[m]
        }
        
        muval <- sum(exp(gamma.vec + X.vec[m+1]*psi.mat[ , m+1]))
        
        # Tval is computed for each interval, but is overwritten until the final interval
        Tval <- m + (muK*T0 - mu.tot)/muval
        mu.tot <- mu.tot + muval
        m <- m + 1
      }
      
      # After exiting the loop, the survival time has been generated as Tval
      # Now need to generate the failure type.
      if (m > N) {
        # In the case of censoring at the Nth interval, no failure.
        Tval <- m - 1
        Z.vec[i] <- 0
      } else {
        # In the case of failure, use the ratio hazards to define the relevant multinomial distribution on the K causes.
        Z.vec[i] <- sample(c(1:K), 1, prob = exp(gamma.vec + X.vec[m]*psi.mat[ , m]))
      }
      
      # Store the outcomes
      T0.vec[i] <- T0
      T.vec[i] <- Tval
      Y.vec[i] <- m - 1
      ID <- c(ID, rep(i,m)) # Individual
      Int <- c(Int, c(1:m)) # Time point
      X <- c(X, X.vec[1:m]) # Time-updated treatment
      L <- c(L, L.vec[1:m]) # Time-updated covariate
      XLast <- c(XLast, XLast.vec[1:m]) # Treatment at last t point
      LLast <- c(LLast, LLast.vec[1:m]) # Covariate at last t point
      LFirst <- c(LFirst, LFirst.vec[1:m]) # Baseline covariate value
      Z <- c(Z, rep(0, m - 1), Z.vec[i]) # Outcome: Z>0 indicates outcome of some type, determined by value)
      tv <- c(1:m); tv[m] <- Tval
      Tv <- c(Tv, tv) # If event occurs, exact time at which it occurred; o.w. equal to Int)
    }
    
    DeathsK.df <- data.frame(ID, Int, Tv, X, XLast, L, LLast, LFirst, Z)
    
    # Trim off the intervals beyond the Nth (loop goes one too far)
    DeathsK.df <- DeathsK.df[DeathsK.df$Int <= N, ]
    DeathsK.df$Z.f <- factor(DeathsK.df$Z)
    DeathsK.df$D <- DeathsK.df$Z > 0
    DeathsK.df$Int0 <- DeathsK.df$Int - 1
    
    return(DeathsK.df)
  }
  
# Intervene on exposure, allow competing risk

# Exposure = 1
set.seed(123)
sim.dat11 <- simulation(exposure=1, K=2) %>% 
  mutate(last = as.numeric(!duplicated(ID, fromLast=T))) %>% 
  filter(last==1) %>% 
  select(-last)
fit <- cuminc(ftime=sim.dat11$Tv, fstatus=sim.dat11$Z)
  truth$r1[1] <- fit$`1 1`$est[length(fit$`1 1`$est)]
  truth$r1[2] <- fit$`1 2`$est[length(fit$`1 2`$est)]

# Exposure = 0
set.seed(123)
sim.dat01 <- simulation(exposure=0, K=2) %>% 
  mutate(last = as.numeric(!duplicated(ID, fromLast=T))) %>% 
  filter(last==1) %>% 
  select(-last)
fit <- cuminc(ftime=sim.dat01$Tv, fstatus=sim.dat01$Z)
  truth$r0[1] <- fit$`1 1`$est[length(fit$`1 1`$est)]
  truth$r0[2] <- fit$`1 2`$est[length(fit$`1 2`$est)]
    
# Intervene on exposure and remove the competing risk
  
# Exposure = 1
set.seed(123)
sim.dat10 <- simulation(exposure=1, K=1) %>% 
  mutate(last = as.numeric(!duplicated(ID, fromLast=T))) %>% 
  filter(last==1) %>% 
  select(-last)
fit <- survfit(Surv(Tv, Z) ~ 1, data=sim.dat10)
  sfit <- data.frame(time = fit$time, surv = fit$surv)
  truth$r1[3] <- 1 - min(sfit$surv)

# Exposure = 0  
set.seed(123)
sim.dat00 <- simulation(exposure=0, K=1) %>% 
  mutate(last = as.numeric(!duplicated(ID, fromLast=T))) %>% 
  filter(last==1) %>% 
  select(-last)
fit <- survfit(Surv(Tv, Z) ~ 1, data=sim.dat00)
  sfit <- data.frame(time = fit$time, surv = fit$surv)
  truth$r0[3] <- 1 - min(sfit$surv)

# Risk difference
truth$rd <- truth$r1 - truth$r0

# Output results
write.table(truth, file="../results/cr-sim_truth.txt", sep="\t", row.names=FALSE)

