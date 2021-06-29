
###################################################################################################
#
# Purpose: Generate time-varying data with competing events
#          Estimate RDs using g-computation (with bootstrapping)
#
# Author: Jacqueline Rudolph (Credit to Young and Moodie for DGM)
#
# Last Update: 08 Jun 2021
#
##################################################################################################

# Read in packages
lib <- "~/R/x86_64-pc-linux-gnu-library/4.0"
packages <- c("survival", "nnet", "cmprsk", "tidyverse", "data.table", "parallel")
for (package in packages) {
  library(package, character.only=T, lib.loc=lib)
}

# Pull in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Define parameters and functions
n <- 1000                          # Number of subjects
N <- 6                             # Number of time points per subject
K <- 2                             # Number of causes of death
montecarlo <- 8000                 # Size of Monte Carlo resample
sim_l <- as.numeric(args[1])       # Bottom of range of sims to run
sim_u <- as.numeric(args[2])       # Top of range of sims to run
nboot <- 200                       # Number of bootstrap resamples
lambda <- 0.1                      # Baseline rate of the outcome

expit <- function(x) {1/(1+exp(-x))}
cores <- detectCores()

# Prepare data set to hold simulation results
sim.res <- data.frame(
  outcome=c("Y1", "Y2", "Y1.KM", "Y1.noY2", "Y1.noY2.AJ"),
  r1=rep(NA, 5), 
  r0=rep(NA, 5),
  rd=rep(NA, 5),
  sim=rep(NA, 5),
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

simloop <- function(s) {
  sim.res$sim <- s
  set.seed(s)
  
  # This is the matrix of parameters of interest, possibly different at each interval
  psi.mat <- matrix(0, nrow=K, ncol=N+1)

  # Here are the effect sizes for the K=2 causes
  psi.mat[1, ] <- log(2)    # X results in more Y1
  psi.mat[2, ] <- log(3)    # X results in more Y2 
  
  # Here the (untreated) all-cause rate is set to lambda, with lambda/K per cause.
  gamma.vec <- rep(log(lambda/K), K)
  muK <- sum(exp(gamma.vec)) 
  X<-L<-ID<-Y<-Z<-Tv<-Int<-XLast<-LLast<-LFirst <- numeric()
  T0.vec<-T.vec<-Y.vec<-Z.vec <- rep(0, n)

  # Here are the coefficients determining the mediation and treatment assignment mechanisms.
  bevec <- c(log(3/7), 2, log(0.5), log(1.5)) # Used to generate time-varying confounder
  alvec <- c(log(2/7), 0.5, 0.5, log(4))      # Used to generate exposure

  # cval is used to introduce the confounding
  cval <- 30

  # Begin the data-generation loop
  simulation <- function (exposure) {
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
    
  # Now create simulation to be used in analysis steps
  sim.dat <- simulation(exposure=NULL)


##################################################################################################
## Bootstrap resample to get CIs

  # Set up data set to hold bootstrap results  
  boot.res <- data.frame(
    outcome=c("Y1", "Y2", "Y1.KM", "Y1.noY2", "Y1.noY2.AJ"),
    r1=rep(NA, 5), 
    r0=rep(NA, 5),
    rd=rep(NA, 5),
    boot_num=rep(NA, 5),
    stringsAsFactors = FALSE
  ) 
  
  bootrep <- function(r) {
    boot.res$boot_num <- r
    set.seed(r+1)
    firstobs <- sim.dat[sim.dat$Int==1, ] 
    samp <- table(firstobs[sample(1:nrow(firstobs), nrow(firstobs), replace=T), (names(sim.dat)=="ID")])
    
    # The step below pulls in the simulated data for boot=0; otherwise grabs all records for the resampled observations
    boot <- NULL
    if(r==0){
      boot <- sim.dat %>% 
        rename(bid = ID)
    } else{
      for(zzz in 1:max(samp)){ 
        cc <- sim.dat[sim.dat$ID %in% names(samp[samp %in% c(zzz:max(samp))]),]
        cc$bid<-paste0(cc$ID,zzz)
        boot <- rbind(boot, cc)
      }
      boot <- select(boot, -ID)
    }


##################################################################################################
## Classical g-computation
  # Time-ordering: LLast, XLast, L, X, Z

    # Model confounder
    mod.L <- glm(L ~ XLast + LLast + as.factor(Int), family=binomial(link="logit"), data=boot)

    # Model exposure (if we want natural course)
    mod.X <- glm(X ~ L + XLast + LLast + as.factor(Int), family=binomial(link="logit"), data=boot)
  
    # Model outcomes using multinomial pooled logistic regression
    mod.Z <- multinom(Z ~ X + L + XLast + LLast + as.factor(Int), data=boot)
    
    # Take a Monte Carlo (MC) sample
    # Select first obs for each person to obtain joint empirical distribution of baseline covariates
    MC0 <- boot %>% 
      filter(Int==1) %>% 
      select(X, L)
    index <- sample(1:nrow(MC0), montecarlo, replace=T)
    MC <- MC0[index, ]
    MC$id <- 1:montecarlo
    
    # Predict follow-up based on g-formula using PGF function
      # ii: individual
      # mc_data: MC sample data set
      # length: number of time points
      # exposure: what are we setting exposure to (NULL for natural course)
      # rmCR: are we removing competing risk?
    pgf<-function(ii, mc_data, length, exposure=NULL, rmCR=FALSE){
      pFunc <- function (mod, dat) {
        p <- predict(mod, newdata=dat, type="response")
        var <- rbinom(1, 1, p)
        return(var)
      }
      simZ<-function(mod, dat){
        # Predicted probabilities of each outcome from the multinomial model
        probZ <- predict(mod, newdata=dat, type="probs") 
        if (rmCR==FALSE) {
          # If we are letting the competing event be, then prob of Y1 and Y2 matters
          Yp <- rbinom(1, 1, probZ[2] + probZ[3])
          # If an event occurred, determine if it's Y2. If so Zp=2; o.w. Zp=1
          Zp <- ifelse(Yp==0, 0,
                       rbinom(1, 1, probZ[3]/(probZ[2] + probZ[3])) + 1) 
        } else if (rmCR==TRUE){
          # If we are intervening to remove competing event, then only Y1 matters
          Yp <- rbinom(1, 1, probZ[2])
          # If an event occurred, it's Y1
          Zp <- Yp
        }
        return(Zp)
      }
      
      d <- mc_data %>% 
        filter(id==ii)
      lngth <- length
      Lp<-Xp<-Zp<-Int <- numeric()
      id <- d$id

      Lp[1] <- d$L
      Int[1] <- 1
      
      if (is.null(exposure)) {
        Xp[1] <- d$X
      } else{
        Xp[1] <- exposure
      }
      
      dZp <- data.table(L=Lp[1], LLast=0, XLast=0, X=Xp[1], Int=Int[1])
      Zp[1] <- simZ(mod.Z, dZp)
      
      for (j in 2:lngth) {
        if (Zp[j-1]==0){
          XLast <- Xp[j-1]; LLast <- Lp[j-1]
          Int[j] <- j
          
          dLp <- data.table(LLast,XLast,Int=Int[j])
          Lp[j] <- pFunc(mod.L, dLp)
          
          if (is.null(exposure)) {
            dXp <- data.table(L=Lp[j], LLast, XLast, Int=Int[j])
            Xp[j] <- pFunc(mod.X, dXp)
          } else{
            Xp[j] <- exposure
          }
          
          dZp <- data.table(X=Xp[j], L=Lp[j], LLast, XLast, Int=Int[j])
          Zp[j] <- simZ(mod.Z, dZp)
        } else {
          break
        }
      }
      
      gdat <- data.table(id, Int, Xp, Lp, Zp)
      return(gdat)
    }
    
    # Intervene on exposure, allow competing event to occur
    res0 <- lapply(1:montecarlo, function(x) {pgf(x, mc_data=MC, length=N, exposure=0, rmCR=FALSE)})
    res0 <- do.call(rbind,res0) %>% 
      mutate(last = as.numeric(!duplicated(id, fromLast=T))) %>% 
      filter(last==1) %>% 
      select(-last)

    res1<-lapply(1:montecarlo, function(x) {pgf(x, mc_data=MC, length=N, exposure=1, rmCR=FALSE)})
    res1<-do.call(rbind,res1) %>% 
      mutate(last = as.numeric(!duplicated(id, fromLast=T))) %>% 
      filter(last==1) %>% 
      select(-last)
    
    gcomp.cr <- data.table(rbind(res1, res0))
    
    # Intervene on exposure and eliminate competing event
    res0 <- lapply(1:montecarlo, function(x) {pgf(x, mc_data=MC, length=N, exposure=0, rmCR=TRUE)})
    res0 <- do.call(rbind,res0) %>% 
      mutate(last = as.numeric(!duplicated(id, fromLast=T))) %>% 
      filter(last==1) %>% 
      select(-last)
    
    res1<-lapply(1:montecarlo, function(x) {pgf(x, mc_data=MC, length=N, exposure=1, rmCR=TRUE)})
    res1<-do.call(rbind,res1) %>% 
      mutate(last = as.numeric(!duplicated(id, fromLast=T))) %>% 
      filter(last==1) %>% 
      select(-last)
    
    gcomp.nocr <- data.table(rbind(res1, res0))
    

# Run outcome models
    # Cumulative incidence functions for both events
    fit <- cuminc(ftime=gcomp.cr$Int, fstatus=gcomp.cr$Zp, group=gcomp.cr$Xp)
      #Outcome 1 risk
      boot.res$r0[1] <- fit$`0 1`$est[length(fit$`0 1`$est)]
      boot.res$r1[1] <- fit$`1 1`$est[length(fit$`1 1`$est)]

      #Outcome 2 risk  
      boot.res$r0[2] <- fit$`0 2`$est[length(fit$`0 2`$est)]
      boot.res$r1[2] <- fit$`1 2`$est[length(fit$`1 2`$est)]
      
    # Cumulative distribution function for Outcome 1, censoring Outcome 2
      fit <- summary(survfit(Surv(Int, (Zp==1)) ~ Xp, data=gcomp.cr))
      surv <- data.frame(time = fit$time, 
                         surv = fit$surv,
                         exposure = fit$strata)
      boot.res$r0[3] <- ifelse(is.infinite(1 - min(surv$surv[surv$exposure=="Xp=0"])), 0, 
                               1 - min(surv$surv[surv$exposure=="Xp=0"])) # In case, there are 0 events
      boot.res$r1[3] <- 1 - min(surv$surv[surv$exposure=="Xp=1"])
      
    # Cumulative distribution function for Outcome 1, with  Outcome 2 removed (KM estimator)
    fit <- summary(survfit(Surv(Int, Zp) ~ Xp, data=gcomp.nocr))
      surv <- data.frame(time = fit$time, 
                         surv = fit$surv,
                         exposure = fit$strata)
      boot.res$r0[4] <- ifelse(is.infinite(1 - min(surv$surv[surv$exposure=="Xp=0"])), 0, 
                               1 - min(surv$surv[surv$exposure=="Xp=0"])) # In case, there are 0 events
      boot.res$r1[4] <- 1 - min(surv$surv[surv$exposure=="Xp=1"])
      
    # Cumulative distribution function for Outcome 1, with  Outcome 2 removed (AJ estimator)
    fit <- cuminc(ftime=gcomp.nocr$Int, fstatus=gcomp.nocr$Zp, group=gcomp.nocr$Xp)
      boot.res$r0[5] <- fit$`0 1`$est[length(fit$`0 1`$est)]
      boot.res$r1[5] <- fit$`1 1`$est[length(fit$`1 1`$est)]
    
    # Estimate risk differences
    boot.res <- mutate(boot.res, rd = r1-r0)
      
      
    return(boot.res)
  }

  #all.boot <- mclapply(0:nboot, function(tt) {bootrep(tt)}, mc.cores=cores, mc.set.seed=FALSE)
  all.boot <- lapply(0:nboot, function(tt) {bootrep(tt)})
  all.boot <- do.call(rbind, all.boot)

  # For point estimate, pull out results where boot=0
  boot0 <- filter(all.boot, boot_num == 0)
  sim.res$r0 <- boot0$r0
  sim.res$r1 <- boot0$r1
  sim.res$rd <- boot0$rd

  
##################################################################################################
##Aggregate results

  # Summarize over bootstraps
  boot.summ <- all.boot %>% 
    group_by(outcome) %>% 
    summarize(se = sd(rd))
  
  sim.res <- left_join(sim.res, boot.summ, by="outcome")

  return(sim.res)
}

all.res <- mclapply(sim_l:sim_u, function(x) {simloop(x)}, mc.cores=cores, mc.set.seed=FALSE)
  # Use lapply to test the code if there are errors
  #all.res <- lapply(sim_l:sim_u, function(x) {simloop(x)})
all.res <- do.call(rbind, all.res)

# Summarize over simulations
file_name <- as.character(args[3])
write.table(all.res, file=file_name, sep="\t", row.names=FALSE)
  # If running just 1 simulation
  #write.table(all.res, file="../results/crsim_mainest.txt", sep="\t", row.names=FALSE)

