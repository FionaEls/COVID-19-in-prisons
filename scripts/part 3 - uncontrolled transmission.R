##################################################################################
##################################################################################
# an R script to estimate R0 using an SEIR model
#
# Project: COVID in Prisons, ICI3D MMED 2023
#
##################################################################################

# Load packages -----------------------------------------------------------

require(boot); require(deSolve); require(ellipse); require(tidyverse); require(data.table); require(cowplot)

# Parameters --------------------------------------------------------------

## Parameters for SEIR model
## to be passed as default parameters to disease_params()
latentPeriod <- 3        # days (from literature)
infectiousPeriod <- 7    # days (from literature)

# Functions ---------------------------------------------------------------

## Function that makes a list of disease parameters with default values ##
## Default values are defined above in the Parameters section

disease_params <- function(R0,                         # basic reproduction number, must be specified
                           gamma = 1/infectiousPeriod, # rate of removal
                           sigma = 1/latentPeriod)     # rate of progression
{
  return(as.list(environment()))
}

## Function that makes a vector of initial state variable conditions ##
## N0 (initial population value) as argument

init_cond <- function(N0) {
  init <- c(S = N0 - 1,  # susceptible
            E = 1,       # exposed
            I = 0,       # infected
            X = 0,       # removed
            C = 0)       # cumulative cases
  
  return(init)
}

## Function defining the SEIR (SEIXC) model ##
## 'X' is used instead of 'R' to avoid confusion with R, the reproduction number
## 'C' is an added compartment that calculates the cumulative infected cases per
## time step The difference between 'C' at time t and t-1 will be the incidence
##
## Overall, this model is a standard deterministic SEIR model with no births,
## deaths, introductions, replenishment of susceptibles, etc.
seixc <- function(tt, yy, parms) with(c(parms,as.list(yy)), {
  
  ## State variables are: S, E, I, X
  N <- S + E + I + X                ## total population
  beta <- R0*gamma                  ## Infectious contact rate
  forceOfInfection <- beta*I/N   
  
  ## state variable derivatives (ODE system)
  deriv <- rep(NA,5)
  deriv[1] <-	-forceOfInfection*S                    # dSdt
  deriv[2] <-	forceOfInfection*S - sigma*E           # dEdt
  deriv[3] <-	sigma*E - gamma*I                      # dIdt
  deriv[4] <-	gamma*I                                # dXdt
  deriv[5] <- sigma*E                                # cumulative cases
  return(list(deriv))
})

## Function to run the deterministic model simulation, based on the ODE system defined in seixc() ##
## This function outputs a dataframe with the following variables:
## time - timestep
## S, E, I, X, C - number of individuals in each compartment per timestep
## N - total population
## incidence - the number of incident cases, calculated as the difference in C over each time step
## inc_rate - the incidence rate (incidence/N) per time step
simEpidemic <- function(init, tseq, modFunction=seixc, parms = disease_params()) {
  modDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  modDat$N <- rowSums(modDat[, c('S','E','I','X')])  # total outbreak population
  modDat$incidence <- c(0,diff(modDat$C))           # incident cases (difference in S over each time step)
  modDat$inc_rate <- with(modDat, incidence/N)       # incidence rate based on total population N
  return(modDat)
}

## Function to create simulated data using the model ##
## From a simulated epidemic, measure the incidence proportion at several time points by drawing
## cross-sectional samples of individuals at each time, testing them, and then calculating sample
## incidence proportion and associated binomial confidence intervals. The incidence per time
## step will then be calculated from the difference in exposed between each time step.
##
## This function be used to generate simulated data for testing the model and optimization functions.

sampleEpidemic <- function(modDat                                      # Simulated "data" which we treat as real 
                           , sampleDates = seq(1, 365, by = 5)         # Sample every 5 days for 1 year
                           , N = rep(1000, length(sampleDates))     # Number of individuals sampled at each time point
){
  inc_rate_at_sample_times <- modDat[modDat$time %in% sampleDates, 'inc_rate']
  inc <- rbinom(length(N), N, inc_rate_at_sample_times)
  lci <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = inc, n = N)
  uci <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = inc, n = N)    
  return(data.frame(time = sampleDates, incidence = inc, N, samp_inc_rate =  inc/N,
                    lci = lci, uci = uci))
}

## Function to calculate negative log-likelihood ## 
## This will be passed to the objective function, which will be passed to optim()
nllikelihood <- function(parms, obsDat, init, tseq) {
  modDat <- simEpidemic(init, parms=parms, tseq)
  ## What are the rows from our simulation at which we have observed data?
  matchedTimes <- modDat$time %in% obsDat$time
  nlls <- -dbinom(obsDat$incidence,                     # incident cases of observed data
                  obsDat$N,                             # total population at timepoints of observed data
                  prob = modDat$inc_rate[matchedTimes], # incidence rate from model data
                  log = T)
  return(sum(nlls))
}

## Function to combine guess (fit parameters) and fixed parameters ##
## This will be used within the objective function for calculating the least
## squares statistic.
subsParms <- function(fit.params, fixed.params=disease_params())
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]        
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })

## Objective function to be passed to optim() for parameter estimation ##
## Make the objective function a function of fixed and fitted parameters
objFXN <- function(fit.params                        # paramters to fit
                   , fixed.params = disease_params() # fixed paramters
                   , obsDat                          # observed data
                   , init
                   , tseq) {                 
  parms <- subsParms(fit.params, fixed.params)
  nllikelihood(parms, obsDat = obsDat, init=init, tseq = tseq)  # then nllikelihood
}

# Test the functions with simulated data ----------------------------------

## Generate simulated data by sampling from model-generated data
trueParms <- disease_params(R0 = 1.5)                  # Default model parameters
init.pop <- init_cond(N0 = 1000)               # initial population of 1000 people
time.out <- seq(0, 365, 1)                     # range of timepoints for model time series

modDat <- simEpidemic(init.pop, parms = trueParms, tseq = time.out) # Simulated epidemic (underlying process)
ggplot(modDat, aes(x = time, y = incidence)) + geom_col()           # quick look at the model incidence

set.seed(1)                                           # Initiate the random number generator
myDat <- sampleEpidemic(modDat)                       # Simulate data from the sampling process (substitute to REAL DATA)
ggplot(myDat, aes(x=time, y=incidence)) + geom_col()  # quick look at the simulated incidence

## Use optim to estimate R0
## Use only SANN for now, for testing purposes
optim.vals <- optim(par = c(log_R0 = log(1))
                    , objFXN
                    , fixed.params = disease_params()  # defined at top of script
                    , obsDat = myDat
                    , init = init.pop
                    , tseq = time.out
                    , control = list(trace = 3, maxit = 150)
                    , method = "SANN")


## View the R0 estimate
exp(unname(optim.vals$par))  # estimate from simulated data
trueParms['R0']              # true value
## The true R0 and estimated R0 are close in value

## Simulate an epidemic using the R0 estimate (fitDat)
fitDat <- simEpidemic(init = init.pop, tseq = time.out, parms = subsParms(optim.vals$par, trueParms))

## Visually compare incidence rates of fitDat and modDat (the true underlying process)
df <- rbind(modDat = modDat %>% 
              mutate(variable = "modDat"),
            fitDat = fitDat %>% 
              mutate(variable = "fitDat"))

myDat <- myDat %>% 
  mutate(variable = "myDat")

ggplot(df, aes(x = time, y = inc_rate, color = variable))+
  geom_line() +
  geom_point(data = myDat, mapping=aes(x = time, y=samp_inc_rate, color = variable)) +
  #  geom_ribbon(data = myDat, mapping=aes(x = time, y = samp_inc_rate, ymin=lci, ymax=uci), color = NA, fill = "black", alpha = 0.2) +
  labs(title = "Comparing incidence rates of true and fitted models",
       x= "Time [days]",
       y= "Incidence rate") +
  scale_color_manual(labels = c("Fitted","Truth","Data"), values = c("red", "blue", "black"))
#  theme(text = element_text(size = 18))

## Visually compare incidence (cases) of fitDat and modDat (the true underlying process)

ggplot(df, aes(x = time, y = incidence, color = variable))+
  geom_line() +
  geom_point(data = myDat, mapping=aes(x = time, y=incidence, color = variable)) +
  labs(title = "Comparing incidence of true and fitted models",
       x= "Time [days]",
       y= "Cases") +
  scale_color_manual(labels = c("Fitted","Truth","Data"), values = c("red", "blue", "black"))

## As expected, the output of the above two plots are about the same

# Estimate R0 with actual data --------------------------------------------

## Import and wrangle dataset
all.data <- read_csv("new_cases_per_day.csv")

## Extract one outbreak: Avenal State Prison (ASP) Outbreak 1
fixed.N0 <- 4286  # initial population of prison
ASP.1.data <- all.data %>%
  filter(Facility_Outbreak == "Avenal State Prison (ASP) Outbreak 1") %>%
  select(Date, Daily_New_Cases) %>% 
  rename(incidence = Daily_New_Cases) %>% 
  mutate(inc_rate = incidence/fixed.N0) %>% 
  mutate(time = as.integer(-((Date[1]-1) - Date))) %>%
  mutate(N = rep(fixed.N0, nrow(.))) %>% 
  filter(time < 100)

ASP.1.data

## Define time sequence and initial conditions
time.out2 <- seq(0, 99, 1)
init.prison <- init_cond(N0 = fixed.N0)

## Use optim to estimate R0
prison.estim <- optim(par = c(log_R0 = log(1))
                      , objFXN
                      , fixed.params = disease_params()  # defined at top of script
                      , obsDat = ASP.1.data
                      , init = init.prison
                      , tseq = time.out2
                      , control = list(trace = 3, maxit = 150)
                      , method = "SANN")


## View the R0 estimate
(R0.prison <- exp(unname(prison.estim$par)))  # estimate from prison data

## Run the SEIR model with the R0 estimate
prison.ts <- simEpidemic(init.prison, time.out2, seixc, disease_params(R0 = R0.prison))

## Visually compare incidence rates of the model and actual data

df2 <- rbind(ASP.1.data %>% 
              select(-Date) %>% 
              mutate(variable = "ASP.1.data"),
            prison.ts %>% 
              select(incidence, inc_rate, time, N) %>% 
              mutate(variable = "prison.ts"))

ggplot(df2, aes(x = time, y = inc_rate, color = variable))+
  geom_line() +
  geom_point(data = subset(df2, df2$variable == "ASP.1.data"), mapping=aes(x = time, y=inc_rate, color = variable)) +
  labs(title = "Comparing incidence rates of data and model",
       subtitle = "Avenal State Prison Outbreak 1 data",
       x= "Time [days]",
       y= "Incidence rate") +
  scale_color_manual(labels = c("Data","Model"), values = c("red", "blue"))
#  theme(text = element_text(size = 18))

## The model incidence and actual incidence clearly do not fit well based on
## visual inspection. The model cannot match the data incidence with two peaks.
## Let's try to cut the data to only include the initial peak, and then assess
## model performance.

ASP.1.data.cut <- ASP.1.data %>%
  filter(time < 60)  # cutoff at day 60 of outbreak

time.out3 <- seq(0, 60, 1)

prison.estim.cut <- optim(par = c(log_R0 = log(1))
                          , objFXN
                          , fixed.params = disease_params()  # defined at top of script
                          , obsDat = ASP.1.data.cut
                          , init = init.prison
                          , tseq = time.out3
                          , control = list(trace = 3, maxit = 150)
                          , method = "SANN")


## View the R0 estimate
(R0.prison.cut <- exp(unname(prison.estim.cut$par)))  # estimate from prison data

## Run the SEIR model with the R0 estimate

prison.ts.cut <- simEpidemic(init.prison, time.out3, seixc, disease_params(R0 = R0.prison.cut))

## Visually compare incidence rates of the model and actual data

df3 <- rbind(ASP.1.data.cut %>% 
               select(-Date) %>% 
               mutate(variable = "ASP.1.data.cut"),
             prison.ts.cut %>% 
               select(incidence, inc_rate, time, N) %>% 
               mutate(variable = "prison.ts.cut"))

ggplot(df3, aes(x = time, y = inc_rate, color = variable))+
  geom_line() +
  geom_point(data = subset(df3, df3$variable == "ASP.1.data.cut"), mapping=aes(x = time, y=inc_rate, color = variable)) +
  labs(title = "Comparing incidence rates of data and model",
       subtitle = "Avenal State Prison Outbreak 1 data (cut off to one peak)",
       x= "Time [days]",
       y= "Incidence rate") +
  scale_color_manual(labels = c("Data","Model"), values = c("red", "blue"))
#  theme(text = element_text(size = 18))

## The model is somewhat closer to the data in shape, but still not a good fit.
## Estimated R0 is likely inaccurate.
##
## > Model needs additional processes or compartments?
## > Initial conditions must be changed?
## > Errors in the code?
## > ...

# Plots used in the report ------------------------------------------------

## Figure 1
ggplot(df, aes(x = time, y = inc_rate, color = variable))+
  geom_line() +
  geom_point(data = myDat, mapping=aes(x = time, y=samp_inc_rate, color = variable)) +
  #  geom_ribbon(data = myDat, mapping=aes(x = time, y = samp_inc_rate, ymin=lci, ymax=uci), color = NA, fill = "black", alpha = 0.2) +
  labs(x= "Time [days]",
       y= "Incidence rate") +
  scale_color_manual(labels = c("Fitted","Truth","Data"), values = c("red", "blue", "black")) +
  theme(legend.title= element_blank())

## Figure 2
A <- ggplot(df2, aes(x = time, y = inc_rate, color = variable))+
  geom_line() +
  geom_point(data = subset(df2, df2$variable == "ASP.1.data"), mapping=aes(x = time, y=inc_rate, color = variable)) +
  labs(x= "Time [days]",
       y= "Incidence rate") +
  scale_color_manual(labels = c("Data","Model"), values = c("red", "blue")) +
  theme(legend.position = "bottom", legend.title= element_blank())

B <- ggplot(df3, aes(x = time, y = inc_rate, color = variable))+
  geom_line() +
  geom_point(data = subset(df3, df3$variable == "ASP.1.data.cut"), mapping=aes(x = time, y=inc_rate, color = variable)) +
  labs(x= "Time [days]",
       y= "Incidence rate") +
  scale_color_manual(labels = c("Data","Model"), values = c("red", "blue"))

AB <- plot_grid(A+theme(legend.position = "none"), 
                B + theme(legend.position = "none"), 
                labels = "AUTO", ncol=2)

legend <- get_legend(A)

plot_grid(AB, legend, ncol=1,rel_heights = c(0.95, 0.05))
