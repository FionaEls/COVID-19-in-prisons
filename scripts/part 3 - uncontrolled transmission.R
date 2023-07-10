require(boot); require(deSolve); require(ellipse);

# Parameters ----

## Parameters for SEIX model
## to be passed to disease_params()
R0.estim <- 1.5          # estimated basic reproduction number (placeholder)
latentPeriod <- 3        # days (from literature)
infectiousPeriod <- 7    # days (from literature)

## Parameters for initial conditions of state variables
## to be passed to init_cond()
N0.estim <- 1000         # initial population

# Functions ----

## Function that makes a list of disease parameters with default values
disease_params <- function(R0 = R0.estim,   # transmission coefficient
                           gamma = 1/infectiousPeriod, # rate of removal
                           sigma = 1/latentPeriod     # rate of progression
){
  return(as.list(environment()))
}

## Function that makes a vector of initial state variable conditions with N0 
## as argument

init_cond <- function(N0) {
  init <- c(S = N0 - 1,  # susceptible
            E = 1,       # exposed
            I = 0,       # infected
            X = 0)       # removed
  
  return(init)
}

## SEIX model
seix <- function(tt, yy, parms) with(c(parms,as.list(yy)), {
  
  ## State variables are: S, E, I, X
  N <- S + E + I + X                ## total population
  beta <- R0*gamma                  ## Infectious contact rate
  forceOfInfection <- beta*I/N   
  
  ## state variable derivatives (ODE system)
  deriv <- rep(NA,4)
  deriv[1] <-	-forceOfInfection*S                    # dSdt
  deriv[2] <-	forceOfInfection*S - sigma*E           # dEdt
  deriv[3] <-	sigma*E - gamma*I                      # dIdt
  deriv[4] <-	gamma*I                                # dXdt
  return(list(deriv))
})

## Function to run the deterministic model simulation, based on the ODE system defined in seix().
simEpidemic <- function(init, tseq, modFunction=seix, parms = disease_params()) {
  modDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  modDat$N <- rowSums(modDat[, c('S','E','I','X')]) # total outbreak population
  modDat$PI <- with(modDat, (E+I)/N)                # prevalence of infection (E + I)
  modDat$PD <- with(modDat, I/N)                    # prevalence of disease (I only)
  modDat$CI <- cumsum(modDat[,'I'])                 # cumulative I
  return(modDat)
}

## Function to 'sample' the population:
## From a simulated epidemic, measure prevalence at several time points by drawing
## cross-sectional samples of individuals at each time, testing them, and then calculating sample
## prevalence and associated binomial confidence intervals

sampleEpidemic <- function(modDat                                      # Simulated "data" which we treat as real 
                           , sampleDates = seq(1, 365, by = 5)         # Sample every 5 days for 1 year
                           , numN = rep(80, length(sampleDates))       # Number of individuals sampled at each time point
){
  prev_at_sample_times <- modDat[modDat$time %in% sampleDates, 'PI']
  cases <- rbinom(length(numN), numN, prev_at_sample_times)
  lci <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = cases, n = numN)
  uci <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = cases, n = numN)    
  return(data.frame(time = sampleDates, cases, numN, sampPrev =  cases/numN,
                    lci = lci, uci = uci))
}

## negative loglikelihood
## uses prevalence as basis for matching model-generated data to real data
nllikelihood <- function(parms, obsDat, init, tseq) {
  modDat <- simEpidemic(init, parms=parms, tseq)
  ## What are the rows from our simulation at which we have observed data?
  matchedTimes <- modDat$time %in% obsDat$time
  nlls <- -dbinom(obsDat$cases,                   # infected cases of observed data
                  obsDat$numN,                    # total population at timepoints of observed data
                  prob = modDat$PI[matchedTimes], # Prevalence from model data
                  log = T)
  return(sum(nlls))
}

## Combine guess (fit parameters) and fixed parameters
subsParms <- function(fit.params, fixed.params=disease_params())
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]        
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })

## Make likelihood a function of fixed and fitted parameters
objFXN <- function(fit.params                        # paramters to fit
                   , fixed.params = disease_params() # fixed paramters
                   , obsDat                          # observed data
                   , init
                   , tseq) {                 
  parms <- subsParms(fit.params, fixed.params)
  nllikelihood(parms, obsDat = obsDat, init=init, tseq = tseq)  # then call likelihood
}

# Plotting code ----

# Initial plots

## Plot model prevalence through time:
# par(bty='n', lwd = 2)
# with(modDat, plot(time, PI, xlab = '', ylab = 'prevalence', type = 'l', ylim = c(0,0.4), col='black', las = 1))
# points(myDat$time, myDat$sampPrev,
#         col = 'red', pch = 20, cex = .5) # Plot sample prevalence at each time point
# arrows(myDat$time, myDat$uci, myDat$time, myDat$lci, col = 'red', len = .025, angle = 90, code = 3) # Plot 95% CIs around the sample prevalences

# ## Plot SEIX model:
# 
# modDat.long <- melt(as.data.table(modDat), id.vars = 'time')
# 
# ## SEIX plot
# 
# (modDat.long %>%
#     filter(variable %in% c('S', 'E', 'I', 'X')) %>% 
#     ggplot()
#   + aes(x = time, y = value, color = variable, linetype = variable)
#   + geom_line()
#   + labs(title = "SEIX model time series",
#          x = "time [days]",
#          y = "cases"))
# 
# ## Prevalence of Infection, Disease
# 
# (ggplot(modDat.long[variable %in% c('PI', 'PD')])
#   + aes(x = time, y = value)
#   + geom_line(aes(color = variable))
#   + scale_color_discrete(labels = c("Infection","Disease"))
#   + labs(title = "Prevalence of infection and disease",
#          x = "time [days]",
#          y = "prevalence")
#   + geom_point(data = myDat, aes(x = time, y = sampPrev), color = "black")
# )

# ## Plot MLE fit time series
# 
# par(bty='n', lwd = 2, las = 1)
# maxPI <- max(modDat$PI)
# ymax = maxPI + maxPI/(maxPI*100)
# 
# ### Model data
# with(modDat, plot(time, PI, xlab = '', ylab = 'prevalence', type = 'l', ylim = c(0,ymax), col='red'))
# 
# ### Fit data (with parameters estimated from "real" data)
# fitDat <- simEpidemic(init, parms = subsParms(optim.vals$par, trueParms))
# with(fitDat, lines(time, PI, col='blue'))
# 
# ### Observed data (from "real" data)
# points(myDat$time, myDat$sampPrev, col = 'red', pch = 16, cex = 0.5)
# # arrows(myDat$time, myDat$uci, myDat$time, myDat$lci, col = 'red', len = .025, angle = 90, code = 3)
# legend("topleft", c('truth', 'observed', 'fitted'), lty = c(1, NA, 1), pch = c(NA,16,NA),
#        col = c('red', 'red', 'blue'), bty = 'n')

# Test functions ----

## Generate simulated data by sampling from model-generated data
trueParms <- disease_params()                  # Default model parameters
init.pop <- init_cond(N0 = 1000)
time.out <- seq(0, 365, 1)  # range of timepoints for model time series
modDat <- simEpidemic(init.pop, parms = trueParms, tseq = time.out) # Simulated epidemic (underlying process)
set.seed(1)                            # Initiate the random number generator
myDat <- sampleEpidemic(modDat)        # Simulate data from the sampling process (substitute to REAL DATA)

## Optim to estimate R0
## init.pars <- c(log_R0 = log(1))
optim.vals <- optim(par = c(log_R0 = log(1))
                    , objFXN
                    , fixed.params = disease_params()  # defined at top of script
                    , obsDat = myDat
                    , init = init.pop
                    , tseq = time.out
                    , control = list(trace = 3, maxit = 150)
                    , method = "SANN")


## R0 estimate
exp(unname(optim.vals$par))  # estimate from simulated data
trueParms['R0']              # true value

## Plot

fitDat <- simEpidemic(init = init.pop, tseq = time.out, parms = subsParms(optim.vals$par, trueParms))

df <- rbind(modDat = modDat %>% 
              mutate(variable = "modDat"),
            fitDat = fitDat %>% 
              mutate(variable = "fitDat"))

myDat %>% 
  filter(sampPrev < 0.065) -> myDat2

ggplot(df, aes(x = time, y = PI, color = variable))+
  geom_line() +
  geom_point(data = myDat2, mapping=aes(x = time, y=sampPrev), color = "black") +
  labs(title = "Simulated data",
       x= "time [days]",
       y= "prevalence") +
  scale_color_discrete(labels = c("Fitted","Truth")) +
  theme(text = element_text(size =18))

# Estimate R0 from data ----

prison.data <- read_csv("Infectous_timeseries.csv")
fixed.N0 <- 6500  # initial population of prison
prison.data <- prison.data %>%
  select(day, I) %>% 
  rename(time = day, cases = I) %>% 
  mutate(sampPrev = cases/fixed.N0) %>% 
  mutate(numN = rep(fixed.N0, nrow(.))) %>% 
  filter(time < 100)

prison.data

time.out2 <- seq(0, 264, 1)
init.prison <- init_cond(N0 = fixed.N0)

prison.estim <- optim(par = c(log_R0 = log(1))
                    , objFXN
                    , fixed.params = disease_params()  # defined at top of script
                    , obsDat = prison.data
                    , init = init.prison
                    , tseq = time.out2
                    , control = list(trace = 3, maxit = 150)
                    , method = "SANN")


## R0 estimate
R0.prison <- exp(unname(prison.estim$par))  # estimate from prison data
R0.prison

## Plot SEIX

prison.ts <- simEpidemic(init.prison, time.out2, seix, disease_params(R0 = R0.prison))
prison.ts.long <- melt(as.data.table(prison.ts), id.vars = 'time')

# (prison.ts.long %>%
#     filter(variable %in% c('S', 'E', 'I', 'X')) %>%
#     ggplot()
#   + aes(x = time, y = value, color = variable, linetype = variable)
#   + geom_line()
#   + labs(title = "SEIX model time series",
#          x = "time [days]",
#          y = "cases"))

## Plot prevalence

(ggplot(prison.ts.long[variable %in% c('PI')])
  + aes(x = time, y = value)
  + geom_line(color = "red")
  + scale_color_discrete(labels = c("Infection","Disease"))
  + labs(title = "Actual data",
         x = "time [days]",
         y = "prevalence")
  + geom_point(data = prison.data, aes(x = time, y = sampPrev), color = "black")
  + theme(text = element_text(size =18))
  + xlim(0,110)
)

