##################################################################################
##################################################################################
# an R script to estimate R0 and N using an SEIR model
#
# Project: COVID in Prisons, ICI3D MMED 2023
#
##################################################################################

# Load packages -----------------------------------------------------------

require(boot); require(deSolve); require(ellipse); require(tidyverse)

# Parameters --------------------------------------------------------------

## Parameters for SEIR model
## to be passed as default parameters to disease_params()
latentPeriod <- 3        # days
infectiousPeriod <- 7    # days

# Functions (edited from Lab 6) --------------------------------------------

## Function that makes a list of disease parameters ##
##
## > R0 and N are undefined (we treat N as an unknown parameter, because the
##   population size for each outbreak in the prisons are not known. Therefore,
##   we have to estimate it.)
## > Rates gamma and sigma are defined with the parameters stated above.

disease_params <- function(R0,                         # basic reproduction number, unknown
                           N,                          # the total population, unknown
                           gamma = 1/infectiousPeriod, # rate of removal
                           sigma = 1/latentPeriod)     # rate of progression
{
  return(as.list(environment()))
}

## Function that makes a vector of initial state variable conditions ##
##
## initI - initial value of I as argument (defaults to 1)

init_cond <- function(initI = 1) {
  init <- c(E = 0,       # exposed
            I = initI,   # infected
            X = 0,       # removed
            C = 0)       # cumulative cases
  
  return(init)
}

## Function defining the SEIR model##
##
## > We remove the S compartment, since we cannot define its initial conditions
##   due to N being treated as an unknown parameter. S will be calculated later
##   in the next function (simEpidemic).
## > 'X' is used instead of 'R' to avoid confusion with R0, the reproduction 
##   number.
## > 'C' is an added compartment that calculates the cumulative infected cases 
##   per time step. This will be equivalent to the cumulative incidence. We
##   will use this to fit the model to the actual data, since our prison data
##   consists of the daily incidence counts per prison outbreak.
## > Overall, this model is a standard deterministic SEIR model with no births,
##   deaths, introductions, replenishment of susceptibles, etc.
seixc <- function(tt, yy, parms) with(c(parms,as.list(yy)), {
  
  ## Infectious contact rate
  beta <- R0*gamma/N                
  
  ## State variable derivatives (ODE system)
  ## State variables are: E, I, X, C
  ## Since there is no S compartment, we substitute S with N-E-I-X
  deriv <- rep(NA,4)
  deriv[1] <-	beta*(N-E-I-X)*I - sigma*E   # dEdt
  deriv[2] <-	sigma*E - gamma*I            # dIdt
  deriv[3] <-	gamma*I                      # dXdt
  deriv[4] <- sigma*E                      # cumulative new cases C
  return(list(deriv))
})

## Function to run the deterministic model simulation, based on the ODE system 
## defined in seixc() ##
##
## This function outputs a dataframe with the following variables:
## > time - time-step
## > S, E, I, X, C - number of individuals in each compartment per timestep
## > N - total population
## > cinc_rate - the cumulative incidence rate (C/N) per time step (this will be
##   used for creating simulated data with sampleEpidemic)
simEpidemic <- function(init, tseq, modFunction=seixc, parms = disease_params()) {
  modDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  modDat$N <- rep(parms$N, nrow(modDat))             # total population N
  modDat$S <- with(modDat, N-E-I-X)                  # Here we add compartment S. Recall that S = N-E-I-X
  modDat$cinc_rate <- with(modDat, C/N)              # cumulative incidence rate based on total population N
  col_order <- c("time","S", "E", "I", "X", "C", "N", "cinc_rate")  # for rearranging the columns
  return(modDat[,col_order])
}

## Function to create simulated data using the model ##
##
## > From a simulated epidemic, measure the cumulative incidence rate at several
##   time points by drawing cross-sectional samples of individuals at each time,
##   testing them, and then calculating sample cumulative incidence rate and 
##   associated binomial confidence intervals.
##
## > This function will be used to generate simulated data for testing the model 
##   and optimization functions.

sampleEpidemic <- function(modDat,                                   # Simulated "data" which we treat as real 
                           sampleDates = seq(1, 365, by = 3),        # Sample every 3 days for 1 year
                           N = rep(1000, length(sampleDates))        # Number of individuals sampled at each time point
){
  cinc_rate_at_sample_times <- modDat[modDat$time %in% sampleDates, 'cinc_rate']
  C <- rbinom(length(N), N, cinc_rate_at_sample_times)
  lci <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = C, n = N)
  uci <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = C, n = N)  
  return(data.frame(time = sampleDates, C, N, samp_cinc_rate =  C/N,
                    lci = lci, uci = uci))
}

## Function to calculate residual sum of squares ##
##
## > We used residual sum of squares because we do not have the total population
##   N of the actual data, so we can't use negative log-likelihood with dbinom
##   like in Lab 6. Are there alternatives to this?
## > Using the residual sum of squares of the cumulative incidence produces
##   better model fits than using the incidence data directly, hence our use
##   of the cumulative incidence instead of the actual daily incidence.
## > This will be passed to the objective function, which will be passed to 
##   optim().
ressum <- function(parms, obsDat, init, tseq) {
  modDat <- simEpidemic(init, parms=parms, tseq)
  ## What are the rows from our simulation at which we have observed data?
  matchedTimes <- modDat$time %in% obsDat$time
  sqe <- (obsDat$C - modDat$C[matchedTimes])^2
  return(sum(sqe))
}

## Function to combine guess (fit parameters) and fixed parameters ##
##
## This will be used within the objective function.
subsParms <- function(fit.params, fixed.params=disease_params())
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]        
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })

## Objective function to be passed to optim() for parameter estimation ##
##
## Make the objective function a function of fixed and fitted parameters
objFXN <- function(fit.params,                        # paramters to fit
                   fixed.params = disease_params(),   # fixed paramters
                   obsDat,                            # observed data
                   init,
                   tseq) {                 
  parms <- subsParms(fit.params, fixed.params)
  ressum(parms, obsDat = obsDat, init=init, tseq = tseq)  # then ressum
}

# Test the functions with simulated data ----------------------------------

## Generate simulated data ##

## We will create "observed data" by sampling from a model with known R0 and N
trueParms <- disease_params(R0 = 1.50, N=1000)  # underlying model parameters
init.pop <- init_cond()                         # keep the default initial conditions
time.out <- seq(0, 365, 1)                      # range of timepoints for model time series

modDat <- simEpidemic(init.pop, parms = trueParms, tseq = time.out) # Simulated epidemic (underlying process)
ggplot(modDat, aes(x = time, y = C)) +          # quick look at the model cumulative incidence
  geom_col() +
  labs(title = "Model Cumulative Incidence")

set.seed(1)                                     # Initiate the random number generator
myDat <- sampleEpidemic(modDat)                 # Simulate data from the sampling process (acts as our "observed data")
ggplot(myDat, aes(x=time, y=C)) +               # quick look at the simulated cumulative incidence
  geom_col() +
  labs(title = "Simulated Cumulative Incidence")

## Use optim to estimate R0 and N ##

## Initial estimation via SANN
optim.vals <- optim(par = c(log_R0 = log(1), log_N = log(1000))
                    , objFXN
                    , fixed.params = disease_params()  # defined at top of script
                    , obsDat = myDat
                    , init = init.pop
                    , tseq = time.out
                    , control = list(trace = 3, maxit = 150)
                    , method = "SANN")

## See initial estimates
estimParms <- exp(unname(optim.vals$par))
names(estimParms) <- c('R0', 'N')
trueParms[c('R0', 'N')]   # true values
estimParms                # estimated values

## Feed the last parameters of SANN in as the first values of Nelder-Mead
optim.vals <- optim(par = optim.vals$par
                    , objFXN
                    , fixed.params = disease_params()
                    , obsDat = myDat
                    , init = init.pop
                    , tseq = time.out
                    , control = list(trace = 3, maxit = 800, reltol = 10^-7)
                    , method = "Nelder-Mead"
                    , hessian = T)


## Store estimates
MLEfits <- optim.vals$par
log_R0.fit <- MLEfits["log_R0"]
log_N.fit <- MLEfits["log_N"]
estimParms <- exp(unname(MLEfits))
names(estimParms) <- c('R0', 'N')

## View the R0 and N estimates
trueParms[c('R0', 'N')]   # true values
estimParms                # estimated values

## Simulate an epidemic using the R0 and N estimates
fitDat <- simEpidemic(init = init.pop, tseq = time.out, parms = subsParms(optim.vals$par, trueParms))

## Visually compare cumulative incidence of fitDat and modDat (the true underlying process)
df <- rbind(modDat = modDat %>% 
              mutate(variable = "modDat"),
            fitDat = fitDat %>% 
              mutate(variable = "fitDat"))

myDat <- myDat %>% 
  mutate(variable = "myDat")

ggplot(df, aes(x = time, y = C, color = variable))+
  geom_line() +
  geom_point(data = myDat, mapping=aes(x = time, y=C, color = variable)) +
  labs(title = "Comparing the cumulative incidence of true and fitted models",
       subtitle = paste("True parameters: R0 =", trueParms['R0'], ", N =", trueParms['N'],
                        "\nEstimated parameters: R0 =", round(estimParms['R0'],2), ", N =", round(estimParms['N'],2)),
       x= "Time [days]",
       y= "Cases") +
  scale_color_manual(labels = c("Fitted","Truth","Data"), 
                     values = c("red", "blue", "black"))


# Contour plots for simulated data (edited from Lab 6) ---------------------

## Generate contour plots with the Hessian ##
fisherInfMatrix <- solve(optim.vals$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates

## Create a sequence of R0 and N values for the grid
res <- 15
R0.seq <- exp(seq(log_R0.fit-1, log_R0.fit+1, l = res))
N.seq <- exp(seq(log_N.fit-1, log_N.fit + 1, l = res))

## Initialize plot of parameters
plot(1,1, type = 'n', log = 'xy',
     # xlim = range(R0.seq), ylim = range(N.seq),
     xlim = c(1.495,1.505), ylim = c(1000,1004),  # use these limits to zoom-in
     las = 1,
     xlab = expression(R0), ylab = expression(N),
     main = "-log(likelihood) contours", bty = "n")
## Add true parameter values to the plot
with(trueParms, points(R0, N, pch = 16, cex = 2, col = 'red'))
## Add MLE to the plot
points(exp(MLEfits['log_R0']), exp(MLEfits['log_N']), pch = 16, cex = 2, col = 'black')
## Add 95% contour ellipse from Hessian
lines(exp(ellipse(fisherInfMatrix, centre = MLEfits, level = .95)))
legend("topleft", c('truth', 'MLE', '95% Confidence Region'), lty = c(NA, NA, 1), pch = c(16,16, NA),
       col = c('red', 'black', 'black'), bg='white', bty = 'n')

## The 95% confidence region is quite small and doesn't cover the true value?

## Generate contour plots with likelihood profiles ##
objXR0_N <- function(R0, N, fixed.params = disease_params(), browse=F, obsDat,
                     init, tseq) {
  objFXN(fit.params = c(log_R0 = log(R0), log_N = log(N))
         , fixed.params = fixed.params, obsDat, init, tseq)}

objXR0_NVEC <- Vectorize(objXR0_N, list("R0", "N"))

mat <- outer(R0.seq, N.seq, objXR0_NVEC, init = init_cond(), tseq=time.out, obsDat=myDat)

ml.val <- optim.vals$value
conf.cutoff <- ml.val + qchisq(.95,2)/2

## Show likelihood contours
par(cex = 1.2)
plot(1,1, type = 'n', log = 'xy',
     # xlim = range(R0.seq), ylim = range(N.seq),
     xlim = c(1.495,1.505), ylim = c(1000,1004),  # use these limits to zoom-in
     xlab = expression(R0), ylab = expression(N),
     main = "-log(likelihood) contours", bty = "n")
.filled.contour(R0.seq, N.seq, mat, levels = seq(min(mat), max(mat), l=20), col = topo.colors(20))
## Add contour for 95% CI from likelihood ratio
## This comes out very small!
contour(R0.seq, N.seq, mat, levels = c(conf.cutoff),
        col = "black", lwd = 2, labels = "", labcex = .2, add = T)
## Add contour for 95% CI from Hessian
lines(exp(ellipse(fisherInfMatrix, centre = MLEfits, level = .95)), lty = 2)
## Add MLE to the plot
points(exp(log_R0.fit), exp(log_N.fit), pch = 16, cex = 1, col = 'black')
## Add true parameter values to the plot
with(trueParms, points(R0, N, pch = 16, cex = 1, col = 'red'))
legend("topleft",
       c('truth', 'MLE', '95% contour (profile likelihood)', '95% contour (Fisher information matrix)')
       , lty = c(NA, NA, 1, 2), pch = c(16,16, NA, NA),
       col = c('red', rep('black',3)), bg='white', bty = 'n')

## The 95% contours are still quite small and don't cover the true value?

# Estimate R0 and N with actual data ---------------------------------------

## Import and wrangle dataset
## This dataset contains the daily incidence of each outbreak of each prison
all.data <- read_csv("new_cases_per_day.csv")

## Extract one outbreak: Centinela State Prison (CEN) Outbreak 5
CEN.data <- all.data %>%
  filter(Facility_Outbreak == "Centinela State Prison (CEN) Outbreak 5") %>%
  select(Date, Daily_New_Cases) %>% 
  slice(-1) %>%                             # remove first date with 0 new cases
  rename(incidence = Daily_New_Cases) %>%
  mutate(C = cumsum(incidence)) %>%         # calculate cumulative incidence
  mutate(time = seq(0,nrow(.)-1)) %>% 
  mutate(variable = "CEN.data")             # will be used for ggplot later

summary(CEN.data)

## Quick look at the daily incidence of CEN Outbreak 5
## > There are some days when there are suddenly a large numbers of cases 
##   detected
## > This could affect the reliability of our estimates
ggplot(CEN.data, aes(time, incidence)) + 
  geom_col() +
  labs(title = "Daily Incidence of CEN Outbreak 5",
       y = "Incidence (cases)")

## Quick look at the cumulative incidence of CEN Outbreak 5
ggplot(CEN.data, aes(time, C)) + 
  geom_col() +
  labs(title = "Cumulative Incidence of CEN Outbreak 5",
       y = "Cumulative incidence")

## Define time sequence and initial conditions
time.out.CEN <- seq(0, max(CEN.data$time), 1)
init.CEN <- init_cond(initI = CEN.data$C[1])

## Use optim to estimate R0
CEN.estim <- optim(par = c(log_R0 = log(1), log_N = log(1000))
                   , objFXN
                   , fixed.params = disease_params()  # defined at top of script
                   , obsDat = CEN.data
                   , init = init.CEN
                   , tseq = time.out.CEN
                   , control = list(trace = 3, maxit = 150)
                   , method = "SANN")

## Feed the last parameters of SANN in as the first values of Nelder-Mead
CEN.estim <- optim(par = CEN.estim$par
                   , objFXN
                   , fixed.params = disease_params()  # defined at top of script
                   , obsDat = CEN.data
                   , init = init.CEN
                   , tseq = time.out.CEN
                   , control = list(trace = 3, maxit = 800, reltol = 10^-7)
                   , method = "Nelder-Mead"
                   , hessian = T)

## Store estimates
CENfits <- CEN.estim$par
log_R0.fit <- CENfits["log_R0"]
log_N.fit <- CENfits["log_N"]
CENParms <- exp(unname(CENfits))
names(CENParms) <- c('R0', 'N')

## View the R0 and N estimates
CENParms

## Run the SEIR model with the R0 and N estimates
CEN.ts <- simEpidemic(init.CEN, time.out.CEN, seixc, disease_params(R0 = CENParms['R0'], N = CENParms['N']))

## Visually compare incidence of the model and actual data

CEN.ts <- CEN.ts %>% 
  mutate(variable = "CEN.ts")

ggplot(CEN.ts, aes(x = time, y = C, color = variable))+
  geom_line() +
  geom_point(data = CEN.data, mapping=aes(x = time, y=C, color = variable)) +
  labs(x= "Time [days]",
       y= "Cumulative cases",
       title = "Comparing incidence of data and model",
       subtitle = paste("CEN Outbreak 5 data
Estimates: R0 = ", round(CENParms['R0'],2), ", N = ", round(CENParms['N'],2)),
x= "Time [days]",
y= "Cases") +
  scale_color_manual(labels = c("Data","Model"), values = c("black", "red")) +
  theme(legend.title= element_blank())

## Are the estimated R0 and N values reasonable?

max(CEN.data$C)   # end size of the outbreak
CENParms['N']     # the estimated N

# Contour plots for actual data (edited from Lab 6) ------------------------
## We will use the estimated R0 and N from the ASP1 data restricted to one peak

## Generate contour plots with the Hessian ##
fisherInfMatrix <- solve(CEN.estim$hessian)

## Create a sequence of R0 and N values for the grid
res <- 15
R0.seq <- exp(seq(log_R0.fit-1, log_R0.fit+1, l = res))
N.seq <- exp(seq(log_N.fit-1, log_N.fit+1, l = res))

## Initialize plot of parameters
plot(1,1, type = 'n', log = 'xy',
     # xlim = range(R0.seq), ylim = range(N.seq),
     xlim = c(2.93,2.95), ylim = c(841,843),  # use these limits to zoom-in
     las = 1,
     xlab = expression(R0), ylab = expression(N),
     main = "-log(likelihood) contours", bty = "n")
## Add MLE to the plot
points(CENParms['R0'], CENParms['N'], pch = 16, cex = 2, col = 'black')
## Add 95% contour ellipse from Hessian
lines(exp(ellipse(fisherInfMatrix, centre = CENfits, level = .95)))
legend("topleft", c('MLE', '95% Confidence Region'), lty = c(NA, 1), pch = c(16, NA),
       col = c('black', 'black'), bg='white', bty = 'n')

## Generate contour plots with likelihood profiles ##
mat <- outer(R0.seq, N.seq, objXR0_NVEC, init = init.CEN, tseq=time.out.CEN, obsDat=CEN.data)

ml.val <- CEN.estim$value
conf.cutoff <- ml.val + qchisq(.95,2)/2

## Show likelihood contours
par(cex = 1.2)
plot(1,1, type = 'n', log = 'xy',
     # xlim = range(R0.seq), ylim = range(N.seq),
     xlim = c(2.93,2.95), ylim = c(841,843),  # use these limits to zoom into the Fisher information matrix contour
     # xlim = c(2.94280,2.94284), ylim = c(842.160,842.164),  # use these limits to zoom into the profile likelihood contour
     xlab = expression(R0), ylab = expression(N),
     main = "-log(likelihood) contours", bty = "n")
.filled.contour(R0.seq, N.seq, mat, levels = seq(min(mat), max(mat), l=20), col = topo.colors(20))
## Add contour for 95% CI from likelihood ratio
## This comes out extremely small
## Change xlim and ylim to zoom-in further
contour(R0.seq, N.seq, mat, levels = c(conf.cutoff),
        col = "black", lwd = 2, labels = "", labcex = .2, add = T)
## Add contour for 95% CI from Hessian
lines(exp(ellipse(fisherInfMatrix, centre = CENfits, level = .95)), lty = 2)
## Add MLE to the plot
points(exp(log_R0.fit), exp(log_N.fit), pch = 16, cex = 1, col = 'black')
## Add true parameter values to the plot
legend("topleft",
       c('MLE', '95% contour (profile likelihood)', '95% contour (Fisher information matrix)')
       , lty = c(NA, 1, 2), pch = c(16, NA, NA),
       col = c(rep('black',3)), bg='white', bty = 'n')


# Estimating R0 and N with other outbreak data ----------------------------

## Extract another outbreak: San Quentin State Prison (SQ) Outbreak 1
SQ.data <- all.data %>%
  filter(Facility_Outbreak == "San Quentin State Prison (SQ) Outbreak 1") %>%
  select(Date, Daily_New_Cases) %>% 
  slice(-1) %>%                             # remove first date with 0 new cases
  rename(incidence = Daily_New_Cases) %>%
  mutate(C = cumsum(incidence)) %>%         # calculate cumulative incidence
  mutate(time = seq(0,nrow(.)-1)) %>% 
  mutate(variable = "SQ.data")             # will be used for ggplot later

SQ.data

## Quick look at the daily incidence of SQ Outbreak 1
ggplot(SQ.data, aes(time, incidence)) + 
  geom_col() +
  labs(title = "Daily Incidence of SQ Outbreak 1",
       y = "Incidence (cases)")

## Quick look at the cumulative incidence of SQ Outbreak 1
ggplot(SQ.data, aes(time, C)) + 
  geom_col() +
  labs(title = "Cumulative Incidence of SQ Outbreak 1",
       y = "Cumulative incidence")

## Define time sequence and initial conditions
time.out.SQ <- seq(0, max(SQ.data$time), 1)
init.SQ <- init_cond(initI = SQ.data$C[1])

## Use optim to estimate R0
SQ.estim <- optim(par = c(log_R0 = log(1), log_N = log(1000))
                  , objFXN
                  , fixed.params = disease_params()  # defined at top of script
                  , obsDat = SQ.data
                  , init = init.SQ
                  , tseq = time.out.SQ
                  , control = list(trace = 3, maxit = 150)
                  , method = "SANN")

## Feed the last parameters of SANN in as the first values of Nelder-Mead
SQ.estim <- optim(par = SQ.estim$par
                  , objFXN
                  , fixed.params = disease_params()  # defined at top of script
                  , obsDat = SQ.data
                  , init = init.SQ
                  , tseq = time.out.SQ
                  , control = list(trace = 3, maxit = 800, reltol = 10^-7)
                  , method = "Nelder-Mead"
                  , hessian = T)

## Store estimates
SQfits <- SQ.estim$par
estimSQ <- exp(unname(SQfits))
names(estimSQ) <- c('R0', 'N')

## View the R0 and N estimates
estimSQ
max(SQ.data$C)  # max cumulative incidence is greater than estimated N

## Extract another outbreak: Avenal State Prison (ASP) Outbreak 1
ASP.data <- all.data %>%
  filter(Facility_Outbreak == "Avenal State Prison (ASP) Outbreak 1") %>%
  select(Date, Daily_New_Cases) %>% 
  slice(-1) %>%                             # remove first date with 0 new cases
  rename(incidence = Daily_New_Cases) %>%
  mutate(C = cumsum(incidence)) %>%         # calculate cumulative incidence
  mutate(time = seq(0,nrow(.)-1)) %>% 
  filter(time < 60) %>%                     # restrict data to one peak of the outbreak
  mutate(variable = "ASP.data")              # will be used for ggplot later

ASP.data

## Quick look at the daily incidence of SQ Outbreak 1
ggplot(ASP.data, aes(time, incidence)) + 
  geom_col() +
  labs(title = "Daily Incidence of ASP Outbreak 1",
       y = "Incidence (cases)")

## Quick look at the cumulative incidence of SQ Outbreak 1
ggplot(ASP.data, aes(time, C)) + 
  geom_col() +
  labs(title = "Cumulative Incidence of ASP Outbreak 1",
       y = "Cumulative incidence")

## Define time sequence and initial conditions
time.out.ASP <- seq(0, max(ASP.data$time), 1)
init.ASP <- init_cond(initI = ASP.data$C[1])

## Use optim to estimate R0
ASP.estim <- optim(par = c(log_R0 = log(1), log_N = log(1000))
                   , objFXN
                   , fixed.params = disease_params()  # defined at top of script
                   , obsDat = ASP.data
                   , init = init.ASP
                   , tseq = time.out.ASP
                   , control = list(trace = 3, maxit = 150)
                   , method = "SANN")

## Feed the last parameters of SANN in as the first values of Nelder-Mead
ASP.estim <- optim(par = ASP.estim$par
                   , objFXN
                   , fixed.params = disease_params()  # defined at top of script
                   , obsDat = ASP.data
                   , init = init.ASP
                   , tseq = time.out.ASP
                   , control = list(trace = 3, maxit = 800, reltol = 10^-7)
                   , method = "Nelder-Mead"
                   , hessian = T)

## Store estimates
ASPfits <- ASP.estim$par
estimASP <- exp(unname(ASPfits))
names(estimASP) <- c('R0', 'N')

## View the R0 and N estimates
estimASP
max(ASP.data$C)  # max cumulative incidence is less than estimated N

# Estimating R0 and N for all outbreaks -----------------------------------

## Initialize the table of all estimates ##
## maxC is the maximum cumulative incidence, which we can compare with the
## estimated N.
estim.table <- tibble(outbreak = character(),  # outbreak name
                      R0 = numeric(),          # R0 estimate
                      N = numeric(),           # N estimate
                      maxC = integer())        # max cumulative incidence

## Turn the estimation into a function ##
## Returns a new row for the table of all estimates
estimateR0N <- function(data) {
  
  ## Define time sequence and initial conditions
  time.out.data <- seq(0, max(data$time), 1)
  init.data <- init_cond(initI = data$C[1])
  
  ## Use optim to estimate R0
  data.estim <- optim(par = c(log_R0 = log(1), log_N = log(1000))
                      , objFXN
                      , fixed.params = disease_params()
                      , obsDat = data
                      , init = init.data
                      , tseq = time.out.data
                      , control = list(maxit = 150)
                      , method = "SANN")
  
  ## Feed the last parameters of SANN in as the first values of Nelder-Mead
  data.estim <- optim(par = data.estim$par
                      , objFXN
                      , fixed.params = disease_params()
                      , obsDat = data
                      , init = init.data
                      , tseq = time.out.data
                      , control = list(maxit = 800, reltol = 10^-7)
                      , method = "Nelder-Mead"
                      , hessian = T)
  
  R0N <- exp(unname(data.estim$par))
  
  estimData <- tibble_row(outbreak = data$outbreak[1], R0 = R0N[1], N = R0N[2], maxC = max(data$C))
  return(estimData)
}

## Turn data wrangling into a function ##
## > An issue here is that some outbreaks have multiple peaks. Our model cannot
##   fit to outbreaks with multiple peaks.
## > These peaks aren't being separated for now, but we will have to figure out
##   a way to do that.
extractOutbreak <- function(outbreak_name) {
  df <- all.outbreaks %>%
    filter(Facility_Outbreak == outbreak_name) %>%
    select(Date, Daily_New_Cases) %>%
    {if(.$Daily_New_Cases[1] == 0) slice(.,-1) else .} %>%
    rename(incidence = Daily_New_Cases) %>%
    mutate(C = cumsum(incidence)) %>%                         
    mutate(time = seq(0,nrow(.)-1)) %>%
    mutate(outbreak = outbreak_name)
  
  return(df)
}

## Extract outbreaks from all data ##
## > There are many "outbreaks" in the dataset that only have very small
##   outbreak sizes (i.e. some have only 1 case). These are not actual outbreaks 
##   and need to be removed.
## > For the purposes of this analysis, wee will only analyze the outbreaks with 
##   outbreak size greater than 10.
all.outbreaks <- all.data %>% 
  filter(outbreak_size > 10)

## Extract names of all outbreaks ##
outbreaks <- unique(all.outbreaks$Facility_Outbreak)

## Estimate R0 and N for all outbreaks ##
for (ii in outbreaks) {
  df <- extractOutbreak(ii)
  new.row <- estimateR0N(df)
  estim.table <- bind_rows(estim.table, new.row)
  print(paste("Done estimating:", ii))
}

## Explore estimates ##

## There are extreme estimated values for R0 and N.
summary(estim.table)

## There are estimated values for N that are smaller than the maximum cumulative
## case count.
any(estim.table$N < estim.table$maxC)
estim.table$outbreak[which(estim.table$N < estim.table$maxC)]

## We need to reevaluate our data, model, and optimization procedure.

# Plots used in the report ------------------------------------------------

## Figure 1 (simulated data)
ggplot(df, aes(x = time, y = C, color = variable))+
  geom_line() +
  geom_point(data = myDat, mapping=aes(x = time, y=C, color = variable)) +
  #  geom_ribbon(data = myDat, mapping=aes(x = time, y = samp_cinc_rate, ymin=lci, ymax=uci), color = NA, fill = "black", alpha = 0.2) +
  labs(x= "Time [days]",
       y= "Cumulative cases") +
  scale_color_manual(labels = c("Fitted","Truth","Data"), values = c("red", "blue", "black")) +
  theme(legend.title= element_blank())

## Figure 2 (actual data)
ggplot(CEN.ts, aes(x = time, y = C, color = variable))+
  geom_line() +
  geom_point(data = CEN.data, mapping=aes(x = time, y=C, color = variable)) +
  labs(x= "Time [days]",
       y= "Cumulative cases") +
  scale_color_manual(labels = c("Data","Model"), values = c("black", "red")) +
  theme(legend.title= element_blank())
