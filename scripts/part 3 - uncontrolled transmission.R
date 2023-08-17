##################################################################################
##################################################################################
# an R script to estimate R0 and N using an SEIR model
#
# Project: COVID in Prisons, ICI3D MMED 2023
#
##################################################################################

# Load packages -----------------------------------------------------------

require(boot); require(deSolve); require(ellipse); require(tidyverse); require(cowplot)

# Parameters --------------------------------------------------------------

## Parameters for SEIR model
## to be passed as default parameters to disease_params()
latentPeriod <- 3        # days
infectiousPeriod <- 7    # days

# Functions ---------------------------------------------------------------

## Function that makes a list of disease parameters ##
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
## initI - initial value of I as argument (defaults to 1)

init_cond <- function(initI = 1) {
  init <- c(E = 0,       # exposed
            I = initI,   # infected
            X = 0,       # removed
            C = 0)       # cumulative cases
  
  return(init)
}

## Function defining the SEIR model ##
## > We remove the S compartment, since we cannot define its initial conditions
##   due to N being treated as an unknown parameter. S will be calculated later
##   in the next function (simEpidemic).
## > 'X' is used instead of 'R' to avoid confusion with R0, the reproduction 
##   number.
## > 'C' is an added compartment that calculates the cumulative infected cases 
##   per time step. The difference between 'C' at time t and t-1 will be the 
##   incidence.
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
  deriv[4] <- sigma*E                      # cumulative cases C
  return(list(deriv))
})

## Function to run the deterministic model simulation, based on the ODE system defined in seixc() ##
## This function outputs a dataframe with the following variables:
## > time - time-step
## > S, E, I, X, C - number of individuals in each compartment per timestep
## > N - total population
## > incidence - the number of incident cases, calculated as the difference in C 
##   over each time step
## > inc_rate - the incidence rate (incidence/N) per time step
simEpidemic <- function(init, tseq, modFunction=seixc, parms = disease_params()) {
  modDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  modDat$N <- rep(parms$N, nrow(modDat))             # total population N
  modDat$S <- with(modDat, N-E-I-X)                  # Here we add compartment S. Recall that S = N-E-I-X
  modDat$incidence <- c(modDat$C[1],diff(modDat$C))  # incident cases (difference in C over each time step)
  modDat$inc_rate <- with(modDat, incidence/N)       # incidence rate based on total population N
  col_order <- c("time","S", "E", "I", "X", "C", "N", "incidence", "inc_rate")  # for rearranging the columns
  return(modDat[,col_order])
}

## Function to create simulated data using the model ##
## From a simulated epidemic, measure the incidence proportion at several time points by drawing
## cross-sectional samples of individuals at each time, testing them, and then calculating sample
## incidence proportion and associated binomial confidence intervals. The incidence per time
## step will then be calculated from the difference in exposed between each time step.
##
## This function be used to generate simulated data for testing the model and optimization functions.

sampleEpidemic <- function(modDat,                                   # Simulated "data" which we treat as real 
                           sampleDates = seq(1, 365, by = 3),        # Sample every 3 days for 1 year
                           N = rep(1000, length(sampleDates))        # Number of individuals sampled at each time point
){
  inc_rate_at_sample_times <- modDat[modDat$time %in% sampleDates, 'inc_rate']
  inc <- rbinom(length(N), N, inc_rate_at_sample_times)
  lci <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = inc, n = N)
  uci <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = inc, n = N)    
  return(data.frame(time = sampleDates, incidence = inc, N, samp_inc_rate =  inc/N,
                    lci = lci, uci = uci))
}

## Function to calculate residual sum of squares ## 
## This will be passed to the objective function, which will be passed to optim()
ressum <- function(parms, obsDat, init, tseq) {
  modDat <- simEpidemic(init, parms=parms, tseq)
  ## What are the rows from our simulation at which we have observed data?
  matchedTimes <- modDat$time %in% obsDat$time
  sqe <- (obsDat$incidence - modDat$incidence[matchedTimes])^2
  return(sum(sqe))
}

## Function to combine guess (fit parameters) and fixed parameters ##
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
ggplot(modDat, aes(x = time, y = incidence)) + 
  geom_col() +
  labs(title = "Model Incidence")              # quick look at the model incidence

set.seed(1)                                    # Initiate the random number generator
myDat <- sampleEpidemic(modDat)                # Simulate data from the sampling process (acts as our "observed data")
ggplot(myDat, aes(x=time, y=incidence)) + 
  geom_col() +
  labs(title = "Simulated Incidence")          # quick look at the simulated incidence

## Use optim to estimate R0 and N ##
## Intial estimation via SANN
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

## See final estimates
trueParms[c('R0', 'N')]   # true values
estimParms                # estimated values

## Simulate an epidemic using the R0 estimate (fitDat)
fitDat <- simEpidemic(init = init.pop, tseq = time.out, parms = subsParms(optim.vals$par, trueParms))

## Visually compare incidence (cases) of fitDat and modDat (the true underlying process)
df <- rbind(modDat = modDat %>% 
              mutate(variable = "modDat"),
            fitDat = fitDat %>% 
              mutate(variable = "fitDat"))

myDat <- myDat %>% 
  mutate(variable = "myDat")

ggplot(df, aes(x = time, y = incidence, color = variable))+
  geom_line() +
  geom_point(data = myDat, mapping=aes(x = time, y=incidence, color = variable)) +
  labs(title = "Comparing incidence of true and fitted models",
       subtitle = paste("True parameters: R0 =", trueParms['R0'], ", N =", trueParms['N'],
                        "\nEstimated parameters: R0 =", round(estimParms['R0'],2), ", N =", round(estimParms['N'],2)),
       x= "Time [days]",
       y= "Cases") +
  scale_color_manual(labels = c("Fitted","Truth","Data"), values = c("red", "blue", "black"))


# Contour plots for simulated data ----------------------------------------

## Generate contour plots with the Hessian ##
fisherInfMatrix <- solve(optim.vals$hessian) ## invert the Hessian, to estimate the covar-var matrix of parameter estimates

## Create a sequence of R0 and N values for the grid
res <- 15
R0.seq <- exp(seq(log_R0.fit-1, log_R0.fit+1, l = res))
N.seq <- exp(seq(log_N.fit-1, log_N.fit + 1, l = res))

## Initialize plot of parameters
plot(1,1, type = 'n', log = 'xy',
     # xlim = range(R0.seq), ylim = range(N.seq),
     xlim = c(1.4,1.6), ylim = c(950,1100),  # use these limits to zoom-in
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

## Change the xlim and ylim values to zoom-into the plot!

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
     xlim = c(1.48,1.52), ylim = c(950,1100),  # use these limits to zoom-in
     xlab = expression(R0), ylab = expression(N),
     main = "-log(likelihood) contours", bty = "n")
.filled.contour(R0.seq, N.seq, mat, levels = seq(min(mat), max(mat), l=20), col = topo.colors(20))
## Add contour for 95% CI from likelihood ratio
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

## Change the xlim and ylim values to zoom-into the plot!

# Estimate R0 and N with actual data ----------------------------------------

## Import and wrangle dataset
all.data <- read_csv("new_cases_per_day.csv")

## Extract one outbreak: Avenal State Prison (ASP) Outbreak 1
ASP.1.data <- all.data %>%
  filter(Facility_Outbreak == "Avenal State Prison (ASP) Outbreak 1") %>%
  select(Date, Daily_New_Cases) %>% 
  rename(incidence = Daily_New_Cases) %>% 
  mutate(time = seq(0,nrow(.)-1)) %>%
  filter(time < 100)

ASP.1.data

## Define time sequence and initial conditions
time.out2 <- seq(0, max(ASP.1.data$time), 1)
init.prison <- init_cond(initI = 25)

## Use optim to estimate R0
prison.estim <- optim(par = c(log_R0 = log(1), log_N = log(1000))
                      , objFXN
                      , fixed.params = disease_params()  # defined at top of script
                      , obsDat = ASP.1.data
                      , init = init.prison
                      , tseq = time.out2
                      , control = list(trace = 3, maxit = 150)
                      , method = "SANN")

## Feed the last parameters of SANN in as the first values of Nelder-Mead
# prison.estim <- optim(par = optim.vals$par
#                     , objFXN
#                     , fixed.params = disease_params()
#                     , obsDat = ASP.1.data
#                     , init = init.pop
#                     , tseq = time.out2
#                     , control = list(trace = 3, maxit = 800, reltol = 10^-7)
#                     , method = "Nelder-Mead"
#                     , hessian = T)

## View the R0 estimate
estimParms <- exp(unname(prison.estim$par))
names(estimParms) <- c('R0', 'N')
estimParms

## Run the SEIR model with the R0 estimate
prison.ts <- simEpidemic(init.prison, time.out2, seixc, disease_params(R0 = estimParms['R0'], N = estimParms['N']))

## Visually compare incidence of the model and actual data

df2 <- rbind(ASP.1.data %>% 
               select(-Date) %>% 
               mutate(variable = "ASP.1.data"),
             prison.ts %>% 
               select(incidence, time) %>% 
               mutate(variable = "prison.ts"))

ggplot(df2, aes(x = time, y = incidence, color = variable))+
  geom_line() +
  geom_point(data = subset(df2, df2$variable == "ASP.1.data"), mapping=aes(x = time, y=incidence, color = variable)) +
  labs(title = "Comparing incidence of data and model",
       subtitle = paste("Avenal State Prison Outbreak 1 data
Estimates: R0 = ", round(estimParms['R0'],2), ", N = ", round(estimParms['N'],2)),
x= "Time [days]",
y= "Cases") +
  scale_color_manual(labels = c("Data","Model"), values = c("red", "blue"))

## The model incidence and actual incidence clearly do not fit well based on
## visual inspection. The model cannot match the data incidence with two peaks,
## and is only fitting one peak.
## Let's try to cut the data to only include the initial peak, and then assess
## model performance.

# Estimate with one-peak data ---------------------------------------------

ASP.1.data.cut <- ASP.1.data %>%
  filter(time < 60)  # cutoff at day 60 of outbreak

time.out3 <- seq(0, max(ASP.1.data.cut$time), 1)

prison.estim.cut <- optim(par = c(log_R0 = log(1), log_N = log(1000))
                          , objFXN
                          , fixed.params = disease_params()
                          , obsDat = ASP.1.data.cut
                          , init = init.prison
                          , tseq = time.out3
                          , control = list(trace = 3, maxit = 150)
                          , method = "SANN"
                          , hessian = T)

## Feed the last parameters of SANN in as the first values of Nelder-Mead
# prison.estim.cut <- optim(par = optim.vals$par
#                     , objFXN
#                     , fixed.params = disease_params()
#                     , obsDat = ASP.1.data.cut
#                     , init = init.prison
#                     , tseq = time.out3
#                     , control = list(trace = 3, maxit = 800, reltol = 10^-7)
#                     , method = "Nelder-Mead"
#                     , hessian = T)

## Store estimates
MLEfits <- prison.estim.cut$par
log_R0.fit <- MLEfits["log_R0"]
log_N.fit <- MLEfits["log_N"]
estimParms2 <- exp(unname(MLEfits))
names(estimParms2) <- c('R0', 'N')

## See final estimates
estimParms2

## Run the SEIR model with the R0 estimate

prison.ts.cut <- simEpidemic(init.prison, time.out3, seixc, disease_params(R0 = estimParms2['R0'], N = estimParms2['N']))

## Visually compare incidence of the model and actual data

df3 <- rbind(ASP.1.data.cut %>% 
               select(-Date) %>% 
               mutate(variable = "ASP.1.data.cut"),
             prison.ts.cut %>% 
               select(incidence, time) %>% 
               mutate(variable = "prison.ts.cut"))

ggplot(df3, aes(x = time, y = incidence, color = variable))+
  geom_line() +
  geom_point(data = subset(df3, df3$variable == "ASP.1.data.cut"), mapping=aes(x = time, y=incidence, color = variable)) +
  labs(title = "Comparing incidence of data and model",
       subtitle = paste("Avenal State Prison Outbreak 1 data (cut off to one peak)
Estimates: R0 = ", round(estimParms2['R0'],2), ", N = ", round(estimParms2['N'],2)),
x= "Time [days]",
y= "Cases") +
  scale_color_manual(labels = c("Data","Model"), values = c("red", "blue"))

## The model is somewhat closer to the data in shape
## Are the estimated R0 and N values reasonable?

# Contour plots for actual data -------------------------------------------
## We will use the estimated R0 and N from the ASP1 data restricted to one peak

## Generate contour plots with the Hessian ##
fisherInfMatrix <- solve(prison.estim.cut$hessian)

## Create a sequence of R0 and N values for the grid
R0.seq <- exp(seq(log_R0.fit-1, log_R0.fit+1, l = res))
N.seq <- exp(seq(log_N.fit-1, log_N.fit + 1, l = res))

## Initialize plot of parameters
plot(1,1, type = 'n', log = 'xy',
     # xlim = range(R0.seq), ylim = range(N.seq),
     xlim = c(7.50,7.55), ylim = c(7980,8000),  # use these limits to zoom-in
     las = 1,
     xlab = expression(R0), ylab = expression(N),
     main = "-log(likelihood) contours", bty = "n")
## Add MLE to the plot
points(estimParms2['R0'], estimParms2['N'], pch = 16, cex = 2, col = 'black')
## Add 95% contour ellipse from Hessian
lines(exp(ellipse(fisherInfMatrix, centre = MLEfits, level = .95)))
legend("topleft", c('MLE', '95% Confidence Region'), lty = c(NA, 1), pch = c(16, NA),
       col = c('black', 'black'), bg='white', bty = 'n')

## Generate contour plots with likelihood profiles ##
mat <- outer(R0.seq, N.seq, objXR0_NVEC, init = init_cond(), tseq=time.out, obsDat=myDat)

ml.val <- prison.estim.cut$value
conf.cutoff <- ml.val + qchisq(.95,2)/2

## Show likelihood contours
par(cex = 1.2)
plot(1,1, type = 'n', log = 'xy',
     xlim = range(R0.seq), ylim = range(N.seq),
     # xlim = c(7.50,7.55), ylim = c(7980,8000),  # use these limits to zoom-in
     xlab = expression(R0), ylab = expression(N),
     main = "-log(likelihood) contours", bty = "n")
.filled.contour(R0.seq, N.seq, mat, levels = seq(min(mat), max(mat), l=20), col = topo.colors(20))
## Add contour for 95% CI from likelihood ratio
## This comes out really large!
contour(R0.seq, N.seq, mat, levels = c(conf.cutoff),
        col = "black", lwd = 2, labels = "", labcex = .2, add = T)
## Add contour for 95% CI from Hessian
lines(exp(ellipse(fisherInfMatrix, centre = MLEfits, level = .95)), lty = 2)
## Add MLE to the plot
points(exp(log_R0.fit), exp(log_N.fit), pch = 16, cex = 1, col = 'black')
## Add true parameter values to the plot
legend("topleft",
       c('MLE', '95% contour (profile likelihood)', '95% contour (Fisher information matrix)')
       , lty = c(NA, 1, 2), pch = c(16, NA, NA),
       col = c(rep('black',3)), bg='white', bty = 'n')

# Plots used in the report ------------------------------------------------

## Figure 1
ggplot(df, aes(x = time, y = incidence, color = variable))+
  geom_line() +
  geom_point(data = myDat, mapping=aes(x = time, y=incidence, color = variable)) +
  #  geom_ribbon(data = myDat, mapping=aes(x = time, y = samp_inc_rate, ymin=lci, ymax=uci), color = NA, fill = "black", alpha = 0.2) +
  labs(x= "Time [days]",
       y= "Cases") +
  scale_color_manual(labels = c("Fitted","Truth","Data"), values = c("red", "blue", "black")) +
  theme(legend.title= element_blank())

## Figure 2
A <- ggplot(df2, aes(x = time, y = incidence, color = variable))+
  geom_line() +
  geom_point(data = subset(df2, df2$variable == "ASP.1.data"), mapping=aes(x = time, y=incidence, color = variable)) +
  labs(x= "Time [days]",
       y= "Cases") +
  scale_color_manual(labels = c("Data","Model"), values = c("red", "blue")) +
  theme(legend.position = "bottom", legend.title= element_blank())

B <- ggplot(df3, aes(x = time, y = incidence, color = variable))+
  geom_line() +
  geom_point(data = subset(df3, df3$variable == "ASP.1.data.cut"), mapping=aes(x = time, y=incidence, color = variable)) +
  labs(x= "Time [days]",
       y= "Cases") +
  scale_color_manual(labels = c("Data","Model"), values = c("red", "blue"))

AB <- plot_grid(A+theme(legend.position = "none"), 
                B + theme(legend.position = "none"), 
                labels = "AUTO", ncol=2)

legend <- get_legend(A)

plot_grid(AB, legend, ncol=1,rel_heights = c(0.95, 0.05))
