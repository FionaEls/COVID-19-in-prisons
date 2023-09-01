# Updated simulation and inference code

## import required packages
library(tidyverse) 
library(stats4)
############################################################
# create a function called log_g_ij that calculates
# log_prob that i cases cause j cases in a single generation
log_g_ij <- function(i,j,rval,kval) {
  lgamma(j +kval*i) - lgamma(j+1) - lgamma(kval*i) + 
    kval * i * log(kval/(rval+kval)) +
    j * log(rval/(rval+kval))
}
# log_prob that m primaries causes n cases in totality
log_r_mn <- function (m,n,rval,kval) {
  log(m/n) + log_g_ij(i=n, j = n-m, rval = rval, kval = kval)
}

# function to simulate the infected caes in current generation and total cases in chain so far
simple_simulation <- function (R, k) {
  
  total_cases <- 1 # total cases in chain so far
  cur_cases <- 1 # infected cases in current generation
  while (cur_cases > 0 & total_cases < OUTBREAK_THRESHOLD) { # run this code while infected cases in the current generation are greater than zero and within the outbreak threshold
    cur_cases <- sum(rnbinom(n = cur_cases, size = k, mu = R)) # assign the new value of infected cases in current generation  
    total_cases <- total_cases + cur_cases # assign the new value of total cases in the chain so far
  }
  cur_cases <- min(total_cases,OUTBREAK_THRESHOLD) # otherwise infected cases in the current generation are a minimum between the total cases in the chain so far, and the outbreak threshold
}

# function to simualte introduction of cases in the population using the simple_simulation() function created above
simulate_intros <- function (R,k) {
  tibble(trial = 1:NUM_TRIALS) %>% group_by (trial) %>% do({
    temp_res <- simple_simulation(R = R, k = k)
    tibble(sim_result = temp_res)
  })  
}
# function to calcualte the negative log likelihood
calc_nLL <- function(R, k) {
  small_sizes <- 1:(OUTBREAK_THRESHOLD-1) # create a vector the size of the outbreak
  nLL_vec <- log_r_mn(1,small_sizes,R,k) #calculate the negative log likelihood vector using the log_prob fucntion that m primaries causes n cases in totality
  nLL_tibble <- tibble(sim_result = small_sizes, nLL = nLL_vec) # create a dataframe with the outbreak size, and the negative log likelihood for each
  outbreak_prob <- 1- sum(exp(nLL_vec)) # calculate the probability of the outbreak which is 1 minus the small blips
  nLL_tibble <- rbind(nLL_tibble, tibble(sim_result = OUTBREAK_THRESHOLD, nLL = log(outbreak_prob))) # put the results in a dataframe
  results <- left_join(sum_results, nLL_tibble, by = 'sim_result') # join the two datasets containing negative log likeihood and the simulated results
  nLL <- -sum(results$nLL * results$n)
  nLL
}
R_TRUE <- 1.2 # the reproductive number
K_TRUE <- 0.5 # the heterogeneity
NUM_TRIALS <- 2000 # number of trials
OUTBREAK_THRESHOLD <- 10 # outbreak threshold
sim_results <- simulate_intros(R = R_TRUE, k = K_TRUE)
sum_results <- sim_results %>% group_by(sim_result) %>% summarize(n = n())
(fit <- mle(calc_nLL, start = list(R = 1, k = 1))) # Use the MLE function to calculate the value of R and k from the data
confint(fit) # find the confidence intervals of R and k

outbreak_size_data <- read_csv(file.choose())


## use actual data

simulate_intros1 <- function (R,k) {
  outbreak_size_data %>% group_by (out) %>% summarize(n = n())#do({
  # temp_res <- simple_simulation(R = R, k = k)
  #  tibble(sim_result = temp_res)
  #  })  
}

sim_results <- simulate_intros1(R = R_TRUE, k = K_TRUE)

sum_results <- outbrsumeak_size_data %>% group_by(chain_size) %>% summarize(n = n()) %>%rename(sim_result=chain_size)
(fit <- mle(calc_nLL, start = list(R = 1, k = 1))) # Use the MLE function to calculate the value of R and k from the data
confint(fit) # find the confidence intervals of R and k

