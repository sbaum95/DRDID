devtools::install_github("sbaum95/DRDID")
library(tidyverse)
library(fixest)



# Multi-group simulation --------------------------------------------------
# Create fake data 
set.seed(123) # For reproducibility

# Parameters
n_units <- 100
n_periods <- 50
treated_units <- 1
control_units <- n_units - treated_units
treatment_period <- 5

# Simulate population sizes for each unit (e.g., between 1,000 and 10,000)
population <- sample(1000:10000, n_units, replace = TRUE)

# Baseline incidence rate: cases per population (use a reasonable baseline rate)
baseline_rate <- 500 / 100000

# The treatment effect (an additional 1000 cases per 100,000 population)
treatment_effect <- 1000 / 100000

# Empty data frame to store results
simulated_data <- data.frame()

# Simulating data
for (unit in 1:n_units) {
  for (period in 1:n_periods) {
    # Determine if the unit is treated and if treatment has started
    is_treated <- unit <= treated_units
    is_treatment_period <- period >= treatment_period
    
    # Calculate the incidence rate
    incidence_rate <- baseline_rate
    if (is_treated && is_treatment_period) {
      incidence_rate <- incidence_rate + treatment_effect
    }
    
    # Simulate cases using Poisson distribution
    cases <- rpois(1, lambda = incidence_rate * population[unit])
    
    
    # Add to dataset
    simulated_data <- rbind(simulated_data, data.frame(
      unit = unit,
      period = period,
      treated = as.integer(is_treated),
      population = population[unit],
      cases = cases
    ))
  }
}

#'@ param dname The name of the column containing the treatment group (=1 if observation is treated in the post-treatment, =0 otherwise)
simulated_data <- simulated_data %>% 
  mutate(D = if_else(treated == 1 & period >= 5, 1, 0))

# View the first few rows of the simulated dataset
dplyr::glimpse(simulated_data)



# 2 x 2 Simulation --------------------------------------------------------
dat_2x2 = data.frame(
  id = rep(1:2, times = 2), 
  tt = rep(1:2, each = 2), 
  pop = c(10000, 10000, 10000, 10000)
)

## How do we force a model 
dat_2x2 <- dat_2x2 %>% 
  within({
    D     = id == 2 & tt == 2
    btrue = ifelse(D, 0.1, 0)
    
    # log-linear model 
    # Counts - log(mu(x)) = XtBeta (where mu = exp(XtBeta)
    # Rates -  log(mu(x)) = XtBeta + log(pop) 
    log_cases = (id + 0.005*tt + btrue * D) + log(pop)
    
    # Expected cases 
    cases = exp(log_cases)
    incidence = cases/pop
  })

ggplot(dat_2x2, aes(x = tt, y = log(incidence), col = factor(id))) + 
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 1.5, lty = 2) +
  scale_x_continuous(breaks = 1:2, labels = c("Pre", "Post")) +
  labs(x = "Time variable", y = "Outcome variable", col = "ID")



## TWFE
twfe_ols_2x2 = feols(log_cases ~ D | id + tt, data = dat_2x2)
twfe_pois_2x2 = fepois(cases ~ D | id + tt, offset = ~log(pop), data = dat_2x2)


# 2x2 with multiple time periods -------------------------
dat_multi_tt = data.frame(
  id = rep(1:2, times = 10), # Two treatment group
  tt = rep(1:10, each = 2),  # Observe each group 10 times 
  pop = c(10000, 10000, 10000, 10000)
)

dat_multi_tt <- dat_multi_tt %>% 
  within({
    D     = id == 2 & tt >= 5
    btrue = ifelse(D, 0.1, 0)
    
    # log-linear model 
    # Counts - log(mu(x)) = XtBeta (where mu = exp(XtBeta)
    # Rates -  log(mu(x)) = XtBeta + log(pop) 
    log_cases = (id + 0.005*tt + btrue * D) + log(pop)
    
    # Expected cases 
    cases = exp(log_cases)
    incidence = cases/pop
  })

ggplot(dat_multi_tt, aes(x = tt, y = log(incidence), col = factor(id))) + 
  geom_point() +
  geom_line() +
  geom_line(dat_multi_tt %>% 
              mutate(counter_t = ), aes(x = tt, y = log(incidence), col = factor(id))) +
  geom_vline(xintercept = 4.5, lty = 2) +
  labs(x = "Time variable", y = "Outcome variable", col = "ID")


## Multi 
twfe_pois_multi_tt = feols(log(incidence) ~ D | id + tt, data = dat_multi_tt)
twfe_pois_multi_tt = fepois(cases ~ D | id + tt, offset = ~log(pop), data = dat_multi_tt)


# 2x2 with multiple ids  -------------------------
dat_multi_id = data.frame(
  id = rep(1:3000, times = 2), # Ten treatment group
  tt = rep(1:2, each = 3000),  # Observe each group 2 times 
  offset = rep(10000, times = 3000)
)

dat_multi_id <- dat_multi_id %>% 
  within({
    D     = id >= 750 & tt == 2
    btrue = ifelse(D, 0.5, 0)
    
    # log-linear model 
    log_cases = (0.005*tt + btrue * D) + log(offset)
    
    # Expected cases 
    cases = exp(log_cases)
    incidence = cases/offset
  })

ggplot(dat_multi_id, aes(x = tt, y = log(incidence), col = factor(id))) + 
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 1.5, lty = 2) +
  labs(x = "Time variable", y = "Outcome variable", col = "ID")


## Multi 
twfe_pois_multi_id = feols(log(incidence) ~ D | id + tt, data = dat_multi_id)
twfe_pois_multi_id = fepois(cases ~ D | id + tt, offset = ~log(pop), data = dat_multi_id)

# Staggered treatment + multiple time periods -------------------------
dat_multi_gtt = data.frame(
  id = rep(1:4, times = 10), # Two treatment group
  tt = rep(1:10, each = 4),  # Observe each group 10 times 
  pop = rep(10000, each = 4*10)
)

dat_multi_gtt <- dat_multi_gtt %>% 
  within({
    D     = id >= 2 & tt >= 5
    btrue = ifelse(D, 0.1, 0)
    
    # log-linear model 
    # Counts - log(mu(x)) = XtBeta (where mu = exp(XtBeta)
    # Rates -  log(mu(x)) = XtBeta + log(pop) 
    log_cases = (id + 0.005*tt + btrue * D) + log(pop)
    
    # Expected cases 
    cases = exp(log_cases)
    incidence = cases/pop
  })

dat_multi_gtt <- dat_multi_gtt %>% mutate(treated = if_else(id >=25, 1, 0))

ggplot(dat_multi_gtt, aes(x = tt, y = incidence, col = factor(id))) + 
  geom_point() +
  geom_line() +
  geom_line(dat_multi_tt %>% 
              mutate(counter_t = ), aes(x = tt, y = log(incidence), col = factor(id))) +
  geom_vline(xintercept = 4.5, lty = 2) +
  labs(x = "Time variable", y = "Outcome variable", col = "ID")


## Multi 
twfe_pois_multi_gtt = feols(log(incidence) ~ D | id + tt, data = dat_multi_gtt)
twfe_pois_multi_gtt = fepois(cases ~ D | id + tt, offset = ~log(pop), data = dat_multi_gtt)




# DRDID -------------------------------------------------------------------

#' # Implement improved DR locally efficient DiD with panel data
#' drdid(yname="re", tname = "year", idname = "id", dname = "experimental",
#'       xformla= ~ age+ educ+ black+ married+ nodegree+ hisp+ re74,
#'       data = eval_lalonde_cps, panel = TRUE)
#'
#' #Implement "traditional" DR locally efficient DiD with panel data
#' drdid(yname="re", tname = "year", idname = "id", dname = "experimental",
#'       xformla= ~ age+ educ+ black+ married+ nodegree+ hisp+ re74,
#'       data = eval_lalonde_cps, panel = TRUE, estMethod = "trad")
#'


## Actuatl DRDID
remotes::install_github("pedrohcgs/DRDID")
library(DRDID)
data(nsw_long)
## Make D numeric
dat_multi_tt <- dat_multi_tt %>% mutate(treated = if_else(id == 1, 1, 0))

dat_2x2 <- dat_2x2 %>% mutate(treated = if_else(id == 1, 1, 0))

# Form the Lalonde sample with CPS comparison group
eval_lalonde_cps <- subset(nsw_long, nsw_long$treated == 0 | nsw_long$sample == 2)
glimpse(eval_lalonde_cps)
drdid(yname = "re", tname = "year", idname = "id", dname = "experimental", data = eval_lalonde_cps)

dat_multi_id <- dat_multi_id %>% mutate(treated = if_else(id >= 750, 1, 0))
glimpse(dat_multi_id)
drdid(yname = "log_cases", tname = "tt", idname = "id", dname = "treated", data = dat_multi_id)


## Actuatl DRDID
remotes::install_github("sbaum95/DRDID")
#remotes::install_github("sbaum95/did")
library(DRDID)

drdid(yname = "cases", tname = "tt", idname = "id", dname = "treated", 
      offset = "pop", panel = TRUE, data = dat_multi_id)

dat_multi_id <- dat_multi_id %>% mutate(treated = if_else(id >= 750, 1, 0))
dp <- pre_process_drdid(yname = "cases", tname = "tt", idname = "id", dname = "treated", 
                  offset = "offset", panel = TRUE, data = dat_multi_id, estMethod = "imp")










