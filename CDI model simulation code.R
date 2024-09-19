# HA-CDI Simulation Model 
# Michael J Ray, PhD MPH 
# April 2024

# Load necessary library
library(GillespieSSA2)

# Set up RTOOLS path - RTOOLS is required to run this simulation 
write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
Sys.which("make")

# Simulation name
sim_name <- "HA-CDI simulation model"

# Define model parameters
params <- c(
  alpha = 85,
  xiUL = .787,
  xiEL = .045,
  xiUH = .150,
  xiEH = .015,
  xiD = .003,
  thetaL = .2,
  thetaH = .167, 
  kappaL = .004,
  kappaH = .02,
  iota = .034,
  gamma = .12, 
  omega = .03, 
  psi = .1, 
  chi = .72, 
  lambda = .011, 
  beta = 0.009454, # Fitted value using Approximate Bayesian Computation 
  phi = .25,
  zeta = .1,
  mu = .1,
  rho = 1.5,
  tau = 1.33
)

# Define simulation time and initial state
final_time <- 365
initial_state <- c(UL = 379, UH = 90, EL = 20, EH = 10, D = 1, DR = 0)

# Define reactions
reactions <- list(
  reaction("mu * UL", c(UL = -1, UH = +1), name = "abx_unexposed"),
  reaction("mu * EL", c(EL = -1, EH = +1), name = "abx_exposed"), 
  reaction("lambda * UH", c(UH = -1, UL = +1), name = "flora_recovery_unexposed"), 
  reaction("lambda * EH", c(EH = -1, EL = +1), name = "flora_recovery_exposed"), 
  reaction("(beta * UL) * ((EH + EL + (tau * D)) / (UL + UH + EL + EH + D))", c(UL = -1, EL = +1), name = "exposure_low_risk"), 
  reaction("(beta * rho * UH) * ((EH + EL + (tau * D)) / (UL + UH + EL + EH + D))", c(UH = -1, EH = +1), name = "exposure_high_risk"), 
  reaction("kappaL * EL", c(EL = -1, D = +1), name = "dev_symptoms_low_risk"),
  reaction("kappaH * EH", c(EH = -1, D = +1), name = "dev_symptoms_high_risk"),
  reaction("zeta * phi * psi  * D", c(D = -1, EH = +1), name = "post_cdi_return_to_exposed"),
  reaction("(1 - zeta) * phi * psi * D", c(D = -1, UH = +1), name = "post_cdi_return_to_unex"), 
  reaction("chi * psi * D", c(D = -1, DR = +1), name = "discharge_from_CDI"), 
  reaction("gamma * iota * DR", c(DR = -1, D = +1), name = "recurrent_CDI"), 
  reaction("alpha * xiUL", c(UL = +1), name = "admit_to_UL"), 
  reaction("alpha * xiUH", c(UH = +1), name = "admit_to_UH"), 
  reaction("alpha * xiEL", c(EL = +1), name = "admit_to_EL"), 
  reaction("alpha * xiEH", c(EH = +1), name = "admit_to_EH"), 
  reaction("alpha * xiD", c(D = +1), name = "admit_to_D"),  
  reaction("omega * psi * D", c(D = -1), name = "death"), 
  reaction("thetaL * UL", c(UL = -1), name = "discharge_from_UL"), 
  reaction("thetaH * UH", c(UH = -1), name = "discharge_from_UH"), 
  reaction("thetaL * EL", c(EL = -1), name = "discharge_from_EL"), 
  reaction("thetaH * EH", c(EH = -1), name = "discharge_from_EH"), 
  reaction("(1 - gamma) * iota * DR", c(DR = -1), name = "non_recurr_cdi")
)

# Pre-compile reactions to run faster simulations
compiled_reactions <- compile_reactions(
  reactions = reactions,
  state_ids = names(initial_state),
  params = params
)

# Number of repetitions
num_repetitions <- 10

# Initialize lists to store results
output_list <- vector("list", length = num_repetitions)
patient_days_list <- numeric(num_repetitions)
incident_cdi_list <- numeric(num_repetitions)
cdi_incidence_rate_list <- numeric(num_repetitions)



# Run simulations
for (i in 1:num_repetitions) {
  out <- ssa(
    initial_state = initial_state,
    reactions = compiled_reactions,
    params = params,
    final_time = final_time,
    method = ssa_exact(),
    sim_name = sim_name,
    census_interval = 1, 
    log_firings = TRUE,
    verbose = TRUE
  )
  
  output_list[[i]] <- out$state
  
  # Calculate N (total patients)
  N <- rowSums(out$state[, c("UL", "UH", "EL", "EH", "D")])
  # 90-day burn-in period 
  NN <- N[90:366]
  
  # Calculate patient days for this simulation
  patient_days <- sum(NN)
  patient_days_list[i] <- patient_days
  
  # Extract CDI events
  cdi_events <- out$firings[90:366, c("dev_symptoms_high_risk", "dev_symptoms_low_risk")]
  
  # Calculate incident CDI for this simulation
  incident_cdi <- sum(cdi_events)
  incident_cdi_list[i] <- incident_cdi
  
  # Calculate CDI incidence rate for this simulation
  cdi_incidence_rate <- (incident_cdi / patient_days) * 10000
  cdi_incidence_rate_list[i] <- cdi_incidence_rate
}

# Output as data frame
output_df <- data.frame(
  Total_patient_days = patient_days_list, 
  Incident_cdi_count = incident_cdi_list, 
  Cdi_incidence_rate = cdi_incidence_rate_list
)

# Calculate mean CDI incidence rate 
mean_cdi_incidence_rate <- mean(cdi_incidence_rate_list)
mean_cdi_incidence_rate

# Calculate the standard error of the CDI incidence rate
cdi_incidence_se <- sd(output_df$Cdi_incidence_rate) / sqrt(num_repetitions)

# Calculate the margin of error for a 95% confidence interval
margin_of_error <- 1.96 * cdi_incidence_se 

# Calculate the 95% confidence interval
confidence_interval <- c(mean_cdi_incidence_rate - margin_of_error, mean_cdi_incidence_rate + margin_of_error)

# Print the results
cat("Mean HA-CDI Incidence Rate:", mean_cdi_incidence_rate, "\n")
cat("95% Confidence Interval:", confidence_interval, "\n")

