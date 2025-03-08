# Load necessary library
library(dplyr)

# Define function to calculate Lq and P0 for M/M/c queue
calculate_Lq_and_P0 <- function(lambda, mu, c, num_simulations = 10000) {
  rho <- lambda / (c * mu)

  if (rho < 1) {
    # Use standard formula when system is stable
    p0 <- 1 / (sum(sapply(0:(c-1), function(n) (lambda / mu)^n / factorial(n))) +
                 ((lambda / mu)^c / (factorial(c) * (1 - rho))))
    Lq <- (p0 * (lambda / mu)^c * rho) / (factorial(c) * (1 - rho)^2)
  } else {
    # Simulate the queue for cases where rho >= 1
    queue_length <- numeric(num_simulations)
    current_queue <- 0

    for (i in 1:num_simulations) {
      arrivals <- rpois(1, lambda)  # Poisson distributed arrivals
      departures <- min(current_queue + arrivals, c)  # Limited by servers
      current_queue <- max(0, current_queue + arrivals - departures)
      queue_length[i] <- current_queue
    }

    Lq <- mean(queue_length)
    p0 <- mean(queue_length == 0)  # Estimate P0 using the simulation (percentage of time system is empty)
  }

  return(list(Lq = Lq, P0 = p0))
}

# Define peak demand multipliers
peak_multipliers <- c(1.25, 1.5, 1.75, 2.0)

# Define constants
num_standard_servers <- 7  # Standard registers
num_priority_servers <- 1  # Priority Pass register
lambda <- 1 / 0.275606273  # Arrival rate - mean
mu <- 1 / 1.420613077  # Service rate - mean

# Initialize results list
peak_simulation_results <- list()

for (multiplier in peak_multipliers) {
  # Adjust arrival rates for peak times
  lambda_peak <- multiplier * lambda
  lambda_standard_peak <- 0.95 * lambda_peak
  lambda_priority_peak <- 0.05 * lambda_peak

  # Compute new utilizations
  rho_standard_peak <- lambda_standard_peak / (num_standard_servers * mu)
  rho_priority_peak <- lambda_priority_peak / (num_priority_servers * (1.25 * mu))

  # Calculate new queue lengths and P0 values
  results_standard <- calculate_Lq_and_P0(lambda_standard_peak, mu, num_standard_servers)
  Lq_standard_peak <- results_standard$Lq
  P0_standard_peak <- results_standard$P0

  results_priority <- calculate_Lq_and_P0(lambda_priority_peak, 1.25 * mu, num_priority_servers)
  Lq_priority_peak <- results_priority$Lq
  P0_priority_peak <- results_priority$P0

  # Calculate new wait times
  Wq_standard_peak <- Lq_standard_peak / lambda_standard_peak
  Wq_priority_peak <- Lq_priority_peak / lambda_priority_peak

  # Calculate new total time in system
  W_standard_peak <- Wq_standard_peak + (1 / mu)
  W_priority_peak <- Wq_priority_peak + (1 / (1.25 * mu))

  # Store results
  peak_simulation_results <- rbind(peak_simulation_results,
                                   data.frame(Peak_Multiplier = multiplier,
                                              Queue_Type = "Standard Queue",
                                              Utilization = rho_standard_peak,
                                              P0 = P0_standard_peak,
                                              Expected_Queue_Length_Lq = Lq_standard_peak,
                                              Expected_Wait_Time_Wq = Wq_standard_peak,
                                              Expected_Total_Time_W = W_standard_peak))

  peak_simulation_results <- rbind(peak_simulation_results,
                                   data.frame(Peak_Multiplier = multiplier,
                                              Queue_Type = "Priority Pass Queue",
                                              Utilization = rho_priority_peak,
                                              P0 = P0_priority_peak,
                                              Expected_Queue_Length_Lq = Lq_priority_peak,
                                              Expected_Wait_Time_Wq = Wq_priority_peak,
                                              Expected_Total_Time_W = W_priority_peak))
}

# Display results
df_results <- as.data.frame(peak_simulation_results)
print(df_results)

