# Load necessary library
library(dplyr)

# Define function to calculate Lq for M/M/c/priority queue
calculate_Lq_priority <- function(lambda_standard, lambda_priority, mu, c, priority_boost = 1.25, num_simulations = 10000) {
  lambda_total <- lambda_standard + lambda_priority
  rho <- lambda_total / (c * mu)

  if (rho < 1) {
    # Compute P0 (steady-state probability of zero customers in system)
    p0 <- 1 / (sum(sapply(0:(c-1), function(n) (lambda_total / mu)^n / factorial(n))) +
                 ((lambda_total / mu)^c / (factorial(c) * (1 - rho))))

    # Compute Lq for total system
    Lq_total <- (p0 * (lambda_total / mu)^c * rho) / (factorial(c) * (1 - rho)^2)

    # Approximate Lq split using weighted arrival rates
    Lq_priority <- (lambda_priority / lambda_total) * Lq_total * priority_boost
    Lq_standard <- Lq_total - Lq_priority
  } else {
    # Simulation for unstable systems
    queue_length <- numeric(num_simulations)
    priority_queue <- numeric(num_simulations)
    standard_queue <- numeric(num_simulations)
    current_priority <- 0
    current_standard <- 0

    for (i in 1:num_simulations) {
      arrivals_standard <- rpois(1, lambda_standard)
      arrivals_priority <- rpois(1, lambda_priority)

      total_arrivals <- arrivals_standard + arrivals_priority
      departures <- min(current_standard + current_priority + total_arrivals, c)

      # Priority customers get served first
      served_priority <- min(current_priority + arrivals_priority, departures)
      served_standard <- min(departures - served_priority, current_standard + arrivals_standard)

      current_priority <- max(0, current_priority + arrivals_priority - served_priority)
      current_standard <- max(0, current_standard + arrivals_standard - served_standard)

      priority_queue[i] <- current_priority
      standard_queue[i] <- current_standard
      queue_length[i] <- current_priority + current_standard
    }

    Lq_standard <- mean(standard_queue)
    Lq_priority <- mean(priority_queue)
    p0 <- mean(queue_length == 0)  # Probability of no customers (P0) from simulation
  }

  return(list(Lq_standard = Lq_standard, Lq_priority = Lq_priority, P0 = p0))
}

# Define peak demand multipliers
peak_multipliers <- c(1.25, 1.5, 1.75, 2.0)

# Define constants
num_servers <- 8  # Total servers
lambda <- 1 / 0.275606273  # Arrival rate - mean
mu <- 1 / 1.420613077  # Service rate - mean

# Initialize results list
peak_simulation_results <- list()

for (multiplier in peak_multipliers) {
  # Adjust arrival rates for peak times
  lambda_peak <- multiplier * lambda
  lambda_standard_peak <- 0.95 * lambda_peak
  lambda_priority_peak <- 0.05 * lambda_peak

  # Compute utilizations
  rho_total_peak <- (lambda_standard_peak + lambda_priority_peak) / (num_servers * mu)

  # Compute Lq values using priority queuing function
  Lq_results <- calculate_Lq_priority(lambda_standard_peak, lambda_priority_peak, mu, num_servers)
  Lq_standard_peak <- Lq_results$Lq_standard
  Lq_priority_peak <- Lq_results$Lq_priority
  P0_peak <- Lq_results$P0

  # Compute wait times
  Wq_standard_peak <- Lq_standard_peak / lambda_standard_peak
  Wq_priority_peak <- Lq_priority_peak / lambda_priority_peak

  # Compute total time in system
  W_standard_peak <- Wq_standard_peak + (1 / mu)
  W_priority_peak <- Wq_priority_peak + (1 / (1.25 * mu))

  # Store results
  peak_simulation_results <- rbind(peak_simulation_results,
                                   data.frame(Peak_Multiplier = multiplier,
                                              Queue_Type = "Standard Queue",
                                              Utilization = rho_total_peak,
                                              P0 = P0_peak,
                                              Expected_Queue_Length_Lq = Lq_standard_peak,
                                              Expected_Wait_Time_Wq = Wq_standard_peak,
                                              Expected_Total_Time_W = W_standard_peak))

  peak_simulation_results <- rbind(peak_simulation_results,
                                   data.frame(Peak_Multiplier = multiplier,
                                              Queue_Type = "Priority Pass Queue",
                                              Utilization = rho_total_peak,
                                              P0 = P0_peak,
                                              Expected_Queue_Length_Lq = Lq_priority_peak,
                                              Expected_Wait_Time_Wq = Wq_priority_peak,
                                              Expected_Total_Time_W = W_priority_peak))
}

# Display results
df_results <- as.data.frame(peak_simulation_results)
print(df_results)
