# Load reticulate
library(reticulate)

# Import Python jsf package
jsf <- import("jsf")

# Initial state
x0 <- c(10, 0)

# Define reaction rates
rates <- function(state, time) {
  return(c(0.1 * state[1], 0.05 * state[2]))
}

# Stoichiometry
stoich <- list(
  nu = list(c(-1, 1), c(1, -1)),
  nuReactant = list(c(1, 0), c(0, 1))
)

# Simulation time
t_max <- 10

# Configuration
config <- list(
  dt = 0.1,
  EnforceDo = c(1, 1),
  SwitchingThreshold = c(10, 10)
)

# Run simulation using operator splitting
result <- jsf$jsf(
  x0,
  rates,
  stoich,
  t_max,
  method = "operator-splitting",
  config = config
)

# Print result
print(result)
