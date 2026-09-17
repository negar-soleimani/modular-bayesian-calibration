###############################################################################
# Article 2 - Full / Cut / Plug-in
###############################################################################

library(mvtnorm)
library(invgamma)
library(truncnorm)
library(coda)
library(ggplot2)

n_total <- 45
m_phys <- 15

g_true <- 9.8
h0_true <- 46.45
sigma_sq_err_true <- 0.01

delta_amplitude <- 2.5

SCENARIOS <- c(
  "no_discrepancy",
  "weak_smooth_positive",
  "moderate_smooth_positive",
  "delayed_smooth_positive",
  "moderate_local_bump",
  "moderate_quadratic_confounding"
)

SCENARIO_LABELS <- c(
  no_discrepancy = "Absence de discrépance",
  weak_smooth_positive = "Discrépance lisse faible",
  moderate_smooth_positive = "Discrépance lisse modérée",
  delayed_smooth_positive = "Discrépance à apparition retardée",
  moderate_local_bump = "Discrépance locale en forme de bosse",
  moderate_quadratic_confounding = "Discrépance quadratique"
)

SCENARIO_NUMBERS <- c(
  no_discrepancy = 1,
  weak_smooth_positive = 2,
  moderate_smooth_positive = 3,
  delayed_smooth_positive = 4,
  moderate_local_bump = 5,
  moderate_quadratic_confounding = 6
)

METHOD_LEVELS <- c(
  "Full Mixture",
  "Conditional Cut Mixture",
  "Conditional Cut No-Mixture",
  "Plug-in No-Mixture"
)

REPLICATE_ID <- 1

# Seeds for the six scenarios

SEED_DATA <- c(101, 102, 103, 104, 105, 106)

SEED_MODULE_P <- c(201, 202, 203, 204, 205, 206)

SEED_PHYS_IDS <- c(301, 302, 303, 304, 305, 306)

SEED_FULL <- c(401, 402, 403, 404, 405, 406)

SEED_CUT_MIX <- c(501, 502, 503, 504, 505, 506)

SEED_CUT_NOMIX <- c(601, 602, 603, 604, 605, 606)

SEED_PLUGIN <- c(701, 702, 703, 704, 705, 706)

# Module P
N_KEEP_PHYS <- 30000
BURN_PHYS <- 3000

# Conditional Cut
N_OUTER_PHYS <- 200
BURN_D_PER_PHYS <- 1000
N_KEEP_D_PER_PHYS <- 1000

# Plug-in
BURN_PLUGIN <- 3000
N_KEEP_PLUGIN <- 5000

# Full Mixture
BURN_FULL <- 3000
N_KEEP_FULL <- 50000

init_phys <- c(g = 9.8, h0 = 46.45, sigma_sq_err = 0.01)

init_D <- c(psi_delta = 0.5, k = 0.1)

init_mix <- c(psi_delta = 0.5, k = 0.1, alpha = 0.5)

PSI_LOWER <- 0.10
PSI_UPPER <- 0.50
PSI_PROPOSAL_SD <- 0.50

K_LOWER <- 0.02
K_UPPER <- 1.00

set.seed(12345)

t <- seq(0, 1, length.out = n_total)

t <- sort(t + rnorm(n_total, mean = 0, sd = 0.003))

t <- (t - min(t)) / (max(t) - min(t))

t_min <- min(t)
t_range <- max(t) - min(t)

n <- length(t)

I_P <- seq_len(m_phys)
I_D <- seq.int( from = m_phys + 1, to = n)

n_P <- length(I_P)
n_D <- length(I_D)