### Get libraries
library(reticulate)
library(hmer)

### Pull in the source file for HPVsim
## Pre-install
# py_install("hpvsim", pip = TRUE)
source_python("Documents/HPVsim/runFunction.py")

### Testing the reticulate stuff
# file_in <- "Documents/HPVsim/TestV2.csv"
# file_out <- "Documents/HPVsim/TestV2Out.csv"
# main(file_in, file_out, n_reps = 10)

## Helper Function - target_hits
target_hits <- function(result, targets, sum_func = "sum") {
  hits <- purrr::map_lgl(names(targets), function(t) {
    if (is.atomic(targets[[t]]))
      return(result[t] <= targets[[t]][2] && result[t] >= targets[[t]][1])
    result[t] <= targets[[t]]$val + 3*targets[[t]]$sigma && result[t] >= targets[[t]]$val - 3*targets[[t]]$sigma
  })
  sum_function <- get(sum_func)
  return(sum_function(hits))
}

## Setting up parameter ranges, targets, and discrepancy
# Ranges
ranges = list(
  beta = c(0.02, 0.25),
  dur_episomal_16 = c(3, 10),
  dur_episomal_p2_16 = c(5, 15),
  dur_episomal_18 = c(3, 10),
  dur_episomal_p2_18 = c(5, 15),
  dur_episomal_hi5 = c(3, 10),
  dur_episomal_p2_hi5 = c(5, 15),
  sev_rate_16 = c(0.1, 0.5),
  sev_rate_18 = c(0.1, 0.5),
  sev_rate_hi5 = c(0.1, 0.5),
  dur_precin_16 = c(0.25, 4),
  dur_precin_18 = c(0.25, 4),
  dur_precin_hi5 = c(0.25, 4),
  trafo_prob_16 = c(1e-10, 1e-8),
  trafo_prob_18 = c(1e-10, 1e-8),
  trafo_prob_hi5 = c(1e-10, 1e-8),
  debut_female = c(12, 19),
  hpv_control_prob = c(0, 1),
  hpv_reactivation = c(0, 0.15)
)

# Targets
target_data <- read.csv("Documents/HPVsim/targets.csv", header = TRUE,
                        stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
targets <- purrr::map(seq_len(nrow(target_data)), ~c(target_data[.,"lower"], target_data[.,"upper"])) |>
  setNames(target_data$name)
targets <- targets[-c(1, 16)]

# Discrepancy
disc_internal <- read.csv("Documents/HPVsim/int_discrep.csv", header = TRUE,
                          stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
disc_int_usage <- purrr::map(seq_len(nrow(disc_internal))[-c(1,16)], ~disc_internal[.,]) |>
  setNames(names(targets))
ext_int_usage <- purrr::map_dbl(targets, mean)*0.095/3

## Initial setup
lhs_train <- lhs::maximinLHS(10*length(ranges), length(ranges))
initial_points <- data.frame(t(apply(lhs_train, 1, function(x) {
  x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
}))) |> setNames(names(ranges))
lhs_valid <- lhs::maximinLHS(5*length(ranges), length(ranges))
valid_points <- data.frame(t(apply(lhs_valid, 1, function(x) {
  x * purrr::map_dbl(ranges, diff) + purrr::map_dbl(ranges, ~.[[1]])
}))) |> setNames(names(ranges))
write.csv(rbind.data.frame(initial_points, valid_points), "Documents/HPVsim/points/wave0.csv", row.names = FALSE)

# Global parameter definitions
n_reps <- 16
n_waves <- 0
em_list <- list()
wave_list <- list()
wave_list_agg <- list()

## Waves of emulation
for (i in 0:n_waves) {
  # Get the points
  input_file <- paste0("Documents/HPVsim/points/wave", i, ".csv")
  output_file <- paste0("Documents/HPVsim/output/wave", i, ".csv")
  # Run the points, save to lists, and partition
  run_wave <- main(input_file, output_file, n_reps)
  this_wave <- read.csv(output_file, header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  wave_list[[i+1]] <- this_wave
  wave_list_agg[[i+1]] <- this_wave |> dplyr::group_by(across(all_of(names(ranges)))) |>
    dplyr::summarise(across(everything(), mean))
  train <- this_wave[1:(10*n_reps*length(ranges)),]
  valid <- this_wave[(10*n_reps*length(ranges)+1):nrow(this_wave),]
  # Train the emulators
  ems_this_wave <- emulator_from_data(train, names(targets), ranges, emulator_type = 'variance', check.ranges = TRUE, verbose = FALSE)
  # Define internal and external discrepancy
  for (j in seq_along(ems_this_wave$expectation)) {
    em_name <- ems_this_wave$expectation[[j]]$output_name
    ems_this_wave$expectation[[j]]$disc <- list(external = ext_int_usage[[em_name]],
                                                internal = disc_int_usage[[em_name]])
  }
  # Validation diagnostics: first check misclassifications, then check comparison
  for (j in seq_along(ems_this_wave$expectation)) {
    how_many <- nrow(validation_diagnostics(ems_this_wave$expectation[[j]], targets,
                                            valid, which_diag = c('ce'), plt = FALSE))
    while(how_many > 0) {
      ems_this_wave$expectation[[j]] <- ems_this_wave$expectation[[j]]$mult_sigma(1.05)
      how_many <- nrow(validation_diagnostics(ems_this_wave$expectation[[j]], targets,
                                              valid, which_diag = c('ce'), plt = FALSE))
    }
  }
  for (j in seq_along(ems_this_wave$expectation)) {
    em_names <- purrr::map_chr(ems_this_wave$expectation, "output_name")
    how_many <- nrow(validation_diagnostics(ems_this_wave$expectation[[j]], targets,
                                            valid, which_diag = c('cd'), plt = FALSE))
    if (how_many > 0.2 * nrow(valid)) {
      em_name <- ems_this_wave$expectation[[j]]$output_name
      ems_this_wave <- subset_emulators(ems_this_wave, em_names[!em_names %in% c(em_name)])
      j <- j - 1
    }
  }
  # Append to emulator list
  em_list[[i+1]] <- ems_this_wave
  # Generate points
  points <- generate_new_design(
    em_list, length(ranges)*15, targets, nth = 1
  )
  # Save points to file
  write.csv(points, file = paste0("Documents/HPVsim/points/wave", i+1, ".csv"), row.names = FALSE)
}

