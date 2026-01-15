#' Run a full ionospheric tomography scan
#'
#' The `omotomo()` function orchestrates a full workflow for ionospheric
#' electron density tomography using multiple data sources (GNSS, 
#' radio occultation, ionosondes, ISR, DIDBase). It prepares the model
#' domain, processes observational data, builds background and prior fields,
#' solves the posterior electron density model, and outputs results.
#'
#' The function requires that a `config.R` file is present in the working
#' directory, which defines the global `tomo` parameter list including
#' file paths, model parameters, and analysis settings. It also requires
#' that `tomo$utils$parameter_file` exists and is valid.
#'
#' @param t0 [POSIXct]  
#'   Start time of the analysis interval.  
#' @param t1 [POSIXct]  
#'   End time of the analysis interval.  
#' @param CALIBRATE [logical] (default: `FALSE`)  
#'   If `TRUE`, initializes the prior electron density field from the
#'   background model only. If `FALSE`, uses results from the previous
#'   tomography run (if available).  
#' @param verbose [logical] (default: `FALSE`)  
#'   If `TRUE`, prints progress messages and diagnostic information.
#'
#' @details
#' The function performs the following main steps:
#' \enumerate{
#'   \item Validates inputs and loads configuration files.
#'   \item Prepares the 3D spatial grid and results directories.
#'   \item Collects and processes multi-instrument observational data:
#'         \itemize{
#'           \item GNSS (Madrigal STEC data)
#'           \item Radio occultation (STEC and profiles)
#'           \item Ionosonde profiles
#'           \item Incoherent scatter radar (ISR) profiles
#'           \item DIDBase parameters
#'         }
#'   \item Extracts ionospheric peak parameters and smooths time series.
#'   \item Builds 2D background fields using GAM interpolation in MODIP coordinates.
#'   \item Constructs 3D Chapman electron density background and prior fields.
#'   \item Performs data assimilation (Kalman filtering / Bayesian inversion) to
#'         obtain the posterior 3D electron density field.
#'   \item Computes VTEC maps and stores results in `tomo`.
#'   \item Plots results and saves outputs to the configured results directory.
#' }
#'
#' Errors are thrown with class `"omotomo_error"` if required inputs or files
#' are missing, invalid, or incorrectly specified.
#'
#' @return
#' Updates and returns the global `tomo` list, containing fields such as:
#' \itemize{
#'   \item `tomo$domain` – 2D/3D grid definitions
#'   \item `tomo$ne` – background, prior, predictive, and posterior electron density
#'   \item `tomo$vtec` – vertical TEC fields derived from electron density
#'   \item `tomo$obs_indirect` – processed indirect observations (GNSS, RO)
#'   \item `tomo$dcb` – station and satellite DCB priors/posteriors
#' }
#'
#' @seealso
#' Functions used internally, such as \code{\link{get_madrigal_gnss_analysis_data}},
#' \code{\link{filter_gnss}}, \code{\link{get_domain_grid}},
#' \code{\link{solve_posterior}}.
#'
#' @examples
#' \dontrun{
#' t0 <- as.POSIXct("2021-11-10 09:05:00", tz = "UTC")
#' t1 <- as.POSIXct("2021-11-10 09:10:00", tz = "UTC")
#' tomo <- omotomo(t0, t1, CALIBRATE = TRUE, verbose = TRUE)
#' }
#'
#' @export
omotomo <- function(t0, t1, tomo = NULL, CALIBRATE = FALSE, verbose = FALSE) {

  tomo <- tomoscand::tomoscand_init(t0,  t1, tomo, verbose)
  tomo <- omotomo_data(t0,  t1, tomo, verbose)
  tomo <- tomoscand::tomoscand_prior(t0, t1, tomo, CALIBRATE, verbose)
  tomo <- tomoscand::tomoscand_solve(tomo, copy_input = TRUE, verbose)

  tomoscand::tomoscand_save(tomo, t_ref = t1, verbose)
  #save_json_data(tomo, rm_negat = TRUE, status_col = NULL)
  omotomo_plot(tomo = tomo, t_ref = t1, extra_label = tomo$plotting$extra_label, verbose = FALSE)

  # image_plot_lat_alt_slice(
  #   tomo$ne,
  #   z_var = "posterior",
  #   long_select = 25,
  #   alt_lim = c(0, 800e3),
  #   lat_lim = c(50, 80),
  #   z_lim = c(0, 1e12),
  #   verbose = verbose
  # )
}