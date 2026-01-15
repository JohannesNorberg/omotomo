#' Read and process data for tomography scan
#'
#'
#' @param t0 [POSIXct]  
#'   Start time of the analysis interval.  
#' @param t1 [POSIXct]  
#'   End time of the analysis interval.  
#' @param verbose [logical] (default: `FALSE`)  
#'   If `TRUE`, prints progress messages and diagnostic information.
#'
#' @details
#' The function performs the following main steps:
#' \enumerate{
#'   \item Collects and processes multi-instrument observational data:
#'         \itemize{
#'           \item GNSS (Madrigal STEC data)
#'           \item Radio occultation (STEC and profiles)
#'           \item Ionosonde profiles
#'           \item Incoherent scatter radar (ISR) profiles
#'           \item DIDBase parameters
#'         }
#'   \item Extracts ionospheric peak parameters and smooths time series.
#' }
#'
#' Errors are thrown with class `"tomoscand_error"` if required inputs or files
#' are missing, invalid, or incorrectly specified.
#'
#' @return
#' Updates and returns the global `tomo` list, containing fields such as:
#' \itemize{
#'   \item `tomo$obs_ne_peaks - Ionospheric peak parameters
#'   \item `tomo$obs_stec - Slant TEC measurement (ground-based GNSS, Radio occultation)
#'   \item `tomo$obs_vtec - Vertica TEC measurements (Altimeter Jason/Topex)    
#' }
#'
#' @export
omotomo_data <- function(t0, t1, tomo = NULL, verbose = FALSE) {


  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # Data
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Madrigal GNSS data
  # ------------------------------------------------------------------------------

  if (isTRUE(tomo$gnss_param$USE)) {
    gnss_stec_raw <- get_omotomo(
      dir_path = tomo$gnss_param$directory,
      t_start = t0,
      t_stop  = t1,
      sat_modeling_alt = tomo$gnss_param$sat_modelling_alt,      
      file_span_min = 15,
      verbose = TRUE,
      lat_filt  = tomo$domain$lat_lim,
      long_filt = tomo$domain$long_lim,
      el_filt   = c(0, 90),
      rm_n = 5
    )
  } else {gnss_stec_raw <- data.table()}

  gnss_stec_filt    <- filter_gnss(tbl = gnss_stec_raw, t_start = t0, t_stop = t1, elevation_min =  tomo$gnss_param$limit_elevation, rm_negat = tomo$gnss_param$rm_negative, systems = tomo$gnss_param$systems, verbose = verbose)

  gnss_stec_reduced <- thin_gnss_per_arc(gnss_stec_filt, n_per_min = tomo$gnss_param$n_per_min, copy_input = FALSE)
  gnss_stec         <- add_modelling_error(gnss_stec_reduced, systems = tomo$gnss_param$systems, modelling_error = tomo$gnss_param$modelling_error, verbose = verbose)

  gnss_stec <- get_plasma_range(gnss_stec, check_los = FALSE, check_earth_block = FALSE, copy_input = FALSE, sat_alt_map = tomo$gnss_param$sat_alt, earth_radius_m = tomo$domain$earth_radius, verbose = verbose)
  gnss_stec <- get_plasma_pierce_points(gnss_stec, h_ps_m = tomo$plasma_param$plasma_pp_alt, earth_radius_m = tomo$domain$earth_radius,copy_input = FALSE, verbose = verbose)

  # plot_gnss_values(tbl = gnss_stec, z_var = "vtec", z_lim = tomo$plotting$tec_lim, na_rm = TRUE)
  # print(plot_gnss_locations(gnss_stec_raw, sat_posit = FALSE, pp = FALSE))#, lat_lim = c(40, 90), long_lim = c(-20, 55)))
  # print(plot_gnss_locations(gnss_stec_raw, sat_posit = FALSE, pp = TRUE))#, lat_lim = c(40, 90), long_lim = c(-20, 55)))  
  # print(plot_gnss_locations(gnss_stec_raw, sat_posit = TRUE, pp = FALSE))#, lat_lim = c(40, 90), long_lim = c(-20, 55)))    
  # print(plot_gnss_locations(gnss_stec_filt, sat_posit = TRUE, pp = TRUE))#, lat_lim = c(40, 90), long_lim = c(-20, 55)))
  # print(plot_gnss_locations(gnss_stec_reduced, sat_posit = TRUE, pp = TRUE))#, lat_lim = c(40, 90), long_lim = c(-20, 55)))
  # print(plot_gnss_locations(gnss_stec, sat_posit = TRUE, pp = TRUE))#, lat_lim = c(40, 90), long_lim = c(-20, 55)))
  # print(plot_gnss_curve(gnss_stec_raw, x_var = "el", y_var = "stec", plot_key = unique(gnss_stec_raw$dcb_key)[1] ))
  # print(plot_gnss_curve(gnss_stec_filt, x_var = "el", y_var = "stec", plot_key = unique(gnss_stec_filt$dcb_key)[1] ))
  # print(plot_gnss_curve(gnss_stec, x_var = "el", y_var = "stec", plot_key = unique(gnss_stec$dcb_key)[1] ))

  # ------------------------------------------------------------------------------
  # Radio occultation stec data
  # ------------------------------------------------------------------------------
  if (isTRUE(tomo$ro_param$USE)) {
    ro_stec_raw <- get_ro_stec_data(t_start          = t0 - tomo$ro_param$times$back,
                                    t_stop           = t1 + tomo$ro_param$times$forward,
                                    base_dir         = tomo$ro_param$directory,
                                    POD_only         = tomo$ro_param$POD_only, 
                                    sat_modeling_alt = tomo$ro_param$sat_modelling_alt,
                                    independent_dcb_for_each_arc = tomo$ro_param$independent_dcb_for_each_arc,
                                    verbose          = verbose)
  } else {ro_stec_raw <- data.table()}  

  ro_stec_filt <- filter_gnss(ro_stec_raw, t_start = t0, t_stop = t1, elevation_min = tomo$ro_param$limit_elevation, rm_negat = tomo$ro_param$rm_negative, lat_lim  = tomo$domain$lat_lim, long_lim = tomo$domain$long_lim, verbose = verbose)

  ro_stec_w_sd <- estimate_stec_sd_from_snr(ro_stec_filt, copy_input = TRUE, verbose = verbose)

  ro_stec_one_per_arc <- thin_gnss_per_arc(ro_stec_w_sd, n_per_min = tomo$ro_param$n_per_min, copy_input = FALSE)

  ro_stec <-  add_modelling_error_ro(ro_stec_one_per_arc, modelling_error = tomo$ro_param$modelling_error, copy_input = TRUE, verbose = verbose) 

  # print(plot_ro_measurements_map(tbl = ro_stec_raw))
  # print(plot_ro_measurements_map(tbl = ro_stec_filt))  
  # print(plot_ro_measurements_map(tbl = ro_stec_avg))    
  # print(plot_ro_measurements_map(tbl = ro_stec))      
  # plot_gnss_curve(ro_stec, x_var = "el", y_var = "stec", plot_key = unique(ro_stec_avg$dcb_key)[1] )


  # ------------------------------------------------------------------------------
  # Satellite altimeter vtec data
  # ------------------------------------------------------------------------------

  if (isTRUE(tomo$altimeter_param$USE)) {
    altimeter_tec_raw <- get_madrigal_altimeter_data(path = tomo$altimeter_param$directory,
                                                     t_start  = t0,# -12*60*60,
                                                     t_stop   = t1,# +12*60*60,
                                                     lat_lim  = tomo$domain$lat_lim,
                                                     long_lim = tomo$domain$long_lim
#                                                      lat_lim  = c(-90, 90),
#                                                      long_lim = c(-180, 180)
                                                   )
  } else{altimeter_tec_raw <- data.table()}                                                 
  altimeter_tec <- add_altimeter_alt_column(vtec_tbl = altimeter_tec_raw, 
                                            altimeter_param = tomo$altimeter_param,
                                            verbose     = verbose)
  # plot_altimeter_values( altimeter_tec, lat_lim = c(-90, 90), long_lim = c(-180, 180), time_range_location = "subtitle", show_time_range = TRUE, verbose = verbose)
  # plot(x = altimeter_tec$t[altimeter_tec$sat_code == "j3"], y = altimeter_tec$tec[altimeter_tec$sat_code == "j3"], type = 'b')

  # ------------------------------------------------------------------------------
  # Ionosonde data
  # ------------------------------------------------------------------------------
  
  # Read ionosonde meta data table
  tomo$ionosonde_param$meta <- transform(
    as.data.frame(tomoscand::ionosondes),
    station_alt = tomo$ionosonde_param$station_alt
  )

  # Read ionosonde profiles
  ionosonde_profile_raw <- get_ionosonde_profile_data(
    base_dir = tomo$ionosonde_param$directory,
    ionosonde_tbl = tomo$ionosonde_param$meta,
    t_start = t0 - tomo$ionosonde_param$times$back,
    t_stop = t1 + tomo$ionosonde_param$times$forward,
    exclude_station = tomo$ionosonde_param$times$exclude_station,
    verbose = verbose)
      
  ionosonde_profile_filt1 <- filter_profile_values(ionosonde_profile_raw, range_max = tomo$ionosonde_param$range_max, alt_lim = c(80e3, 800e3), ne_na_rm = tomo$ionosonde_param$ne_na_rm, verbose = verbose)
  ionosonde_profile_filt  <- filter_profile_outliers_by_sd(ionosonde_profile_filt1,  ne_sigma = tomo$ionosonde_param$ne_sigma, sd_sigma = tomo$ionosonde_param$sd_sigma, verbose = verbose)
  ionosonde_profile_peaks <- find_profile_peak_parameters(ionosonde_profile_filt, F_max_alt = tomo$ionosonde_param$F_max_alt, EF_boundary_alt = tomo$ionosonde_param$EF_boundary_alt, ne_max_limit = tomo$ionosonde_param$ne_max_limit, verbose = verbose)

  # i <- 1
  # i<- i + 1
  # key <- unique(ionosonde_profile_peaks$profile_key)[i]
  # plot_profiles_time_alt(ionosonde_profile_peaks, plot_key = key, hm_plot = TRUE, alt_res = 5000, alt_lim = c(80e3, 450e3), z_lim = c(0, 20e11), t_res = 60 * 15)#, alt_ax = tomo$domain$alt_ax, alt_grid = tomo$domain$alt_grid)
  # 
  # t_p <- t0
  # t_p <- t_p + 60 * 15 + 1
  # plot_profile(ionosonde_profile_peaks, z_var = "ne", plot_key = key, t_plot = t_p, alt_lim = c(80e3, 800e3), z_lim = c(0, 15e11), verbose = verbose)
  
  # ------------------------------------------------------------------------------
  # ISR data
  # ------------------------------------------------------------------------------
  isr_profile_raw <- get_madrigal_isr_data(
    base_dir = tomo$isr_param$directory,
    t_start = t0 - tomo$isr_param$times$back,
    t_stop = t1 + tomo$isr_param$times$forward,
    el_precision = tomo$isr_param$el_precision,
    az_precision = tomo$isr_param$az_precision,
    exclude_station = tomo$isr_param$exclude_station,
    verbose = verbose)
  
  isr_profile_filt1 <- filter_profile_outliers_by_sd(isr_profile_raw, verbose = verbose,  ne_sigma = 5, sd_sigma = 2)
  isr_profile_filt  <- filter_profile_values(tbl = isr_profile_filt1, range_max = tomo$isr_param$range_max, sd_max = tomo$isr_param$sd_max, alt_lim = tomo$isr_param$alt_lim, ne_na_rm = tomo$isr_param$ne_na_rm, sd_na_rm = tomo$isr_param$sd_na_rm, ne_lim = tomo$isr_param$ne_lim, verbose = verbose)
  isr_profile_peaks <- find_profile_peak_parameters(isr_profile_filt, F_max_alt = tomo$isr_param$F_max_alt, verbose = verbose)

  # key <- unique(isr_profile_peaks$profile_key)[1]
  # t_p <- t1
  # t_p <- t_p + 60 * 5
  # plot_profile(isr_profile_peaks, z_var = "ne", plot_key = key, t_plot = t_p, alt_lim = c(80e3, 800e3))
  # plot_profiles_time_alt(tbl = isr_profile_filt, hm_plot = TRUE, plot_key = key, z_var = "ne", alt_res = 5000, alt_lim = c(0e3, 800e3), t_res = 120)#, alt_ax = tomo$domain$alt_ax, alt_grid = tomo$domain$alt_grid)

  # print(plot_profile_locations(isr_profile_peaks, lat_lim = NULL, long_lim = NULL, verbose = verbose))    
  # profile_peaks <- rbindlist(list(isr_profile_peaks, ionosonde_profile_peaks), fill = TRUE)
  # print(plot_profile_locations(profile_peaks, lat_lim = NULL, long_lim = NULL, verbose = verbose))
    
  # ------------------------------------------------------------------------------
  # Radio occultation Ne profile data
  # ------------------------------------------------------------------------------
  ro_profile_raw  <- get_ro_profile_data(t_start     = t0 - tomo$ro_ionden_param$times$back,
                                         t_stop      = t1 + tomo$ro_ionden_param$times$forward,
                                         base_dir    = tomo$ro_ionden_param$directory,
                                         lat_lim     = tomo$ro_ionden_paramlat_lim,  
                                         long_lim    = tomo$ro_ionden_paramlong_lim,
                                         RO_use_only = tomo$ro_ionden_paramRO_use_only,
                                         verbose     = verbose)

  ro_profile_filt   <- filter_profile_values(ro_profile_raw, range_max = tomo$ro_ionden_param$range_max, alt_lim = tomo$ro_ionden_param$alt_lim, ne_lim = tomo$ro_ionden_param$ne_lim, ne_na_rm = tomo$ro_ionden_param$ne_na_rm, verbose = verbose)
  ro_profile_peaks0 <- find_profile_peak_parameters(profiles = ro_profile_filt, F_max_alt = tomo$isr_param$F_max_alt, verbose = verbose)
  ro_profile_peaks  <- filter_profile_outliers_by_sd(ro_profile_peaks0,  ne_sigma = tomo$ro_ionden_param$ne_sigma, verbose = verbose)

  # i <- 0
  # i <- 1 + i
  # key <- unique(ro_profile_peaks$profile_key)[i]
  # t_p <- t1
  # t_p <- t_p + 60 * 60
  # print(plot_profile(ro_profile_peaks, plot_key = key, t_plot = t_p, alt_lim = c(0, 500e3), z_lim = c(0, 2e12)))
  # print(plot_profile_locations(ro_profile_peaks, lat_lim = NULL, long_lim = NULL))
  
  # ------------------------------------------------------------------------------
  # DIDBase parameter data
  # ------------------------------------------------------------------------------

  if (!identical(tomo$didbase_param$USE, FALSE)) {
    didbase_peaks <- get_didbase_data(
      path            = tomo$didbase_param$path,
      t_start         = t0 - tomo$didbase_param$times$back,
      t_stop          = t1 + tomo$didbase_param$times$forward,
      exclude_station = tomo$didbase_param$exclude_station, 
      verbose         = verbose
    )
  } else {didbase_peaks <- data.table()}
  # print(plot_profile_locations(didbase_peaks, lat_lim = NULL, long_lim = NULL))
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # Ionospheric profile parameters
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Parameters from groundbased instruments
  # ------------------------------------------------------------------------------

  # Extract peak parameters only
  if (!identical(tomo$isr_param$USE, FALSE)) {
    isr_peaks <- reduce_profiles_to_peaks_only(isr_profile_peaks)
  } else {isr_peaks <- data.table()}

  if (!identical(tomo$ionosonde_param$USE, FALSE)) {  
    ionosonde_peaks <- reduce_profiles_to_peaks_only(ionosonde_profile_peaks)  
  } else {ionosonde_peaks <- data.table()}

  # Combine all groundbased 
  groundbased_peaks <- rbindlist(list(didbase_peaks, isr_peaks, ionosonde_peaks), use.names = TRUE, fill = TRUE)

  # Smooth interpolation for parameter time series from constant groundbased locations
  groundbased_peaks_smooth <- smooth_peak_parameters_gam(tbl = groundbased_peaks, param_list = tomo$gam1d_param, verbose = verbose)

  # print(plot_profile_locations(groundbased_peaks, lat_lim = NULL, long_lim = NULL))    
  #
  # i <- 1 + i
  # plot_temporal_peak_evolution(
  #   raw_tbl = groundbased_peaks,
  #   smooth_tbl = groundbased_peaks_smooth,
  #   t_lim = c(min(groundbased_peaks$t), max(groundbased_peaks$t)),
  #   t_ref = t1,
  #   plot_key = unique(groundbased_peaks$profile_key)[i],
  # )
  #
  # plot_temporal_peak_evolution_all(
  #   raw_tbl = groundbased_peaks,
  #   smooth_tbl = groundbased_peaks_smooth,
  #   output_dir = "/Users/norberg/Desktop/tmp",
  #   t_lim = c(min(groundbased_peaks$t), max(groundbased_peaks$t)),
  #   t_ref = t1,
  #   ne_lim = c(0, 30),
  #   ne_lim_E = c(0, 10),  
  #   altlim = c(100, 600),
  #   dpi = 120
  # )  

  # ------------------------------------------------------------------------------
  # Parameters from radio occultation profiles
  # ------------------------------------------------------------------------------

  # Extract peak parameters only
  if (isTRUE(tomo$ro_ionden_param$USE)) {
    ro_peaks <- reduce_profiles_to_peaks_only(ro_profile_peaks)
  } else {ro_peaks <- data.table()}

  # Outlier removal for parameters with arbitrary locations
  ro_peaks_filt <- filter_ungrouped_peak_outliers(ro_peaks, param_list = tomo$ro_peak_param, verbose = verbose)

  # plot_iono_param(tbl = ro_peaks_filt, z_var = "NmF2",
  #   hmE_lim = NULL,
  #   NmE_lim = NULL,
  #   hmF2_lim = c(180, 500),
  #   NmF2_lim = c(0, 25),
  #   B0_lim = NULL,
  #   B1_lim = NULL,
  #   lat_lim = NULL,
  #   long_lim = NULL)

  # ------------------------------------------------------------------------------
  # Combine all profile parameters
  # ------------------------------------------------------------------------------

  peaks_all_filtered <- rbindlist(list(groundbased_peaks_smooth, ro_peaks_filt), fill = TRUE)

  # Remove values from start and end of peak profile parameter timeseries
  # For possible boundary effects
  peaks_filtered <- filter_parameter_timeseries_boundaries(
    tbl       = peaks_all_filtered,
    t_start   = t0,
    t_stop    = t1,
    cut_start = 60 * 60,
    cut_end   = 60 * 60,
    verbose   = verbose
  )
  
  # plot_profile_locations(peaks_all, lat_lim = NULL, long_lim = NULL)
  # 
  # plot_iono_param(tbl = peaks_all, z_var = "NmF2",
  #   hmE_lim = NULL,
  #   NmE_lim = NULL,
  #   hmF2_lim = c(180, 500),
  #   NmF2_lim = c(0, 25),
  #   B0_lim = NULL,
  #   B1_lim = NULL,
  #   lat_lim = NULL,
  #   long_lim = NULL,
  #   verbose = verbose)
  
  tomo$profiles_raw <- rbindlist(list(isr_profile_raw, ionosonde_profile_raw, ro_profile_raw), fill = TRUE)
  tomo$obs_ne_peaks_raw <- rbindlist(list(groundbased_peaks, ro_peaks), fill = TRUE)
  tomo$obs_ne_peaks <- peaks_filtered  
  tomo$obs_stec     <- rbindlist(list(gnss_stec, ro_stec), fill = TRUE)
  tomo$obs_vtec     <- altimeter_tec

  return(tomo)
}