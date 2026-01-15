#' Save numeric reconstruction results as JSON (new tomo$ne format)
#'
#' Writes a JSON file containing voxel-wise electron density values on the
#' reconstruction grid, using the table format:
#'
#' tomo$ne: data.table with columns:
#'   lat, alt (m), long, background, posterior
#'
#' The JSON output uses:
#'   ne            = posterior / 1e11       (10^11 / m^3)
#'   background_ne = background / 1e11      (10^11 / m^3)
#'
#' @param tomo A list with at least:
#'   - tomo$ne (data.table)
#'   - tomo$domain$t0, tomo$domain$t1 (POSIXct or character)
#'   - tomo$style$latlim, tomo$style$altlim (km), tomo$style$longlim
#'   - tomo$utils$results_directory, tomo$utils$label_t
#' @param rm_negat Logical. If TRUE, negative posterior/background are set to zero.
#' @param status_col Optional quality indicator stored to JSON.
#'
#' @export
save_json_data <- function(tomo, rm_negat = TRUE, status_col = NULL) {

  if (is.null(tomo$ne)) {
    rlang::abort("tomo$ne is missing.", class = "tomoscand_error")
  }

  ne_dt <- data.table::as.data.table(tomo$ne)

  required_cols <- c("lat", "alt", "long", "posterior", "background")
  missing_cols <- setdiff(required_cols, names(ne_dt))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      message = sprintf(
        "tomo$ne is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      class = "tomoscand_error"
    )
  }

  # ------------------------------------------------------------------
  # Optional negative clipping
  # ------------------------------------------------------------------
  if (isTRUE(rm_negat)) {
    ne_dt[is.finite(posterior) & posterior < 0,  posterior  := 0]
    ne_dt[is.finite(background) & background < 0, background := 0]
  }

  # ------------------------------------------------------------------
  # Domain filtering (style limits)
  # style$altlim is in km in the old code
  # ------------------------------------------------------------------
  if (!is.null(tomo$plotting$lat_lim) || !is.null(tomo$plotting$alt_lim) || !is.null(tomo$plotting$long_lim)) {

    latlim  <- tomo$plotting$lat_lim
    longlim <- tomo$plotting$long_lim
    altlim_m <- tomo$plotting$alt_lim * 1e3
  
    V_indices <- (ne_dt$lat  >= latlim[1]) &
                 (ne_dt$lat  <= latlim[2]) &
                 (ne_dt$alt  >= altlim_m[1]) &
                 (ne_dt$alt  <= altlim_m[2]) &
                 (ne_dt$long >= longlim[1]) &
                 (ne_dt$long <= longlim[2])
  
    ne_sub <- ne_dt[V_indices]
    
  } else {ne_sub <- ne_dt}
  # ------------------------------------------------------------------
  # Build JSON payload
  # Posterior/background assumed to be in m^-3
  # Convert to 10^11 / m^3
  # ------------------------------------------------------------------
  result_file <- list()

  result_file$status  <- status_col
  result_file$t_start <- tomo$domain$t0
  result_file$t_end   <- tomo$domain$t1

  result_file$lat  <- ne_sub$lat
  result_file$alt  <- ne_sub$alt / 1e3
  result_file$long <- ne_sub$long

  result_file$ne            <- ne_sub$posterior  / 1e11
  result_file$background_ne <- ne_sub$background / 1e11

  result_file$README <- paste(
    "t_start & t_end (UTC): Start and end times of the epoch used for the reconstruction",
    "status: Optional quality indicator",
    "lat & long (WGS84 decimal degrees): Latitude and longitude coordinates",
    "alt (kilometres): Altitude",
    "ne (10^11/m^3): Reconstructed electron density (posterior)",
    "background_ne (10^11/m^3): Background electron density",
    sep = "; "
  )

  # ------------------------------------------------------------------
  # Write JSON
  # ------------------------------------------------------------------

  t_stamp <- tomoscand::resolve_timestamp(t_ref = tomo$t1, tomo)$stamp

  if (is.null(tomo$utils$results_directory)) {
    rlang::abort("tomo$utils$results_directory and tomo$utils$label_t must be defined.", class = "tomoscand_error")
  }

  out_dir <- file.path(tomo$utils$results_directory, "numeric_result")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_path <- file.path(out_dir, paste0("numeric_result_", t_stamp, ".json"))

  jsonlite::write_json(
    result_file,
    path = out_path,
    auto_unbox = TRUE,
    digits = NA
  )

  invisible(out_path)
}