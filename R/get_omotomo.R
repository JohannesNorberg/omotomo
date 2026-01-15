#' Parse Omotomo file end time from filename
#'
#' Omotomo filenames are expected to follow the pattern
#' `omotomo<yydoyHHMM>_par.h5`, where:
#' \itemize{
#'   \item `yy` is the 2 digit year (assumed 2000 + yy)
#'   \item `doy` is day of year (001 to 366)
#'   \item `HHMM` is hour and minute in UTC
#' }
#'
#' The parsed time corresponds to the last sample time in the file.
#'
#' @param fname Character. Filename or full path.
#' @param tz Character. Timezone used for output, default "UTC".
#'
#' @return POSIXct. Parsed file end time.
#' @export
parse_omotomo_filename_time <- function(fname, tz = "UTC") {
  m <- regexec("^omotomo(\\d{9})_par\\.h5$", basename(fname))
  reg <- regmatches(basename(fname), m)[[1]]
  if (length(reg) != 2L) {
    rlang::abort(
      message = sprintf("Filename does not match expected pattern: %s", fname),
      class = "tomoscand_error"
    )
  }

  stamp <- reg[2]
  yy   <- as.integer(substr(stamp, 1, 2))
  doy  <- as.integer(substr(stamp, 3, 5))
  HH   <- as.integer(substr(stamp, 6, 7))
  MM   <- as.integer(substr(stamp, 8, 9))

  if (is.na(yy) || is.na(doy) || is.na(HH) || is.na(MM)) {
    rlang::abort(
      message = sprintf("Failed to parse time stamp from filename: %s", fname),
      class = "tomoscand_error"
    )
  }

  year <- 2000L + yy

  if (doy < 1L || doy > 366L) {
    rlang::abort(
      message = sprintf("Invalid DOY in filename: %s", fname),
      class = "tomoscand_error"
    )
  }
  if (HH < 0L || HH > 23L || MM < 0L || MM > 59L) {
    rlang::abort(
      message = sprintf("Invalid HHMM in filename: %s", fname),
      class = "tomoscand_error"
    )
  }

  t0 <- as.POSIXct(sprintf("%04d-01-01 00:00:00", year), tz = tz)
  t0 + (doy - 1L) * 86400 + HH * 3600 + MM * 60
}


#' Read GNSS LOS data from FMI Omotomo HDF5 format
#'
#' Reads per station groups from an Omotomo HDF5 file, applies station domain
#' filtering, optional day and hour filters, and returns a combined table.
#'
#' This implementation aligns the key output concepts with
#' `read_madrigal_gnss_stec` where reasonable:
#' \itemize{
#'   \item `t` is returned as POSIXct in UTC
#'   \item `dcb_key` is used as the grouping key name (renamed from older `bias_subject`)
#'   \item `pp_lat`, `pp_long`, `pp_alt` are computed for a 350 km pierce point
#'   \item `lat`, `long`, `alt` are computed for a modeling altitude `sat_modeling_alt`
#' }
#'
#' Notes:
#' \itemize{
#'   \item Station longitude is stored as `station_lon` in attributes and is renamed to `station_long`.
#'   \item Receiver altitude is stored as `station_h` in attributes and is renamed to `station_alt`.
#'   \item Elevation angles outside `el_filt` are filtered when using the `day_filt + h_filt` branch.
#'     If you want the same elevation filter always applied, add it before returning.
#' }
#'
#' @param h5_path Character. Path to the Omotomo HDF5 file.
#' @param earth_rad Numeric. Earth radius in meters. Default 6371e3.
#' @param return_as_datatable Logical. If TRUE returns a data.table, else a data.frame.
#' @param lat_filt Numeric vector length 2. Data latitude filter.
#' @param long_filt Numeric vector length 2. Data longitude filter.
#' @param el_filt Numeric vector length 2. Elevation filter in degrees. Used in the hour filter branch.
#' @param day_filt Integer vector or NULL. If not NULL, keeps only measurements from these DOYs.
#' @param h_filt Integer or NULL. If not NULL, keeps only measurements around this hour (on day_filt).
#' @param h_look_back Numeric. Hours to include before h_filt. Default 0.5.
#' @param rm_n Integer. Remove measurement arcs with fewer than rm_n samples. Default 5.
#' @param verbose Logical. Print progress messages.
#' @param sat_modeling_alt Numeric. Modeling altitude in meters for `lat,long,alt` pierce coordinates. Default 1200e3.
#'
#' @return data.table or data.frame of GNSS observations.
#' @export
#'
#' @importFrom data.table :=
#' @importFrom data.table .N
read_omotomo <- function(
  h5_path,
  earth_rad = 6371e3,
  return_as_datatable = FALSE,
  lat_filt = c(-90, 90),
  long_filt = c(-180, 180),
  el_filt = c(0, 90),
  day_filt = NULL,
  h_filt = NULL,
  h_look_back = 0.5,
  rm_n = 5,
  verbose = FALSE,
  sat_modeling_alt
) {
  h5_file <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(h5_file$close_all(), add = TRUE)

  h5_ls <- h5_file$ls(recursive = TRUE)
  stations_w_measurements <- h5_ls$name[h5_ls$dataset.maxdims > 0]

  take_variables <- c(
    "t", "stec", "stec_std", "az", "el", "satellite", "station", "x", "y", "z",
    "code", "arcid", "year", "rangetec", "dcb"
  )

  take_variables_station <- c(
    "station_lat", "station_lon", "station_h",
    "station_x", "station_y", "station_z"
  )

  h5data <- list()

  pb <- NULL
  if (!isTRUE(verbose)) pb <- utils::txtProgressBar(style = 3)

  for (i in seq_along(stations_w_measurements)) {
    station_name <- stations_w_measurements[i]

    if (isTRUE(verbose)) {
      cat("\n  Reading hdf5 data from station:", station_name, " - ",
          i, "/", length(stations_w_measurements), "\n")
    }

    station_param <- hdf5r::h5attributes(h5_file[[station_name]])[take_variables_station]
    station_lat <- station_param$station_lat
    station_long <- station_param$station_lon

    in_domain <- (lat_filt[1] <= station_lat) & (lat_filt[2] >= station_lat) &
      (long_filt[1] <= station_long) & (long_filt[2] >= station_long)

    if (!in_domain) {
      if (isTRUE(verbose)) cat("\n    station out of domain\n")
      if (!isTRUE(verbose)) utils::setTxtProgressBar(pb, i / length(stations_w_measurements))
      next
    }

    file_i <- h5_file[[station_name]][][, take_variables]

    if (is.null(day_filt)) {
      h5data[[length(h5data) + 1]] <- cbind(file_i, station_param)

    } else if (is.null(h_filt)) {
      if (isTRUE(verbose)) cat("\n    checking for doy(s)", day_filt, "\n")
      day_filt_indices <- trunc(file_i$t) %in% day_filt
      if (any(day_filt_indices)) {
        h5data[[length(h5data) + 1]] <- cbind(file_i[day_filt_indices, ], station_param)
      }

    } else {
      if (isTRUE(verbose)) cat("\n    checking for hour", h_filt, "\n")

      h_filt_indices <- (file_i$t > (as.numeric(day_filt) + (h_filt - h_look_back) / 24)) &
        (file_i$t < (as.numeric(day_filt) + (h_filt + 1) / 24)) &
        (file_i$el > el_filt[1]) &
        (file_i$el < el_filt[2])

      if (any(h_filt_indices)) {
        h5data[[length(h5data) + 1]] <- cbind(file_i[h_filt_indices, ], station_param)
      }
    }

    if (!isTRUE(verbose)) utils::setTxtProgressBar(pb, i / length(stations_w_measurements))
  }

  if (!isTRUE(verbose) && !is.null(pb)) close(pb)

  if (length(h5data) == 0L) {
    rlang::abort(
      message = "Empty h5 file (no stations passed filters).",
      class = "tomoscand_error"
    )
  }

  if (isTRUE(verbose)) cat("\n  Combine data to one data.table\n")
  h5data <- data.table::rbindlist(h5data, use.names = TRUE, fill = TRUE)

  if (nrow(h5data) == 0L) {
    rlang::abort(
      message = "Empty h5 file after filtering.",
      class = "tomoscand_error"
    )
  }

  # Standardize station column names
  data.table::setnames(h5data, "station_h", "station_alt")
  data.table::setnames(h5data, "station_lon", "station_long")
  if ("stec_std" %in% names(h5data) && !"stec_sd" %in% names(h5data)) {
    data.table::setnames(h5data, old = "stec_std", new = "stec_sd")
  }

  # Satellite type + dcb_key (compatible idea with Madrigal naming)
  if (isTRUE(verbose)) cat("\n  Construct DCB grouping key\n")
  sat_code <- substring(h5data$satellite, 1, 1)
  satellite_type <- character(length(sat_code))
  satellite_type[sat_code == "G"] <- "GPS"
  satellite_type[sat_code == "E"] <- "GALILEO"
  satellite_type[sat_code == "R"] <- "GLONASS"
  satellite_type[satellite_type == ""] <- NA_character_
  h5data[, satellite_type := satellite_type]

  # Station type (tomoscand compatible)
  h5data[, station_type := paste0("ground_", satellite_type)]
  
  h5data[, dcb_key := paste(station, satellite_type, code, sep = "__")]

  # Remove short arcs
  if (rm_n > 0) {
    if (isTRUE(verbose)) cat("\n  Removing arcs with less than rm_n =", rm_n, "\n")
    h5data[, dcb_key_count := .N, by = c("dcb_key", "arcid")]
    h5data <- h5data[dcb_key_count >= rm_n]
    h5data[, dcb_key_count := NULL]
  }

  # Combine satellite id and code if desired
  h5data[, satellite := paste(satellite, "__", code, sep = "")]

  # Time conversion: Omotomo stores `t` as DOY with fraction, with `year` column.
  if (isTRUE(verbose)) cat("\n  Time conversions\n")

  h5data[, t_num := as.numeric(
    as.POSIXct(paste(year, "01-01 00:00:00"), format = "%Y %m-%d %H:%M:%S", tz = "UTC") +
      (t - 1) * 86400
  )]
  h5data[, t := as.POSIXct(t_num, origin = "1970-01-01", tz = "UTC")]
  h5data[, t_num := NULL]

  # Coordinate conversions + ranges to modeling altitude
  if (isTRUE(verbose)) cat("\n  Coordinate conversions\n")

  xyz_at_alt <- get_xyz_at_alt(
    x0 = h5data$station_x,
    y0 = h5data$station_y,
    z0 = h5data$station_z,
    x1 = h5data$x,
    y1 = h5data$y,
    z1 = h5data$z,
    alt = earth_rad + sat_modeling_alt
  )

  h5data[, range := sqrt((x - station_x)^2 + (y - station_y)^2 + (z - station_z)^2)]

  h5data[, plasma_range := range - sqrt(
    (xyz_at_alt[, 1] - station_x)^2 +
      (xyz_at_alt[, 2] - station_y)^2 +
      (xyz_at_alt[, 3] - station_z)^2
  )]

  h5data[, c("x", "y", "z", "station_x", "station_y", "station_z") := NULL]

  latlonalt <- xyzToLatLongAltCpp(t(xyz_at_alt))
  h5data[, `:=`(
    lat = latlonalt[1, ],
    long = latlonalt[2, ],
    alt = latlonalt[3, ]
  )]

  # -------------------------------------------------------------
  # Domain filtering in lat/long/alt (sat_modeling_alt locations)
  # -------------------------------------------------------------
  if (verbose) cat("\n  Domain filtering (lat/long/alt)\n")
  
  h5data <- h5data[
    is.finite(lat) & is.finite(long) & is.finite(alt) &
      lat  >= lat_filt[1]  & lat  <= lat_filt[2] &
      long >= long_filt[1] & long <= long_filt[2]
  ]

  # Pierce points at 350 km
  if (isTRUE(verbose)) cat("\n  Piercepoint computations\n")
  pp350 <- compute_pierce_point(
    lat  = h5data$station_lat,
    long = h5data$station_long,
    el   = h5data$el,
    az   = h5data$az,
    h    = 350e3,
    earth_rad = earth_rad,
    wrap_lon = TRUE
  )
  h5data[, `:=`(
    pp_lat = pp350$lat_pp,
    pp_long = pp350$long_pp,
    pp_alt = 350e3
  )]

  # Defensive station altitude
  if (any(h5data$station_alt < 0, na.rm = TRUE)) {
    if (isTRUE(verbose)) cat("\n  Setting negative station altitudes to 1 m\n")
    h5data$station_alt[h5data$station_alt < 0] <- 1
  }

  if (isTRUE(verbose)) cat("\n  Close h5 and return\n")

  if (isTRUE(return_as_datatable)) return(h5data)
  as.data.frame(h5data)
}


#' Load Omotomo GNSS data from a directory using file time windows
#'
#' Reads one or more Omotomo HDF5 files from a directory such that the returned
#' measurements cover the requested time interval.
#'
#' The Omotomo filename time stamp (`yydoyHHMM`) is interpreted as the file end time.
#' Each file is assumed to cover a fixed interval of length `file_span_min`
#' ending at that time.
#'
#' This wrapper:
#' \itemize{
#'   \item lists Omotomo files under `dir_path`
#'   \item parses end times from filenames
#'   \item optionally checks that file end times have constant spacing
#'   \item selects files whose coverage overlaps `[t_start, t_stop]`
#'   \item reads selected files with \code{\link{read_omotomo}}
#'   \item filters rows to exactly `[t_start, t_stop]`
#' }
#'
#' @param dir_path Character. Directory containing Omotomo files.
#' @param t_start POSIXct. Start time (UTC recommended).
#' @param t_stop POSIXct. Stop time (UTC recommended).
#' @param file_span_min Numeric. File coverage length in minutes. Default 15.
#' @param lat_filt Numeric vector length 2. Data latitude filter.
#' @param long_filt Numeric vector length 2. Data longitude filter.
#' @param sat_modeling_alt Numeric. Modeling altitude in meters for `lat,long,alt` pierce coordinates. Default 1200e3.
#' @param tz Character. Timezone for filename parsing. Default "UTC".
#' @param pattern Character. Regex for selecting filenames. Default matches `omotomo<9 digits>_par.h5`.
#' @param check_constant_step Logical. If TRUE, checks if end time spacing is constant.
#' @param verbose Logical. Print progress.
#' @param ... Additional arguments forwarded to \code{\link{read_omotomo}}.
#'
#' @return data.table. Combined observations from selected files.
#'   Attribute `files_selected` contains a data.table of selected file coverage windows.
#' @export
get_omotomo <- function(
  dir_path,
  t_start,
  t_stop,
  file_span_min = 15,
  lat_filt = c(-90, 90),
  long_filt = c(-180, 180),  
  sat_modeling_alt = 1200e3,  
  tz = "UTC",
  pattern = "^omotomo\\d{9}_par\\.h5$",
  check_constant_step = TRUE,
  verbose = FALSE,
  ...
) {
  if (!inherits(t_start, "POSIXct") || !inherits(t_stop, "POSIXct")) {
    rlang::abort(
      message = "t_start and t_stop must be POSIXct.",
      class = "tomoscand_error"
    )
  }
  if (t_stop < t_start) {
    rlang::abort(
      message = "t_stop must be >= t_start.",
      class = "tomoscand_error"
    )
  }

  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    rlang::abort(
      message = sprintf("No omotomo files found in directory: %s", dir_path),
      class = "tomoscand_error"
    )
  }
  file_span_s <- as.numeric(file_span_min) * 60

  t_end_num <- vapply(files, parse_omotomo_filename_time, as.POSIXct(NA), tz = tz) 
  t_end <- as.POSIXct(t_end_num, origin = "1970-01-01", tz = tz) - file_span_s

  o <- order(t_end)
  files <- files[o]
  t_end <- t_end[o]

  t_begin <- t_end - file_span_s

#   if (isTRUE(check_constant_step) && length(t_end) >= 2L) {
#     dt <- diff(as.numeric(t_end))
#     dt_u <- sort(unique(dt))
#     if (length(dt_u) != 1L) {
#       warning("Omotomo file end times are not on a constant step.", call. = FALSE)
#       if (isTRUE(verbose)) {
#         message(sprintf("Unique step sizes (s): %s", paste(dt_u, collapse = ", ")))
#       }
#     } else {
#       if (isTRUE(verbose)) {
#         message(sprintf("Detected constant file step: %g seconds", dt_u))
#       }
#     }
#   }

  keep <- (t_end >= t_start) & (t_begin < t_stop)
  if (!any(keep)) {
    rlang::abort(
      message = "No omotomo files overlap the requested time window.",
      class = "tomoscand_error"
    )
  }

  files_sel <- files[keep]
  t_end_sel <- t_end[keep]
  t_begin_sel <- t_begin[keep]

  if (isTRUE(verbose)) {
    message(sprintf(
      "Selected %d file(s) covering %s to %s",
      length(files_sel),
      format(min(t_begin_sel), "%Y-%m-%d %H:%M:%S %Z"),
      format(max(t_end_sel), "%Y-%m-%d %H:%M:%S %Z")
    ))
  }

  out_list <- vector("list", length(files_sel))
  for (i in seq_along(files_sel)) {
    if (isTRUE(verbose)) {
      message(sprintf("Reading %d/%d: %s", i, length(files_sel), basename(files_sel[i])))
    }
    out_list[[i]] <- read_omotomo(
      h5_path = files_sel[i],
      sat_modeling_alt = sat_modeling_alt,
      return_as_datatable = TRUE,
      lat_filt = lat_filt,
      long_filt = long_filt,
      ...
    )
  }

  out <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)

  if (!("t" %in% names(out))) {
    rlang::abort(
      message = "read_omotomo output does not contain column 't'.",
      class = "tomoscand_error"
    )
  }
  if (!inherits(out$t, "POSIXct")) {
    rlang::abort(
      message = "read_omotomo output column 't' is not POSIXct. Convert it inside read_omotomo.",
      class = "tomoscand_error"
    )
  }

  out <- out[out$t >= t_start & out$t <= t_stop, ]

  attr(out, "files_selected") <- data.table::data.table(
    file = basename(files_sel),
    path = files_sel,
    t_begin = t_begin_sel,
    t_end = t_end_sel
  )

  out
}