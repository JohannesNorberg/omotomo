#' @export
omotomo_plot <- function(tomo, t_ref, extra_label = NULL, verbose = FALSE) {
  plot_ne_panels(tomo = tomo, t_ref = t_ref, extra_label = extra_label, verbose = verbose)
  plot_vtec_panels(tomo = tomo, t_ref = t_ref, extra_label = extra_label, verbose = verbose)
  plot_residual_panels(tomo = tomo, t_ref = t_ref, extra_label = extra_label, verbose = verbose)    
  plot_gnss_panel(tomo = tomo, t_ref = t_ref, extra_label = extra_label, verbose = verbose)
#  plot_gnss_vtec_panel(tomo = tomo, t_ref = t_ref, extra_label = extra_label, verbose = verbose)
  plot_parameter_panels(tomo = tomo, t_ref = t_ref, extra_label = extra_label, verbose = verbose)
  plot_altimeter_panel(tomo = tomo, t_ref = t_ref, extra_label = extra_label, verbose = verbose)
  plot_profile_panels(tomo = tomo, t_ref = t_ref, extra_label = extra_label, verbose = verbose)  
  plot_parameter_value_panels(tomo, t_ref = t_ref, extra_label = extra_label, verbose = verbose)  
}

plot_parameter_value_panels <- function(
  tomo, t_ref, extra_label = NULL, verbose = FALSE,
  params = c("NmF2", "hmF2", "NmE", "hmE", "B0", "B1", "scaleF2"),
  param_limits = list(
    NmF2_pred     = tomo$plotting$ne_lim,
    hmF2_pred     = c(1.8e5,   4.50e5),
    NmE_pred      = c(0.00e11, 10e11),
    hmE_pred      = c(0.8e5,   1.50e5),
    B0_pred       = c(40e3,    200e3),
    B1_pred       = c(0,       4.0),
    scaleF2_pred  = c(10e3,    100e3),
    NmE_oval      = c(0.00e11, 5e11),
    NmE_sd        = c(0.00e11, 10e11)
  )
) {
  verbose_message("[plot_parameter_value_panels]", verbose)

  peaks <- add_inv_modip_coords_auto_layers(tomo$obs_ne_peaks, verbose = verbose)
  params_in_use <- params[params %in% names(peaks)]

  # Helper: pick limits for a raw param or its *_pred fallback
  get_limits <- function(p) {
    if (!is.null(param_limits[[p]])) return(param_limits[[p]])
    pred_name <- paste0(p, "_pred")
    if (!is.null(param_limits[[pred_name]])) return(param_limits[[pred_name]])
    NULL
  }

  # Call the plotting function ONCE per param, passing limits to the right argument.
  plots <- lapply(params_in_use, function(p) {
    lim <- get_limits(p)
    lim_args <- switch(p,
      "NmE"     = list(NmE_lim = lim),
      "hmE"     = list(hmE_lim = lim),
      "NmF2"    = list(NmF2_lim = lim),
      "hmF2"    = list(hmF2_lim = lim),
      "B0"      = list(B0_lim = lim),
      "B1"      = list(B1_lim = lim),
      "scaleF2" = list(scaleF2_lim = lim),
      list()
    )

    do.call(
      plot_iono_param_solar_invmodip,
      c(list(
          tbl   = peaks,
          z_var = p,
          long_lim = tomo$plotting$long_lim,
          lat_lim  = tomo$plotting$lat_lim
        ),
        lim_args
      )
    )
  })

  base_dir <- tomo$utils$results_directory
  out_dir  <- file.path(base_dir, "plots", "parameter_value_panels")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  t_stamp <- resolve_timestamp(t_ref = t_ref, tomo)$stamp
  lbl <- if (is.null(extra_label)) "" else as.character(extra_label)
  outfile <- file.path(out_dir, paste0("parameter_value_panels_", t_stamp, lbl, ".png"))

  ncol <- 3
  nrow <- ceiling(length(plots) / ncol)
  combined <- patchwork::wrap_plots(plots, ncol = ncol)


  grDevices::png(
    filename = outfile,
    width = ncol * 4, height = nrow * 4, units = "in", res = 200,
    type = "cairo"
  )
  print(combined)
  grDevices::dev.off()

}

plot_altimeter_panel <- function(tomo, t_ref, extra_label = NULL, z_var = "tec", verbose = FALSE) {
  vtbl <- tomo$obs_vtec

  # Early exit if table is NULL or empty
  if (is.null(vtbl) || nrow(vtbl) == 0L) {
    verbose_message("[plot_altimeter_panel] No data to plot (NULL or empty table).", verbose)
    return(invisible(NULL))
  }

  verbose_message("[plot_altimeter_panel]", verbose)

  # Base results directory from tomo
  base_dir <- tomo$utils$results_directory

  # First-level folder: "plots"
  plots_dir <- file.path(base_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
    verbose_message(glue::glue("[init] Created directory: {plots_dir}"), verbose)
  }

  # Second-level folder: "altimeter_panels"
  alt_panel_dir <- file.path(plots_dir, "altimeter_panels")
  if (!dir.exists(alt_panel_dir)) {
    dir.create(alt_panel_dir, recursive = TRUE)
    verbose_message(glue::glue("[init] Created directory: {alt_panel_dir}"), verbose)
  }

  # Timestamp label
  t_stamp <- resolve_timestamp(t_ref = t_ref, tomo)$stamp
  extra   <- if (is.null(extra_label)) "" else as.character(extra_label)

  # Build plot
  p <- plot_altimeter_values(
    tbl    = vtbl,
    z_var  = z_var,
    lat_lim  = tomo$plotting$lat_lim,
    long_lim = tomo$plotting$long_lim,
    z_lim    = tomo$plotting$tec_lim,
    verbose  = verbose
  )

  # Save PNG
  outfile <- file.path(alt_panel_dir, paste0("altimeter_panel_", t_stamp, extra, ".png"))
  grDevices::png(outfile, width = 1200, height = 1200, res = 200)
  print(p)
  grDevices::dev.off()

}

plot_parameter_panels <- function(
  tomo, t_ref, extra_label = NULL, verbose = FALSE,
  params = c("NmF2_pred", "hmF2_pred", "NmE_pred", "hmE_pred",
             "B0_pred", "B1_pred", "NmE_sd", "scaleF2_pred"),
  # --- fixed default limits per parameter (edit these as you like) ---
  param_limits = list(
    NmF2_pred     = tomo$plotting$ne_lim,
    hmF2_pred     = c(1.8e5,   4.50e5),
    NmE_pred      = c(0.00e11, 10e11),
    hmE_pred      = c(0.8e5,   1.50e5),
    B0_pred       = c(40e3,    200e3),
    B1_pred       = c(0,     4.0),
    scaleF2_pred  = c(10e3,    100e3),    
    NmE_oval      = c(0.00e11, 5e11),
    NmE_sd        = c(0.00e11, 10e11)
  ),
  ncol = 4,
  panel_size_px = 900,
  res = 200,
  # --- legend tuning (effective if image_plot_field(...) → fields::image.plot) ---
  legend_horizontal = TRUE,  # horizontal colorbar at the bottom
  legend_shrink = 0.9,       # larger → longer bar
  legend_width  = 1.1,       # larger → thicker bar
  legend_mar    = 2,         # space around legend
  legend_axis_cex = 1.0,     # legend tick label size
  # --- base-graphics cosmetics (margins/labels/ticks) ---
  par_mar = c(5.2, 5.2, 3.2, 2.8),  # bottom, left, top, right
  par_mgp = c(2.2, 0.7, 0),         # axis title, labels, line
  par_tcl = -0.35,                  # tick length (negative = ticks inside)
  axis_cex = 1.0, lab_cex = 1.05, main_cex = 1.05
) {

  verbose_message("[plot_parameter_panels]", verbose)

  # Source grid
  grid <- tomo$background_fields$prediction_grid_all
  if (is.null(grid)) stop("prediction_grid_all not found in tomo$background_fields.")

  # Keep requested order; warn for missing columns
  have <- params[params %in% names(grid)]
  missing <- setdiff(params, names(grid))
  if (!length(have)) stop("None of the requested parameters were found in the grid.")

  # Output paths
  base_dir <- tomo$utils$results_directory
  plots_dir <- file.path(base_dir, "plots"); if (!dir.exists(plots_dir)) dir.create(plots_dir, TRUE)
  out_dir <- file.path(plots_dir, "parameter_panels"); if (!dir.exists(out_dir)) dir.create(out_dir, TRUE)

  # Timestamp + filename
  t_stamp <- resolve_timestamp(t_ref = t_ref, tomo)$stamp
  lbl <- if (is.null(extra_label)) "" else as.character(extra_label)
  outfile <- file.path(out_dir, paste0("parameter_panels_", t_stamp, lbl, ".png"))

  # Helper: fetch (min, max) for a parameter; if missing, fall back to data range
  get_lim <- function(pn) {
    if (!is.null(param_limits[[pn]]) && is.numeric(param_limits[[pn]]) && length(param_limits[[pn]]) == 2) {
      return(param_limits[[pn]])
    } else {
      v <- as.numeric(grid[[pn]])
      v <- v[is.finite(v)]
      if (!length(v)) return(NULL)
      lo <- min(v, na.rm = TRUE); hi <- max(v, na.rm = TRUE)
      if (lo == hi) { pad <- max(1e-12, abs(lo) * 1e-6); lo <- lo - pad; hi <- hi + pad }
      return(c(lo, hi))
    }
  }

  # Device size computed from grid geometry
  k <- length(have)
  if (is.null(ncol) || ncol < 1) ncol <- ceiling(sqrt(k))
  nrow <- ceiling(k / ncol)
  width_px  <- ncol * panel_size_px
  height_px <- nrow * panel_size_px

  grDevices::png(outfile, width = width_px, height = height_px, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Save and restore par settings
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  # VTEC-style layout: simple mfrow grid
  graphics::par(mfrow = c(nrow, ncol),
                mar = par_mar, oma = c(0, 0, 0, 6),
                mgp = par_mgp, tcl = par_tcl,
                cex.axis = axis_cex, cex.lab = lab_cex, cex.main = main_cex)

  # Legend arguments passed via ... (if image_plot_field forwards them)
  legend_args <- list(
    horizontal    = legend_horizontal,
    legend.shrink = legend_shrink,
    legend.width  = legend_width,
    legend.mar    = legend_mar,
    axis.args     = list(cex.axis = legend_axis_cex)
  )

  # Draw each panel
  for (pn in have) {
    args <- c(
      list(grid = grid, param = pn, limits = get_lim(pn), verbose = verbose),
      legend_args
    )
    ok <- try(do.call(image_plot_field, args), silent = TRUE)
    if (inherits(ok, "try-error")) {
      # Fallback if your wrapper does not accept ...
      image_plot_field(grid = grid, param = pn, limits = get_lim(pn), verbose = verbose)
    }
    title(main = paste("Predicted", pn), cex.main = main_cex)
  }
}


plot_gnss_panel <- function(tomo, t_ref, extra_label = NULL, verbose = FALSE) {

  tbl <- tomo$obs_stec
  
  # Early exit if tbl is NULL or has no rows
  if (is.null(tbl) || nrow(tbl) == 0L) {
    verbose_message("[plot_gnss_panel] No data to plot (NULL or empty table).", verbose)
    return(invisible(NULL))
  }
  
  verbose_message("[plot_gnss_panel]", verbose)

  # Base results directory from tomo
  base_dir <- tomo$utils$results_directory

  # First-level folder: "plots"
  plots_dir <- file.path(base_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
    verbose_message(glue::glue("[init] Created directory: {plots_dir}"), verbose)
  }

  # Second-level folder: "residual_panels"
  gnss_panel_dir <- file.path(plots_dir, "gnss_panels")
  if (!dir.exists(gnss_panel_dir)) {
    dir.create(gnss_panel_dir, recursive = TRUE)
    verbose_message(glue::glue("[init] Created directory: {gnss_panel_dir}"), verbose)
  }

  t_stamp <- resolve_timestamp(t_ref = t_ref, tomo)$stamp
  
  p <- plot_gnss_locations(tomo$obs_stec, sat_posit = TRUE, pp = FALSE, lat_lim = tomo$plotting$lat_lim, long_lim = tomo$plotting$long_lim)

  grDevices::png(paste0(gnss_panel_dir, "/gnss_panel_", t_stamp, extra_label, ".png"), width = 1200, height = 1200, res = 200)
  print(p) 
  grDevices::dev.off()
}


#' Plot residual panels with a fixed middle band (±band) and variable outer bands.
#'
#' Each panel (El, Lat, Long) is built from three vertically stacked strips:
#' bottom [ymin, -band], middle [-band, band], and top [band, ymax].
#' The strips have fixed relative heights, so the ±band lines land at the same
#' pixel positions across runs. The gap between strips is controlled in pixels.
#' Only the middle strip carries the Y-axis title; the X title is drawn outside
#' the strips to keep the top/bottom strips exactly the same height.
plot_residual_panels <- function(tomo, t_ref, extra_label = NULL, verbose = FALSE,
                                 band = 10,                      # half-width of fixed middle band
                                 rel_heights = c(1, 2, 1),       # bottom / middle / top strip heights
                                 gap_px = 6,                     # gap between strips IN PIXELS
                                 x_ticks = TRUE,                 # keep x tick labels in bottom strip?
                                 xlab_rel_height = 0.08,         # external x-label band height (relative)
                                 outer_y_pad_frac = 0.03,        # pad for top/bottom strips (fraction of their range)
                                 outer_y_pad_min  = 0.5,         # min pad in TECU for top/bottom strips
                                 outer_top_pad_px = 8,           # extra blank space above the WHOLE figure (pixels)
                                 point_size = 0.5,
                                 panel_height_px = 600,          # height per panel in pixels
                                 base_size = 10) {               # theme base size

  # Data
  tbl <- tomo$obs_stec

  # Early exit if tbl is NULL or has no rows
  if (is.null(tbl) || nrow(tbl) == 0L) {
    verbose_message("[plot_residual_panels] No data to plot (NULL or empty table).", verbose)
    return(invisible(NULL))
  }
  verbose_message("[plot_residual_panels]", verbose)

  # Output directories
  base_dir <- tomo$utils$results_directory
  plots_dir <- file.path(base_dir, "plots"); if (!dir.exists(plots_dir)) dir.create(plots_dir, TRUE)
  residual_panels_dir <- file.path(plots_dir, "residual_panels"); if (!dir.exists(residual_panels_dir)) dir.create(residual_panels_dir, TRUE)

  # Timestamp and filename label
  t_stamp <- resolve_timestamp(t_ref = t_ref, tomo)$stamp
  lbl <- if (is.null(extra_label)) "" else as.character(extra_label)

  # Residuals
  tbl$resid <- tbl$stec_pred - tbl$stec

  # Units
  el_unit   <- "\u00B0"; lat_unit <- "\u00B0"; long_unit <- "\u00B0"; y_unit <- "TECU"

  # Base theme: remove vertical plot margins so gap_px fully controls spacing
  base_theme <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 0, r = 6, b = 0, l = 6)
    )

  # Convert desired pixel gap to cowplot relative height.
  core_px <- panel_height_px / (1 + xlab_rel_height)
  sum_rel <- sum(rel_heights)
  max_gap <- floor((core_px - 1) / 2)
  if (gap_px > max_gap) gap_px <- max_gap
  g_rel <- (gap_px / core_px) * sum_rel / (1 - 2 * gap_px / core_px)

  # --- UPDATED: resolve x-limits per panel using domain boundaries -------------
  resolve_xlim <- function(xvar, df) {
    if (xvar == "el") {
      return(c(0, 90))
    } else if (xvar == "lat") {
      vec <- try(tomo$domain$lat_boundary_deg, silent = TRUE)
      if (is.numeric(vec) && any(is.finite(vec))) {
        vec <- vec[is.finite(vec)]
        if (length(vec)) return(range(vec))
      }
    } else if (xvar == "long") {
      vec <- try(tomo$domain$long_boundary_deg, silent = TRUE)
      if (is.numeric(vec) && any(is.finite(vec))) {
        vec <- vec[is.finite(vec)]
        if (length(vec)) return(range(vec))
      }
    }
    # Fallback: data range
    range(df[[xvar]], na.rm = TRUE)
  }

  # Helper to build one panel from three fixed-height strips
  broken_panel <- function(df, xvar, xlab, color) {
    xlim <- resolve_xlim(xvar, df)

    y_min <- min(df$resid, na.rm = TRUE)
    y_max <- max(df$resid, na.rm = TRUE)
    y_bottom <- min(y_min, -band)
    y_top    <- max(y_max,  band)

    # Padding so points/lines don't touch the strip borders
    pad_top <- max(outer_y_pad_min,  (y_top - band)     * outer_y_pad_frac)
    pad_bot <- max(outer_y_pad_min,  (-band - y_bottom) * outer_y_pad_frac)

    d_bot <- df[df$resid <  -band, ]
    d_mid <- df[df$resid >= -band & df$resid <= band, ]
    d_top <- df[df$resid >   band, ]

    xscale <- ggplot2::scale_x_continuous(limits = xlim, expand = ggplot2::expansion(mult = 0.02))

    # Top strip (no x labels or titles)
    p_top <- ggplot2::ggplot(d_top, ggplot2::aes(x = .data[[xvar]], y = resid)) +
      ggplot2::geom_point(size = point_size, colour = color, alpha = 0.9) +
      xscale +
      ggplot2::coord_cartesian(ylim = c(band, y_top + pad_top), expand = FALSE) +
      ggplot2::geom_hline(yintercept = band, linetype = "dashed", linewidth = 0.3, alpha = 0.8) +
      base_theme +
      ggplot2::theme(
        axis.title  = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x= ggplot2::element_blank()
      )

    # Middle strip (fixed y-scale; carries the Y title)
    p_mid <- ggplot2::ggplot(d_mid, ggplot2::aes(x = .data[[xvar]], y = resid)) +
      ggplot2::geom_point(size = point_size, colour = color, alpha = 0.9) +
      xscale +
      ggplot2::coord_cartesian(ylim = c(-band, band), expand = FALSE) +
      ggplot2::geom_hline(yintercept = c(-band, band), linetype = "dashed", linewidth = 0.3, alpha = 0.8) +
      ggplot2::labs(y = paste0("STEC residuals (", y_unit, ")")) +
      base_theme +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )

    # Bottom strip (optionally shows x ticks)
    p_bot <- ggplot2::ggplot(d_bot, ggplot2::aes(x = .data[[xvar]], y = resid)) +
      ggplot2::geom_point(size = point_size, colour = color, alpha = 0.9) +
      xscale +
      ggplot2::coord_cartesian(ylim = c(y_bottom - pad_bot, -band), expand = FALSE) +
      ggplot2::geom_hline(yintercept = -band, linetype = "dashed", linewidth = 0.3, alpha = 0.8) +
      base_theme +
      ggplot2::theme(
        axis.title  = ggplot2::element_blank(),
        axis.text.x = if (x_ticks) ggplot2::element_text() else ggplot2::element_blank(),
        axis.ticks.x= if (x_ticks) ggplot2::element_line()  else ggplot2::element_blank()
      )

    # Three strips + two gaps (gaps have exact pixel height via g_rel weight)
    core <- cowplot::plot_grid(
      p_top, NULL, p_mid, NULL, p_bot,
      ncol = 1,
      rel_heights = c(rel_heights[3], g_rel, rel_heights[2], g_rel, rel_heights[1]),
      align = "v", axis = "l"
    )

    # External x-axis label band (fixed relative height; keeps strip heights intact)
    xlab_grob <- cowplot::ggdraw() + cowplot::draw_label(xlab, x = 0.5, y = 0.5)
    cowplot::plot_grid(core, xlab_grob, ncol = 1, rel_heights = c(1, xlab_rel_height))
  }

  # Build the three panels
  p_el   <- broken_panel(tbl, "el",   paste0("El (",   el_unit,   ")"), "purple")
  p_lat  <- broken_panel(tbl, "lat",  paste0("Lat (",  lat_unit,  ")"), "pink")
  p_long <- broken_panel(tbl, "long", paste0("Long (", long_unit, ")"), "cyan")

  # Align left axes and stack panels equally
  aligned <- cowplot::align_plots(p_el, p_lat, p_long, align = "v", axis = "l")
  p_all <- cowplot::plot_grid(plotlist = aligned, ncol = 1, rel_heights = c(1, 1, 1))

  # Add a small blank spacer ABOVE the whole figure to prevent clipping at the very top
  top_pad_rel <- outer_top_pad_px / (3 * panel_height_px)
  p_all_padded <- cowplot::plot_grid(NULL, p_all, ncol = 1, rel_heights = c(top_pad_rel, 1))

  # Save
  outfile <- file.path(residual_panels_dir, paste0("residual_panels_", t_stamp, lbl, ".png"))
  grDevices::png(outfile, width = 1200, height = 3L * panel_height_px, res = 200)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p_all_padded)

}

plot_profile_panels <- function(tomo, t_ref, extra_label = NULL, verbose = FALSE) {

  profiles <- tomo$profiles_raw[tomo$profiles_raw$instrument %in% c("ISR", "ionosonde"),]

  if (is.null(profiles) | nrow(profiles) == 0) {
    verbose_message("[plot_profile_panels] No profile data available", verbose)
    return(invisible(NULL))
  }
  
  verbose_message("[plot_profile_panels]", verbose)

  # Base results directory from tomo
  base_dir <- tomo$utils$results_directory

  # First-level folder: "plots"
  plots_dir <- file.path(base_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
    verbose_message(glue::glue("[init] Created directory: {plots_dir}"), verbose)
  }

  # Second-level folder: "profile_panels"
  profile_panels_dir <- file.path(plots_dir, "profile_panels")
  if (!dir.exists(profile_panels_dir)) {
    dir.create(profile_panels_dir, recursive = TRUE)
    verbose_message(glue::glue("[init] Created directory: {profile_panels_dir}"), verbose)
  }

  t_stamp <- resolve_timestamp(t_ref = t_ref, tomo)$stamp
  
  
  el_precision <- tomo$isr_param$el_precision
  az_precision <- tomo$isr_param$az_precision

  profiles[, el := el - el %% el_precision]
  profiles[, az := az - az %% az_precision]
   
  fixed_cols <- c("profile_key", "station", "instrument", "station_lat", "station_long", "station_alt", "el", "az")  
  unique_profiles <- unique(profiles[, lapply(.SD, identity), .SDcols = fixed_cols])
  fixed_stations <- unique_profiles

  rows <- ceiling(sqrt(nrow(fixed_stations)))

  p <- plot_validation(
    profiles = profiles,
    tomo_ne = tomo$ne,
    include_profiles = c("posterior","background"),#,"predictive", "ci_minus", "ci_plus"),
    t_start = t0,
    t_stop = t1,
    fixed_stations = fixed_stations,
    alt_lim = tomo$plotting$alt_lim / 1e3,
    ne_lim = tomo$plotting$ne_lim / 1e11,
    lwd = 0.5)

  grDevices::png(paste0(profile_panels_dir, "/profile_panels_", t_stamp, extra_label, ".png"), width = rows * 220, height = rows * 200, res = 150)
  print(p) 
  grDevices::dev.off()
}

plot_ne_panels <- function(tomo, t_ref, extra_label = NULL, verbose = FALSE) {

  verbose_message("[plot_ne_panels]", verbose)
  # Base results directory from tomo
  base_dir <- tomo$utils$results_directory
  
  # First-level folder: "plots"
  plots_dir <- file.path(base_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
    verbose_message(glue::glue("[init] Created directory: {plots_dir}"), verbose)
  }

  # Second-level folder: "ne_panels"
  ne_panels_dir <- file.path(plots_dir, "ne_panels")
  if (!dir.exists(ne_panels_dir)) {
    dir.create(ne_panels_dir, recursive = TRUE)
    verbose_message(glue::glue("[init] Created directory: {ne_panels_dir}"), verbose)
  }
  
  # Now ne_panels_dir is guaranteed to exist

  t_stamp <- resolve_timestamp(t_ref = t_ref, tomo)$stamp
  
  ## --- common limits and selections ---
  ne_lim    <- tomo$plotting$ne_lim  
  alt_lim   <- tomo$plotting$alt_lim 
  lat_lim   <- tomo$plotting$lat_lim 
  long_lim  <- tomo$plotting$long_lim
  long_plot <- tomo$plotting$long_plot
  lat_plot  <- tomo$plotting$lat_plot  
  
  ## save to PNG
  grDevices::png(paste0(ne_panels_dir, "/ne_panels_", t_stamp, extra_label, ".png"), width = 3600, height = 2200, res = 200)
  
  ## 2 x 4 layout
  old_par <- graphics::par(no.readonly = TRUE)
  graphics::layout(matrix(1:8, nrow = 2, byrow = TRUE),
         widths = c(1, 1, 1, 1), heights = c(1, 1))
  
  ## panel-specific margins + outer margin on the right for colorbar
  graphics::par(mar = c(4, 4, 2, 2),    # bottom, left, top, right (per panel)
      oma = c(0, 0, 0, 8))    # outer margins: extra space on the right
  
  ## --- top row: lat-alt (@ long_select = long_plot) ---
  image_plot_lat_alt_slice(
    tomo$ne, z_var = "background",
    long_select = long_plot, alt_lim = alt_lim, lat_lim = lat_lim, z_lim = ne_lim,
    legend = FALSE, xlab = "", t_ref = t_ref, verbose = verbose
  )
  
  image_plot_lat_alt_slice(
    tomo$ne, z_var = "prior",
    long_select = long_plot, alt_lim = alt_lim, lat_lim = lat_lim, z_lim = ne_lim,
    legend = FALSE, xlab = "", ylab = "", t_ref = t_ref, verbose = verbose
  )
  
  image_plot_lat_alt_slice(
    tomo$ne, z_var = "predictive",
    long_select = long_plot, alt_lim = alt_lim, lat_lim = lat_lim, z_lim = ne_lim,
    legend = FALSE, xlab = "", ylab = "", t_ref = t_ref, verbose = verbose
  )
  
  image_plot_lat_alt_slice(
    tomo$ne, z_var = "posterior",
    long_select = long_plot, alt_lim = alt_lim, lat_lim = lat_lim, z_lim = ne_lim,
    legend = TRUE, xlab = "", ylab = "",  t_ref = t_ref, verbose = verbose
  )
  
  ## --- bottom row: long-alt (@ lat_select = lat_plot) ---
  image_plot_long_alt_slice(
    tomo$ne, z_var = "background",
    lat_select = lat_plot, alt_lim = alt_lim, long_lim = long_lim, z_lim = ne_lim,
    legend = FALSE, t_ref = t_ref, verbose = verbose
  )
  
  image_plot_long_alt_slice(
    tomo$ne, z_var = "prior",
    lat_select = lat_plot, alt_lim = alt_lim, long_lim = long_lim, z_lim = ne_lim,
    legend = FALSE, ylab = "", t_ref = t_ref, verbose = verbose
  )
  
  image_plot_long_alt_slice(
    tomo$ne, z_var = "predictive",
    lat_select = lat_plot, alt_lim = alt_lim, long_lim = long_lim, z_lim = ne_lim,
    legend = FALSE, ylab = "", t_ref = t_ref, verbose = verbose
  )
  
  image_plot_long_alt_slice(
    tomo$ne, z_var = "posterior",
    lat_select = lat_plot, alt_lim = alt_lim, long_lim = long_lim, z_lim = ne_lim,
    legend = TRUE, ylab = "",  t_ref = t_ref, verbose = verbose
  )
  
  ## restore par settings
  graphics::par(old_par)
  
  ## close PNG device
  grDevices::dev.off()
}

#' Plot VTEC maps for all parameters (background, prior, predictive, posterior)
#'
#' Creates a 2x2 panel plot of VTEC maps and saves it to results directory.
#'
#' @param tomo List object containing at least:
#'   - utils$results_directory
#'   - vtec data.table with columns: lat, long, background, prior, predictive, posterior
#' @param t_ref Optional reference time (POSIXct or string) for title and filename.
#' @param verbose Logical. Print progress messages.
#'
#' @return Invisibly returns path to saved PNG file.
#' @export
plot_vtec_panels <- function(tomo, t_ref, extra_label = NULL, verbose = FALSE) {

  verbose_message("[plot_vtec_panels]", verbose)  
  # Base results directory from tomo
  base_dir <- tomo$utils$results_directory
  
  # Ensure plots/vtec_panels directories exist
  plots_dir <- file.path(base_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  vtec_panels_dir <- file.path(plots_dir, "vtec_panels")
  if (!dir.exists(vtec_panels_dir)) dir.create(vtec_panels_dir, recursive = TRUE)
  
  # Resolve timestamp
  t_stamp <- resolve_timestamp(t_ref = t_ref, tomo)$stamp
  
  # Output filename
  outfile <- file.path(vtec_panels_dir, paste0("vtec_panels_", t_stamp, extra_label, ".png"))
  
  # --- common plotting limits ---
  lat_lim  <- tomo$plotting$lat_lim
  long_lim <- tomo$plotting$long_lim
  tec_lim  <- tomo$plotting$tec_lim
  
  # Open PNG device
  grDevices::png(outfile, width = 1600, height = 1400, res = 200)
  
  # Layout 2x2
  old_par <- graphics::par(no.readonly = TRUE)
  graphics::layout(matrix(1:4, nrow = 2, byrow = TRUE),
         widths = c(1, 1), heights = c(1, 1))
  
  graphics::par(mar = c(4, 4, 2, 2), oma = c(0, 0, 0, 6))
  
  # --- four panels ---
  plot_vtec_map(tomo$vtec, z_var = "background",
                lat_lim = lat_lim, long_lim = long_lim, z_lim = tec_lim,
                legend = TRUE, t_ref = t_ref, verbose = verbose,
                xlab = "", ylab = "")
  
  plot_vtec_map(tomo$vtec, z_var = "prior",
                lat_lim = lat_lim, long_lim = long_lim, z_lim = tec_lim,
                legend = TRUE, t_ref = t_ref, verbose = verbose,
                xlab = "", ylab = "")
  
  plot_vtec_map(tomo$vtec, z_var = "predictive",
                lat_lim = lat_lim, long_lim = long_lim, z_lim = tec_lim,
                legend = TRUE, t_ref = t_ref, verbose = verbose,
                xlab = "", ylab = "")
  
  plot_vtec_map(tomo$vtec, z_var = "posterior",
                lat_lim = lat_lim, long_lim = long_lim, z_lim = tec_lim,
                legend = TRUE, t_ref = t_ref, verbose = verbose,
                xlab = "", ylab = "")
                  
  # Restore par
  graphics::par(old_par)
  
  grDevices::dev.off()
  verbose_message(glue::glue("[plot_results_vtec] Saved: {outfile}"), verbose)
}



plot_gnss_vtec_panel <- function(tomo, t_ref, extra_label = NULL, verbose = FALSE) {
  
  tbl <- tomo$obs_stec
  
  # Early exit if tbl is NULL or has no rows
  if (is.null(tbl) || nrow(tbl) == 0L) {
    verbose_message("[plot_gnss_vtec_panel] No data to plot (NULL or empty table).", verbose)
    return(invisible(NULL))
  }
  
  verbose_message("[plot_gnss_vtec_panel]", verbose)
  
  # Base results directory from tomo object
  base_dir <- tomo$utils$results_directory
  
  # First-level folder: "plots"
  plots_dir <- file.path(base_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
    verbose_message(glue::glue("[init] Created directory: {plots_dir}"), verbose)
  }
  
  # Second-level folder: "gnss_vtec_panels"
  vtec_panel_dir <- file.path(plots_dir, "gnss_vtec_panels")
  if (!dir.exists(vtec_panel_dir)) {
    dir.create(vtec_panel_dir, recursive = TRUE)
    verbose_message(glue::glue("[init] Created directory: {vtec_panel_dir}"), verbose)
  }
  
  # Resolve timestamp for filename
  t_stamp <- resolve_timestamp(t_ref = t_ref, tomo)$stamp

  lat_lim  <- tomo$plotting$lat_lim
  long_lim <- tomo$plotting$long_lim
  tec_lim  <- tomo$plotting$tec_lim

  # Generate the plot using plot_gnss_values instead of plot_gnss_locations
  p <- plot_gnss_values(
    tbl          = tbl,
    z_var        = "vtec",
    point_size   = 1.5,
    alpha        = 0.8,
    map_fill     = "gray95",
    map_border   = "gray70",
    expand_frac  = 0.05,
    lat_lim      = lat_lim,
    long_lim     = long_lim,
    z_lim        = tec_lim, # Use limits from tomo settings
    legend_title = NULL,
    na_rm        = TRUE
  )
  
  # Construct the full filename
  filename <- paste0(vtec_panel_dir, "/gnss_vtec_panel_", t_stamp, extra_label, ".png")
  
  # Save the plot
  grDevices::png(filename, width = 1200, height = 1200, res = 200)
  print(p) 
  grDevices::dev.off()
}