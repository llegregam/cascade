#' Improve signal
#'
#' @include cascade_defaults.R
#' @include filter_fft.R
#' @include signal_sharpening.R
#'
#' @param df Dataframe with columns 'rtime' and 'intensity'
#' @param fourier_components Fraction of Fourier components to keep for
#'   filtering. Default is 0.01 (1%). Lower values provide more smoothing.
#' @param frequency Acquisition frequency in Hz. Default is 2.
#' @param resample Resampling factor. Default is 1.
#' @param time_min Time min in minutes. Default is 0.
#' @param time_max Time max in minutes. Default is Inf.
#' @param intensity_floor Small positive value to ensure all intensities are
#'   strictly positive after shifting. Default is 0.001.
#' @param k2 K2 parameter for signal sharpening. Default is 250.
#' @param k4 K4 parameter for signal sharpening. Default is 1250000.
#' @param sigma Sigma parameter for signal sharpening. Default is 0.05.
#' @param smoothing_width Smoothing width for signal sharpening. Default is 8.
#'
#' @return A dataframe with improved signal
#'
#' @examples NULL
improve_signal <-
  function(
    df,
    fourier_components = cascade_defaults$fourier_components,
    frequency = cascade_defaults$frequency,
    resample = cascade_defaults$resample,
    time_min = cascade_defaults$time_min,
    time_max = cascade_defaults$time_max,
    intensity_floor = cascade_defaults$intensity_floor,
    k2 = cascade_defaults$k2,
    k4 = cascade_defaults$k4,
    sigma = cascade_defaults$sigma,
    smoothing_width = cascade_defaults$smoothing_width
  ) {
    ## Ensure intensity values are strictly positive
    ## Shift by the absolute minimum value plus a small floor
    min_intensity <- min(df$intensity, na.rm = TRUE)
    shift_amount <- if (min_intensity <= 0) {
      abs(min_intensity) + intensity_floor
    } else {
      0
    }

    df_fourier <- df |>
      tidytable::mutate(intensity = intensity + shift_amount) |>
      tidytable::mutate(
        intensity_fourier = filter_fft(
          x = intensity,
          components = fourier_components
        )
      )

    f <- stats::approxfun(
      x = df_fourier$rtime,
      y = df_fourier$intensity_fourier
    )

    time <- seq(
      from = time_min,
      to = min(max(df_fourier$rtime), time_max),
      by = 1 / (frequency * 60 * resample)
    )

    intensity <- f(seq(
      from = time_min,
      to = min(max(df_fourier$rtime), time_max),
      by = 1 / (frequency * 60 * resample)
    ))

    intensity_sharpened <- signal_sharpening(
      time = time,
      intensity = intensity,
      k2 = k2,
      k4 = k4,
      sigma = sigma,
      smoothing_width = smoothing_width
    )

    ## The signal sharpening function removes the first 4 points due to
    ## derivative calculations and smoothing operations
    trim_start <- 5
    return(data.frame(
      "rtime" = time[trim_start:length(time)],
      "intensity" = intensity_sharpened
    ))
  }
